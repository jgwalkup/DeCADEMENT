"""
Full model: sweep the enzyme Vmax-Km slope across a temperature gradient for MMRT (mmrt_Vmax_Km_slope) and Arrhenius (Vmax_Km)
kinetics. The slope scales Km only (Vmax is unchanged): MMRT Km = MMRT(Km params) x substrate-specific factor x slope + intercept;
Arrhenius Km0 = Vmax0 x slope + intercept (then scaled by its own Arrhenius term). Flat climate (Psi -0.5, FM 50), one pulse of
365 days, HSP production as in the input folder, seed-replicated communities.

Usage:
    python full_km_sweep.py ../input_e1500_enz2 ../analysis_output/km_sweep [--slopes 7 3 1 0.3 0.1] [--seeds 0 1 2]
        [--temps 4 8 ...] [--side 20] [--days 365] [--plot-only]
Writes <out_dir>/timeseries.csv (daily, per run, with the remaining C of every substrate) and the figures km_sweep_summary.png and
km_sweep_substrates.png.
"""
import argparse
import os
import shutil
import sys
import tempfile
from multiprocessing import Pool

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "src"))
sys.path.insert(0, HERE)
from grid import Grid  # noqa: E402
from initialization import initialize_data  # noqa: E402
from sweep_minimal import set_params, set_runtime  # noqa: E402

plt.rcParams.update({"font.size": 10, "axes.titlesize": 11, "axes.labelsize": 10, "legend.fontsize": 9, "figure.dpi": 150})
TEMPS = [4, 8, 12, 16, 20, 24, 28, 32, 36, 40]
KIN = {"MMRT": (1, "mmrt_Vmax_Km_slope"), "Arrhenius": (0, "Vmax_Km")}
UPTAKE_PARAM = {"MMRT": "Uptake_mmrt_Vmax_Km_slope", "Arrhenius": "Uptake_Vmax_Km"}  # uptake Km = uptake Vmax x slope


def run(job):
    src, kin, slope, temp, seed, side, days = job[:7]
    uslope = job[7] if len(job) > 7 else None  # optional uptake Vmax-Km slope
    model, pname = KIN[kin]
    tmp = tempfile.mkdtemp()
    w = os.path.join(tmp, "in")
    shutil.copytree(src, w)
    set_runtime(os.path.join(w, "runtime.txt"), dict(x=side, y=side, gridsize=side * side, pulse=1, end_time=days,
                                                     Kinetics_Model=model))
    changes = {pname: float(slope)}
    if uslope is not None:
        changes[UPTAKE_PARAM[kin]] = float(uslope)
    set_params(os.path.join(w, "parameters.csv"), changes)
    pd.DataFrame({"Temp": float(temp), "FM": 50.0, "Psi": -0.5, "Date": "1/1/00"}, index=range(days + 5)).to_csv(
        os.path.join(w, "climate.csv"))
    cwd = os.getcwd()
    os.chdir(w)
    rows = []
    try:
        np.random.seed(seed)
        rt = pd.read_csv("runtime.txt", header=None, index_col=0, sep="\t").drop("output_format").astype(float)
        data = initialize_data(rt)
        eco = Grid(rt, data)
        G = side * side
        subs = [s for s in data["Substrate_index"] if s not in ("DeadMic", "DeadEnz")]
        sub_c = lambda s: float(np.asarray(eco.Substrates.loc[s, "C"]).sum())
        sub0 = {s: sub_c(s) for s in subs}
        tot0 = sum(sub0.values())
        for i in range(days):
            eco.degradation(0, i)
            eco.uptake(0, i)
            eco.metabolism(i)
            eco.mortality(i)
            eco.reproduction(i)
            C = np.asarray(eco.Microbes["C"])
            E = np.asarray(eco.Microbes["E"])
            alive = C > 0
            n = int(alive.sum())
            rem = {s: sub_c(s) for s in subs}
            mon_C = float(np.asarray(eco.Monomers["C"]).sum()) / G
            rows.append(dict(kin=kin, slope=slope, uslope=uslope if uslope is not None else np.nan, mon_C=mon_C, Temp=temp, seed=seed, day=i + 1, biomass_C=C.sum() / G, alive_rows=n,
                             energy_to_C=E[alive].sum() / C[alive].sum() if n else np.nan, CUE=float(eco.CUE_system),
                             remaining_pct=100 * sum(rem.values()) / tot0,
                             **{f"rem_{s}": 100 * rem[s] / sub0[s] if sub0[s] > 0 else np.nan for s in subs}))
    finally:
        os.chdir(cwd)
        shutil.rmtree(tmp, ignore_errors=True)
    return pd.DataFrame(rows)


def drop_nan_runs(ts):
    """Runs in which a substrate pool became NaN (0/0 in the cellulose-lignin coupling once a cell has neither left) are dropped:
    from that day the cellulose monomer pool is NaN in every cell, so the rest of the run is not trustworthy."""
    key = ["kin", "slope", "Temp", "seed"]
    bad = ts.groupby(key)["remaining_pct"].apply(lambda x: x.isna().any())
    bad = bad[bad].index
    return ts[~ts.set_index(key).index.isin(bad)], len(bad)


def final(ts, day):
    g = ts[ts["day"] == day].groupby(["kin", "slope", "Temp"])
    a = g.mean(numeric_only=True)
    a = a[g.size() >= 2]  # at least 2 clean communities per point
    a["degraded"] = 100 - a["remaining_pct"]
    c = ts[(ts["day"] > day - 30) & (ts["day"] <= day)].groupby(["kin", "slope", "Temp"])["CUE"].mean()
    a["CUE30"] = c
    return a.reset_index()


def summary_plot(ts, out_dir, n_seeds, day):
    f = final(ts, day)
    slopes = sorted(f["slope"].unique(), reverse=True)
    cm = plt.cm.viridis(np.linspace(0.05, 0.9, len(slopes)))
    fig, axes = plt.subplots(4, 2, figsize=(11, 13), sharex=True)
    fig.suptitle(f"Full model, enzyme Vmax-Km slope sweep (Km only), day {day}, mean of {n_seeds} communities", fontsize=11,
                 fontweight="bold")
    panels = [("degraded", "Substrate degraded (% of initial C)", (0, 100)), ("biomass_C", "Biomass C per grid cell", (0, None)),
              ("CUE30", "CUE (last 30 days)", (0, 0.7)), ("energy_to_C", "Energy : C", (0, None))]
    for r, (col, lab, yl) in enumerate(panels):
        for c, kin in enumerate(KIN):
            ax = axes[r][c]
            for s, colr in zip(slopes, cm):
                d = f[(f["kin"] == kin) & (f["slope"] == s)].sort_values("Temp")
                ax.plot(d["Temp"], d[col], marker="o", ms=4, lw=1.5, color=colr, label=f"slope {s:g}")
            ax.set_ylim(*yl)
            if r == 0:
                ax.set_title(f"{kin} ({KIN[kin][1]})")
            if c == 0:
                ax.set_ylabel(lab)
            if r == len(panels) - 1:
                ax.set_xlabel("Temperature (°C)")
            ax.spines[["top", "right"]].set_visible(False)
    axes[0][0].axhline(50, color="gray", lw=0.8, ls="--")
    axes[0][1].axhline(50, color="gray", lw=0.8, ls="--")
    h, l = axes[0][0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncols=len(slopes), framealpha=0.85, bbox_to_anchor=(0.5, -0.01))
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "km_sweep_summary.png"), dpi=150, bbox_inches="tight")
    plt.close(fig)


def substrate_plot(ts, out_dir, day):
    f = final(ts, day)
    subs = [c for c in f.columns if c.startswith("rem_")]
    slopes = sorted(f["slope"].unique(), reverse=True)
    best = {k: f[(f["kin"] == k) & (f["slope"] == slopes[0])].sort_values("degraded").iloc[-1]["Temp"] for k in KIN}
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.8), sharey=True)
    fig.suptitle(f"Substrate degraded by day {day} at the temperature of maximum degradation of the baseline slope",
                 fontsize=11, fontweight="bold")
    w = 0.8 / len(slopes)
    for ax, kin in zip(axes, KIN):
        t = best[kin]
        for k, s in enumerate(slopes):
            d = f[(f["kin"] == kin) & (f["slope"] == s) & (f["Temp"] == t)]
            if d.empty:  # no clean run left at this slope and temperature
                continue
            d = d.iloc[0]
            ax.bar(np.arange(len(subs)) + k * w, [100 - d[c] for c in subs], width=w, label=f"slope {s:g}")
        ax.set_xticks(np.arange(len(subs)) + 0.4 - w / 2)
        ax.set_xticklabels([c[4:] for c in subs], rotation=45, ha="right")
        ax.set_title(f"{kin}, {t:g} °C")
        ax.spines[["top", "right"]].set_visible(False)
    axes[0].set_ylabel("% of initial substrate C degraded")
    axes[0].legend(frameon=False, ncols=2)
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "km_sweep_substrates.png"), dpi=150, bbox_inches="tight")
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("input_folder")
    ap.add_argument("out_dir")
    ap.add_argument("--slopes", type=float, nargs="+", default=[7, 3, 1, 0.3, 0.1])
    ap.add_argument("--seeds", type=int, nargs="+", default=[0, 1, 2])
    ap.add_argument("--temps", type=float, nargs="+", default=TEMPS)
    ap.add_argument("--side", type=int, default=20)
    ap.add_argument("--days", type=int, default=365)
    ap.add_argument("--plot-only", action="store_true")
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    csv = os.path.join(a.out_dir, "timeseries.csv")
    if not a.plot_only:
        src = os.path.abspath(a.input_folder)
        jobs = [(src, k, s, t, sd, a.side, a.days) for k in KIN for s in a.slopes for t in a.temps for sd in a.seeds]
        print(f"{len(jobs)} runs of {a.days} days", file=sys.stderr, flush=True)
        with Pool(min(len(jobs), os.cpu_count())) as p:
            parts = p.map(run, jobs, chunksize=1)
        pd.concat(parts, ignore_index=True).to_csv(csv, index=False)
    ts = pd.read_csv(csv)
    ts, n_bad = drop_nan_runs(ts)
    print(f"dropped {n_bad} runs with NaN substrate pools", file=sys.stderr)
    day = int(ts["day"].max())
    summary_plot(ts, a.out_dir, ts["seed"].nunique(), day)
    ts["n"] = 1
    print(ts[ts["day"] == day].groupby(["kin", "slope", "Temp"])["n"].sum().unstack("Temp").to_string(), file=sys.stderr)
    substrate_plot(ts, a.out_dir, day)
    print("wrote", a.out_dir)


if __name__ == "__main__":
    main()
