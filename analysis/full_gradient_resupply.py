"""
Temperature gradient of the full multi-taxon, multi-substrate model with 365-day substrate resupply: consecutive pulses of 365 days
at a constant temperature, with and without HSP production (no-HSP control), seed-replicated (each seed = a different random
community). Flat climate (Psi -0.5, FM 50). At the end of each pulse Grid.repopulation() resets substrate, monomers and enzymes and
re-seeds the community from the taxa's cumulative abundance (a stub provides the per-taxon counts that Output normally tracks), as
in src/dementpy.py.

Usage:
    python full_gradient_resupply.py ../input_e1500_enz2 ../analysis_output/full_gradient_resupply [--pulses 3] [--seeds 0 1 2 3 4]
        [--side 20] [--temps 4 8 ...] [--kinetics 0|1] [--plot-only]
Writes <out_dir>/timeseries.csv (one row per run and day) and the figures
  full_gradient_resupply_summary.png     end-of-pulse outcomes vs temperature, HSP (solid) vs control (dashed), line shade = pulse
  full_gradient_resupply_timeseries.png  CUE, energy:C, biomass C trajectories over all pulses at selected temperatures
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
from full_climate_run import CYCLE, _StubOutput  # noqa: E402
from full_gradient_730 import make_folders  # noqa: E402
from grid import Grid  # noqa: E402
from initialization import initialize_data  # noqa: E402
from sweep_minimal import set_runtime  # noqa: E402

plt.rcParams.update({"font.size": 10, "axes.titlesize": 11, "axes.labelsize": 10, "legend.fontsize": 9, "figure.dpi": 150})
TEMPS = [4, 8, 12, 16, 20, 24, 28, 32, 36, 40]
SHOW = [4, 12, 20, 28, 36, 40]
PCOL = ["#9ecae1", "#4292c6", "#08306b", "#6a51a3", "#54278f"]


def run(job):
    src, temp, side, pulses, seed = job
    tmp = tempfile.mkdtemp()
    w = os.path.join(tmp, "in")
    shutil.copytree(src, w)
    set_runtime(os.path.join(w, "runtime.txt"), dict(x=side, y=side, gridsize=side * side, pulse=pulses, end_time=CYCLE))
    n = pulses * CYCLE + 5
    pd.DataFrame({"Temp": float(temp), "FM": 50.0, "Psi": -0.5, "Date": "1/1/00"}, index=range(n)).to_csv(
        os.path.join(w, "climate.csv"))
    cwd = os.getcwd()
    os.chdir(w)
    rows = []
    try:
        np.random.seed(seed)
        rt = pd.read_csv("runtime.txt", header=None, index_col=0, sep="\t").drop("output_format").astype(float)
        data = initialize_data(rt)
        eco = Grid(rt, data)
        G, T = side * side, int(rt.loc["n_taxa", 1])
        stub = _StubOutput(T, pulses * CYCLE, pd.unique(eco.Microbes.index))
        subs = [s for s in data["Substrate_index"] if s not in ("DeadMic", "DeadEnz")]
        sub_total = lambda: float(sum(np.asarray(eco.Substrates.loc[s, "C"]).sum() for s in subs))
        sub0 = sub_total()
        ps = lambda x: float(np.asarray(x).sum())
        for pl in range(pulses):
            for i in range(pl * CYCLE, (pl + 1) * CYCLE):
                eco.degradation(pl, i)
                eco.uptake(pl, i)
                eco.metabolism(i)
                eco.mortality(i)
                eco.reproduction(i)
                C = np.asarray(eco.Microbes["C"])
                E = np.asarray(eco.Microbes["E"])
                alive = C > 0
                nal = int(alive.sum())
                cnt = alive.reshape(G, T).sum(axis=0)
                stub.Taxon_count_repop.iloc[:, i + 1] = cnt
                rows.append(dict(
                    input=os.path.basename(src), Temp=temp, seed=seed, pulse=pl + 1, day=i + 1, doy=i - pl * CYCLE + 1,
                    biomass_C=C.sum() / G, alive_rows=nal, taxa_alive=int((cnt > 0).sum()),
                    energy_to_C=E[alive].sum() / C[alive].sum() if nal else np.nan, CUE=float(eco.CUE_system),
                    substrate_remaining_pct=100 * sub_total() / sub0,
                    HSP_prop=float(np.asarray(eco.HSP_prop)[alive].mean()) if nal else np.nan,
                    net_energy=ps(eco.Microbe_E_Gain) / max(nal, 1)))
                if i == (pl + 1) * CYCLE - 1 and pl < pulses - 1:
                    eco.repopulation(stub, i, True)
    finally:
        os.chdir(cwd)
        shutil.rmtree(tmp, ignore_errors=True)
    return pd.DataFrame(rows)


def end_values(ts):
    """Values on day 365 of each pulse (CUE = mean of the last 30 days of the pulse), seed-mean."""
    a = ts[ts["doy"] == CYCLE].groupby(["hsp", "pulse", "Temp"])[
        ["biomass_C", "energy_to_C", "alive_rows", "taxa_alive", "substrate_remaining_pct"]].mean()
    c = ts[ts["doy"] > CYCLE - 30].groupby(["hsp", "pulse", "Temp"])["CUE"].mean()
    return a.join(c).reset_index()


def summary_plot(ts, out_dir):
    e = end_values(ts)
    n_seeds, pulses = ts["seed"].nunique(), int(ts["pulse"].max())
    fig, axes = plt.subplots(2, 3, figsize=(14, 8), sharex=True)
    fig.suptitle(f"Full model, 365-day resupply: HSP (solid) vs no-HSP control (dashed) on day 365 of each pulse, mean of {n_seeds} "
                 f"communities", fontsize=11, fontweight="bold")
    panels = [("biomass_C", "Biomass C", "C per grid cell"), ("rel", "Change in biomass vs control", "%"),
              ("CUE", "CUE (mean of last 30 days)", "CUE"), ("energy_to_C", "Energy : C", "living biomass"),
              ("alive_rows", "Cell-taxa alive", "cell-taxa"), ("deg", "Substrate degraded in the pulse", "% of initial C")]
    for ax, (col, title, yl) in zip(axes.flat, panels):
        for p in range(1, pulses + 1):
            colr = PCOL[p - 1]
            if col == "rel":
                h = e[e["hsp"] & (e["pulse"] == p)].set_index("Temp")["biomass_C"]
                c = e[~e["hsp"] & (e["pulse"] == p)].set_index("Temp")["biomass_C"]
                ax.plot(h.index, 100 * (h - c) / c, color=colr, marker="o", ms=4, lw=1.6, label=f"pulse {p}")
                continue
            for hsp, ls in ((True, "-"), (False, "--")):
                d = e[(e["hsp"] == hsp) & (e["pulse"] == p)].sort_values("Temp")
                y = 100 - d["substrate_remaining_pct"] if col == "deg" else d[col]
                ax.plot(d["Temp"], y, color=colr, ls=ls, marker="o", ms=4, lw=1.5, label=f"pulse {p}" if hsp else None)
        if col == "rel":
            ax.axhline(0, color="gray", lw=0.8, ls="--")
        if col == "CUE":
            ax.set_ylim(0, 1)
        ax.set_title(title, fontsize=10)
        ax.set_ylabel(yl, fontsize=9)
        ax.spines[["top", "right"]].set_visible(False)
    for ax in axes[1]:
        ax.set_xlabel("Temperature (°C)")
    h, l = axes[0][0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncols=pulses, framealpha=0.85, bbox_to_anchor=(0.5, -0.02))
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "full_gradient_resupply_summary.png"), dpi=150, bbox_inches="tight")
    plt.close(fig)


def timeseries_plot(ts, out_dir):
    n_seeds = ts["seed"].nunique()
    m = ts.groupby(["hsp", "Temp", "day"], as_index=False)[["CUE", "energy_to_C", "biomass_C"]].mean()
    rows = [("CUE", "CUE", (0, 1)), ("energy_to_C", "Energy : C", (0, None)), ("biomass_C", "Biomass C", (0, None))]
    fig, axes = plt.subplots(3, len(SHOW), figsize=(3.1 * len(SHOW), 8.5), sharex=True, squeeze=False)
    fig.suptitle(f"Full model, 365-day resupply: HSP (solid) vs no-HSP control (dashed), mean of {n_seeds} communities, 7-day mean",
                 fontsize=11, fontweight="bold")
    for i, (col, lab, yl) in enumerate(rows):
        for j, t in enumerate(SHOW):
            ax = axes[i][j]
            for hsp, ls, c in ((False, "--", "gray"), (True, "-", "#1f77b4")):
                d = m[(m["hsp"] == hsp) & (m["Temp"] == t)].sort_values("day")
                ax.plot(d["day"], d[col].rolling(7, center=True, min_periods=1).mean(), color=c, ls=ls, lw=1.4,
                        label="HSP" if hsp else "no HSP (control)")
            for x in range(CYCLE, int(m["day"].max()), CYCLE):
                ax.axvline(x, color="lightgray", lw=0.8)
            ax.set_ylim(*yl) if col != "biomass_C" else ax.set_ylim(0, None)
            if i == 0:
                ax.set_title(f"{t}°C")
            if j == 0:
                ax.set_ylabel(lab)
            if i == 2:
                ax.set_xlabel("Day")
            ax.spines[["top", "right"]].set_visible(False)
    h, l = axes[0][0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncols=2, framealpha=0.85, bbox_to_anchor=(0.5, -0.02))
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "full_gradient_resupply_timeseries.png"), dpi=150, bbox_inches="tight")
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("input_folder")
    ap.add_argument("out_dir")
    ap.add_argument("--pulses", type=int, default=3)
    ap.add_argument("--seeds", type=int, nargs="+", default=[0, 1, 2, 3, 4])
    ap.add_argument("--side", type=int, default=20)
    ap.add_argument("--temps", type=float, nargs="+", default=TEMPS)
    ap.add_argument("--kinetics", type=int, choices=[0, 1], default=None,
                    help="override Kinetics_Model in runtime.txt (0 = Arrhenius, 1 = MMRT); default: as in the input folder")
    ap.add_argument("--plot-only", action="store_true")
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    csv = os.path.join(a.out_dir, "timeseries.csv")
    if not a.plot_only:
        tmp = tempfile.mkdtemp()
        on, off = make_folders(os.path.abspath(a.input_folder), tmp)
        if a.kinetics is not None:
            for f in (on, off):
                set_runtime(os.path.join(f, "runtime.txt"), {"Kinetics_Model": a.kinetics})
        jobs = [(f, t, a.side, a.pulses, s) for t in a.temps for s in a.seeds for f in (on, off)]
        print(f"{len(jobs)} runs of {a.pulses} x {CYCLE} days", file=sys.stderr, flush=True)
        with Pool(min(len(jobs), os.cpu_count())) as p:
            parts = p.map(run, jobs, chunksize=1)
        shutil.rmtree(tmp, ignore_errors=True)
        ts = pd.concat(parts, ignore_index=True)
        ts["hsp"] = ts["input"] == "hsp_on"
        ts.to_csv(csv, index=False)
    ts = pd.read_csv(csv)
    summary_plot(ts, a.out_dir)
    timeseries_plot(ts, a.out_dir)
    print("wrote", a.out_dir)


if __name__ == "__main__":
    main()
