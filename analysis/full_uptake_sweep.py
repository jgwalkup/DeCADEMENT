"""
Full model: uptake Vmax-Km slope (Uptake_Vmax_Km, shared by MMRT and Arrhenius) crossed with the enzyme Vmax-Km
slope, across a temperature gradient. Uptake Km = uptake Vmax x slope, so for monomer pools well below Km (the case here) the
uptake of a monomer per unit transporter is about pool / slope. Flat climate (Psi -0.5, FM 50), one 365-day pulse, HSP production
as in the input folder, seed-replicated communities.

Usage:
    python full_uptake_sweep.py ../input_e1500_enz2 ../analysis_output/uptake_sweep [--enz-slopes 2 3 4]
        [--uptake-slopes 10 5 2.5 1] [--temps 4 12 20 24 28 40] [--seeds 0 1 2] [--plot-only]
Writes <out_dir>/timeseries.csv, uptake_sweep_summary.png (outcomes vs temperature; columns = kinetics x enzyme slope, lines =
uptake slope) and uptake_sweep_enzyme_slopes.png (the same data with the enzyme slope as the line colour).
"""
import argparse
import os
import sys
from multiprocessing import Pool

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from full_km_sweep import KIN, drop_nan_runs, run  # noqa: E402

plt.rcParams.update({"font.size": 9, "axes.titlesize": 10, "axes.labelsize": 9, "legend.fontsize": 8, "figure.dpi": 150})
PANELS = [("degraded", "Substrate degraded (% of initial C)", (0, 100)), ("biomass_C", "Biomass C per grid cell", (0, None)),
          ("CUE30", "CUE (last 30 days)", (0, 0.7)), ("energy_to_C", "Energy : C", (0, None)),
          ("mon_C", "Monomer C pool per cell", (0, None))]


def final(ts, day):
    g = ts[ts["day"] == day].groupby(["kin", "slope", "uslope", "Temp"])
    a = g.mean(numeric_only=True)
    a = a[g.size() >= 2]
    a["degraded"] = 100 - a["remaining_pct"]
    a["CUE30"] = ts[(ts["day"] > day - 30) & (ts["day"] <= day)].groupby(["kin", "slope", "uslope", "Temp"])["CUE"].mean()
    return a.reset_index()


def summary_plot(f, out_dir, n, day, enz):
    uptakes = sorted(f["uslope"].unique(), reverse=True)
    cm = plt.cm.plasma(np.linspace(0.05, 0.85, len(uptakes)))
    cols = [(k, e) for k in KIN for e in enz]
    fig, axes = plt.subplots(len(PANELS), len(cols), figsize=(2.9 * len(cols), 2.5 * len(PANELS)), sharex=True, squeeze=False)
    fig.suptitle(f"Full model: uptake Vmax-Km slope (line) x enzyme Vmax-Km slope (column), day {day}, mean of {n} communities",
                 fontsize=11, fontweight="bold")
    for r, (col, lab, yl) in enumerate(PANELS):
        for c, (kin, e) in enumerate(cols):
            ax = axes[r][c]
            for u, colr in zip(uptakes, cm):
                d = f[(f["kin"] == kin) & (f["slope"] == e) & (f["uslope"] == u)].sort_values("Temp")
                ax.plot(d["Temp"], d[col], marker="o", ms=3, lw=1.3, color=colr, label=f"uptake slope {u:g}")
            ax.set_ylim(*yl)
            if r == 0:
                ax.set_title(f"{kin}, enzyme slope {e:g}")
            if c == 0:
                ax.set_ylabel(lab)
            if r == len(PANELS) - 1:
                ax.set_xlabel("Temperature (°C)")
            ax.spines[["top", "right"]].set_visible(False)
    h, l = axes[0][0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncols=len(uptakes), framealpha=0.85, bbox_to_anchor=(0.5, -0.005))
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "uptake_sweep_summary.png"), dpi=150, bbox_inches="tight")
    plt.close(fig)


def enzyme_plot(f, out_dir, n, day, enz):
    uptakes = sorted(f["uslope"].unique(), reverse=True)
    cm = plt.cm.viridis(np.linspace(0.1, 0.85, len(enz)))
    cols = [(k, u) for k in KIN for u in uptakes]
    fig, axes = plt.subplots(3, len(cols), figsize=(2.5 * len(cols), 7.2), sharex=True, squeeze=False)
    fig.suptitle(f"Same data, enzyme slope as line colour (columns = kinetics x uptake slope), day {day}, {n} communities",
                 fontsize=11, fontweight="bold")
    for r, (col, lab, yl) in enumerate(PANELS[:2] + PANELS[3:4]):
        for c, (kin, u) in enumerate(cols):
            ax = axes[r][c]
            for e, colr in zip(enz, cm):
                d = f[(f["kin"] == kin) & (f["slope"] == e) & (f["uslope"] == u)].sort_values("Temp")
                ax.plot(d["Temp"], d[col], marker="o", ms=3, lw=1.3, color=colr, label=f"enzyme slope {e:g}")
            ax.set_ylim(*yl)
            if r == 0:
                ax.set_title(f"{kin}, uptake {u:g}")
            if c == 0:
                ax.set_ylabel(lab)
            if r == 2:
                ax.set_xlabel("Temperature (°C)")
            ax.spines[["top", "right"]].set_visible(False)
    h, l = axes[0][0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncols=len(enz), framealpha=0.85, bbox_to_anchor=(0.5, -0.01))
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "uptake_sweep_enzyme_slopes.png"), dpi=150, bbox_inches="tight")
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("input_folder")
    ap.add_argument("out_dir")
    ap.add_argument("--enz-slopes", type=float, nargs="+", default=[2, 3, 4])
    ap.add_argument("--uptake-slopes", type=float, nargs="+", default=[10, 5, 2.5, 1])
    ap.add_argument("--temps", type=float, nargs="+", default=[4, 12, 20, 24, 28, 40])
    ap.add_argument("--seeds", type=int, nargs="+", default=[0, 1, 2])
    ap.add_argument("--side", type=int, default=20)
    ap.add_argument("--days", type=int, default=365)
    ap.add_argument("--plot-only", action="store_true")
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    csv = os.path.join(a.out_dir, "timeseries.csv")
    if not a.plot_only:
        src = os.path.abspath(a.input_folder)
        jobs = [(src, k, e, t, sd, a.side, a.days, u) for k in KIN for e in a.enz_slopes for u in a.uptake_slopes
                for t in a.temps for sd in a.seeds]
        print(f"{len(jobs)} runs of {a.days} days", file=sys.stderr, flush=True)
        with Pool(min(len(jobs), os.cpu_count())) as p:
            parts = p.map(run, jobs, chunksize=1)
        pd.concat(parts, ignore_index=True).to_csv(csv, index=False)
    ts = pd.read_csv(csv)
    ts, n_bad = drop_nan_runs_u(ts)
    print(f"dropped {n_bad} runs with NaN substrate pools", file=sys.stderr)
    day = int(ts["day"].max())
    f = final(ts, day)
    enz = sorted(f["slope"].unique())
    summary_plot(f, a.out_dir, ts["seed"].nunique(), day, enz)
    enzyme_plot(f, a.out_dir, ts["seed"].nunique(), day, enz)
    print("wrote", a.out_dir)


def drop_nan_runs_u(ts):
    key = ["kin", "slope", "uslope", "Temp", "seed"]
    bad = ts.groupby(key)["remaining_pct"].apply(lambda x: x.isna().any())
    bad = bad[bad].index
    return ts[~ts.set_index(key).index.isin(bad)], len(bad)


if __name__ == "__main__":
    main()
