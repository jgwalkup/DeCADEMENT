"""
Daily substrate degradation and CUE at one multiple of baseline (default 1x), from two maint_endo_sweep.py runs
(e.g. HSP effect on maintenance Ea on vs off), endogenous respiration (R_endo) off and on: four lines per panel, one column per temperature.

Usage:
    python plot_baseline_cue_degradation.py DIR_A DIR_B OUT.png [--labels "HSP effect on" "HSP effect off"] [--mult 1]
Colours = run (A blue, B orange); wide pale = R_endo off, thin dashed = R_endo on. CUE is left blank on days with no living cells.
"""
import argparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd

plt.rcParams.update({"font.size": 9, "axes.titlesize": 10, "axes.labelsize": 9, "legend.fontsize": 8, "figure.dpi": 150,
                     "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True, "grid.color": "#e6e5e0",
                     "grid.linewidth": 0.6, "axes.edgecolor": "#8a8981", "xtick.color": "#5f5e57", "ytick.color": "#5f5e57"})
COLORS = ["#2a78d6", "#eb6834"]
ENDO = r"R$_{\mathrm{endo}}$"  # endogenous respiration


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("dir_a")
    ap.add_argument("dir_b")
    ap.add_argument("out")
    ap.add_argument("--labels", nargs=2, default=["HSP effect on (slope -90,000)", "HSP effect off (slope 0)"])
    ap.add_argument("--mult", type=float, default=1.0)
    a = ap.parse_args()
    s = pd.read_csv(os.path.join(a.dir_a, "summary.csv"), nrows=1)
    param, base = s["param"].iloc[0], float(s["value"].iloc[0] / s["mult"].iloc[0])
    parts = []
    for d, lab in zip((a.dir_a, a.dir_b), a.labels):
        ts = pd.read_csv(os.path.join(d, "timeseries.csv.gz"),
                         usecols=["mult", "Temp", "R_endo", "seed", "day", "CUE", "substrate_remaining_pct", "alive_cells"])
        parts.append(ts[ts["mult"].round(6) == round(a.mult, 6)].assign(run=lab))
    ts = pd.concat(parts, ignore_index=True)
    ts["degraded"] = 100 - ts["substrate_remaining_pct"]
    ts.loc[ts["alive_cells"] == 0, "CUE"] = np.nan
    m = ts.groupby(["run", "R_endo", "Temp", "day"])[["degraded", "CUE"]].mean()
    temps = sorted(ts["Temp"].unique())

    fig, axes = plt.subplots(2, len(temps), figsize=(3.1 * len(temps), 5.4), sharex=True, sharey="row", squeeze=False)
    for c, t in enumerate(temps):
        for r, (v, lab) in enumerate((("degraded", "Substrate degraded (%)"), ("CUE", "CUE"))):
            ax = axes[r][c]
            for i, run in enumerate(a.labels):
                # R_endo off: wide pale line underneath; R_endo on: thin dashed line on top (they coincide until
                # cells run out of energy, so both stay visible where they overlap)
                for e, ls, lw, alpha in ((0, "-", 4.0, 0.35), (1, (0, (3, 2)), 1.3, 1.0)):
                    y = m.loc[(run, e, t)][v]
                    ax.plot(y.index, y.values, color=COLORS[i], ls=ls, lw=lw, alpha=alpha, solid_capstyle="butt",
                            label=f"{run}, {ENDO} {'off' if e == 0 else 'on'}")
            ax.set_ylim(0, 100 if v == "degraded" else 0.55)
            if r == 0:
                ax.set_title(f"{t:g} °C")
            if c == 0:
                ax.set_ylabel(lab)
            if r == 1:
                ax.set_xlabel("Day")
    h, l = axes[0][0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncol=2, frameon=False, bbox_to_anchor=(0.5, -0.06))
    fig.suptitle(f"Daily substrate degradation and CUE at {param} = {base * a.mult:g} ({a.mult:g}× baseline; mean of seeds)\n"
                 f"Endogenous respiration ({ENDO}) off (wide pale) and on (thin dashed) coincide until cells run out of energy")
    fig.tight_layout()
    fig.savefig(a.out, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
