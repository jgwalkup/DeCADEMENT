"""
CUE, biomass C and remaining substrate over time from full_gradient_resupply.py output (HSP vs no-HSP control).

Usage:
    python plot_resupply_timeseries.py ../analysis_output/full_gradient_resupply [--temps 4 12 20 28 36 40] [--out NAME.png]
Lines are the mean of the communities (7-day rolling mean); shading is the min-max across communities; vertical lines mark
the start of each 365-day pulse (substrate reset and reseeding).
"""
import argparse
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

plt.rcParams.update({"font.size": 10, "axes.titlesize": 11, "axes.labelsize": 10, "legend.fontsize": 9, "figure.dpi": 150})
ROWS = [("CUE", "CUE", (0, 0.7)), ("biomass_C", "Biomass C per grid cell", (0, None)),
        ("substrate_remaining_pct", "Substrate C remaining (% of initial)", (0, 100))]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("folder")
    ap.add_argument("--temps", type=float, nargs="+", default=[4, 12, 20, 28, 36, 40])
    ap.add_argument("--out", default="full_gradient_resupply_cue_biomass_substrate.png")
    a = ap.parse_args()
    ts = pd.read_csv(os.path.join(a.folder, "timeseries.csv"))
    n = ts["seed"].nunique()
    last = int(ts["day"].max())
    fig, axes = plt.subplots(3, len(a.temps), figsize=(3.2 * len(a.temps), 8.8), sharex=True, squeeze=False)
    fig.suptitle(f"Full model, 365-day resupply: HSP vs no-HSP control, mean (line) and range (shading) of {n} communities, "
                 f"7-day mean", fontsize=11, fontweight="bold")
    for i, (col, lab, yl) in enumerate(ROWS):
        for j, t in enumerate(a.temps):
            ax = axes[i][j]
            for hsp, c, ls, name in ((False, "#7f7f7f", "--", "no HSP (control)"), (True, "#1f77b4", "-", "HSP")):
                d = ts[(ts["Temp"] == t) & (ts["hsp"] == hsp)].pivot(index="day", columns="seed", values=col)
                d = d.rolling(7, center=True, min_periods=1).mean()
                ax.fill_between(d.index, d.min(axis=1), d.max(axis=1), color=c, alpha=0.15, lw=0)
                ax.plot(d.index, d.mean(axis=1), color=c, ls=ls, lw=1.4, label=name)
            for x in range(365, last, 365):
                ax.axvline(x, color="lightgray", lw=0.8)
            ax.set_ylim(*yl)
            if i == 0:
                ax.set_title(f"{t:g}°C")
            if j == 0:
                ax.set_ylabel(lab)
            if i == len(ROWS) - 1:
                ax.set_xlabel("Day")
            ax.spines[["top", "right"]].set_visible(False)
    h, l = axes[0][0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncols=2, framealpha=0.85, bbox_to_anchor=(0.5, -0.01))
    fig.tight_layout()
    out = os.path.join(a.folder, a.out)
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print("wrote", out)


if __name__ == "__main__":
    main()
