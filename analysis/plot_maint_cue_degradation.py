"""
Daily substrate degradation and CUE for every run of maint_erosion_sweep.py, one column per temperature, one line per
multiple of baseline (mean of seeds).

Usage:
    python plot_maint_cue_degradation.py ../analysis_output/maint_erosion_sweep
Writes <dir>/maint_cue_degradation.png. Rows: substrate degraded (%; identical with endogenous respiration (R_endo)
off and on), CUE with R_endo off, CUE with R_endo on. CUE is the grid value of each day (after rebalancing losses, floored at 0); days with
no living cells are left blank.
"""
import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm, LinearSegmentedColormap
from matplotlib.cm import ScalarMappable
import numpy as np
import pandas as pd

plt.rcParams.update({"font.size": 9, "axes.titlesize": 10, "axes.labelsize": 9, "legend.fontsize": 8, "figure.dpi": 150,
                     "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True, "grid.color": "#e6e5e0",
                     "grid.linewidth": 0.6, "axes.edgecolor": "#8a8981", "xtick.color": "#5f5e57", "ytick.color": "#5f5e57"})
# one-hue blue ramp, light (low maintenance) -> dark (high maintenance); ordinal steps 250-700
RAMP = ["#86b6ef", "#6da7ec", "#5598e7", "#3987e5", "#2a78d6", "#256abf", "#1c5cab", "#184f95", "#104281", "#0d366b"]
BASE_COLOR = "#eb6834"  # the 1x baseline, highlighted


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("dir")
    ap.add_argument("--note", default="", help="text appended to the figure title, e.g. the fixed setting of the run")
    a = ap.parse_args()
    s = pd.read_csv(os.path.join(a.dir, "summary.csv"), nrows=1)
    param, base = s["param"].iloc[0], float(s["value"].iloc[0] / s["mult"].iloc[0])
    ts = pd.read_csv(os.path.join(a.dir, "timeseries.csv.gz"),
                     usecols=["mult", "Temp", "erosion", "seed", "day", "CUE", "substrate_remaining_pct", "alive_cells"])
    ts["degraded"] = 100 - ts["substrate_remaining_pct"]
    ts.loc[ts["alive_cells"] == 0, "CUE"] = np.nan
    m = ts.groupby(["Temp", "erosion", "mult", "day"])[["degraded", "CUE"]].mean()
    temps = sorted(ts["Temp"].unique())
    mults = sorted(ts["mult"].unique())
    cmap = LinearSegmentedColormap.from_list("ramp", RAMP, N=len(mults))
    norm = BoundaryNorm(np.arange(len(mults) + 1) - 0.5, cmap.N)   # one equal band per multiple
    idx = {mu: i for i, mu in enumerate(mults)}

    ENDO = r"R$_{\mathrm{endo}}$"  # endogenous respiration
    rows = [("degraded", 0, f"Substrate degraded (%)\n({ENDO} off = on)"), ("CUE", 0, f"CUE, {ENDO} off"),
            ("CUE", 1, f"CUE, {ENDO} on")]
    fig, axes = plt.subplots(len(rows), len(temps), figsize=(3.1 * len(temps), 2.5 * len(rows)), sharex=True,
                             sharey="row", squeeze=False)
    for c, t in enumerate(temps):
        for r, (var, e, lab) in enumerate(rows):
            ax = axes[r][c]
            for mu in mults:
                if mu == 1.0:
                    continue
                y = m.loc[(t, e, mu)][var]
                ax.plot(y.index, y.values, color=cmap(norm(idx[mu])), lw=1.0)
            y = m.loc[(t, e, 1.0)][var]
            ax.plot(y.index, y.values, color=BASE_COLOR, lw=2.2, label=f"1× ({base:g})")
            if var == "CUE":
                ax.set_ylim(0, 0.6)
            else:
                ax.set_ylim(0, 100)
            if r == 0:
                ax.set_title(f"{t:g} °C")
            if c == 0:
                ax.set_ylabel(lab)
            if r == len(rows) - 1:
                ax.set_xlabel("Day")
    axes[0][0].legend(loc="upper left", frameon=False)
    fig.tight_layout(rect=(0, 0, 0.93, 0.93 if a.note else 0.97))
    cax = fig.add_axes([0.945, 0.15, 0.012, 0.7])
    cb = fig.colorbar(ScalarMappable(norm=norm, cmap=cmap), cax=cax, ticks=range(len(mults)))
    cb.ax.set_yticklabels([f"{x:g}" for x in mults], fontsize=7)
    cb.set_label(f"{param} (× baseline {base:g})")
    fig.suptitle(f"Substrate degradation and CUE over time by {param} (mean of seeds)" + (f"\n{a.note}" if a.note else ""))
    fig.savefig(os.path.join(a.dir, "maint_cue_degradation.png"), bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
