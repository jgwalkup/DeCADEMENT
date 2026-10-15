"""
Plot the HSP-vs-no-HSP effect from sweep_hsp_effect.py: % change in biomass C relative to the
no-HSP control against temperature, one panel per (row param value, column param value), one line
per value of the line parameter. Cold (<10 C) and warm (>25 C) target zones are shaded.

Usage:
    python plot_hsp_window.py ../analysis_output/hsp_window_screen2 \\
        --rows Maint_E_HSP_Ea_Slope --cols Maint_E_Rate_Const --lines HSP_prod_E_cost
Reads <dir>/summary.csv, writes <dir>/hsp_window.png.
"""
import argparse
import os
import sys

import matplotlib.pyplot as plt
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from scenario_utils import scenario_color_map  # noqa: E402

plt.rcParams.update({
    "font.size": 10, "axes.titlesize": 11, "axes.labelsize": 10,
    "legend.fontsize": 9, "figure.dpi": 150,
})


def vfmt(v):
    return "0" if v == 0 else (f"{v:.0f}" if abs(v) >= 100 else f"{v:.3g}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("sweep_dir")
    ap.add_argument("--rows", required=True)
    ap.add_argument("--cols", required=True)
    ap.add_argument("--lines", required=True)
    ap.add_argument("--cold", type=float, default=10.0, help="cold zone upper bound (C)")
    ap.add_argument("--warm", type=float, default=25.0, help="warm zone lower bound (C)")
    ap.add_argument("--ymin", type=float, default=-40.0)
    ap.add_argument("--ymax", type=float, default=100.0)
    ap.add_argument("--meta", default="20x20 grid, 365 days, seed 0; control = no HSP production")
    a = ap.parse_args()

    s = pd.read_csv(os.path.join(a.sweep_dir, "summary.csv"))
    rows, cols, lines = (sorted(s[c].unique()) for c in (a.rows, a.cols, a.lines))
    colors = scenario_color_map([vfmt(v) for v in lines], palette="discrete")

    fig, axes = plt.subplots(len(rows), len(cols), figsize=(4.2 * len(cols), 3.6 * len(rows)),
                             sharex=True, sharey=True, squeeze=False)
    fig.suptitle("Effect of HSP production on biomass C vs temperature (% change vs no-HSP control)",
                 fontsize=11, fontweight="bold")
    tmin, tmax = s["Temp"].min(), s["Temp"].max()
    for i, r in enumerate(rows):
        for j, c in enumerate(cols):
            ax = axes[i][j]
            d = s[(s[a.rows] == r) & (s[a.cols] == c)]
            ax.axvspan(tmin - 1, a.cold, color="#2166ac", alpha=0.07, lw=0)
            ax.axvspan(a.warm, tmax + 1, color="#d6172a", alpha=0.07, lw=0)
            ax.axhline(0, color="gray", lw=0.8, ls="--")
            for v in lines:
                g = d[d[a.lines] == v].sort_values("Temp")
                ax.plot(g["Temp"], g["rel_biomC"].clip(a.ymin, a.ymax), color=colors[vfmt(v)], lw=1.5,
                        marker="o", ms=3.5, label=vfmt(v))
            ax.set_title(f"{a.rows.replace('Maint_E_HSP_Ea_', '')} = {vfmt(r)}\n{a.cols.replace('Maint_E_', '')} = {vfmt(c)}",
                         fontsize=10)
            ax.set_ylim(a.ymin, a.ymax)
            ax.spines[["top", "right"]].set_visible(False)
            if i == len(rows) - 1:
                ax.set_xlabel("Temperature (°C)")
            if j == 0:
                ax.set_ylabel("Δ biomass C vs no HSP (%)")
    handles, labels = axes[0][0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncols=len(lines), title=a.lines, framealpha=0.85,
               bbox_to_anchor=(0.5, -0.07))
    fig.text(0.5, -0.125, f"Shaded: cold (<{a.cold:g}°C, blue) and warm (>{a.warm:g}°C, red) target zones. {a.meta}",
             fontsize=8, color="dimgray", ha="center")
    fig.tight_layout()
    out = os.path.join(a.sweep_dir, "hsp_window.png")
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print("wrote", out)


if __name__ == "__main__":
    main()
