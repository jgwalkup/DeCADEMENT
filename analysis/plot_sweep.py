"""
Two figures for a one-parameter sweep saved by sweep_minimal.py (--csv and --timeseries):
  sweep_end_of_run.png  end-of-run outcomes against the swept value
  sweep_timeseries.png  daily trajectories, one line per swept value

Styled to match this analysis folder's other figure scripts (see cue_figures.py,
energy_allocation_figures.py, etc.): scenario_utils' blue->orange->red gradient
for the swept parameter, plain white background, top/right spines off, dpi=150 PNG.

Usage:
    python plot_sweep.py ../analysis_output/maint_sweep Maint_E_Rate_Const --baseline 2.0

Reads <sweep_dir>/summary.csv and <sweep_dir>/timeseries.csv; writes the figures into <sweep_dir>.
"""
import argparse
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from scenario_utils import scenario_color_map  # noqa: E402 - registers 'blue_orange_red'

BLUE = "#1f77b4"  # scenario_utils._PALETTE[0]; used for single-series panels

plt.rcParams.update({
    "font.size": 10, "axes.titlesize": 11, "axes.labelsize": 10,
    "legend.fontsize": 9, "figure.dpi": 150,
})


def fmt(v):
    return f"{v:.3g}"


def vfmt(v):
    return f"{v:.1f}" if abs(v * 10 - round(v * 10)) < 1e-6 else f"{v:.3g}"


def end_of_run(summary, param, baseline, n_cells, meta, out):
    panels = [
        ("deg_pct", "Substrate degraded", "% of initial C, end of run"),
        ("biomC", "Biomass C", "C per grid cell, end of run"),
        ("CUE", "Carbon use efficiency", "mean of last 30 days"),
        ("resp_cum", "Cumulative respiration", "C per grid cell, whole run"),
        ("EC_end", "Energy : C ratio", "living biomass, end of run"),
        ("cells_lost", "Cells lost", f"of {n_cells} cells, end of run"),
    ]
    s = summary.copy()
    s["cells_lost"] = n_cells - s["alive"]
    x = s.index.values

    fig, axes = plt.subplots(2, 3, figsize=(13, 7.5), sharex=True)
    fig.suptitle(f"End-of-Run Outcomes vs {param}", fontsize=11, fontweight="bold")
    for ax, (col, title, ylabel) in zip(axes.flat, panels):
        y = s[col].values
        counts = col == "cells_lost"
        if counts:
            ax.plot(x, y, "o", color=BLUE, ms=5)
        else:
            ax.plot(x, y, color=BLUE, lw=1.6, marker="o", ms=5)
        if baseline in x:
            ax.plot([baseline], [s.loc[baseline, col]], marker="o", ms=12,
                    mfc="none", mec="black", mew=1.3, zorder=4)
        ax.set_xlabel(param)
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.set_ylim(0, 1) if col == "CUE" else ax.set_ylim(bottom=0)
        ax.spines[["top", "right"]].set_visible(False)
    fig.text(0.5, -0.01, f"Ring = baseline ({vfmt(baseline)}). {meta}",
              fontsize=8, color="dimgray", ha="center")
    fig.tight_layout()
    fig.savefig(f"{out}.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


def timeseries(ts, param, baseline, meta, out):
    panels = [
        ("energy_to_C", "Energy : C ratio", "living biomass"),
        ("net_energy", "Net Energy", "per living cell per day, 7-day mean"),
        ("enzyme_pool", "Enzyme Pool", "C per grid cell"),
        ("CUE", "Carbon Use Efficiency", "7-day mean"),
    ]
    vals = np.sort(ts[param].unique())
    labels = [vfmt(v) for v in vals]
    color = scenario_color_map(labels, palette="temperature")

    fig, axes = plt.subplots(2, 2, figsize=(13, 9), sharex=True)
    fig.suptitle(f"Trajectories Across {param}\n(each line one run; blue = low, red = high)",
                 fontsize=11, fontweight="bold")
    for ax, (col, title, ylabel) in zip(axes.flat, panels):
        piv = ts.pivot(index="day", columns=param, values=col)
        piv = piv.rolling(7, center=True, min_periods=1).mean() if col in ("net_energy", "CUE") else piv
        for v, lab in zip(vals, labels):
            is_base = v == baseline
            ax.plot(piv.index, piv[v], color=color[lab], lw=2.4 if is_base else 1.4,
                    label=lab + (" (baseline)" if is_base else ""),
                    zorder=4 if is_base else 2)
        if col == "net_energy":
            ax.axhline(0, ls="--", lw=0.7, color="gray", alpha=0.5)
        ax.set_xlabel("Day")
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        if col == "CUE":
            ax.set_ylim(0, 1)
        ax.spines[["top", "right"]].set_visible(False)
        ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.18), ncols=5,
                  framealpha=0.85, fontsize=7, title=param, title_fontsize=7.5)
    fig.text(0.5, -0.01, meta, fontsize=8, color="dimgray", ha="center")
    fig.tight_layout()
    fig.savefig(f"{out}.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("sweep_dir")
    ap.add_argument("param")
    ap.add_argument("--baseline", type=float, required=True, help="swept value to highlight as the baseline run")
    ap.add_argument("--meta", default="20x20 grid, 365 days, seed 0", help="footnote describing the runs")
    a = ap.parse_args()
    summary = pd.read_csv(os.path.join(a.sweep_dir, "summary.csv"), index_col=0)
    summary.index = summary.index.round(6)
    ts = pd.read_csv(os.path.join(a.sweep_dir, "timeseries.csv"))
    ts[a.param] = ts[a.param].round(6)
    n_cells = int(ts["alive_cells"].max())
    end_of_run(summary, a.param, a.baseline, n_cells, a.meta, os.path.join(a.sweep_dir, "sweep_end_of_run"))
    timeseries(ts, a.param, a.baseline, a.meta, os.path.join(a.sweep_dir, "sweep_timeseries"))
    print("wrote sweep_end_of_run.png and sweep_timeseries.png in", a.sweep_dir)


if __name__ == "__main__":
    main()
