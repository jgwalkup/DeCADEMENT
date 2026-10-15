"""
Figures for a temperature-axis sweep saved by sweep_temperature.py (--csv and --timeseries),
for the HSP / maintenance thermal-response work. Same house style as plot_sweep.py (plain white
background, dpi=150 PNG, scenario_utils palettes) with panels specific to this mechanism:
effective maintenance Ea, HSP_prop, the raw per-capita maintenance rate, plus outcome variables.

1D (temperature only):
    python plot_temperature.py ../analysis_output/maint_thermal
    -> temp_end_of_run.png, temp_timeseries.png (one line per temperature, blue->orange->red)

2D (temperature x one swept parameter; pass the parameter column name and its baseline value):
    python plot_temperature.py ../analysis_output/maint_hsp_slope_thermal \\
        --param Maint_E_HSP_Ea_Slope --baseline -30000 [--ts-temps 28 36 40]
    -> temp_end_of_run.png   x = temperature, one line per parameter value (discrete palette)
       temp_heatmaps.png     temperature x parameter grids of end-of-run outcomes
       temp_timeseries.png   trajectories at selected temperatures, one line per parameter value
"""
import argparse
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from scenario_utils import scenario_color_map  # noqa: E402 - registers 'blue_orange_red'

BLUE = "#1f77b4"

plt.rcParams.update({
    "font.size": 10, "axes.titlesize": 11, "axes.labelsize": 10,
    "legend.fontsize": 9, "figure.dpi": 150,
})

# (summary column, title, y label, clamp y-axis at zero)
PANELS_1D = [
    ("deg_pct", "Substrate degraded", "% of initial C, end of run", True),
    ("biomC", "Biomass C", "C per grid cell, end of run", True),
    ("net_E", "Net energy", "mean of last 30 days, per living cell", False),
    ("maint_Ea", "Effective maintenance Ea", "J/mol, mean of last 30 days", False),
    ("HSP_prop", "HSP proportion of biomass", "mean of last 30 days", True),
    ("maint_rate_percap", "Raw maintenance rate", "per-unit-biomass, mean of last 30 days", True),
]
PANELS_2D = [
    ("deg_pct", "Substrate degraded", "% of initial C, end of run", True),
    ("biomC", "Biomass C", "C per grid cell, end of run", True),
    ("CUE", "Carbon use efficiency", "mean of last 30 days", True),
    ("net_E", "Net energy", "mean of last 30 days, per living cell", False),
    ("EC_end", "Energy : C ratio", "living biomass, end of run", True),
    ("maint_Ea", "Effective maintenance Ea", "J/mol, mean of last 30 days", False),
    ("HSP_prop", "HSP proportion of biomass", "mean of last 30 days", True),
    ("maint_rate_percap", "Raw maintenance rate", "per-unit-biomass, mean of last 30 days", True),
]


def vfmt(v):
    if v == 0:
        return "0"
    return f"{v:.1f}" if abs(v * 10 - round(v * 10)) < 1e-6 and abs(v) < 1000 else f"{v:g}"


# ── 1D: temperature only ──────────────────────────────────────────────────────
def end_of_run_1d(summary, out):
    x = summary.index.values
    fig, axes = plt.subplots(2, 3, figsize=(13, 7.5), sharex=True)
    fig.suptitle("Maintenance / HSP Thermal Response vs Temperature", fontsize=11, fontweight="bold")
    for ax, (col, title, ylabel, zero) in zip(axes.flat, PANELS_1D):
        ax.plot(x, summary[col].values, color=BLUE, lw=1.6, marker="o", ms=5)
        ax.set_xlabel("Temperature (°C)")
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        if zero:
            ax.set_ylim(bottom=0)
        ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(f"{out}.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


def timeseries_1d(ts, out):
    panels = [
        ("maint_rate_percap", "Raw Maintenance Rate", "per-unit-biomass"),
        ("maint_Ea_effective", "Effective Maintenance Ea", "J/mol (Intercept=30000)"),
        ("HSP_prop", "HSP Proportion of Biomass", ""),
        ("energy_to_C", "Energy : C Ratio", "living biomass"),
    ]
    vals = np.sort(ts["Temp"].unique())
    labels = [vfmt(v) for v in vals]
    color = scenario_color_map(labels, palette="temperature")
    fig, axes = plt.subplots(2, 2, figsize=(13, 9), sharex=True)
    fig.suptitle("Trajectories Across Temperature\n(each line one run; blue = low, red = high temperature)",
                 fontsize=11, fontweight="bold")
    for ax, (col, title, ylabel) in zip(axes.flat, panels):
        piv = ts.pivot_table(index="day", columns="Temp", values=col, aggfunc="mean")
        for v, lab in zip(vals, labels):
            ax.plot(piv.index, piv[v], color=color[lab], lw=1.4, label=lab)
        ax.set_xlabel("Day")
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.spines[["top", "right"]].set_visible(False)
        ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.18), ncols=5,
                  framealpha=0.85, fontsize=7, title="Temp (°C)", title_fontsize=7.5)
    fig.tight_layout()
    fig.savefig(f"{out}.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


# ── 2D: temperature x parameter ───────────────────────────────────────────────
def _param_colors(vals):
    labels = [vfmt(v) for v in vals]
    cmap = scenario_color_map(labels, palette="discrete")
    return {v: cmap[lab] for v, lab in zip(vals, labels)}


def end_of_run_2d(summary, param, baseline, out):
    vals = sorted(summary[param].unique())
    colors = _param_colors(vals)
    fig, axes = plt.subplots(2, 4, figsize=(17, 8), sharex=True)
    fig.suptitle(f"Thermal Response vs Temperature, by {param}", fontsize=11, fontweight="bold")
    for ax, (col, title, ylabel, zero) in zip(axes.flat, PANELS_2D):
        for v in vals:
            d = summary[summary[param] == v].sort_values("Temp")
            is_base = baseline is not None and np.isclose(v, baseline)
            ax.plot(d["Temp"], d[col], color=colors[v], lw=2.6 if is_base else 1.4, marker="o", ms=3.5,
                    label=vfmt(v) + (" (baseline)" if is_base else ""), zorder=4 if is_base else 2)
        ax.set_xlabel("Temperature (°C)")
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        if col == "CUE":
            ax.set_ylim(0, 1)
        elif zero:
            ax.set_ylim(bottom=0)
        ax.spines[["top", "right"]].set_visible(False)
    handles, labels = axes.flat[0].get_legend_handles_labels()
    fig.tight_layout(rect=(0, 0.07, 1, 1))
    fig.legend(handles, labels, loc="lower center", ncols=len(vals), framealpha=0.85, fontsize=8,
               title=param, title_fontsize=8.5)
    fig.savefig(f"{out}.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


def heatmaps_2d(summary, param, baseline, out):
    vals = sorted(summary[param].unique())            # most negative (strongest mitigation) first -> top row
    temps = sorted(summary["Temp"].unique())
    specs = [
        ("deg_pct", "Substrate degraded (%)", "YlOrBr", "{:.0f}"),
        ("biomC", "Biomass C per grid cell", "YlOrBr", "{:.1f}"),
        ("EC_end", "Energy : C ratio, end of run", "YlOrBr", "{:.0f}"),
        ("zero_enz_prod_days", "Days with zero enzyme production", "Oranges", "{:.0f}"),
    ]
    fig, axes = plt.subplots(2, 2, figsize=(12, 9))
    fig.suptitle(f"Temperature x {param}: End-of-Run Outcomes", fontsize=11, fontweight="bold")
    for ax, (col, title, cmap, f) in zip(axes.flat, specs):
        mat = summary.pivot(index=param, columns="Temp", values=col).loc[vals, temps].values
        vmax = np.nanmax(mat) if np.nanmax(mat) > 0 else 1
        im = ax.imshow(mat, aspect="auto", cmap=cmap, vmin=0, vmax=vmax)
        ax.set_xticks(range(len(temps)))
        ax.set_xticklabels([f"{t:g}" for t in temps])
        ax.set_yticks(range(len(vals)))
        ax.set_yticklabels([vfmt(v) + ("*" if baseline is not None and np.isclose(v, baseline) else "") for v in vals])
        ax.set_xlabel("Temperature (°C)")
        ax.set_ylabel(param)
        ax.set_title(title)
        for i in range(mat.shape[0]):
            for j in range(mat.shape[1]):
                if np.isfinite(mat[i, j]):
                    ax.text(j, i, f.format(mat[i, j]), ha="center", va="center", fontsize=7,
                            color="white" if mat[i, j] > 0.6 * vmax else "black")
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
    fig.text(0.5, 0.005, "* = baseline value", fontsize=8, color="dimgray", ha="center")
    fig.tight_layout(rect=(0, 0.02, 1, 0.97))
    fig.savefig(f"{out}.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


def timeseries_2d(ts, param, baseline, ts_temps, out):
    rows = [
        ("maint_Ea_effective", "Effective maintenance Ea", "J/mol"),
        ("HSP_prop", "HSP proportion of biomass", ""),
        ("energy_to_C", "Energy : C ratio", "living biomass"),
    ]
    vals = sorted(ts[param].unique())
    colors = _param_colors(vals)
    fig, axes = plt.subplots(len(rows), len(ts_temps), figsize=(5.2 * len(ts_temps), 3.6 * len(rows)),
                             sharex=True, squeeze=False)
    fig.suptitle(f"Trajectories at Selected Temperatures, by {param}", fontsize=11, fontweight="bold")
    for j, T in enumerate(ts_temps):
        sub = ts[np.isclose(ts["Temp"], T)]
        for i, (col, title, ylabel) in enumerate(rows):
            ax = axes[i][j]
            piv = sub.pivot_table(index="day", columns=param, values=col, aggfunc="mean")
            for v in vals:
                is_base = baseline is not None and np.isclose(v, baseline)
                ax.plot(piv.index, piv[v], color=colors[v], lw=2.4 if is_base else 1.3,
                        label=vfmt(v) + (" (baseline)" if is_base else ""), zorder=4 if is_base else 2)
            ax.set_title(f"{title}  |  {T:g}°C")
            ax.set_ylabel(ylabel)
            if i == len(rows) - 1:
                ax.set_xlabel("Day")
            ax.spines[["top", "right"]].set_visible(False)
    handles, labels = axes[0][0].get_legend_handles_labels()
    fig.tight_layout(rect=(0, 0.06, 1, 0.97))
    fig.legend(handles, labels, loc="lower center", ncols=len(vals), framealpha=0.85, fontsize=8,
               title=param, title_fontsize=8.5)
    fig.savefig(f"{out}.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("sweep_dir")
    ap.add_argument("--param", default=None, help="second swept parameter column (2D mode)")
    ap.add_argument("--baseline", type=float, default=None, help="baseline value of --param, highlighted in 2D mode")
    ap.add_argument("--ts-temps", type=float, nargs="+", default=[28, 36, 40],
                    help="temperatures shown in the 2D trajectory figure (default: 28 36 40)")
    a = ap.parse_args()
    ts = pd.read_csv(os.path.join(a.sweep_dir, "timeseries.csv"))
    if a.param is None:
        summary = pd.read_csv(os.path.join(a.sweep_dir, "summary.csv"), index_col=0)
        summary.index.name = "Temp"
        end_of_run_1d(summary, os.path.join(a.sweep_dir, "temp_end_of_run"))
        timeseries_1d(ts, os.path.join(a.sweep_dir, "temp_timeseries"))
        print("wrote temp_end_of_run.png and temp_timeseries.png in", a.sweep_dir)
    else:
        summary = pd.read_csv(os.path.join(a.sweep_dir, "summary.csv"))
        end_of_run_2d(summary, a.param, a.baseline, os.path.join(a.sweep_dir, "temp_end_of_run"))
        heatmaps_2d(summary, a.param, a.baseline, os.path.join(a.sweep_dir, "temp_heatmaps"))
        timeseries_2d(ts, a.param, a.baseline, a.ts_temps, os.path.join(a.sweep_dir, "temp_timeseries"))
        print("wrote temp_end_of_run.png, temp_heatmaps.png and temp_timeseries.png in", a.sweep_dir)


if __name__ == "__main__":
    main()
