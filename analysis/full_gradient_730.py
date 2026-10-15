"""
730-day temperature gradient of the full multi-taxon, multi-substrate model, with and without HSP production (no-HSP control),
seed-replicated (each seed is a different random community). Flat climate (Psi -0.5, FM 50) as in full_model_test.py; one pulse
of 730 days (no substrate resupply), like the minimal-model gradient in gradient_730.py.

Usage:
    python full_gradient_730.py ../input_e1500_enz2 ../analysis_output/full_gradient_730 [--seeds 0 1 2 3 4] [--side 20]
        [--temps 4 8 ...] [--days 730] [--minimal ../analysis_output/gradient_730_v2/timeseries.csv] [--plot-only]
Writes <out_dir>/timeseries.csv (seed-level, long format) and the figures
  full_gradient_730_summary.png     end-of-run outcomes vs temperature, HSP (solid) vs control (dashed)
  full_gradient_730_timeseries.png  CUE, energy:C, biomass C trajectories at selected temperatures
  full_vs_minimal_hsp_effect.png    HSP effect (% change vs control) in the full model vs the minimal model
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
from full_model_test import run  # noqa: E402
from sweep_minimal import set_params  # noqa: E402

plt.rcParams.update({"font.size": 10, "axes.titlesize": 11, "axes.labelsize": 10, "legend.fontsize": 9, "figure.dpi": 150})
TEMPS = [4, 8, 12, 16, 20, 24, 28, 32, 36, 40]
SHOW = [4, 12, 20, 28, 36, 40]
KEEP = ["biomass_C", "alive_rows", "taxa_alive", "energy_to_C", "CUE", "substrate_remaining_pct", "net_energy", "HSP_prop", "maint_Ea"]
COLOR = "#1f77b4"


def make_folders(src, tmp):
    on, off = os.path.join(tmp, "hsp_on"), os.path.join(tmp, "hsp_off")
    shutil.copytree(src, on)
    shutil.copytree(src, off)
    set_params(os.path.join(off, "parameters.csv"), {"HSP_Induci_Prod": 0.0, "HSP_Consti_Prod": 0.0})
    return on, off


def end_values(ts, day):
    cols = [c for c in ("biomass_C", "energy_to_C", "alive_rows", "taxa_alive", "substrate_remaining_pct", "HSP_prop")
            if c in ts.columns]
    a = ts[ts["day"] == day].groupby(["hsp", "Temp"])[cols].mean()
    c = ts[(ts["day"] > day - 30) & (ts["day"] <= day)].groupby(["hsp", "Temp"])["CUE"].mean()
    return a.join(c).reset_index()


def summary_plot(ts, out_dir, n_seeds, last):
    e365, elast = end_values(ts, 365), end_values(ts, last)
    fig, axes = plt.subplots(2, 3, figsize=(14, 8), sharex=True)
    fig.suptitle(f"Full model, {last}-day temperature gradient: HSP (solid) vs no-HSP control (dashed), mean of {n_seeds} communities",
                 fontsize=11, fontweight="bold")
    panels = [("biomass_C", f"Biomass C (day {last})", "C per grid cell"),
              ("rel", "Change in biomass vs control", f"% (solid = day {last}, dotted = day 365)"),
              ("CUE", f"CUE (mean of last 30 days, day {last})", "CUE"),
              ("energy_to_C", f"Energy : C (day {last})", "living biomass"),
              ("alive_rows", f"Cell-taxa alive (day {last})", "cell-taxa"),
              ("deg", f"Substrate degraded (day {last})", "% of initial C")]
    for ax, (col, title, yl) in zip(axes.flat, panels):
        if col == "rel":
            for ee, ls in ((elast, "-"), (e365, ":")):
                h = ee[ee["hsp"]].set_index("Temp")["biomass_C"]
                c = ee[~ee["hsp"]].set_index("Temp")["biomass_C"]
                ax.plot(h.index, 100 * (h - c) / c, color=COLOR, ls=ls, marker="o" if ls == "-" else None, ms=4, lw=1.6)
            ax.axhline(0, color="gray", lw=0.8, ls="--")
        else:
            for hsp, ls in ((True, "-"), (False, "--")):
                d = elast[elast["hsp"] == hsp].sort_values("Temp")
                y = 100 - d["substrate_remaining_pct"] if col == "deg" else d[col]
                ax.plot(d["Temp"], y, color=COLOR, ls=ls, marker="o", ms=4, lw=1.6, label="HSP" if hsp else "no-HSP control")
        ax.set_title(title, fontsize=10)
        ax.set_ylabel(yl, fontsize=9)
        if col == "CUE":
            ax.set_ylim(0, 1)
        ax.spines[["top", "right"]].set_visible(False)
    for ax in axes[1]:
        ax.set_xlabel("Temperature (°C)")
    h, l = axes[0][0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncols=2, framealpha=0.85, bbox_to_anchor=(0.5, -0.02))
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "full_gradient_730_summary.png"), dpi=150, bbox_inches="tight")
    plt.close(fig)


def timeseries_plot(ts, out_dir, n_seeds, last):
    m = ts.groupby(["hsp", "Temp", "day"], as_index=False)[["CUE", "energy_to_C", "biomass_C"]].mean()
    rows = [("CUE", "CUE", (0, 1)), ("energy_to_C", "Energy : C", (0, None)), ("biomass_C", "Biomass C", (0, None))]
    fig, axes = plt.subplots(3, len(SHOW), figsize=(3.1 * len(SHOW), 8.5), sharex=True, squeeze=False)
    fig.suptitle(f"Full model: HSP (solid) vs no-HSP control (dashed), {last} days, mean of {n_seeds} communities, 7-day mean",
                 fontsize=11, fontweight="bold")
    for i, (col, lab, yl) in enumerate(rows):
        for j, t in enumerate(SHOW):
            ax = axes[i][j]
            for hsp, ls, c in ((False, "--", "gray"), (True, "-", COLOR)):
                d = m[(m["hsp"] == hsp) & (m["Temp"] == t)].sort_values("day")
                ax.plot(d["day"], d[col].rolling(7, center=True, min_periods=1).mean(), color=c, ls=ls, lw=1.5,
                        label="HSP" if hsp else "no HSP (control)")
            ax.axvline(365, color="lightgray", lw=0.8)
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
    fig.savefig(os.path.join(out_dir, "full_gradient_730_timeseries.png"), dpi=150, bbox_inches="tight")
    plt.close(fig)


def comparison_plot(ts, minimal_csv, out_dir, last):
    mn = pd.read_csv(minimal_csv)
    mn = mn[mn["set"] == "E_resp 1500"].rename(columns={"alive_cells": "alive_rows"})
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.6), sharex=True)
    fig.suptitle(f"HSP effect relative to the no-HSP control: full model ({ts['seed'].nunique()} communities) vs minimal model "
                 f"({mn['seed'].nunique()} seeds); solid = day {last}, dotted = day 365", fontsize=11, fontweight="bold")
    for lab, d, c in (("full model", ts, "#1f77b4"), ("minimal model", mn, "#ff7f0e")):
        for day, ls in ((last, "-"), (365, ":")):
            e = end_values(d, day)
            for ax, col in zip(axes, ("biomass_C", "deg", "CUE")):
                h, k = e[e["hsp"]].set_index("Temp"), e[~e["hsp"]].set_index("Temp")
                if col == "deg":
                    y = (k["substrate_remaining_pct"] - h["substrate_remaining_pct"])  # extra % of initial C degraded
                elif col == "CUE":
                    y = h["CUE"] - k["CUE"]
                else:
                    y = 100 * (h[col] - k[col]) / k[col]
                ax.plot(y.index, y.values, color=c, ls=ls, marker="o" if ls == "-" else None, ms=4, lw=1.6,
                        label=lab if ls == "-" else None)
    titles = ["Biomass C: % change vs control", "Extra substrate degraded (percentage points of initial C)", "CUE: difference vs control"]
    for ax, t in zip(axes, titles):
        ax.axhline(0, color="gray", lw=0.8, ls="--")
        ax.set_title(t, fontsize=10)
        ax.set_xlabel("Temperature (°C)")
        ax.spines[["top", "right"]].set_visible(False)
    axes[0].legend(frameon=False)
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "full_vs_minimal_hsp_effect.png"), dpi=150, bbox_inches="tight")
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("input_folder")
    ap.add_argument("out_dir")
    ap.add_argument("--seeds", type=int, nargs="+", default=[0, 1, 2, 3, 4])
    ap.add_argument("--side", type=int, default=20)
    ap.add_argument("--temps", type=float, nargs="+", default=TEMPS)
    ap.add_argument("--days", type=int, default=730)
    ap.add_argument("--minimal", default=os.path.join(HERE, "..", "analysis_output", "gradient_730_v2", "timeseries.csv"))
    ap.add_argument("--plot-only", action="store_true")
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    csv = os.path.join(a.out_dir, "timeseries.csv")
    if not a.plot_only:
        tmp = tempfile.mkdtemp()
        on, off = make_folders(os.path.abspath(a.input_folder), tmp)
        jobs = [(f, t, a.side, a.days, s, 1.0, 1.0) for t in a.temps for s in a.seeds for f in (on, off)]
        print(f"{len(jobs)} runs of {a.days} days", file=sys.stderr, flush=True)
        with Pool(min(len(jobs), os.cpu_count())) as p:
            parts = p.map(run, jobs, chunksize=1)
        shutil.rmtree(tmp, ignore_errors=True)
        ts = pd.concat(parts, ignore_index=True)
        ts["hsp"] = ts["input"] == "hsp_on"
        ts = ts.rename(columns={"Temp": "Temp"})
        ts.to_csv(csv, index=False)
    ts = pd.read_csv(csv)
    n, last = ts["seed"].nunique(), int(ts["day"].max())
    summary_plot(ts, a.out_dir, n, last)
    timeseries_plot(ts, a.out_dir, n, last)
    if os.path.exists(a.minimal):
        comparison_plot(ts, a.minimal, a.out_dir, last)
    print("wrote", a.out_dir)


if __name__ == "__main__":
    main()
