"""
730-day temperature gradient for named parameter sets, with and without HSP production (no-HSP control),
seed-replicated. Saves the daily trajectories and plots:
  gradient_730_summary.png           end-of-run (day 365 vs day 730) outcomes vs temperature, HSP vs control
  gradient_730_timeseries_<set>.png  CUE, energy:C, biomass C trajectories (rows) at selected temperatures (columns)

Usage:
    python gradient_730.py ../input_minimal ../analysis_output/gradient_730 [--seeds 0 1 2] [--plot-only]
Parameter sets are defined in SETS below. Writes <out_dir>/timeseries.csv (long format, seed-level).
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
sys.path.insert(0, HERE)
from sweep_minimal import run_one  # noqa: E402

plt.rcParams.update({
    "font.size": 10, "axes.titlesize": 11, "axes.labelsize": 10,
    "legend.fontsize": 9, "figure.dpi": 150,
})

BASE = {"HSP_Induci_Prod": 0.3, "HSP_Consti_Prod": 0.0003, "Maint_E_Rate_Const": 2.5,
        "Maint_E_HSP_Ea_Intercept": 45000, "Maint_E_HSP_Ea_Slope": -90000, "HSP_prod_E_cost": 300}
OFF = {"HSP_Induci_Prod": 0.0, "HSP_Consti_Prod": 0.0}
SETS = {  # name -> (overrides, E_resp, colour)
    "E_resp 1500": (BASE, 1500, "#1f77b4"),
    "E_resp 1700": (BASE, 1700, "#ff7f0e"),
}
TEMPS = [4, 8, 12, 16, 20, 24, 28, 32, 36, 40]
SHOW = [4, 12, 20, 28, 36, 40]
DAYS = 730
KEEP = ["day", "CUE", "energy_to_C", "biomass_C", "alive_cells", "substrate_remaining_pct", "net_energy", "HSP_prop"]


def _job(args):
    folder, changes, seed, temp, name, hsp, side = args
    _, series = run_one(folder, changes, side, DAYS, seed, temp=temp)
    d = pd.DataFrame(series)[KEEP]
    d.insert(0, "hsp", hsp)
    d.insert(0, "set", name)
    d.insert(0, "Temp", temp)
    d.insert(0, "seed", seed)
    return d


def make_folder(src, e_resp, tmp, tag):
    w = os.path.join(tmp, tag)
    shutil.copytree(src, w)
    p = os.path.join(w, "monomer_energy.csv")
    out = []
    for l in open(p).read().split("\n"):
        c = l.split(",")
        out.append(f"{c[0]},{e_resp},{c[2]},{c[3]}" if c[0] == "Mon_IdealSub" else l)
    open(p, "w").write("\n".join(out))
    return w


def end_values(ts, day):
    """Seed-mean values at a given day (CUE = mean of the 30 days ending there)."""
    a = ts[ts["day"] == day].groupby(["set", "hsp", "Temp"])[["biomass_C", "energy_to_C", "alive_cells",
                                                              "substrate_remaining_pct"]].mean()
    c = ts[(ts["day"] > day - 30) & (ts["day"] <= day)].groupby(["set", "hsp", "Temp"])["CUE"].mean()
    return a.join(c).reset_index()


def summary_plot(ts, out_dir, n_seeds):
    e365, e730 = end_values(ts, 365), end_values(ts, DAYS)
    names = list(SETS)
    fig, axes = plt.subplots(2, 3, figsize=(14, 8), sharex=True)
    fig.suptitle(f"730-day temperature gradient: HSP (solid) vs no-HSP control (dashed), mean of {n_seeds} seeds",
                 fontsize=11, fontweight="bold")
    panels = [("biomass_C", "Biomass C (day 730)", "C per grid cell", e730),
              ("rel", "Change in biomass vs control", "% (solid = day 730, dotted = day 365)", None),
              ("CUE", "CUE (mean of last 30 days, day 730)", "CUE", e730),
              ("energy_to_C", "Energy : C (day 730)", "living biomass", e730),
              ("alive_cells", "Cells alive (day 730)", "cells alive", e730),
              ("deg", "Substrate degraded (day 730)", "% of initial C", e730)]
    for ax, (col, title, yl, e) in zip(axes.flat, panels):
        for nm in names:
            colr = SETS[nm][2]
            if col == "rel":
                for ee, ls in ((e730, "-"), (e365, ":")):
                    h = ee[(ee["set"] == nm) & ee["hsp"]].set_index("Temp")["biomass_C"]
                    c = ee[(ee["set"] == nm) & ~ee["hsp"]].set_index("Temp")["biomass_C"]
                    ax.plot(h.index, 100 * (h - c) / c, color=colr, ls=ls, marker="o" if ls == "-" else None, ms=4,
                            lw=1.6, label=nm if ls == "-" else None)
                ax.axhline(0, color="gray", lw=0.8, ls="--")
                continue
            for hsp, ls in ((True, "-"), (False, "--")):
                d = e[(e["set"] == nm) & (e["hsp"] == hsp)].sort_values("Temp")
                y = 100 - d["substrate_remaining_pct"] if col == "deg" else d[col]
                ax.plot(d["Temp"], y, color=colr, ls=ls, marker="o", ms=4, lw=1.6, label=nm if hsp else None)
        ax.set_title(title, fontsize=10)
        ax.set_ylabel(yl, fontsize=9)
        if col == "CUE":
            ax.set_ylim(0, 1)
        ax.spines[["top", "right"]].set_visible(False)
    for ax in axes[1]:
        ax.set_xlabel("Temperature (°C)")
    h, l = axes[0][0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncols=len(names), framealpha=0.85, bbox_to_anchor=(0.5, -0.03))
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "gradient_730_summary.png"), dpi=150, bbox_inches="tight")
    plt.close(fig)


def timeseries_plots(ts, out_dir, n_seeds):
    m = ts.groupby(["set", "hsp", "Temp", "day"], as_index=False)[["CUE", "energy_to_C", "biomass_C"]].mean()
    rows = [("CUE", "CUE", (0, 1)), ("energy_to_C", "Energy : C", (0, None)), ("biomass_C", "Biomass C", (0, None))]
    for nm in SETS:
        colr = SETS[nm][2]
        fig, axes = plt.subplots(3, len(SHOW), figsize=(3.1 * len(SHOW), 8.5), sharex=True, squeeze=False)
        fig.suptitle(f"{nm}: HSP (solid) vs no-HSP control (dashed), 730 days, mean of {n_seeds} seeds, 7-day mean",
                     fontsize=11, fontweight="bold")
        for i, (col, lab, yl) in enumerate(rows):
            for j, t in enumerate(SHOW):
                ax = axes[i][j]
                for hsp, ls, c in ((False, "--", "gray"), (True, "-", colr)):
                    d = m[(m["set"] == nm) & (m["hsp"] == hsp) & (m["Temp"] == t)].sort_values("day")
                    ax.plot(d["day"], d[col].rolling(7, center=True, min_periods=1).mean(), color=c, ls=ls, lw=1.5,
                            label="HSP" if hsp else "no HSP (control)")
                ax.axvline(365, color="lightgray", lw=0.8)
                if col == "biomass_C":
                    ax.set_ylim(0, None)
                else:
                    ax.set_ylim(*yl)
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
        fig.savefig(os.path.join(out_dir, f"gradient_730_timeseries_{nm.replace(' ', '_')}.png"), dpi=150,
                    bbox_inches="tight")
        plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("input_folder")
    ap.add_argument("out_dir")
    ap.add_argument("--seeds", type=int, nargs="+", default=[0, 1, 2])
    ap.add_argument("--side", type=int, default=20, help="grid side (side x side cells)")
    ap.add_argument("--set", action="append", default=[], metavar="NAME=VALUE",
                    help="override a parameter of the HSP/control parameter set (repeatable)")
    ap.add_argument("--plot-only", action="store_true")
    a = ap.parse_args()
    BASE.update({k: float(v) for k, v in (x.split("=") for x in a.set)})
    os.makedirs(a.out_dir, exist_ok=True)
    csv = os.path.join(a.out_dir, "timeseries.csv")
    if not a.plot_only:
        src = os.path.abspath(a.input_folder)
        tmp = tempfile.mkdtemp()
        jobs = []
        for i, (nm, (ch, er, _)) in enumerate(SETS.items()):
            folder = make_folder(src, er, tmp, f"in{i}")
            for t in TEMPS:
                for s in a.seeds:
                    jobs.append((folder, dict(ch), s, t, nm, True, a.side))
                    jobs.append((folder, {**ch, **OFF}, s, t, nm, False, a.side))
        print(f"{len(jobs)} runs of {DAYS} days", file=sys.stderr, flush=True)
        with Pool(os.cpu_count()) as pool:
            parts = pool.map(_job, jobs, chunksize=1)
        pd.concat(parts, ignore_index=True).to_csv(csv, index=False)
        shutil.rmtree(tmp, ignore_errors=True)
    ts = pd.read_csv(csv)
    n = ts["seed"].nunique()
    summary_plot(ts, a.out_dir, n)
    timeseries_plots(ts, a.out_dir, n)
    print("wrote figures in", a.out_dir)


if __name__ == "__main__":
    main()
