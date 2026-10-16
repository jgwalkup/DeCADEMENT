"""
Daily CUE (and energy:C) trajectories for a set of named parameter sets, with and without HSP production,
across temperature. Runs the minimal scenario in-process (run_one) on a process pool, then plots:
  {cue,energy_to_C,biomass_C}_timeseries_by_set.png   one panel per parameter set, one line per temperature (HSP run)
  {cue,energy_to_C,biomass_C}_timeseries_by_temp.png  one panel per temperature, one line per parameter set (solid = HSP, dashed = no-HSP control)

Usage:
    python cue_timeseries.py ../input_minimal ../analysis_output/cue_timeseries [--seeds 0 1 2] [--plot-only]
The parameter sets are defined in SETS below. Writes <out_dir>/timeseries.csv (long format, seed-level).
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
from scenario_utils import scenario_color_map  # noqa: E402

plt.rcParams.update({
    "font.size": 10, "axes.titlesize": 11, "axes.labelsize": 10,
    "legend.fontsize": 9, "figure.dpi": 150,
})

INDUCI = {"HSP_Induci_Prod": 0.3, "HSP_Consti_Prod": 0.0003}
OFF = {"HSP_Induci_Prod": 0.0, "HSP_Consti_Prod": 0.0}
# name -> (parameters.csv overrides, E_resp for Mon_IdealSub or None)
SETS = {
    "Original baseline": ({"HSP_Induci_Prod": 0.003, "HSP_Consti_Prod": 0.0015}, None),
    "Inducible baseline": ({**INDUCI}, None),
    "Window: Const 2, slope -60k, cost 300": ({**INDUCI, "Maint_E_HSP_Ea_Slope": -60000, "HSP_prod_E_cost": 300}, None),
    "E_max route: Const 1.5, E_max 1050, slope -90k, cost 200": (
        {**INDUCI, "Maint_E_Rate_Const": 1.5, "Maint_E_HSP_Ea_Slope": -90000, "HSP_prod_E_cost": 200,
         "Erat_b": 350, "E_max_ratio_b": 1050, "E_max_ratio_f": 1050}, None),
    "E_resp 1000, Const 2.5, slope -90k, cost 300": (
        {**INDUCI, "Maint_E_Rate_Const": 2.5, "Maint_E_HSP_Ea_Slope": -90000, "HSP_prod_E_cost": 300}, 1000),
    "E_resp 900, Const 2, slope -90k, cost 200": (
        {**INDUCI, "Maint_E_HSP_Ea_Slope": -90000, "HSP_prod_E_cost": 200}, 900),
}
TEMPS = [4, 8, 12, 16, 20, 24, 28, 32, 36]
FOCUS = [4, 16, 28, 36]
KEEP = ["day", "CUE", "energy_to_C", "biomass_C"]


def _job(args):
    folder, changes, seed, temp, name, hsp = args
    _, series = run_one(folder, changes, 20, 365, seed, temp=temp)
    d = pd.DataFrame(series)[KEEP]
    d.insert(0, "hsp", hsp)
    d.insert(0, "set", name)
    d.insert(0, "Temp", temp)
    d.insert(0, "seed", seed)
    return d


def make_folder(src, e_resp, tmp, tag):
    if e_resp is None:
        return src
    w = os.path.join(tmp, tag)
    shutil.copytree(src, w)
    p = os.path.join(w, "monomer_energy.csv")
    out = []
    for l in open(p).read().split("\n"):
        c = l.split(",")
        out.append(f"{c[0]},{e_resp},{c[2]},{c[3]}" if c[0] == "Mon_IdealSub" else l)
    open(p, "w").write("\n".join(out))
    return w


def smooth(d, col):
    return d.sort_values("day")[col].rolling(7, center=True, min_periods=1).mean()


VARS = {
    # column: (file prefix, title name, y label, ymin, ymax)
    "CUE": ("cue", "CUE", "CUE", 0, 1),
    "energy_to_C": ("energy_to_C", "Energy : C ratio", "Energy : C (living biomass)", 0, None),
    "biomass_C": ("biomass_C", "Biomass C", "Biomass C (per grid cell)", 0, None),
}


def plots(ts, out_dir, n_seeds):
    names = list(SETS)
    labels = [str(t) for t in TEMPS]
    tcolor = scenario_color_map(labels, palette="temperature")
    scol = scenario_color_map([str(i) for i in range(len(names))], palette="discrete")
    for var, (prefix, title, ylab, ymin, ymax) in VARS.items():
        m = ts.groupby(["set", "hsp", "Temp", "day"], as_index=False)[[var]].mean()
        top = ymax if ymax is not None else float(m[var].max()) * 1.05

        # ---- one panel per set, lines per temperature ----
        fig, axes = plt.subplots(2, 3, figsize=(14, 8), sharex=True, sharey=True)
        fig.suptitle(f"{title} over time by parameter set (HSP runs, mean of {n_seeds} seeds, 7-day mean)",
                     fontsize=11, fontweight="bold")
        for ax, nm in zip(axes.flat, names):
            for t in TEMPS:
                d = m[(m["set"] == nm) & (m["hsp"]) & (m["Temp"] == t)]
                ax.plot(d.sort_values("day")["day"], smooth(d, var), color=tcolor[str(t)], lw=1.4, label=f"{t}")
            ax.set_title(nm, fontsize=9)
            ax.set_ylim(ymin, top)
            ax.spines[["top", "right"]].set_visible(False)
        for ax in axes[1]:
            ax.set_xlabel("Day")
        for ax in axes[:, 0]:
            ax.set_ylabel(ylab)
        h, l = axes[0][0].get_legend_handles_labels()
        fig.legend(h, l, loc="lower center", ncols=len(TEMPS), title="Temperature (°C)", framealpha=0.85,
                   bbox_to_anchor=(0.5, -0.07))
        fig.tight_layout()
        fig.savefig(os.path.join(out_dir, f"{prefix}_timeseries_by_set.png"), dpi=150, bbox_inches="tight")
        plt.close(fig)

        # ---- one panel per temperature, lines per set ----
        fig, axes = plt.subplots(2, 2, figsize=(13, 8.5), sharex=True, sharey=var != "biomass_C")
        fig.suptitle(f"{title} over time: HSP (solid) vs no-HSP control (dashed), mean of {n_seeds} seeds, 7-day mean",
                     fontsize=11, fontweight="bold")
        for ax, t in zip(axes.flat, FOCUS):
            for i, nm in enumerate(names):
                for hsp, ls, lw in ((True, "-", 1.7), (False, "--", 1.0)):
                    d = m[(m["set"] == nm) & (m["hsp"] == hsp) & (m["Temp"] == t)]
                    ax.plot(d.sort_values("day")["day"], smooth(d, var), color=scol[str(i)], ls=ls, lw=lw,
                            alpha=1 if hsp else 0.7, label=nm if hsp else None)
            ax.set_title(f"{t}°C")
            ax.set_ylim(ymin, None if var == "biomass_C" else top)
            ax.spines[["top", "right"]].set_visible(False)
        for ax in axes[1]:
            ax.set_xlabel("Day")
        for ax in axes[:, 0]:
            ax.set_ylabel(ylab)
        h, l = axes[0][0].get_legend_handles_labels()
        fig.legend(h, l, loc="lower center", ncols=2, framealpha=0.85, bbox_to_anchor=(0.5, -0.07))
        fig.tight_layout()
        fig.savefig(os.path.join(out_dir, f"{prefix}_timeseries_by_temp.png"), dpi=150, bbox_inches="tight")
        plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("input_folder")
    ap.add_argument("out_dir")
    ap.add_argument("--seeds", type=int, nargs="+", default=[0, 1, 2])
    ap.add_argument("--plot-only", action="store_true")
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    csv = os.path.join(a.out_dir, "timeseries.csv")
    if not a.plot_only:
        src = os.path.abspath(a.input_folder)
        tmp = tempfile.mkdtemp()
        jobs = []
        for i, (nm, (ch, er)) in enumerate(SETS.items()):
            folder = make_folder(src, er, tmp, f"in{i}")
            for t in TEMPS:
                for s in a.seeds:
                    jobs.append((folder, dict(ch), s, t, nm, True))
                    jobs.append((folder, {**ch, **OFF}, s, t, nm, False))
        print(f"{len(jobs)} runs", file=sys.stderr, flush=True)
        with Pool(os.cpu_count()) as pool:
            parts = pool.map(_job, jobs, chunksize=1)
        pd.concat(parts, ignore_index=True).to_csv(csv, index=False)
        shutil.rmtree(tmp, ignore_errors=True)
    ts = pd.read_csv(csv)
    plots(ts, a.out_dir, ts["seed"].nunique())
    print("wrote figures in", a.out_dir)


if __name__ == "__main__":
    main()
