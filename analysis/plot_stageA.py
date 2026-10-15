"""
Figures and summary for the stage-A one-parameter sweeps re-run at the settled set (stageA_rerun/<param>/seed*.csv from
sweep_hsp_effect.py): outcomes against the swept value (as a multiple of baseline), one line per temperature,
HSP runs solid and no-HSP controls dashed.

Usage:
    python plot_stageA.py ../analysis_output/stageA_rerun
Writes <dir>/<param>/stageA_<param>.png and <dir>/stageA_summary.csv (collapse multipliers, tolerance ranges).
"""
import argparse
import glob
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from scenario_utils import scenario_color_map  # noqa: E402

plt.rcParams.update({
    "font.size": 10, "axes.titlesize": 11, "axes.labelsize": 10,
    "legend.fontsize": 9, "figure.dpi": 150,
})

BASE = {"Maint_E_Rate_Const": 2.5, "Enz_Prod": 0.002, "HSP_Consti_Prod": 0.0015, "Uptake_C_cost": 0.0375,
        "HSP_prod_E_cost": 300.0, "Enz_Prod_E_cost": 600.0, "Uptake_Prod_E_cost": 600.0}
OLD_COLLAPSE_15C = {"Maint_E_Rate_Const": 2.2, "HSP_Consti_Prod": 3.7, "HSP_prod_E_cost": 3.3, "Uptake_Prod_E_cost": 27,
                    "Enz_Prod_E_cost": 47, "Enz_Prod": 19, "Uptake_C_cost": np.nan}


def load(d, name):
    s = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(os.path.join(d, name, "seed*.csv")))], ignore_index=True)
    s["mult"] = (s[name] / BASE[name]).round(6)
    return s


def fmt(m):
    return "0" if m == 0 else f"{m:g}"


def plot_param(s, name, outdir):
    temps = sorted(s["Temp"].unique())
    color = scenario_color_map([f"{t:g}" for t in temps], palette="temperature")
    g = s.groupby(["mult", "Temp"]).mean(numeric_only=True).reset_index()
    mults = sorted(g["mult"].unique())
    xpos = {m: i for i, m in enumerate(mults)}
    g["x"] = g["mult"].map(xpos)
    g["lost"] = 400 - g["alive"]
    g["ref_lost"] = 400 - g["ref_alive"]
    panels = [("deg_pct", "ref_deg_pct", "Substrate degraded", "% of initial C"),
              ("biomC", "ref_biomC", "Biomass C", "C per grid cell"),
              ("CUE", "ref_CUE", "CUE", "last 30 days"),
              ("EC_end", "ref_EC_end", "Energy : C ratio", "living biomass, end of run"),
              ("lost", "ref_lost", "Cells lost", "of 400 cells"),
              ("rel_biomC", None, "Biomass change with HSP vs control", "%")]
    fig, axes = plt.subplots(2, 3, figsize=(14, 8), sharex=True)
    fig.suptitle(f"{name}: sweep as a multiple of baseline ({BASE[name]:g}); HSP solid, no-HSP control dashed; mean of "
                 f"{s['seed'].nunique() if 'seed' in s else 3} seeds", fontsize=11, fontweight="bold")
    for ax, (col, ref, title, yl) in zip(axes.flat, panels):
        for t in temps:
            d = g[g["Temp"] == t].sort_values("x")
            ax.plot(d["x"], d[col], color=color[f"{t:g}"], lw=1.6, marker="o", ms=3.5, label=f"{t:g}")
            if ref:
                ax.plot(d["x"], d[ref], color=color[f"{t:g}"], lw=1.0, ls="--", alpha=0.7)
        if col == "rel_biomC":
            ax.axhline(0, color="gray", lw=0.8, ls="--")
        if 1.0 in xpos:
            ax.axvline(xpos[1.0], color="lightgray", lw=1.0, zorder=0)
        ax.set_title(title, fontsize=10)
        ax.set_ylabel(yl)
        if col == "CUE":
            ax.set_ylim(0, 1)
        elif col not in ("rel_biomC",):
            ax.set_ylim(bottom=0)
        ax.spines[["top", "right"]].set_visible(False)
    for ax in axes[1]:
        ax.set_xticks(range(len(mults)), [fmt(m) for m in mults], rotation=60, ha="right", fontsize=8)
        ax.set_xlabel(f"{name} (x baseline)")
    h, l = axes[0][0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncols=len(temps), title="Temperature (°C)", framealpha=0.85,
               bbox_to_anchor=(0.5, -0.075))
    fig.text(0.5, -0.125, "Grey vertical line = baseline. 20x20 grid, 365 days, input_minimal_e1500.", fontsize=8,
             color="dimgray", ha="center")
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, name, f"stageA_{name}.png"), dpi=150, bbox_inches="tight")
    plt.close(fig)
    return g


def summarise(g, name):
    rows = []
    for t in sorted(g["Temp"].unique()):
        d = g[g["Temp"] == t].sort_values("mult")
        for arm, ec, deg, lost in (("HSP", "EC_end", "deg_pct", "lost"), ("control", "ref_EC_end", "ref_deg_pct", "ref_lost")):
            base = d[d["mult"] == 1.0]
            b_deg = float(base[deg].iloc[0]) if len(base) else np.nan
            above = d[d["mult"] >= 1.0]
            below = d[d["mult"] <= 1.0]
            hit = above[above[ec] < 1.0]
            collapse = float(hit["mult"].min()) if len(hit) else np.nan
            hit_lo = below[below[ec] < 1.0]
            collapse_lo = float(hit_lo["mult"].max()) if len(hit_lo) else np.nan
            ok = d[(d[deg] >= 0.9 * b_deg) & (d[deg] <= 1.1 * b_deg) & (d[lost] < 20)]
            rows.append(dict(param=name, Temp=t, arm=arm, baseline_deg=b_deg,
                             collapse_above_x=collapse, collapse_below_x=collapse_lo,
                             within10pct_min_x=float(ok["mult"].min()) if len(ok) else np.nan,
                             within10pct_max_x=float(ok["mult"].max()) if len(ok) else np.nan,
                             deg_at_1p5x_rel=100 * (float(d[d["mult"] == 1.5][deg].iloc[0]) / b_deg - 1)
                             if len(d[d["mult"] == 1.5]) and b_deg else np.nan))
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("sweep_dir")
    a = ap.parse_args()
    rows = []
    for name in BASE:
        if not glob.glob(os.path.join(a.sweep_dir, name, "seed*.csv")):
            continue
        s = load(a.sweep_dir, name)
        g = plot_param(s, name, a.sweep_dir)
        rows += summarise(g, name)
        g.to_csv(os.path.join(a.sweep_dir, name, "pooled_means.csv"), index=False)
    out = pd.DataFrame(rows)
    out.to_csv(os.path.join(a.sweep_dir, "stageA_summary.csv"), index=False)
    print("wrote figures and stageA_summary.csv in", a.sweep_dir)


if __name__ == "__main__":
    main()
