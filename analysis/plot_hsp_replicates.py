"""
Pool seed replicates from sweep_hsp_effect.py (seed*.csv) and
  1. screen parameter sets for the "cold-detrimental / warm-beneficial" window,
  2. plot the HSP effect on biomass C (mean +/- SD over seeds) for every grid point, and
  3. plot carbon use efficiency vs temperature (HSP run vs no-HSP control) for the shortlisted sets.

Window criteria (all on the seed means, % change in biomass C vs the no-HSP control):
  cold:  every temperature <= --cold-max has a mean change in [--cold-lo, --cold-hi]  (default -10..-5 %... see args)
  warm:  every temperature >= --warm-min has a mean change > 0 and mean - SD > 0

Usage:
    python plot_hsp_replicates.py ../analysis_output/hsp_window_replicates \\
        --rows Maint_E_HSP_Ea_Slope --cols Maint_E_Rate_Const --lines HSP_prod_E_cost
Writes into the same directory: pooled.csv, window_candidates.csv, replicates_biomass.png, replicates_cue.png.
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

PARAMS = ["Maint_E_Rate_Const", "Maint_E_HSP_Ea_Slope", "HSP_prod_E_cost"]


def vfmt(v):
    return "0" if v == 0 else (f"{v:.0f}" if abs(v) >= 100 else f"{v:.3g}")


def pool(d):
    parts = []
    for f in sorted(glob.glob(os.path.join(d, "seed*.csv"))):
        s = pd.read_csv(f)
        s["seed"] = int(os.path.basename(f)[4:-4])
        parts.append(s)
    return pd.concat(parts, ignore_index=True)


def main():
    global PARAMS
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("sweep_dir")
    ap.add_argument("--params", nargs="+", default=PARAMS,
                    help="all grid parameters that define a parameter set (default: Const, slope, cost)")
    ap.add_argument("--rows", nargs="+", default=["Maint_E_HSP_Ea_Slope"])
    ap.add_argument("--cols", nargs="*", default=["Maint_E_Rate_Const"])
    ap.add_argument("--lines", default="HSP_prod_E_cost")
    ap.add_argument("--cue-min", type=float, default=0.35, help="target mean CUE (HSP run) over --cue-tmax and below")
    ap.add_argument("--cue-tmax", type=float, default=32.0)
    ap.add_argument("--cold-max", type=float, default=10.0, help="'cold' = temperatures at or below this (C)")
    ap.add_argument("--cold-lo", type=float, default=-10.0)
    ap.add_argument("--cold-hi", type=float, default=-5.0)
    ap.add_argument("--warm-min", type=float, default=25.0, help="'warm' = temperatures at or above this (C)")
    ap.add_argument("--ymin", type=float, default=-40.0)
    ap.add_argument("--ymax", type=float, default=100.0)
    ap.add_argument("--top", type=int, default=6, help="number of shortlisted sets to show in the CUE figure")
    a = ap.parse_args()

    PARAMS = a.params
    s = pool(a.sweep_dir)
    n_seeds = s["seed"].nunique()
    s.to_csv(os.path.join(a.sweep_dir, "pooled.csv"), index=False)
    g = s.groupby(PARAMS + ["Temp"])
    m = g[["rel_biomC", "CUE", "ref_CUE", "d_CUE", "d_deg_pct", "HSP_prop"]].mean()
    sd = g[["rel_biomC", "CUE", "ref_CUE", "d_CUE"]].std()
    st = m.join(sd, rsuffix="_sd").reset_index()

    # ---- window screen -------------------------------------------------------------------
    rows = []
    for key, d in st.groupby(PARAMS):
        cold = d[d["Temp"] <= a.cold_max]
        warm = d[d["Temp"] >= a.warm_min]
        cold_ok = ((cold["rel_biomC"] >= a.cold_lo) & (cold["rel_biomC"] <= a.cold_hi)).all()
        warm_ok = ((warm["rel_biomC"] > 0) & (warm["rel_biomC"] - warm["rel_biomC_sd"] > 0)).all()
        d = d.sort_values("Temp")
        cu = d[d["Temp"] <= a.cue_tmax]["CUE"]
        pos = d[d["rel_biomC"] > 0]["Temp"]
        rows.append({**dict(zip(PARAMS, key)),
                     "cold_mean": cold["rel_biomC"].mean(), "cold_worst": cold["rel_biomC"].min(),
                     "warm_mean": warm["rel_biomC"].mean(), "warm_min_minus_sd": (warm["rel_biomC"] - warm["rel_biomC_sd"]).min(),
                     "first_T_positive": pos.min() if len(pos) else np.nan,
                     "cue_mean": cu.mean(), "cue_min": cu.min(), "cue_at_4": d[d["Temp"] == d["Temp"].min()]["CUE"].iloc[0],
                     "cold_ok": cold_ok, "warm_ok": warm_ok, "window": cold_ok and warm_ok,
                     "cue_ok": cu.mean() >= a.cue_min,
                     "warm_pos": bool((d[d["Temp"].isin([28.0, 32.0])]["rel_biomC"] > 0).all())})
    cand = pd.DataFrame(rows)
    cand["both"] = cand["window"] & cand["cue_ok"]
    cand = cand.sort_values(["both", "window", "warm_pos", "cue_mean"], ascending=[False, False, False, False])
    cand.to_csv(os.path.join(a.sweep_dir, "window_candidates.csv"), index=False)
    print(f"{n_seeds} seeds, {len(cand)} parameter sets: {int(cand['window'].sum())} meet the window, "
          f"{int(cand['cue_ok'].sum())} meet CUE >= {a.cue_min:g}, {int(cand['both'].sum())} meet both")
    print(cand.round(3).to_string(index=False))

    # ---- biomass figure ------------------------------------------------------------------
    def keys(cols):
        if not cols:
            return [()]
        return sorted(set(map(tuple, s[cols].drop_duplicates().values.tolist())))

    def lab(cols, k):
        short = lambda n: n.replace("Maint_E_HSP_Ea_", "").replace("Maint_E_", "").replace("HSP_prod_E_", "")
        return ", ".join(f"{short(n)} {vfmt(v)}" for n, v in zip(cols, k))

    def sel(df, cols, k):
        for n, v in zip(cols, k):
            df = df[df[n] == v]
        return df
    rv, cv = keys(a.rows), keys(a.cols)
    lv = sorted(s[a.lines].unique())
    colors = scenario_color_map([vfmt(v) for v in lv], palette="discrete")
    tmin, tmax = s["Temp"].min(), s["Temp"].max()
    fig, axes = plt.subplots(len(rv), len(cv), figsize=(4.4 * len(cv), 3.4 * len(rv)), sharex=True, sharey=True,
                             squeeze=False)
    fig.suptitle(f"Effect of HSP production on biomass C vs temperature (mean ± SD, {n_seeds} seeds)",
                 fontsize=11, fontweight="bold")
    for i, r in enumerate(rv):
        for j, c in enumerate(cv):
            ax = axes[i][j]
            ax.axvspan(tmin - 1, a.cold_max, color="#2166ac", alpha=0.07, lw=0)
            ax.axvspan(a.warm_min, tmax + 1, color="#d6172a", alpha=0.07, lw=0)
            ax.axhspan(a.cold_lo, a.cold_hi, xmin=0, xmax=(a.cold_max - tmin + 1) / (tmax - tmin + 2),
                       color="#2166ac", alpha=0.12, lw=0)
            ax.axhline(0, color="gray", lw=0.8, ls="--")
            for v in lv:
                d = sel(sel(st[st[a.lines] == v], a.rows, r), a.cols, c).sort_values("Temp")
                col = colors[vfmt(v)]
                ax.plot(d["Temp"], d["rel_biomC"].clip(a.ymin, a.ymax), color=col, lw=1.5, marker="o", ms=3, label=vfmt(v))
                ax.fill_between(d["Temp"], (d["rel_biomC"] - d["rel_biomC_sd"]).clip(a.ymin, a.ymax),
                                (d["rel_biomC"] + d["rel_biomC_sd"]).clip(a.ymin, a.ymax), color=col, alpha=0.12, lw=0)
            ax.set_title(f"{lab(a.rows, r)}\n{lab(a.cols, c)}", fontsize=9)
            ax.set_ylim(a.ymin, a.ymax)
            ax.spines[["top", "right"]].set_visible(False)
            if i == len(rv) - 1:
                ax.set_xlabel("Temperature (°C)")
            if j == 0:
                ax.set_ylabel("Δ biomass C vs no HSP (%)")
    h, l = axes[0][0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncols=len(lv), title=a.lines, framealpha=0.85, bbox_to_anchor=(0.5, -0.04))
    fig.text(0.5, -0.075, f"Dark blue band: target cold penalty ({a.cold_lo:g} to {a.cold_hi:g}%). "
             "20x20 grid, 365 days; control = no HSP production.", fontsize=8, color="dimgray", ha="center")
    fig.tight_layout()
    fig.savefig(os.path.join(a.sweep_dir, "replicates_biomass.png"), dpi=150, bbox_inches="tight")
    plt.close(fig)

    # ---- CUE figure for shortlisted sets ---------------------------------------------------
    short = cand.head(a.top)
    ncol = min(3, len(short))
    nrow = int(np.ceil(len(short) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(4.6 * ncol, 3.6 * nrow), sharex=True, sharey=True, squeeze=False)
    fig.suptitle(f"Carbon use efficiency, HSP vs no-HSP control (mean ± SD, {n_seeds} seeds)",
                 fontsize=11, fontweight="bold")
    for ax in axes.flat[len(short):]:
        ax.set_visible(False)
    for ax, (_, r) in zip(axes.flat, short.iterrows()):
        d = st
        for p in PARAMS:
            d = d[d[p] == r[p]]
        d = d.sort_values("Temp")
        ax.axvspan(tmin - 1, a.cold_max, color="#2166ac", alpha=0.07, lw=0)
        ax.axvspan(a.warm_min, tmax + 1, color="#d6172a", alpha=0.07, lw=0)
        ax.plot(d["Temp"], d["ref_CUE"], color="gray", lw=1.5, ls="--", marker="s", ms=3, label="no HSP (control)")
        ax.fill_between(d["Temp"], d["ref_CUE"] - d["ref_CUE_sd"], d["ref_CUE"] + d["ref_CUE_sd"], color="gray", alpha=0.15, lw=0)
        ax.plot(d["Temp"], d["CUE"], color="#1f77b4", lw=1.8, marker="o", ms=4, label="with HSP")
        ax.fill_between(d["Temp"], d["CUE"] - d["CUE_sd"], d["CUE"] + d["CUE_sd"], color="#1f77b4", alpha=0.15, lw=0)
        tag = "window+CUE" if r["both"] else ("window" if r["window"] else ("CUE ok" if r["cue_ok"] else ""))
        ax.set_title(", ".join(f"{n.replace('Maint_E_HSP_Ea_', '').replace('Maint_E_Rate_', '').replace('HSP_prod_E_', '')} {vfmt(r[n])}"
                               for n in PARAMS) + (f"  [{tag}]" if tag else ""), fontsize=8)
        ax.set_ylim(0, 1)
        ax.spines[["top", "right"]].set_visible(False)
        ax.set_xlabel("Temperature (°C)")
        ax.set_ylabel("CUE (mean of last 30 days)")
    h, l = axes[0][0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncols=2, framealpha=0.85, bbox_to_anchor=(0.5, -0.04))
    fig.tight_layout()
    fig.savefig(os.path.join(a.sweep_dir, "replicates_cue.png"), dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("wrote pooled.csv, window_candidates.csv, replicates_biomass.png, replicates_cue.png")


if __name__ == "__main__":
    main()
