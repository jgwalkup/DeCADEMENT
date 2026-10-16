"""
Figures and summary table for maint_erosion_sweep.py output.

Usage:
    python plot_maint_erosion.py ../analysis_output/maint_erosion_sweep
Writes <dir>/maint_erosion.png (end-of-run outcomes vs the multiple of baseline, one column per temperature, erosion off
vs on), <dir>/maint_erosion_timeseries.png (daily biomass and energy:C for selected multiples) and
<dir>/maint_erosion_summary.csv (per temperature and erosion setting: collapse multiple, values at baseline, unpaid energy).
"""
import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

plt.rcParams.update({"font.size": 9, "axes.titlesize": 10, "axes.labelsize": 9, "legend.fontsize": 9, "figure.dpi": 150,
                     "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True, "grid.color": "#e6e5e0",
                     "grid.linewidth": 0.6, "axes.edgecolor": "#8a8981", "xtick.color": "#5f5e57", "ytick.color": "#5f5e57"})
STYLE = {0: dict(color="#2a78d6", ls="-", marker="o", label="erosion off"),
         1: dict(color="#eb6834", ls="--", marker="s", label="erosion on")}
METRICS = [("deg_pct", "Substrate degraded (%)"), ("biomC", "Biomass C, day 365"), ("CUE", "CUE (last 30 days)"),
           ("alive", "Cells alive, day 365"), ("EC_end", "Energy:C, day 365"),
           ("E_unpaid_cum", "Unpaid energy\n(per cell, year)"), ("erosion_C_cum", "C eroded\n(per cell, year)")]
LOST = 100  # collapse: >= 100 of 400 cells lost (stage-A convention)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("dir")
    ap.add_argument("--ts-mult", type=float, nargs="+", default=[1, 2, 3, 4])
    a = ap.parse_args()
    d = pd.read_csv(os.path.join(a.dir, "summary.csv"))
    param = d["param"].iloc[0]
    base = float((d["value"] / d["mult"]).iloc[0])
    temps = sorted(d["Temp"].unique())
    g = d.groupby(["Temp", "erosion", "mult"])
    mean, lo, hi = g.mean(numeric_only=True), g.min(numeric_only=True), g.max(numeric_only=True)

    # ── End-of-run outcomes ───────────────────────────────────────────────────────────────────────────────────
    fig, axes = plt.subplots(len(METRICS), len(temps), figsize=(3.0 * len(temps), 1.9 * len(METRICS)), sharex=True,
                             sharey="row", squeeze=False)
    for c, t in enumerate(temps):
        for r, (m, lab) in enumerate(METRICS):
            ax = axes[r][c]
            for e in sorted(d["erosion"].unique()):
                if (t, e) not in mean.index.droplevel(2):
                    continue
                x = mean.loc[(t, e)].index.values
                ax.fill_between(x, lo.loc[(t, e)][m], hi.loc[(t, e)][m], color=STYLE[e]["color"], alpha=0.15, lw=0)
                ax.plot(x, mean.loc[(t, e)][m], color=STYLE[e]["color"], ls=STYLE[e]["ls"], marker=STYLE[e]["marker"],
                        ms=3.5, lw=1.6, label=STYLE[e]["label"])
            ax.axvline(1, color="#8a8981", lw=0.8, ls=":")
            if m == "CUE":
                ax.set_ylim(0, 0.5)  # CUE is floored at 0 by the model; blank where the population is extinct
            if r == 0:
                ax.set_title(f"{t:g} °C")
            if c == 0:
                ax.set_ylabel(lab)
            if r == len(METRICS) - 1:
                ax.set_xlabel(f"{param} (× baseline {base:g})")
    axes[0][0].legend(loc="best", frameon=False)
    fig.suptitle(f"{param} sensitivity, energy-deficit erosion off vs on (mean of seeds; band = seed range)", y=1.0)
    fig.tight_layout()
    fig.savefig(os.path.join(a.dir, "maint_erosion.png"), bbox_inches="tight")
    plt.close(fig)

    # ── Summary table ─────────────────────────────────────────────────────────────────────────────────────────
    rows = []
    for (t, e), sub in mean.groupby(level=[0, 1]):
        sub = sub.droplevel([0, 1])
        lost = sub.index[(400 - sub["alive"]) >= LOST]
        b = sub.loc[1.0] if 1.0 in sub.index else sub.iloc[(np.abs(sub.index - 1)).argmin()]
        rows.append(dict(Temp=t, erosion=e, collapse_mult=float(lost.min()) if len(lost) else np.nan,
                         deg_pct_at_1x=b["deg_pct"], biomC_at_1x=b["biomC"], CUE_at_1x=b["CUE"], EC_end_at_1x=b["EC_end"],
                         first_mult_with_deficit=float(sub.index[sub["deficit_days"] > 0].min()) if (sub["deficit_days"] > 0).any() else np.nan,
                         unpaid_E_at_max_mult=sub["E_unpaid_cum"].iloc[-1], eroded_C_at_max_mult=sub["erosion_C_cum"].iloc[-1]))
    s = pd.DataFrame(rows)
    s.to_csv(os.path.join(a.dir, "maint_erosion_summary.csv"), index=False)
    print(s.to_string(index=False, float_format=lambda v: f"{v:.4g}"))

    # ── Time series for selected multiples ────────────────────────────────────────────────────────────────────
    tp = os.path.join(a.dir, "timeseries.csv.gz")
    if os.path.exists(tp):
        ts = pd.read_csv(tp, usecols=["mult", "Temp", "erosion", "seed", "day", "biomass_C", "energy_to_C", "alive_cells"])
        ts = ts[ts["mult"].round(6).isin([round(m, 6) for m in a.ts_mult])]
        tm = ts.groupby(["Temp", "erosion", "mult", "day"]).mean(numeric_only=True)
        vars_ = [("biomass_C", "Biomass C"), ("energy_to_C", "Energy:C"), ("alive_cells", "Cells alive")]
        fig, axes = plt.subplots(len(vars_), len(temps), figsize=(3.0 * len(temps), 2.1 * len(vars_)), sharex=True,
                                 sharey="row", squeeze=False)
        shades = dict(zip(sorted(a.ts_mult), np.linspace(1.0, 0.35, len(a.ts_mult))))
        for c, t in enumerate(temps):
            for r, (v, lab) in enumerate(vars_):
                ax = axes[r][c]
                for e in (0, 1):
                    for m in sorted(a.ts_mult):
                        if (t, e, m) not in tm.index.droplevel(3):
                            continue
                        y = tm.loc[(t, e, m)][v]
                        ax.plot(y.index, y.values, color=STYLE[e]["color"], ls=STYLE[e]["ls"], lw=1.2, alpha=shades[m],
                                label=f"{STYLE[e]['label']}, {m:g}×" if (r == 0 and c == 0) else None)
                if r == 0:
                    ax.set_title(f"{t:g} °C")
                if c == 0:
                    ax.set_ylabel(lab)
                if r == len(vars_) - 1:
                    ax.set_xlabel("Day")
        fig.legend(loc="lower center", ncol=len(a.ts_mult), frameon=False, bbox_to_anchor=(0.5, -0.06))
        fig.suptitle(f"Daily trajectories at {', '.join(f'{m:g}' for m in sorted(a.ts_mult))} × baseline {param} (mean of seeds; lighter = higher multiple)", y=1.0)
        fig.tight_layout()
        fig.savefig(os.path.join(a.dir, "maint_erosion_timeseries.png"), bbox_inches="tight")
        plt.close(fig)


if __name__ == "__main__":
    main()
