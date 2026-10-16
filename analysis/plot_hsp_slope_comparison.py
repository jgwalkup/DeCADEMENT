"""
Compare two maint_erosion_sweep.py runs that differ in one setting, e.g. the HSP effect on the maintenance Ea on
(Maint_E_HSP_Ea_Slope -90000) vs off (slope 0).

Usage:
    python plot_hsp_slope_comparison.py DIR_A DIR_B OUT_DIR [--labels "HSP effect on" "HSP effect off"] [--ts-mult 1 2]
Writes OUT_DIR/hsp_slope_endpoints.png (end-of-run outcomes vs multiple of baseline; solid = endogenous respiration (R_endo) off,
dashed = on),
OUT_DIR/hsp_slope_timeseries.png (daily substrate degradation and CUE at the --ts-mult multiples, R_endo off) and
OUT_DIR/hsp_slope_comparison.csv (seed means per temperature, multiple, R_endo setting and run).
"""
import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

plt.rcParams.update({"font.size": 9, "axes.titlesize": 10, "axes.labelsize": 9, "legend.fontsize": 8, "figure.dpi": 150,
                     "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True, "grid.color": "#e6e5e0",
                     "grid.linewidth": 0.6, "axes.edgecolor": "#8a8981", "xtick.color": "#5f5e57", "ytick.color": "#5f5e57"})
COLORS = ["#2a78d6", "#eb6834"]
ENDO = r"R$_{\mathrm{endo}}$"  # endogenous respiration
MARKERS = ["o", "s"]
METRICS = [("deg_pct", "Substrate degraded (%)"), ("biomC", "Biomass C, day 365"), ("CUE", "CUE (last 30 days)"),
           ("alive", "Cells alive, day 365"), ("EC_end", "Energy:C, day 365"), ("HSP_prop", "HSP share of biomass\n(last 30 days)"),
           ("maint_Ea", "Maintenance Ea\n(last 30 days)")]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("dir_a")
    ap.add_argument("dir_b")
    ap.add_argument("out_dir")
    ap.add_argument("--labels", nargs=2, default=["HSP effect on (slope -90,000)", "HSP effect off (slope 0)"])
    ap.add_argument("--ts-mult", type=float, nargs="+", default=[1, 2])
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    runs = [pd.read_csv(os.path.join(d, "summary.csv")).assign(run=lab) for d, lab in zip((a.dir_a, a.dir_b), a.labels)]
    d = pd.concat(runs, ignore_index=True)
    param = d["param"].iloc[0]
    base = float((d["value"] / d["mult"]).iloc[0])
    temps = sorted(d["Temp"].unique())
    mean = d.groupby(["run", "Temp", "erosion", "mult"]).mean(numeric_only=True)
    mean.to_csv(os.path.join(a.out_dir, "hsp_slope_comparison.csv"))

    # ── End-of-run outcomes vs multiple ───────────────────────────────────────────────────────────────────────
    fig, axes = plt.subplots(len(METRICS), len(temps), figsize=(3.0 * len(temps), 1.85 * len(METRICS)), sharex=True,
                             sharey="row", squeeze=False)
    for c, t in enumerate(temps):
        for r, (m, lab) in enumerate(METRICS):
            ax = axes[r][c]
            for i, run in enumerate(a.labels):
                for e, ls in ((0, "-"), (1, "--")):
                    if (run, t, e) not in mean.index.droplevel(3):
                        continue
                    y = mean.loc[(run, t, e)][m]
                    ax.plot(y.index, y.values, color=COLORS[i], ls=ls, marker=MARKERS[i] if e == 0 else None, ms=3,
                            lw=1.5 if e == 0 else 1.1, label=f"{run}, {ENDO} {'off' if e == 0 else 'on'}")
            ax.axvline(1, color="#8a8981", lw=0.8, ls=":")
            if m == "CUE":
                ax.set_ylim(0, 0.55)
            if r == 0:
                ax.set_title(f"{t:g} °C")
            if c == 0:
                ax.set_ylabel(lab)
            if r == len(METRICS) - 1:
                ax.set_xlabel(f"{param} (× baseline {base:g})")
    h, l = axes[0][0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncol=2, frameon=False, bbox_to_anchor=(0.5, -0.035))
    fig.suptitle(f"{param} sweep: {a.labels[0]} vs {a.labels[1]} (mean of seeds)", y=1.0)
    fig.tight_layout()
    fig.savefig(os.path.join(a.out_dir, "hsp_slope_endpoints.png"), bbox_inches="tight")
    plt.close(fig)

    # ── Time series at selected multiples (R_endo off) ────────────────────────────────────────────────────────
    tss = []
    for dd, lab in zip((a.dir_a, a.dir_b), a.labels):
        ts = pd.read_csv(os.path.join(dd, "timeseries.csv.gz"),
                         usecols=["mult", "Temp", "erosion", "seed", "day", "CUE", "substrate_remaining_pct", "alive_cells"])
        ts = ts[(ts["erosion"] == 0) & ts["mult"].round(6).isin([round(x, 6) for x in a.ts_mult])]
        tss.append(ts.assign(run=lab))
    ts = pd.concat(tss, ignore_index=True)
    ts["degraded"] = 100 - ts["substrate_remaining_pct"]
    ts.loc[ts["alive_cells"] == 0, "CUE"] = np.nan
    tm = ts.groupby(["run", "Temp", "mult", "day"])[["degraded", "CUE"]].mean()
    rows = [(v, mu) for mu in sorted(a.ts_mult) for v in ("degraded", "CUE")]
    fig, axes = plt.subplots(len(rows), len(temps), figsize=(3.0 * len(temps), 2.0 * len(rows)), sharex=True,
                             sharey="row", squeeze=False)
    for c, t in enumerate(temps):
        for r, (v, mu) in enumerate(rows):
            ax = axes[r][c]
            for i, run in enumerate(a.labels):
                if (run, t, mu) not in tm.index.droplevel(3):
                    continue
                y = tm.loc[(run, t, mu)][v]
                ax.plot(y.index, y.values, color=COLORS[i], lw=1.6, ls="-" if i == 0 else "--", label=run)
            ax.set_ylim(0, 100 if v == "degraded" else 0.55)
            if r == 0:
                ax.set_title(f"{t:g} °C")
            if c == 0:
                ax.set_ylabel(f"{'Substrate degraded (%)' if v == 'degraded' else 'CUE'}\n{mu:g}× ({base * mu:g})")
            if r == len(rows) - 1:
                ax.set_xlabel("Day")
    h, l = axes[0][0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncol=2, frameon=False, bbox_to_anchor=(0.5, -0.03))
    fig.suptitle(f"Daily substrate degradation and CUE, endogenous respiration ({ENDO}) off (mean of seeds)", y=1.0)
    fig.tight_layout()
    fig.savefig(os.path.join(a.out_dir, "hsp_slope_timeseries.png"), bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
