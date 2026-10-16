"""Figures and table for analysis/param_sweep.py output.

    python analysis/plot_param_sweep.py OUT_DIR [--title TEXT]

For each varied parameter writes OUT_DIR/<param>_timeseries.png (rows: substrate degraded, biomass C, CUE, dead
biomass C; columns: climates; 15-day running means for CUE and dead C), plus OUT_DIR/param_sweep_annual.png (year
totals or means against the parameter value, one line per climate) and OUT_DIR/param_sweep_annual.csv (mean of seeds).
Annual CUE is the uptake-weighted mean of daily system CUE.
"""
import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

CLIM = [("ambient_temp", "Ambient"), ("moderate_warm", "+5 °C"), ("extreme_warm", "+10 °C")]
CLIM_COL = {"ambient_temp": "#3b6fb6", "moderate_warm": "#e0a030", "extreme_warm": "#c0392b"}


def disp(p):
    """Display name: 'X_min+X_max' -> 'X'."""
    return p.split("+")[0].replace("_min", "") if "+" in p else p


def load(out):
    t = pd.read_csv(os.path.join(out, "timeseries.csv.gz"))
    base = pd.read_csv(os.path.join(out, "baseline_values.csv"), header=None, index_col=0)[1]
    params = [p for p in t.param.unique() if p != "baseline"]
    b = t[t.param == "baseline"]
    t = pd.concat([t[t.param != "baseline"]] + [b.assign(param=p, value=base[p]) for p in params])
    t["degraded"] = 100 - t["substrate_remaining_pct"]
    return t, base, params


def annual(t):
    k = ["param", "value", "climate", "seed"]
    def f(g):
        w = g["C_uptake"].clip(lower=0)
        return pd.Series(dict(
            degraded_365=g.loc[g.doy.idxmax(), "degraded"],
            day_50pct=g.loc[g.degraded >= 50, "doy"].min() if (g.degraded >= 50).any() else np.nan,
            biomass_mean=g["biomass_C"].mean(),
            biomass_min_dry=g.loc[g.Psi < -4, "biomass_C"].min(),
            CUE_annual=np.average(g["CUE"].fillna(0), weights=w) if w.sum() > 0 else np.nan,
            dead_C=g["dead_C"].sum(),
            drought_kills=g["killed"].sum(),
            starved=g["starved"].sum(),
            enzyme_pool_mean=g["enzyme_pool"].mean(),
            enzyme_saturation=g["d_deg_saturation"].mean() if "d_deg_saturation" in g else np.nan,
            uptake_capacity_ratio=g["d_upt_capacity_ratio"].mean() if "d_upt_capacity_ratio" in g else np.nan,
            dry_season_biomass_change=(g.loc[g.Psi < -4, "biomass_C"].iloc[-1] - g.loc[g.Psi < -4, "biomass_C"].iloc[0])
            if (g.Psi < -4).any() else np.nan,
            min_E_to_C=g["energy_to_C"].min()))
    return t.groupby(k).apply(f, include_groups=False).reset_index()


def single(t, base, params, rows, clim, dry, a):
    """One climate: one figure per parameter, the variables as panels in two columns."""
    cl, cn = clim
    for p in params:
        vals = sorted(t.loc[t.param == p, "value"].unique())
        cols = plt.cm.viridis(np.linspace(0, 0.9, len(vals)))
        nr = int(np.ceil(len(rows) / 2))
        fig, axs = plt.subplots(nr, 2, figsize=(12, 3.0 * nr), sharex=True)
        axs = axs.ravel()
        for v, c in zip(vals, cols):
            g = t[(t.param == p) & (t.value == v) & (t.climate == cl)].groupby("doy").mean(numeric_only=True)
            isb = np.isclose(v, base[p])
            for ax, (var, lab, win) in zip(axs, rows):
                y = g[var].rolling(win, min_periods=1, center=True).mean() if win > 1 else g[var]
                ax.plot(g.index, y, color="0.15" if isb else c, lw=2.6 if isb else 1.5, ls="--" if isb else "-",
                        label=f"{v:g}" + (" (baseline)" if isb else ""))
        for ax, (var, lab, win) in zip(axs, rows):
            if len(dry):
                ax.axvspan(dry.min(), dry.max(), color="0.92", zorder=0, lw=0)
            if var == "d_upt_capacity_ratio":
                ax.set_yscale("log")
                ax.axhline(1, color="0.5", lw=0.8, ls=":")
            ax.set_ylabel(lab)
            ax.grid(alpha=0.3, lw=0.6)
            ax.spines[["top", "right"]].set_visible(False)
        for ax in axs[len(rows):]:
            ax.set_visible(False)
        for ax in axs[-2:]:
            ax.set_xlabel("Day of year")
        axs[0].legend(title=disp(p), fontsize=8, title_fontsize=9, frameon=False, loc="upper left")
        fig.suptitle(f"{disp(p)}, {cn} climate (shaded: Psi below −4 MPa)\n{a.title}", fontsize=11)
        fig.tight_layout(rect=(0, 0, 1, 0.97))
        fig.savefig(os.path.join(a.out, f"{disp(p)}_timeseries.png"), dpi=140)
        plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("out")
    ap.add_argument("--title", default="Bacterial minimal model, first year, R_endo on, beta 1; mean of 2 seeds")
    ap.add_argument("--logx", nargs="*", default=[], help="parameters to plot on a log x axis in the annual figure")
    a = ap.parse_args()
    t, base, params = load(a.out)
    clims = [c for c in CLIM if c[0] in set(t.climate)]
    dry = t[t.Psi < -4].doy
    rows = [("degraded", "Substrate degraded (%)", 1), ("biomass_C", "Biomass C per grid cell", 1),
            ("CUE", "CUE (15-day mean)", 15), ("dead_C", "Dead biomass C per grid cell\nper day (15-day mean)", 15)]
    if len(clims) == 1:
        extra = [("enzyme_pool", "Enzyme pool per grid cell", 1), ("d_deg_saturation", "Enzyme saturation S / (Km + S)", 1),
                 ("monomer_pool", "Organic monomer pool per grid cell", 1),
                 ("d_upt_capacity_ratio", "Uptake capacity / monomer pool", 1)]
        rows = rows[:2] + extra[:2] + rows[2:] + extra[2:]
        rows = [r for r in rows if r[0] in t.columns]
        single(t, base, params, rows, clims[0], dry, a)
        clims = []
    for p in (params if clims else []):
        vals = sorted(t.loc[t.param == p, "value"].unique())
        cols = plt.cm.viridis(np.linspace(0, 0.9, len(vals)))
        fig, ax = plt.subplots(len(rows), len(clims), figsize=(4.4 * len(clims), 12), sharex=True, sharey="row",
                               squeeze=False)
        for j, (cl, cn) in enumerate(clims):
            for v, c in zip(vals, cols):
                g = t[(t.param == p) & (t.value == v) & (t.climate == cl)].groupby("doy").mean(numeric_only=True)
                isb = np.isclose(v, base[p])
                for r, (var, lab, win) in enumerate(rows):
                    y = g[var].rolling(win, min_periods=1, center=True).mean() if win > 1 else g[var]
                    ax[r, j].plot(g.index, y, color="0.15" if isb else c, lw=2.6 if isb else 1.5,
                                  ls="--" if isb else "-", label=f"{v:g}" + (" (baseline)" if isb else ""))
            ax[0, j].set_title(cn)
            for r, (var, lab, win) in enumerate(rows):
                if len(dry):
                    ax[r, j].axvspan(dry.min(), dry.max(), color="0.92", zorder=0, lw=0)
                ax[r, j].grid(alpha=0.3, lw=0.6)
                ax[r, j].spines[["top", "right"]].set_visible(False)
                if j == 0:
                    ax[r, j].set_ylabel(lab)
            ax[-1, j].set_xlabel("Day of year")
        ax[2, 0].set_ylim(bottom=0)
        ax[0, 0].legend(title=disp(p), fontsize=8, title_fontsize=9, frameon=False, loc="upper left")
        ax[0, -1].text(0.98, 0.03, "shaded: Psi below −4 MPa", transform=ax[0, 2].transAxes, ha="right", fontsize=8, color="0.4")
        fig.suptitle(f"{disp(p)}: {a.title}", fontsize=12)
        fig.tight_layout(rect=(0, 0, 1, 0.97))
        fig.savefig(os.path.join(a.out, f"{disp(p)}_timeseries.png"), dpi=140)
        plt.close(fig)

    ann = annual(t)
    m = ann.groupby(["param", "value", "climate"]).mean(numeric_only=True).drop(columns="seed")
    m.round(4).to_csv(os.path.join(a.out, "param_sweep_annual.csv"))
    mets = [("degraded_365", "Degraded, day 365 (%)"), ("day_50pct", "Day 50 % degraded\n(missing: not reached)"), ("biomass_mean", "Mean biomass C per grid cell"),
            ("CUE_annual", "Annual CUE (uptake-weighted)"), ("dead_C", "Dead biomass C per grid cell (year)"),
            ("drought_kills", "Drought kills (cells, year)"), ("enzyme_pool_mean", "Mean enzyme pool per grid cell"),
            ("enzyme_saturation", "Mean enzyme saturation")]
    mets = [mm for mm in mets if m[mm[0]].notna().any()]
    fig, ax = plt.subplots(len(mets), len(params), figsize=(5.2 * len(params), 2.6 * len(mets)), sharey="row", squeeze=False)
    for j, p in enumerate(params):
        for i, (mv, lab) in enumerate(mets):
            for cl, cn in [c for c in CLIM if c[0] in set(t.climate)]:
                d = m.loc[(p, slice(None), cl), mv].droplevel([0, 2])
                ax[i, j].plot(d.index, d.values, marker="o", ms=6, lw=1.8, color=CLIM_COL[cl], label=cn)
            ax[i, j].axvline(base[p], color="0.6", lw=1, ls=":")
            if p in getattr(a, "logx", []):
                ax[i, j].set_xscale("log")
            ax[i, j].grid(alpha=0.3, lw=0.6)
            ax[i, j].spines[["top", "right"]].set_visible(False)
            if j == 0:
                ax[i, j].set_ylabel(lab)
        ax[0, j].set_title(disp(p))
        ax[-1, j].set_xlabel(f"{disp(p)} (dotted: baseline)")
    ax[0, 0].legend(frameon=False, fontsize=8)
    fig.suptitle(a.title, fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(os.path.join(a.out, "param_sweep_annual.png"), dpi=140)
    print(m.round(3).to_string())


if __name__ == "__main__":
    main()
