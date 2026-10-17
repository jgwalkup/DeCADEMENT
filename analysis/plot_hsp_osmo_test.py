"""Figure and table for analysis/hsp_osmo_test.py output.

    python analysis/plot_hsp_osmo_test.py OUT_DIR

Writes OUT_DIR/hsp_osmo_test.png (rows: substrate remaining, biomass C, E:C, cumulative drought kills, HSP share of
biomass, effective maintenance Ea; columns: climates) and OUT_DIR/hsp_osmo_test.csv (means over seeds).
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

CLIM = [("ambient_temp", "Ambient"), ("moderate_warm", "+5 °C"), ("extreme_warm", "+10 °C")]
STYLE = {"baseline": ("0.1", "--", 2.6),
         "HSP effect off": ("#6baed6", "-", 1.6), "osmolyte effect off": ("#2171b5", "-", 1.6),
         "both effects off": ("#08306b", "-", 1.6),
         "HSP production off": ("#fc9272", "-", 1.6), "osmolyte production off": ("#de2d26", "-", 1.6),
         "both production off": ("#67000d", "-", 1.6)}


def main(out):
    t = pd.read_csv(os.path.join(out, "timeseries.csv.gz"))
    t["degraded"] = 100 - t.substrate_remaining_pct
    t = t.sort_values(["variant", "climate", "seed", "doy"])
    t["kills_cum"] = t.groupby(["variant", "climate", "seed"]).killed.cumsum()

    def f(g):
        g = g.sort_values("doy")
        dg = g.degraded.to_numpy(); d = g.doy.to_numpy()
        S = g.substrate_remaining_pct.to_numpy() > 5
        lim = ((g.R_endo_C > 0) | (g.d_r_consti < 0.999) | (g.d_r_induci < 0.999)).to_numpy()
        w = g.C_uptake.clip(lower=0)
        return pd.Series(dict(
            degraded_365=dg[-1], day_50pct=d[np.argmax(dg >= 50)] if (dg >= 50).any() else np.nan,
            day_95pct=d[np.argmax(dg >= 95)] if (dg >= 95).any() else np.nan,
            limited_days_substrate=int((lim & S).sum()), first_limited=d[lim].min() if lim.any() else np.nan,
            min_E_to_C_substrate=g.energy_to_C.to_numpy()[S].min() if S.any() else np.nan,
            biomass_max=g.biomass_C.max(), drought_kills=g.killed.sum(),
            mean_drought_tol_dry=g.loc[g.Psi < -1, "Drought_tol"].mean(), HSP_share_mean=g.HSP_prop.mean(),
            maint_Ea_mean=g.maint_Ea_effective.mean(),
            CUE_annual=np.average(g.CUE.fillna(0), weights=w) if w.sum() > 0 else np.nan))
    m = t.groupby(["variant", "climate", "seed"]).apply(f, include_groups=False).groupby(["variant", "climate"]).mean()
    m = m.reindex(pd.MultiIndex.from_product([list(STYLE), [c for c, _ in CLIM]]), level=None).dropna(how="all")
    m.round(3).to_csv(os.path.join(out, "hsp_osmo_test.csv"))

    rows = [("substrate_remaining_pct", "Substrate remaining (%)"), ("biomass_C", "Biomass C per grid cell"),
            ("energy_to_C", "E:C"), ("kills_cum", "Drought kills (cells, cumulative)"),
            ("HSP_prop", "HSP share of biomass"), ("maint_Ea_effective", "Maintenance Ea (J/mol)")]
    fig, ax = plt.subplots(len(rows), 3, figsize=(14, 17), sharex=True, sharey="row")
    dry = t[t.Psi < -4].doy
    for j, (cl, cn) in enumerate(CLIM):
        for v, (c, ls, lw) in STYLE.items():
            g = t[(t.variant == v) & (t.climate == cl)].groupby("doy").mean(numeric_only=True)
            if not len(g):
                continue
            for i, (var, lab) in enumerate(rows):
                ax[i, j].plot(g.index, g[var], color=c, ls=ls, lw=lw, label=v)
        ax[0, j].set_title(cn)
        for i, (var, lab) in enumerate(rows):
            if len(dry):
                ax[i, j].axvspan(dry.min(), dry.max(), color="0.92", zorder=0, lw=0)
            ax[i, j].grid(alpha=0.3, lw=0.6)
            ax[i, j].spines[["top", "right"]].set_visible(False)
            if j == 0:
                ax[i, j].set_ylabel(lab)
        ax[-1, j].set_xlabel("Day of year")
    ax[0, 0].legend(fontsize=8, frameon=False, loc="lower left")
    fig.suptitle("HSP and osmolyte effects and production: bacterial minimal model, first year, current defaults "
                 "(MMRT, R_endo on; mean of 2 seeds; shaded: Psi below −4 MPa)\nblue: effect removed, production and "
                 "costs kept; red: production removed (no cost, no effect)", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(os.path.join(out, "hsp_osmo_test.png"), dpi=130)
    print(m.round(2).to_string())


if __name__ == "__main__":
    main(sys.argv[1])
