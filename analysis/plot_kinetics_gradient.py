"""Figures and table for analysis/kinetics_gradient.py output.

    python analysis/plot_kinetics_gradient.py OUT_DIR

Writes OUT_DIR/kinetics_gradient.png (against temperature: enzyme and transporter Vmax on day 1, day-365 % degraded,
day of 50 % degraded, mean biomass C, annual CUE, minimum E:C) and OUT_DIR/kinetics_gradient_timeseries.png
(substrate remaining over the year at each temperature), plus OUT_DIR/kinetics_gradient.csv (means over seeds).
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

KIN = [("MMRT", "#1f77b4", "-", "o"), ("Arrhenius", "#d62728", "--", "s")]


def main(out):
    t = pd.read_csv(os.path.join(out, "timeseries.csv.gz"))
    t["degraded"] = 100 - t.substrate_remaining_pct

    def f(g):
        g = g.sort_values("doy")
        dg = g.degraded.to_numpy(); d = g.doy.to_numpy()
        w = g.C_uptake.clip(lower=0)
        return pd.Series(dict(
            enzyme_Vmax=g.d_deg_Vmax.iloc[0], uptake_Vmax=g.d_upt_Vmax.iloc[0],
            enzyme_Km=g.d_deg_Km.iloc[0], degraded_365=dg[-1],
            day_50pct=d[np.argmax(dg >= 50)] if (dg >= 50).any() else np.nan,
            day_95pct=d[np.argmax(dg >= 95)] if (dg >= 95).any() else np.nan,
            biomass_max=g.biomass_C.max(), biomass_mean=g.biomass_C.mean(), biomass_end=g.biomass_C.iloc[-1],
            CUE_annual=np.average(g.CUE.fillna(0), weights=w) if w.sum() > 0 else np.nan,
            min_E_to_C=g.energy_to_C.min(), starved=g.starved.sum(),
            maint_rate=g.maint_rate_percap.mean()))
    m = t.groupby(["kinetics", "Temp_set", "seed"]).apply(f, include_groups=False).groupby(["kinetics", "Temp_set"]).mean()
    m.round(4).to_csv(os.path.join(out, "kinetics_gradient.csv"))

    panels = [("enzyme_Vmax", "Enzyme Vmax (day 1)"), ("uptake_Vmax", "Transporter Vmax (day 1)"),
              ("degraded_365", "Degraded, day 365 (%)"), ("day_50pct", "Day 50 % degraded\n(missing: not reached)"),
              ("biomass_max", "Maximum biomass C per grid cell"), ("CUE_annual", "Annual CUE (uptake-weighted)"),
              ("biomass_end", "Biomass C at day 365"), ("min_E_to_C", "Minimum E:C")]
    fig, ax = plt.subplots(2, 4, figsize=(17, 8))
    for a, (var, lab) in zip(ax.ravel(), panels):
        for kin, c, ls, mk in KIN:
            d = m.loc[kin, var]
            a.plot(d.index, d.values, color=c, ls=ls, marker=mk, ms=5, lw=1.7, label=kin)
        a.set_ylabel(lab)
        a.set_xlabel("Temperature (°C)")
        a.grid(alpha=0.3, lw=0.6)
        a.spines[["top", "right"]].set_visible(False)
    ax[0, 0].legend(frameon=False)
    fig.suptitle("MMRT vs Arrhenius across constant temperatures (Psi -0.5 MPa): bacterial minimal model, first year, "
                 "production x1.4, Vmax_Km 13, R_endo on; mean of 2 seeds", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(os.path.join(out, "kinetics_gradient.png"), dpi=140)
    plt.close(fig)

    temps = sorted(t.Temp_set.unique())
    n = len(temps); nc = 5; nr = int(np.ceil(n / nc))
    fig, ax = plt.subplots(nr, nc, figsize=(3.4 * nc, 2.9 * nr), sharex=True, sharey=True)
    ax = ax.ravel()
    for a, T in zip(ax, temps):
        for kin, c, ls, mk in KIN:
            g = t[(t.kinetics == kin) & (t.Temp_set == T)].groupby("doy").substrate_remaining_pct.mean()
            a.plot(g.index, g.values, color=c, ls=ls, lw=1.6, label=kin)
        a.set_title(f"{T:g} °C", fontsize=10)
        a.axhline(20, color="0.7", lw=0.7, ls=":")
        a.grid(alpha=0.3, lw=0.6)
        a.spines[["top", "right"]].set_visible(False)
    for a in ax[n:]:
        a.set_visible(False)
    for a in ax[::nc]:
        a.set_ylabel("Substrate remaining (%)")
    for a in ax[(nr - 1) * nc:]:
        a.set_xlabel("Day")
    ax[0].legend(frameon=False, fontsize=8)
    fig.suptitle("Substrate remaining at constant temperatures: MMRT (solid) vs Arrhenius (dashed)", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(os.path.join(out, "kinetics_gradient_timeseries.png"), dpi=140)
    print(m.round(2).to_string())


if __name__ == "__main__":
    main(sys.argv[1])
