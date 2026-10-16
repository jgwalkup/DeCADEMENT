"""Figures and table for analysis/drought_sweep.py.

    python analysis/plot_drought_sweep.py OUT_DIR

Writes OUT_DIR/drought_sweep.png (day-365 % degraded by climate; % degraded and live cells over the year, +10 C) and
OUT_DIR/drought_sweep_summary.csv (means over seeds).
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

CLIM = [("ambient_temp", "Ambient"), ("moderate_warm", "+5 °C"), ("extreme_warm", "+10 °C")]
PARAMS = [("beta", "beta_bac (drought-mortality slope)", 10.0), ("Osmo_tol_K", "Osmo_tol_K (tolerance half-saturation)", 0.05)]


def main(out):
    s = pd.read_csv(os.path.join(out, "summary.csv"))
    t = pd.read_csv(os.path.join(out, "timeseries.csv.gz"))
    base = s[(s.param == "beta") & (s.value == 10)].assign(param="Osmo_tol_K", value=0.05)
    s = pd.concat([s, base])
    tb = t[(t.param == "beta") & (t.value == 10)].assign(param="Osmo_tol_K", value=0.05)
    t = pd.concat([t, tb])
    t["degraded"] = 100 - t["substrate_remaining_pct"]
    k = ["param", "value", "climate"]
    dry = t[t.Psi < -1].groupby(k + ["seed"]).agg(drought_tol=("Drought_tol", "mean"), min_cells=("alive_cells", "min"))
    kil = t.groupby(k + ["seed"]).agg(killed=("killed", "sum"), max_biomass_C=("biomass_C", "max"))
    m = s.merge(dry.reset_index(), on=k + ["seed"]).merge(kil.reset_index(), on=k + ["seed"])
    summ = m.groupby(k)[["deg_pct", "biomC", "max_biomass_C", "killed", "min_cells", "drought_tol", "zero_enz_prod_days"]].mean()
    summ.round(3).to_csv(os.path.join(out, "drought_sweep_summary.csv"))

    fig, ax = plt.subplots(3, 2, figsize=(11, 11.5))
    for j, (p, label, b) in enumerate(PARAMS):
        vals = sorted(s.loc[s.param == p, "value"].unique(), reverse=True)   # baseline (largest) first
        cols = plt.cm.Blues(np.linspace(0.35, 1.0, len(vals)))[::-1] if p == "beta" else plt.cm.Blues(np.linspace(0.35, 1.0, len(vals)))[::-1]
        for v, c in zip(vals, cols):
            lw, ls = (2.6, "-") if v == b else (1.8, "-")
            name = f"{v:g}" + (" (baseline)" if v == b else "")
            y = [summ.loc[(p, v, cl), "deg_pct"] for cl, _ in CLIM]
            ax[0, j].plot(range(3), y, marker="o", ms=7, lw=lw, color=c if v != b else "0.25", ls=ls, label=name)
            for r, var in ((1, "degraded"), (2, "alive_cells")):
                g = t[(t.param == p) & (t.value == v) & (t.climate == "extreme_warm")].groupby("doy")[var].mean()
                ax[r, j].plot(g.index, g.values, lw=lw if v == b else 1.6, color=c if v != b else "0.25", label=name)
        ax[0, j].set_xticks(range(3), [n for _, n in CLIM])
        ax[0, j].set_ylabel("Substrate degraded, day 365 (%)")
        ax[0, j].set_title(label, fontsize=11)
        ax[0, j].legend(title=p, fontsize=8, title_fontsize=9, frameon=False)
        ax[1, j].set_ylabel("Substrate degraded (%), +10 °C")
        ax[2, j].set_ylabel("Live cells (of 400), +10 °C")
        ax[2, j].set_xlabel("Day of year")
        for r in (1, 2):
            dry_days = t[(t.climate == "extreme_warm") & (t.Psi < -4)].doy.unique()
            if len(dry_days):
                ax[r, j].axvspan(dry_days.min(), dry_days.max(), color="0.92", zorder=0, lw=0)
        for a in ax[:, j]:
            a.grid(alpha=0.3, lw=0.6)
            a.spines[["top", "right"]].set_visible(False)
    ax[1, 0].text(0.02, 0.95, "shaded: Psi below −4 MPa", transform=ax[1, 0].transAxes, fontsize=8, color="0.4", va="top")
    fig.suptitle("Drought-mortality parameters, bacterial minimal model, first year (R_endo on; mean of 2 seeds)", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(os.path.join(out, "drought_sweep.png"), dpi=150)
    print(summ.round(2).to_string())


if __name__ == "__main__":
    main(sys.argv[1])
