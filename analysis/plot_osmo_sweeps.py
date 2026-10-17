"""Figure for the osmolyte sweeps: beta (osmolytes on vs production off), Osmo_prod_E_cost, Osmo_Induci_Prod.

    python analysis/plot_osmo_sweeps.py OUT_PNG BETA_ON_DIR BETA_OFF_DIR COST_INDUCI_DIR

Rows: day 95 % degraded, energy-limited days while > 5 % substrate is left, minimum E:C then, drought kills per year;
columns: the three sweeps; colours: climates (beta column: solid osmolytes on, dashed osmolyte production off). Writes
OUT_PNG and OUT_PNG with .csv.
"""
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

CLIM = [("ambient_temp", "Ambient", "#3b6fb6"), ("moderate_warm", "+5 °C", "#e0a030"), ("extreme_warm", "+10 °C", "#c0392b")]


def summ(d):
    t = pd.read_csv(f"{d}/timeseries.csv.gz")
    b = pd.read_csv(f"{d}/baseline_values.csv", header=None, index_col=0)[1]
    out = []
    for (p, v, c, s), g in t.groupby(["param", "value", "climate", "seed"]):
        keys = [(q, b[q]) for q in b.index] if p == "baseline" else [(p, v)]
        g = g.sort_values("doy")
        dg = 100 - g.substrate_remaining_pct.to_numpy(); dd = g.doy.to_numpy()
        S = g.substrate_remaining_pct.to_numpy() > 5
        lim = ((g.R_endo_C > 0) | (g.d_r_consti < 0.999) | (g.d_r_induci < 0.999)).to_numpy()
        for q, val in keys:
            out.append(dict(param=q, value=val, climate=c, seed=s, degraded_365=dg[-1],
                            day_95pct=dd[np.argmax(dg >= 95)] if (dg >= 95).any() else np.nan,
                            limited_days=int((lim & S).sum()), min_E_to_C=g.energy_to_C.to_numpy()[S].min(),
                            drought_kills=g.killed.sum(), drought_tol_dry=g.loc[g.Psi < -1, "Drought_tol"].mean()))
    return pd.DataFrame(out).groupby(["param", "value", "climate"]).mean(numeric_only=True).drop(columns="seed")


def main(out, d_on, d_off, d_ci):
    on, off, ci = summ(d_on), summ(d_off), summ(d_ci)
    pd.concat([on.assign(osmolytes="on"), off.assign(osmolytes="production off"), ci.assign(osmolytes="on")]
              ).round(3).to_csv(out.rsplit(".", 1)[0] + ".csv")
    cols = [("beta_bac+beta_fun", "beta (drought-mortality slope)", on, off),
            ("Osmo_prod_E_cost", "Osmo_prod_E_cost", ci, None),
            ("Osmo_Induci_Prod_min+Osmo_Induci_Prod_max", "Osmo_Induci_Prod", ci, None)]
    rows = [("day_95pct", "Day 95 % degraded"), ("limited_days", "Energy-limited days\n(substrate > 5 %)"),
            ("min_E_to_C", "Minimum E:C (substrate > 5 %)"), ("drought_kills", "Drought kills (cells, year)")]
    fig, ax = plt.subplots(len(rows), 3, figsize=(15, 13), sharey="row")
    for j, (p, lab, d1, d2) in enumerate(cols):
        for cl, cn, c in CLIM:
            for d, ls, suf in ((d1, "-", ""), (d2, "--", ", osmolyte production off")):
                if d is None:
                    continue
                s = d.xs((p, cl), level=["param", "climate"]) if p in d.index.get_level_values(0) else None
                if s is None:
                    continue
                for i, (var, yl) in enumerate(rows):
                    ax[i, j].plot(s.index, s[var], color=c, ls=ls, marker="o", ms=5, lw=1.6, label=cn + suf)
        for i, (var, yl) in enumerate(rows):
            ax[i, j].grid(alpha=0.3, lw=0.6)
            ax[i, j].spines[["top", "right"]].set_visible(False)
            if j == 0:
                ax[i, j].set_ylabel(yl)
            if j > 0:
                ax[i, j].set_xscale("log")
        ax[0, j].set_title(lab)
        ax[-1, j].set_xlabel(lab)
    ax[0, 0].legend(fontsize=7.5, frameon=False)
    fig.suptitle("Osmolyte sweeps from the current defaults (bacterial minimal model, first year, MMRT, R_endo on; mean of "
                 "2 seeds). Baselines: beta 1, Osmo_prod_E_cost 600, Osmo_Induci_Prod 0.05", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(out, dpi=130)


if __name__ == "__main__":
    main(*sys.argv[1:5])
