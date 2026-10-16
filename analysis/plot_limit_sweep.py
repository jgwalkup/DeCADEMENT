"""Figure and table for the degradation-limitation runs (analysis/param_sweep.py output).

    python analysis/plot_limit_sweep.py OUT_PNG RUN_DIR [RUN_DIR ...] [--label DIR=TEXT ...]

Each RUN_DIR is a param_sweep output; a run is labelled "<prefix><param> = <value>" (prefix from --label, for runs made
with extra --set changes). Rows: substrate degraded, biomass C, enzyme pool, uptake capacity / monomer pool (log; >= 1
means the community can take up the whole organic monomer pool each day), organic monomer pool; columns: climates.
Writes OUT_PNG and OUT_PNG with .csv (annual summaries, mean of seeds).
"""
import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

CLIM = [("ambient_temp", "Ambient"), ("moderate_warm", "+5 °C"), ("extreme_warm", "+10 °C")]
SHORT = {"Uptake_Vmax_Km": "Uptake Km slope", "Enz_Prod_min+Enz_Prod_max": "Enz_Prod", "Vmax_Km": "Enzyme Km slope"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("out_png")
    ap.add_argument("dirs", nargs="+")
    ap.add_argument("--label", nargs="*", default=[], metavar="DIR=TEXT")
    ap.add_argument("--keep", nargs="*", default=None, help="labels to plot (default all)")
    a = ap.parse_args()
    pref = dict(x.split("=", 1) for x in a.label)
    parts = []
    for d in a.dirs:
        t = pd.read_csv(os.path.join(d, "timeseries.csv.gz"))
        p = pref.get(d, "")
        t["run"] = np.where(t.param == "baseline", (p + "baseline") if p else "baseline",
                            p + t.param.map(lambda s: SHORT.get(s, s)) + " = " + t.value.map(lambda v: f"{v:g}"))
        parts.append(t)
    t = pd.concat(parts)
    t["degraded"] = 100 - t["substrate_remaining_pct"]

    def f(g):
        dry = g[g.Psi < -4]
        return pd.Series(dict(degraded_365=g.degraded.iloc[-1], biomass_mean=g.biomass_C.mean(),
                              dry_season_biomass_change=dry.biomass_C.iloc[-1] - dry.biomass_C.iloc[0],
                              enzyme_pool_mean=g.enzyme_pool.mean(), enzyme_saturation=g.d_deg_saturation.mean(),
                              uptake_capacity_ratio=g.d_upt_capacity_ratio.mean(), monomer_pool=g.monomer_pool.mean(),
                              min_E_to_C=g.energy_to_C.min(), starved=g.starved.sum()))
    summ = t.groupby(["run", "climate", "seed"]).apply(f, include_groups=False).groupby(["run", "climate"]).mean()
    summ.round(4).to_csv(os.path.splitext(a.out_png)[0] + ".csv")

    runs = a.keep or list(dict.fromkeys(t.run))
    cols = plt.cm.tab10(np.arange(len(runs)) % 10)
    rows = [("degraded", "Substrate degraded (%)", False), ("biomass_C", "Biomass C per grid cell", False),
            ("enzyme_pool", "Enzyme pool per grid cell", False),
            ("d_upt_capacity_ratio", "Uptake capacity / monomer pool\n(≥ 1: uptake not limiting)", True),
            ("monomer_pool", "Organic monomer pool per grid cell", False)]
    fig, ax = plt.subplots(len(rows), 3, figsize=(14, 14), sharex=True, sharey="row")
    dry = t[(t.climate == "extreme_warm") & (t.Psi < -4)].doy
    for j, (cl, cn) in enumerate(CLIM):
        for r, c in zip(runs, cols):
            g = t[(t.run == r) & (t.climate == cl)].groupby("doy").mean(numeric_only=True)
            isb = r == "baseline"
            for i, (var, lab, lg) in enumerate(rows):
                ax[i, j].plot(g.index, g[var], color="0.1" if isb else c, lw=2.6 if isb else 1.5, ls="--" if isb else "-",
                              label=r)
        ax[0, j].set_title(cn)
        for i, (var, lab, lg) in enumerate(rows):
            if lg:
                ax[i, j].set_yscale("log")
                ax[i, j].axhline(1, color="0.5", lw=0.8, ls=":")
            if len(dry):
                ax[i, j].axvspan(dry.min(), dry.max(), color="0.92", zorder=0, lw=0)
            ax[i, j].grid(alpha=0.3, lw=0.6)
            ax[i, j].spines[["top", "right"]].set_visible(False)
            if j == 0:
                ax[i, j].set_ylabel(lab)
        ax[-1, j].set_xlabel("Day of year")
    ax[0, 0].legend(fontsize=7.5, frameon=False, loc="upper left")
    fig.suptitle("What limits degradation: bacterial minimal model, first year (slopes 0.1, beta 1, R_endo on; "
                 "mean of 2 seeds; shaded: Psi below −4 MPa)", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(a.out_png, dpi=140)
    print(summ.round(3).to_string())


if __name__ == "__main__":
    main()
