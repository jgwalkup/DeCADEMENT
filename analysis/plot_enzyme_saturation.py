"""Forward Michaelis-Menten decay (Decay_MM_Form 0) vs decay with enzyme saturation (Decay_MM_Form 1).

    python analysis/plot_enzyme_saturation.py OUT_PNG FORM0_DIR FORM1_DIR

FORM0_DIR / FORM1_DIR are param_sweep outputs over the three climates (FORM1_DIR varying Enz_Sat_K). Rows: substrate
remaining, biomass C, enzyme pool, enzyme-saturation factor K_E S / (K_E S + E); columns: climates.
"""
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

CLIM = [("ambient_temp", "Ambient"), ("moderate_warm", "+5 °C"), ("extreme_warm", "+10 °C")]


def main(out, d0, d1):
    t0 = pd.read_csv(f"{d0}/timeseries.csv.gz").assign(run="Forward MM (form 0), defaults")
    t1 = pd.read_csv(f"{d1}/timeseries.csv.gz")
    base = pd.read_csv(f"{d1}/baseline_values.csv", header=None, index_col=0)[1]["Enz_Sat_K"]
    t1.loc[t1.param == "baseline", "value"] = base
    t1["run"] = "Enzyme saturation (form 1), Enz_Sat_K " + t1["value"].map(lambda v: f"{v:g}")
    t = pd.concat([t0, t1])
    runs = list(dict.fromkeys(t.run))
    cols = ["0.15", "#1f77b4", "#d62728", "#2ca02c", "#9467bd"]
    rows = [("substrate_remaining_pct", "Substrate remaining (%)"), ("biomass_C", "Biomass C per grid cell"),
            ("enzyme_pool", "Enzyme pool per grid cell"), ("d_enz_sat_factor", "Enzyme-saturation factor\nK_E S / (K_E S + E)")]
    fig, ax = plt.subplots(len(rows), 3, figsize=(13.5, 11.5), sharex=True, sharey="row")
    dry = t[t.Psi < -4].doy
    for j, (cl, cn) in enumerate(CLIM):
        for r, c in zip(runs, cols):
            g = t[(t.run == r) & (t.climate == cl)].groupby("doy").mean(numeric_only=True)
            for i, (var, lab) in enumerate(rows):
                ax[i, j].plot(g.index, g[var], color=c, lw=2.4 if r.startswith("Forward") else 1.7,
                              ls="--" if r.startswith("Forward") else "-", label=r)
        ax[0, j].set_title(cn)
        for i, (var, lab) in enumerate(rows):
            if len(dry):
                ax[i, j].axvspan(dry.min(), dry.max(), color="0.92", zorder=0, lw=0)
            ax[i, j].grid(alpha=0.3, lw=0.6)
            ax[i, j].spines[["top", "right"]].set_visible(False)
            if j == 0:
                ax[i, j].set_ylabel(lab)
        ax[-1, j].set_xlabel("Day of year")
    ax[0, 0].axhline(20, color="0.6", lw=0.8, ls=":")
    ax[0, 0].legend(fontsize=8, frameon=False, loc="lower left")
    fig.suptitle("Decay rate law: bacterial minimal model, first year (form 1 with Constit_Prod 0.001; Enz_Prod 0.005, "
                 "Vmax_Km 2.5, Uptake_Vmax_Km 0.25, R_endo on; mean of 2 seeds; shaded: Psi below −4 MPa)", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(out, dpi=140)


if __name__ == "__main__":
    main(*sys.argv[1:4])
