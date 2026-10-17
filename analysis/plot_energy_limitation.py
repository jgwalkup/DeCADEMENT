"""Energy limitation of a parameter set under the scenario climates (param_sweep output).

    python analysis/plot_energy_limitation.py OUT_PNG RUN_DIR [VALUE]

Plots the run of RUN_DIR whose varied value is VALUE (default: the baseline): substrate remaining, biomass C, E:C, and
the energy-limited days (shortfall paid by endogenous respiration, or production capped) while substrate remains.
"""
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

CLIM = [("ambient_temp", "Ambient", "#3b6fb6"), ("moderate_warm", "+5 °C", "#e0a030"), ("extreme_warm", "+10 °C", "#c0392b")]


def main(out, d, value=None):
    t = pd.read_csv(f"{d}/timeseries.csv.gz")
    t = t[t.param == "baseline"] if value is None else t[np.isclose(t.value, float(value)) & (t.param != "baseline")]
    t = t.assign(limited=((t.R_endo_C > 0) | (t.d_r_consti < 0.999) | (t.d_r_induci < 0.999))
                 & (t.substrate_remaining_pct > 5))
    rows = [("substrate_remaining_pct", "Substrate remaining (%)"), ("biomass_C", "Biomass C per grid cell"),
            ("energy_to_C", "E:C (energy per unit biomass C)"), ("limited", "Share of seeds energy-limited\n(substrate > 5 %)")]
    fig, ax = plt.subplots(len(rows), 1, figsize=(10, 11), sharex=True)
    dry = t[t.Psi < -4].doy
    for cl, cn, c in CLIM:
        g = t[t.climate == cl].groupby("doy").mean(numeric_only=True)
        if not len(g):
            continue
        for a, (var, lab) in zip(ax, rows):
            a.plot(g.index, g[var].astype(float), color=c, lw=1.7, label=cn)
    for a, (var, lab) in zip(ax, rows):
        if len(dry):
            a.axvspan(dry.min(), dry.max(), color="0.92", zorder=0, lw=0)
        a.set_ylabel(lab)
        a.grid(alpha=0.3, lw=0.6)
        a.spines[["top", "right"]].set_visible(False)
    ax[0].legend(frameon=False)
    ax[-1].set_xlabel("Day of year")
    fig.suptitle(f"Energy limitation: {d.split('/')[-1]}" + (f", value {value}" if value else " (baseline)") +
                 "\nshaded: Psi below −4 MPa; mean of seeds", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(out, dpi=140)


if __name__ == "__main__":
    main(*sys.argv[1:4])
