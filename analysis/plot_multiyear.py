"""Multi-year runs of the minimal model under MMRT and Arrhenius kinetics.

    python analysis/plot_multiyear.py OUT_DIR SIX_MMRT SIX_ARRH CONT_MMRT CONT_ARRH

SIX_*: climate_endo_run.py output with --pulses 6 (litter reset and community re-placed each year); CONT_*: one
continuous pulse over several years (no reset). Writes OUT_DIR/multiyear_six_pulse.png, multiyear_continuous.png and
multiyear_pulses.csv (per-pulse day-365 degradation, day of 50 %, start and maximum biomass; mean of seeds).
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

CLIM = [("ambient_temp", "Ambient"), ("moderate_warm", "+5 °C"), ("extreme_warm", "+10 °C")]
KIN = [("MMRT", "#1f77b4", "-"), ("Arrhenius", "#d62728", "--")]


def load(d, kin):
    return pd.read_csv(os.path.join(d, "timeseries.csv.gz")).assign(kinetics=kin)


def figure(t, rows, out, title, pulse_lines):
    fig, ax = plt.subplots(len(rows), 3, figsize=(15, 3.0 * len(rows) + 1), sharex=True, sharey="row")
    for j, (cl, cn) in enumerate(CLIM):
        for kin, c, ls in KIN:
            g = t[(t.kinetics == kin) & (t.climate == cl)].groupby("day").mean(numeric_only=True)
            for i, (var, lab, win) in enumerate(rows):
                y = g[var].rolling(win, min_periods=1, center=True).mean() if win > 1 else g[var]
                ax[i, j].plot(g.index, y, color=c, ls=ls, lw=1.5, label=kin)
        ax[0, j].set_title(cn)
        for i, (var, lab, win) in enumerate(rows):
            for x in pulse_lines:
                ax[i, j].axvline(x, color="0.75", lw=0.7)
            ax[i, j].grid(alpha=0.3, lw=0.6)
            ax[i, j].spines[["top", "right"]].set_visible(False)
            if j == 0:
                ax[i, j].set_ylabel(lab)
        ax[-1, j].set_xlabel("Day of run")
    ax[0, 0].legend(frameon=False, fontsize=8)
    fig.suptitle(title, fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(out, dpi=140)
    plt.close(fig)


def main(out, s1, s0, c1, c0):
    six = pd.concat([load(s1, "MMRT"), load(s0, "Arrhenius")])
    cont = pd.concat([load(c1, "MMRT"), load(c0, "Arrhenius")])
    six["degraded"] = 100 - six.substrate_remaining_pct
    cont["degraded"] = 100 - cont.substrate_remaining_pct
    rows = [("substrate_remaining_pct", "Substrate remaining (%)", 1), ("biomass_C", "Biomass C per grid cell", 1),
            ("enzyme_pool", "Enzyme pool per grid cell", 1), ("CUE", "CUE (15-day mean)", 15)]
    base = ("bacterial minimal model: production x1.4 (Constit_Prod 0.00035, Enz_Prod 0.007), Vmax_Km 13, " + os.environ.get("MULTIYEAR_NOTE", "") +
            "Uptake_Vmax_Km 0.25, slopes 0.1, beta 1, R_endo on")
    figure(six, rows, os.path.join(out, "multiyear_six_pulse.png"),
           f"Six yearly pulses (litter reset, community re-placed each year; mean of 2 seeds)\n{base}",
           [365 * k for k in range(1, 6)])
    figure(cont, rows, os.path.join(out, "multiyear_continuous.png"),
           f"One continuous 3-year pulse (no reset; 1 seed)\n{base}", [365, 730])

    def per_pulse(g):
        g = g.sort_values("doy")
        dg = g.degraded.to_numpy(); d = g.doy.to_numpy()
        return pd.Series(dict(degraded_365=dg[-1], day_50pct=d[np.argmax(dg >= 50)] if (dg >= 50).any() else np.nan,
                              biomass_start=g.biomass_C.iloc[0], biomass_max=g.biomass_C.max(),
                              biomass_end=g.biomass_C.iloc[-1], min_E_to_C=g.energy_to_C.min()))
    pp = six.groupby(["kinetics", "climate", "pulse", "seed"]).apply(per_pulse, include_groups=False)
    pp = pp.groupby(["kinetics", "climate", "pulse"]).mean()
    pp.round(3).to_csv(os.path.join(out, "multiyear_pulses.csv"))
    yr = cont.assign(year=(cont.day - 1) // 365 + 1).groupby(["kinetics", "climate", "year"]).agg(
        remaining_end=("substrate_remaining_pct", "last"), biomass_end=("biomass_C", "last"),
        biomass_max=("biomass_C", "max"), min_E_to_C=("energy_to_C", "min"), starved=("starved", "sum"))
    yr.round(3).to_csv(os.path.join(out, "multiyear_continuous_years.csv"))
    print(pp.round(1).to_string()); print(yr.round(2).to_string())


if __name__ == "__main__":
    main(*sys.argv[1:6])
