"""MMRT vs Arrhenius under the scenario climates: daily trajectories, one column per climate.

    python analysis/plot_kinetics_climate.py OUT_PNG MMRT_DIR ARRHENIUS_DIR [SUBTITLE]

Each directory is a param_sweep.py (or climate_endo_run.py) run; only the baseline rows (param == "baseline", or all rows
when there is no param column) with R_endo 1 are used. Rows: daily temperature, substrate remaining, biomass C, energy:C
of living biomass, daily CUE (7-day mean, shown while > 1 % substrate is left) and endogenous respiration C (7-day
mean). Mean of seeds. Writes OUT_PNG and OUT_PNG with .csv (daily means).
"""
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402

CLIM = [("ambient_temp", "Ambient"), ("moderate_warm", "+5 °C"), ("extreme_warm", "+10 °C")]
KIN = [("MMRT", "#1f77b4", "-"), ("Arrhenius", "#d62728", "--")]
COLS = ["Temp", "substrate_remaining_pct", "biomass_C", "energy_to_C", "CUE", "R_endo_C"]


def daily(d):
    t = pd.read_csv(f"{d}/timeseries.csv.gz")
    if "param" in t.columns:
        t = t[t.param == "baseline"]
    if "R_endo" in t.columns:
        t = t[t.R_endo == 1]
    if "pulse" in t.columns:
        t = t[t.pulse == 1]
    return t.groupby(["climate", "doy"])[COLS].mean().reset_index()


def main(out, d_mmrt, d_arr, subtitle=""):
    data = {"MMRT": daily(d_mmrt), "Arrhenius": daily(d_arr)}
    pd.concat([v.assign(kinetics=k) for k, v in data.items()]).round(4).to_csv(out.rsplit(".", 1)[0] + ".csv", index=False)
    rows = [("Temp", "Temperature (°C)"), ("substrate_remaining_pct", "Substrate remaining (%)"),
            ("biomass_C", "Biomass C per grid cell"), ("energy_to_C", "Energy:C of living biomass"),
            ("CUE", "Daily CUE (7-day mean)"), ("R_endo_C", "Endogenous respiration C\n(per grid cell per day, 7-day mean)")]
    fig, ax = plt.subplots(len(rows), 3, figsize=(15, 16), sharex=True, sharey="row")
    for j, (cl, cn) in enumerate(CLIM):
        for k, col, ls in KIN:
            g = data[k][data[k].climate == cl].sort_values("doy")
            for i, (var, _) in enumerate(rows):
                if var == "Temp" and k == "Arrhenius":
                    continue
                y = g[var]
                if var in ("CUE", "R_endo_C"):
                    y = y.rolling(7, center=True, min_periods=1).mean()
                if var == "CUE":
                    y = y.where(g.substrate_remaining_pct > 1)
                c, l = ("#555555", "-") if var == "Temp" else (col, ls)
                ax[i, j].plot(g.doy, y, color=c, ls=l, lw=1.2 if var == "Temp" else 1.8, label=k)
        ax[0, j].set_title(cn)
        ax[-1, j].set_xlabel("Day of year")
        for i, (_, yl) in enumerate(rows):
            ax[i, j].grid(alpha=0.3, lw=0.6)
            ax[i, j].spines[["top", "right"]].set_visible(False)
            if j == 0:
                ax[i, j].set_ylabel(yl)
    ax[1, 0].legend(fontsize=9, frameon=False)
    ax[0, 0].set_xlim(0, 365)
    fig.suptitle("MMRT vs Arrhenius under the scenario climates (bacterial minimal model, first year, R_endo on; mean of seeds)"
                 + (f"\n{subtitle}" if subtitle else ""), fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(out, dpi=120)


if __name__ == "__main__":
    main(*sys.argv[1:5])
