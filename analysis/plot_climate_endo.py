"""
Figure and table for climate_endo_run.py output: daily time courses under each climate scenario, endogenous respiration
(R_endo) off (wide pale line) and on (thin dashed line), mean of seeds.

Usage:
    python plot_climate_endo.py ../analysis_output/climate_endo [--note TEXT]
Writes <dir>/climate_endo.png and <dir>/climate_endo_summary.csv (day-365 seed means and SDs).
"""
import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

plt.rcParams.update({"font.size": 9, "axes.titlesize": 10, "axes.labelsize": 9, "legend.fontsize": 8, "figure.dpi": 150,
                     "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True, "grid.color": "#e6e5e0",
                     "grid.linewidth": 0.6, "axes.edgecolor": "#8a8981", "xtick.color": "#5f5e57", "ytick.color": "#5f5e57"})
ENDO = r"R$_{\mathrm{endo}}$"
ORDER = ["ambient_temp", "moderate_warm", "extreme_warm"]
NAMES = {"ambient_temp": "Ambient", "moderate_warm": "Moderate warming (+5 °C)", "extreme_warm": "Extreme warming (+10 °C)"}
COLORS = {"ambient_temp": "#2a78d6", "moderate_warm": "#eb6834", "extreme_warm": "#1baf7a"}
PANELS = [("Temp", "Temperature (°C)"), ("Psi", "Water potential (MPa)"),
          ("degraded", "Substrate degraded (%)"), ("biomass_C", "Biomass C (per grid cell)"),
          ("alive_cells", "Cells alive (of 400)"), ("CUE", "CUE"),
          ("energy_to_C", "Energy:C of living biomass"), ("Drought_tol", "Drought tolerance (mean of living cells)")]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("dir")
    ap.add_argument("--note", default="")
    a = ap.parse_args()
    ts = pd.read_csv(os.path.join(a.dir, "timeseries.csv.gz"))
    ts["degraded"] = 100 - ts["substrate_remaining_pct"]
    ts.loc[ts["alive_cells"] == 0, ["CUE", "energy_to_C", "Drought_tol"]] = np.nan
    m = ts.groupby(["climate", "endo", "day"]).mean(numeric_only=True)
    clims = [c for c in ORDER if c in ts["climate"].unique()] + sorted(set(ts["climate"]) - set(ORDER))

    fig, axes = plt.subplots(4, 2, figsize=(10, 10.5), sharex=True)
    for ax, (v, lab) in zip(axes.ravel(), PANELS):
        for c in clims:
            col = COLORS.get(c, "#5f5e57")
            for e, ls, lw, alpha in ((0, "-", 4.0, 0.35), (1, (0, (3, 2)), 1.3, 1.0)):
                if (c, e) not in m.index.droplevel(2):
                    continue
                y = m.loc[(c, e)][v]
                if v == "Psi" and c != clims[0]:
                    continue  # same moisture record in every scenario
                ax.plot(y.index, y.values, color=col if v != "Psi" else "#5f5e57", ls=ls, lw=lw, alpha=alpha,
                        solid_capstyle="butt", label=f"{NAMES.get(c, c)}, {ENDO} {'off' if e == 0 else 'on'}")
        ax.set_ylabel(lab)
        if v == "Psi":
            ax.axhline(-1.0, color="#8a8981", lw=0.8, ls=":")
            ax.text(5, -1.0, " wp_fc", va="bottom", fontsize=7, color="#5f5e57")
        if v in ("CUE", "Drought_tol"):
            ax.set_ylim(bottom=0)
    for ax in axes[-1]:
        ax.set_xlabel("Day")
    h, l = axes[0][0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncol=3, frameon=False, bbox_to_anchor=(0.5, -0.045))
    fig.suptitle("Minimal model under the climate scenarios (first year), endogenous respiration "
                 f"({ENDO}) off (wide pale) vs on (thin dashed); mean of seeds" + (f"\n{a.note}" if a.note else ""))
    fig.tight_layout()
    fig.savefig(os.path.join(a.dir, "climate_endo.png"), bbox_inches="tight")
    plt.close(fig)

    s = pd.read_csv(os.path.join(a.dir, "summary.csv"))
    cols = ["deg_pct", "biomC", "CUE", "alive", "EC_end", "E_unpaid_cum", "erosion_C_cum", "deficit_days", "HSP_prop"]
    out = s.groupby(["climate", "endo"])[cols].agg(["mean", "std"])
    out.columns = [f"{x}_{y}" for x, y in out.columns]
    out.to_csv(os.path.join(a.dir, "climate_endo_summary.csv"))
    print(s.groupby(["climate", "endo"])[cols].mean().round(3).to_string())


if __name__ == "__main__":
    main()
