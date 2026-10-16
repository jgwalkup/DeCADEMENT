"""
Moisture test of the minimal model under the ambient climate: which moisture effects limit decomposition.

Usage:
    python plot_moisture_test.py OUT_DIR  (reads the run folders listed in RUNS, relative to OUT_DIR)
Writes OUT_DIR/moisture_test.png and OUT_DIR/moisture_test.csv (day-365 seed means and SDs).
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

plt.rcParams.update({"font.size": 9, "axes.titlesize": 10, "axes.labelsize": 9, "legend.fontsize": 8, "figure.dpi": 150,
                     "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True, "grid.color": "#e6e5e0",
                     "grid.linewidth": 0.6, "axes.edgecolor": "#8a8981", "xtick.color": "#5f5e57", "ytick.color": "#5f5e57"})
# label, folder (relative to OUT_DIR), climate name, colour, line style
RUNS = [("Real climate (Psi down to -6.5 MPa)", "real", "ambient_temp", "#2a78d6", "-"),
        ("Real climate, no moisture effect on enzymes", "enz0", "ambient_temp", "#eb6834", "-"),
        ("Real climate, no moisture effect on uptake", "upt0", "ambient_temp", "#1baf7a", (0, (3, 2))),
        ("Real climate, no moisture effect on enzymes or uptake", "both0", "ambient_temp", "#eda100", "-"),
        ("... and no drought mortality (tolerance ~1)", "both0_tol1", "ambient_temp", "#e87ba4", "-"),
        ("Psi held at -0.5 MPa (real temperature)", "wetpsi", "ambient_wetpsi", "#008300", (0, (4, 2)))]


def main():
    out = sys.argv[1]
    rows, series = [], {}
    for lab, d, clim, col, ls in RUNS:
        s = pd.read_csv(os.path.join(out, d, "summary.csv"))
        t = pd.read_csv(os.path.join(out, d, "timeseries.csv.gz"))
        s, t = s[s["climate"] == clim], t[t["climate"] == clim]
        g = t.groupby("day").mean(numeric_only=True)
        series[lab] = (g, col, ls)
        rows.append(dict(run=lab, degraded_d365=s["deg_pct"].mean(), degraded_sd=s["deg_pct"].std(),
                         biomass_C_max=g["biomass_C"].max(), biomass_C_d365=s["biomC"].mean(),
                         cells_killed=t.groupby("seed")["killed"].sum().mean(), min_cells_alive=g["alive_cells"].min()))
    pd.DataFrame(rows).to_csv(os.path.join(out, "moisture_test.csv"), index=False)
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.2), sharex=True)
    for ax, (v, lab) in zip(axes, (("degraded", "Substrate degraded (%)"), ("biomass_C", "Biomass C (per grid cell)"),
                                   ("alive_cells", "Cells alive (of 400)"))):
        for name, (g, col, ls) in series.items():
            y = 100 - g["substrate_remaining_pct"] if v == "degraded" else g[v]
            ax.plot(y.index, y.values, color=col, ls=ls, lw=1.6, label=name)
        ax.set_ylabel(lab)
        ax.set_xlabel("Day")
        ax.axvspan(200, 295, color="#f0efec", zorder=0)
    axes[0].text(247, axes[0].get_ylim()[1] * 0.97, "dry season\n(Psi -4 to -6.5)", ha="center", va="top", fontsize=7,
                 color="#5f5e57")
    h, l = axes[0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncol=3, frameon=False, bbox_to_anchor=(0.5, -0.12))
    fig.suptitle("Bacterial minimal model, ambient climate, first year: which moisture effects limit decomposition "
                 "(R$_{\\mathrm{endo}}$ on; mean of 3 seeds)")
    fig.tight_layout()
    fig.savefig(os.path.join(out, "moisture_test.png"), bbox_inches="tight")
    plt.close(fig)
    print(pd.DataFrame(rows).round(2).to_string(index=False))


if __name__ == "__main__":
    main()
