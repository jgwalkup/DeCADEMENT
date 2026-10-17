"""Beta sweep at two osmolyte energy costs against osmolyte production off.

    python analysis/plot_osmo_cost_beta.py OUT_PNG COST600_DIR COST150_DIR OFF_DIR

Each directory is a param_sweep.py run varying beta_bac+beta_fun. Rows: % degraded at day 365, day 95 % degraded,
energy-limited days while > 5 % substrate is left, minimum E:C then, drought kills per year; columns: climates; lines:
Osmo_prod_E_cost 600, Osmo_prod_E_cost 150, osmolyte production off. Writes OUT_PNG and OUT_PNG with .csv.
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from plot_osmo_sweeps import CLIM, summ  # noqa: E402

P = "beta_bac+beta_fun"
SETS = [("Osmo_prod_E_cost 600", "#7f7f7f", "-"), ("Osmo_prod_E_cost 150", "#2a9d8f", "-"),
        ("Osmolyte production off", "#555555", "--")]
ROWS = [("degraded_365", "% degraded, day 365"), ("day_95pct", "Day 95 % degraded"),
        ("limited_days", "Energy-limited days\n(substrate > 5 %)"), ("min_E_to_C", "Minimum E:C (substrate > 5 %)"),
        ("drought_kills", "Drought kills (cells, year)")]


def main(out, d600, d150, doff):
    data = [summ(d).xs(P, level="param") for d in (d600, d150, doff)]
    pd.concat([d.assign(case=s[0]) for d, s in zip(data, SETS)]).round(3).to_csv(out.rsplit(".", 1)[0] + ".csv")
    fig, ax = plt.subplots(len(ROWS), 3, figsize=(14, 15), sharey="row", sharex=True)
    for j, (cl, cn, _) in enumerate(CLIM):
        for d, (lab, c, ls) in zip(data, SETS):
            s = d.xs(cl, level="climate")
            for i, (var, _) in enumerate(ROWS):
                ax[i, j].plot(s.index, s[var], color=c, ls=ls, marker="o", ms=5, lw=1.7, label=lab)
        ax[0, j].set_title(cn)
        ax[-1, j].set_xlabel("beta (drought-mortality slope)")
        for i, (_, yl) in enumerate(ROWS):
            ax[i, j].set_xscale("log")
            ax[i, j].set_xticks([1, 2, 3, 5, 10], ["1", "2", "3", "5", "10"])
            ax[i, j].minorticks_off()
            ax[i, j].grid(alpha=0.3, lw=0.6)
            ax[i, j].spines[["top", "right"]].set_visible(False)
            if j == 0:
                ax[i, j].set_ylabel(yl)
    ax[0, 0].legend(fontsize=8, frameon=False)
    fig.suptitle("Beta sweep at Osmo_prod_E_cost 600 and 150 vs osmolyte production off (bacterial minimal model, first "
                 "year, MMRT, R_endo on; mean of 2 seeds)", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(out, dpi=130)


if __name__ == "__main__":
    main(*sys.argv[1:5])
