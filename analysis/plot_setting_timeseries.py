"""Daily substrate remaining, biomass C and CUE under the scenario climates for a few parameter settings.

    python analysis/plot_setting_timeseries.py OUT_PNG TITLE LABEL=DIR[:PARAM=VALUE] [LABEL=DIR[:PARAM=VALUE] ...]

Each DIR is a param_sweep.py (or climate_endo_run.py) run; PARAM=VALUE picks the rows of that varied parameter at that
value (PARAM=baseline picks the baseline rows; no selector takes the baseline rows, or all rows when the run has no
param column). Only R_endo 1 and the first pulse are used. Lines in the order given use light-to-dark blue; a label
starting with "ref" is drawn grey and dashed. Rows: substrate remaining, biomass C per grid cell, daily CUE (7-day mean,
shown while > 1 % substrate is left); columns: climates. Mean of seeds. Writes OUT_PNG and OUT_PNG with .csv.
"""
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

CLIM = [("ambient_temp", "Ambient"), ("moderate_warm", "+5 °C"), ("extreme_warm", "+10 °C")]
BLUES = ["#86b6ef", "#3987e5", "#1c5cab", "#0d366b"]
ROWS = [("substrate_remaining_pct", "Substrate remaining (%)"), ("biomass_C", "Biomass C per grid cell"),
        ("CUE", "Daily CUE (7-day mean)")]


def load(spec):
    lab, rest = spec.split("=", 1)
    d, sel = (rest.split(":", 1) + [None])[:2]
    t = pd.read_csv(f"{d}/timeseries.csv.gz")
    if "param" in t.columns:
        if sel is None or sel.startswith("baseline"):
            t = t[t.param == "baseline"]
        else:
            p, v = sel.split("=")
            t = t[(t.param == p) & np.isclose(t.value, float(v))]
    if "R_endo" in t.columns:
        t = t[t.R_endo == 1]
    if "pulse" in t.columns:
        t = t[t.pulse == 1]
    return lab, t.groupby(["climate", "doy"])[[r for r, _ in ROWS]].mean().reset_index()


def main(out, title, specs):
    data = [load(s) for s in specs]
    pd.concat([d.assign(setting=l) for l, d in data]).round(4).to_csv(out.rsplit(".", 1)[0] + ".csv", index=False)
    fig, ax = plt.subplots(len(ROWS), 3, figsize=(14, 10), sharex=True, sharey="row")
    nb = 0
    for lab, d in data:
        ref = lab.lower().startswith("ref")
        col, ls, lw = ("#8a8a85", "--", 1.5) if ref else (BLUES[nb % len(BLUES)], "-", 1.9)
        nb += 0 if ref else 1
        for j, (cl, _) in enumerate(CLIM):
            g = d[d.climate == cl].sort_values("doy")
            for i, (var, _) in enumerate(ROWS):
                y = g[var]
                if var == "CUE":
                    y = y.rolling(7, center=True, min_periods=1).mean().where(g.substrate_remaining_pct > 1)
                ax[i, j].plot(g.doy, y, color=col, ls=ls, lw=lw, label=lab)
    for j, (_, cn) in enumerate(CLIM):
        ax[0, j].set_title(cn)
        ax[-1, j].set_xlabel("Day of year")
        for i, (_, yl) in enumerate(ROWS):
            ax[i, j].grid(alpha=0.3, lw=0.6)
            ax[i, j].spines[["top", "right"]].set_visible(False)
            if j == 0:
                ax[i, j].set_ylabel(yl)
    ax[0, 0].set_xlim(0, 365)
    ax[0, 0].legend(fontsize=8.5, frameon=False, loc="lower left")
    fig.suptitle(title, fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(out, dpi=120)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3:])
