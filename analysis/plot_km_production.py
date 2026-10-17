"""Enzyme Km slope x enzyme production options (forward Michaelis-Menten decay).

    python analysis/plot_km_production.py OUT_PNG GRID_DIR CLIM_DIR

GRID_DIR holds param_sweep outputs c<Constit>_e<EnzProd>/ (ambient, varying Vmax_Km); CLIM_DIR holds
c<Constit>_e<EnzProd>_km<Vmax_Km>/ runs over the three climates. Top: day-365 % degraded at ambient over the grid
(production as a multiple of Constit_Prod 0.00025 / Enz_Prod 0.005); below: substrate remaining, biomass C and enzyme
saturation for the climate runs. Writes OUT_PNG and OUT_PNG with .csv (the grid table).
"""
import glob
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

CLIM = [("ambient_temp", "Ambient"), ("moderate_warm", "+5 °C"), ("extreme_warm", "+10 °C")]


def main(out, gdir, cdir):
    rows = []
    for d in sorted(glob.glob(os.path.join(gdir, "c*"))):
        c = float(os.path.basename(d)[1:].split("_e")[0])
        t = pd.read_csv(os.path.join(d, "timeseries.csv.gz"))
        t.loc[t.param == "baseline", "value"] = 2.5
        for (v, s), g in t.groupby(["value", "seed"]):
            g = g.sort_values("doy")
            dg = 100 - g.substrate_remaining_pct.to_numpy()
            rows.append(dict(production=round(c / 0.00025, 2), Vmax_Km=v, seed=s, degraded_365=dg[-1],
                             last_quarter_share=(dg[-1] - dg[272]) / max(dg[-1], 1e-9)))
    grid = pd.DataFrame(rows).groupby(["production", "Vmax_Km"]).mean().drop(columns="seed")
    grid.round(3).to_csv(os.path.splitext(out)[0] + ".csv")

    fig = plt.figure(figsize=(14, 15))
    gs = fig.add_gridspec(4, 3, height_ratios=[1.15, 1, 1, 1])
    ax0 = fig.add_subplot(gs[0, :2])
    prods = sorted(grid.index.get_level_values(0).unique())
    cols = plt.cm.viridis(np.linspace(0, 0.85, len(prods)))
    for p, c in zip(prods, cols):
        d = grid.loc[p]
        ax0.plot(d.index, d["degraded_365"], marker="o", ms=6, color=c, lw=1.8, label=f"production x{p:g}")
    ax0.axhspan(75, 85, color="0.92", zorder=0, lw=0)
    ax0.set_xlabel("Vmax_Km (enzyme Km slope)")
    ax0.set_ylabel("Ambient: degraded, day 365 (%)")
    ax0.legend(frameon=False, fontsize=8, title="Constit_Prod / Enz_Prod\nx 0.00025 / 0.005", title_fontsize=8)
    ax0.grid(alpha=0.3, lw=0.6)
    ax0.spines[["top", "right"]].set_visible(False)
    ax0.text(0.99, 0.03, "shaded: 75-85 %", transform=ax0.transAxes, ha="right", fontsize=8, color="0.4")

    runs = sorted(glob.glob(os.path.join(cdir, "c*")),
                  key=lambda d: float(os.path.basename(d).split("_km")[1]))
    rc = plt.cm.plasma(np.linspace(0, 0.85, len(runs)))
    rows = [("substrate_remaining_pct", "Substrate remaining (%)"), ("biomass_C", "Biomass C per grid cell"),
            ("d_deg_saturation", "Enzyme saturation S / (Km + S)")]
    axes = [[fig.add_subplot(gs[i + 1, j]) for j in range(3)] for i in range(3)]
    dry = None
    for d, c in zip(runs, rc):
        b = os.path.basename(d)
        cp = float(b[1:].split("_e")[0]); km = float(b.split("_km")[1])
        lab = f"x{cp / 0.00025:g} production, Vmax_Km {km:g}"
        t = pd.read_csv(os.path.join(d, "timeseries.csv.gz"))
        dry = t[t.Psi < -4].doy
        for j, (cl, cn) in enumerate(CLIM):
            g = t[t.climate == cl].groupby("doy").mean(numeric_only=True)
            for i, (var, yl) in enumerate(rows):
                axes[i][j].plot(g.index, g[var], color=c, lw=1.7, label=lab)
    for i, (var, yl) in enumerate(rows):
        for j, (cl, cn) in enumerate(CLIM):
            a = axes[i][j]
            if dry is not None and len(dry):
                a.axvspan(dry.min(), dry.max(), color="0.92", zorder=0, lw=0)
            a.grid(alpha=0.3, lw=0.6)
            a.spines[["top", "right"]].set_visible(False)
            if i == 0:
                a.set_title(cn)
            if j == 0:
                a.set_ylabel(yl)
            if i == 2:
                a.set_xlabel("Day of year")
        ys = [a.get_ylim() for a in axes[i]]
        for a in axes[i]:
            a.set_ylim(min(y[0] for y in ys), max(y[1] for y in ys))
    axes[0][0].axhline(20, color="0.6", lw=0.8, ls=":")
    axes[0][0].legend(fontsize=7.5, frameon=False, loc="lower left")
    fig.suptitle("Enzyme Km slope and enzyme production (forward M-M, Decay_MM_Form 0): bacterial minimal model, first year,\n"
                 "Uptake_Vmax_Km 0.25, slopes 0.1, beta 1, R_endo on; mean of 2 seeds; shaded: Psi below −4 MPa", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(out, dpi=140)
    print(grid.round(2).to_string())


if __name__ == "__main__":
    main(*sys.argv[1:4])
