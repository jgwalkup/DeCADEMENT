"""
Figures and annual N budgets from n_tracker.py output.

Usage: python plot_n_budget.py OUT_DIR LABEL1=PATH1/n_daily.csv LABEL2=PATH2/n_daily.csv ...
  e.g. python plot_n_budget.py ../analysis_output/n_tracking "N x1=../analysis_output/n_tracking/x1/n_daily.csv" "N x4=../analysis_output/n_tracking/x4/n_daily.csv"
Writes OUT_DIR/n_pools.png (stacked N pools over the year), OUT_DIR/n_fluxes.png (uptake by group, overflow C share, N fates) and prints
annual budgets for pulses 1, 3 and 6.
"""
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

plt.rcParams.update({"font.size": 9.5, "axes.titlesize": 10, "axes.labelsize": 9.5, "legend.fontsize": 8.5, "figure.dpi": 150})
out = sys.argv[1]
runs = {}
for arg in sys.argv[2:]:
    lab, path = arg.split("=", 1)
    runs[lab] = pd.read_csv(path)
labs = list(runs)
PULSES = [1, 3, 6]
POOLS = [("sub_N_litter", "litter substrates", "#8c6d31"), ("sub_N_dead", "dead microbes and enzymes (substrate)", "#bcbd22"),
         ("mon_N_org", "organic monomers", "#17becf"), ("mon_N_nh4", "NH4", "#9edae5"), ("enz_N", "enzymes", "#7f7f7f"),
         ("mic_N_bac", "bacteria", "#1f77b4"), ("mic_N_fun", "fungi", "#d6172a")]
BAC, FUN = "#1f77b4", "#d6172a"
pd.set_option("display.width", 220)


def sm(x, w=7):
    return x.rolling(w, center=True, min_periods=1).mean()


fig, axes = plt.subplots(len(labs), len(PULSES), figsize=(5.2 * len(PULSES), 3.6 * len(labs)), sharex=True, squeeze=False)
for i, lab in enumerate(labs):
    d = runs[lab]
    for j, p in enumerate(PULSES):
        m = d[d.pulse == p].groupby("doy").mean(numeric_only=True)
        ax = axes[i][j]
        ax.stackplot(m.index, *[m[c] for c, _, _ in POOLS], colors=[c for _, _, c in POOLS], labels=[l for _, l, _ in POOLS])
        ax.set_title(f"{lab}, pulse {p}")
        if j == 0:
            ax.set_ylabel("N per grid cell")
        if i == len(labs) - 1:
            ax.set_xlabel("Day of pulse")
        ax.spines[["top", "right"]].set_visible(False)
h, l = axes[0][0].get_legend_handles_labels()
fig.legend(h, l, loc="lower center", ncols=4, framealpha=0.9, bbox_to_anchor=(0.5, -0.02))
fig.suptitle("Where the N is: pools over the year (mean of communities; the total declines only by leaching)", fontsize=11, fontweight="bold")
fig.tight_layout(rect=(0, 0.04, 1, 0.97))
fig.savefig(os.path.join(out, "n_pools.png"), dpi=150, bbox_inches="tight")
plt.close(fig)

fig, axes = plt.subplots(len(labs), 4, figsize=(19, 3.6 * len(labs)), sharex=True, squeeze=False)
for i, lab in enumerate(labs):
    d = runs[lab]
    for j, p in enumerate((1, 6)):
        m = d[d.pulse == p].groupby("doy").mean(numeric_only=True)
        ax = axes[i][j]
        ax.plot(m.index, sm(m["uptake_N_bac"]), color=BAC, label="bacteria: N uptake")
        ax.plot(m.index, sm(m["uptake_N_fun"]), color=FUN, label="fungi: N uptake")
        ax.plot(m.index, sm(m["decay_N"]), color="k", lw=1, ls=":", label="N released by decay")
        ax.set_title(f"{lab}, pulse {p}: N uptake per day")
        if j == 0:
            ax.set_ylabel("N per grid cell per day")
        ax.spines[["top", "right"]].set_visible(False)
    ax = axes[i][2]
    for p, ls in ((1, "-"), (6, "--")):
        m = d[d.pulse == p].groupby("doy").mean(numeric_only=True)
        for g, c in (("bac", BAC), ("fun", FUN)):
            share = sm(m[f"overflow_C_{g}"]) / sm(m[f"uptake_C_{g}"]).replace(0, np.nan)
            ax.plot(m.index, 100 * share, color=c, ls=ls, label=f"{'bacteria' if g == 'bac' else 'fungi'}, pulse {p}")
    ax.set_title(f"{lab}: C released to respiration by cells with too little N, % of C uptake")
    ax.set_ylim(0, None)
    ax.spines[["top", "right"]].set_visible(False)
    ax = axes[i][3]
    m = d[d.pulse == 1].groupby("doy").mean(numeric_only=True)
    ax.plot(m.index, sm(m["death_N_bac"] + m["death_N_fun"]), color="#8c564b", label="N in cells killed (to dead pool)")
    ax.plot(m.index, sm(m["release_N_bac"] + m["release_N_fun"]), color="#2ca02c", label="N released to NH4 by living cells")
    ax.plot(m.index, sm(m["leach_N"]), color="#ff7f0e", label="N leached (lost)")
    ax.plot(m.index, sm(m["enz_prod_N"]), color="#7f7f7f", label="N into enzymes")
    ax.set_title(f"{lab}, pulse 1: other N fluxes")
    ax.spines[["top", "right"]].set_visible(False)
    if i == len(labs) - 1:
        for a in axes[i]:
            a.set_xlabel("Day of pulse")
for a in axes[0]:
    a.legend(frameon=False, fontsize=7.5)
fig.suptitle("N fluxes and the C lost when N is short", fontsize=11, fontweight="bold")
fig.tight_layout(rect=(0, 0, 1, 0.97))
fig.savefig(os.path.join(out, "n_fluxes.png"), dpi=150, bbox_inches="tight")
plt.close(fig)

for lab in labs:
    d = runs[lab]
    print(f"\n=== {lab}: max |daily N residual| {d.resid_N.abs().max():.1e}; N lost to leaching per year (pulse 1) {d[d.pulse == 1].leach_N.sum():.3f}")
    rows = {}
    for p in PULSES:
        x = d[d.pulse == p].groupby("seed").agg(
            decay=("decay_N", "sum"), up_bac=("uptake_N_bac", "sum"), up_fun=("uptake_N_fun", "sum"), enz=("enz_prod_N", "sum"),
            deadenz=("deadenz_N", "sum"), death_bac=("death_N_bac", "sum"), death_fun=("death_N_fun", "sum"),
            rel_bac=("release_N_bac", "sum"), rel_fun=("release_N_fun", "sum"), leach=("leach_N", "sum"),
            ovC_bac=("overflow_C_bac", "sum"), ovC_fun=("overflow_C_fun", "sum"), upC_bac=("uptake_C_bac", "sum"), upC_fun=("uptake_C_fun", "sum")).mean()
        x["overflow_%_of_C_uptake_bac"] = 100 * x.ovC_bac / x.upC_bac
        x["overflow_%_of_C_uptake_fun"] = 100 * x.ovC_fun / x.upC_fun
        end = d[(d.pulse == p) & (d.doy == 365)].mean(numeric_only=True)
        x["end_N_litter"], x["end_N_dead"], x["end_N_mic_bac"], x["end_N_mic_fun"] = end.sub_N_litter, end.sub_N_dead, end.mic_N_bac, end.mic_N_fun
        rows[f"pulse {p}"] = x
    print(pd.DataFrame(rows).round(3).to_string())
