"""
Why energy:C falls across pulses: reseeded starting biomass (fungal share) and the single-pulse start-biomass test.
Usage: python plot_energy_decline.py ../analysis_output/full_climate_km25_osmorange2_50x50_12p ../analysis_output/start_biomass_test
Needs fungal_share.csv (fungal_share.py) and timeseries.csv in the first folder, and timeseries.csv of start_biomass_test.py in the second.
Writes <first folder>/energy_decline.png.
"""
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

plt.rcParams.update({"font.size": 10, "axes.titlesize": 10, "axes.labelsize": 10, "legend.fontsize": 8.5, "figure.dpi": 150})
a, b = sys.argv[1], sys.argv[2]
ts = pd.read_csv(os.path.join(a, "timeseries.csv"))
fs = pd.read_csv(os.path.join(a, "fungal_share.csv"))
sb = pd.read_csv(os.path.join(b, "timeseries.csv"))
fig, axes = plt.subplots(1, 4, figsize=(19, 4.6))

ax = axes[0]
m = fs.groupby("pulse")[["init_biomass", "fungal_share_rows_day1"]].mean()
ax.plot(m.index, m["init_biomass"], "-o", ms=4, color="#1f77b4", label="starting biomass C per cell")
ax.set_xlabel("Pulse")
ax.set_ylabel("Starting biomass C per cell (day 1)", color="#1f77b4")
ax2 = ax.twinx()
ax2.plot(m.index, 100 * m["fungal_share_rows_day1"], "-s", ms=4, color="#d6172a", label="fungal share of cells")
ax2.set_ylabel("Fungal share of reseeded cells (%)", color="#d6172a")
ax.set_title("Reseeded start gets larger as fungi take over")
for x in (ax, ax2):
    x.spines[["top"]].set_visible(False)

ax = axes[1]
for s, g in ts[ts["doy"] == 365].groupby("seed"):
    ax.plot(g["pulse"], g["energy_to_C"], "-", lw=1, alpha=0.6, label=f"community {s}")
ax.plot(ts[ts["doy"] == 365].groupby("pulse")["energy_to_C"].mean(), "k-o", ms=4, lw=2, label="mean")
ax.set_xlabel("Pulse")
ax.set_ylabel("Energy : C on day 365")
ax.set_title("Energy:C at the end of each pulse")
ax.legend(frameon=False, ncols=2)
ax.spines[["top", "right"]].set_visible(False)

ax = axes[2]
r = fs.merge(ts[ts["doy"] == 365][["pulse", "seed", "energy_to_C"]], on=["pulse", "seed"])
sc = ax.scatter(r["init_biomass"], r["energy_to_C"], c=r["pulse"], cmap="viridis", s=22)
ax.set_xscale("log")
ax.set_xlabel("Starting biomass C per cell (log scale)")
ax.set_ylabel("Energy : C on day 365")
ax.set_title("12 pulses x 5 communities")
fig.colorbar(sc, ax=ax, label="pulse")
ax.spines[["top", "right"]].set_visible(False)

ax = axes[3]
order = ["bac1", "base", "x3", "bac25", "fun25", "x10"]
lab = {"bac1": "1 bacterial\nrow, C=1", "base": "pulse-1\nstart", "x3": "pulse-1\nstart x3", "bac25": "bacterial\nC=25", "fun25": "fungal\nC=25", "x10": "pulse-1\nstart x10"}
e = sb[sb["day"] == 365].groupby(["cond", "seed"])["energy_to_C"].mean().unstack("seed")
i0 = sb[sb["day"] == 1].groupby("cond")["biomass_C"].mean()
xs = np.arange(len(order))
for s in e.columns:
    ax.plot(xs, [e.loc[c, s] for c in order], ".", color="gray", alpha=0.6)
ax.bar(xs, [e.loc[c].mean() for c in order], color=["#9ecae1", "#6baed6", "#4292c6", "#fd8d3c", "#d6172a", "#08519c"], alpha=0.8)
ax.set_xticks(xs)
ax.set_xticklabels([f"{lab[c]}\n({i0[c]:.0f} C)" for c in order], fontsize=8)
ax.set_ylabel("Energy : C on day 365")
ax.set_title("Single pulse, only the start changed (20x20, 5 communities)")
ax.spines[["top", "right"]].set_visible(False)
fig.tight_layout()
out = os.path.join(a, "energy_decline.png")
fig.savefig(out, dpi=150, bbox_inches="tight")
print("wrote", out)
