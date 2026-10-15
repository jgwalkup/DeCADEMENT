"""
Figure for necromass_check.py output: first-order decay rates by substrate, the necromass pool and its decay over the year, and its monomer pool.
Usage: python plot_necromass.py ../analysis_output/necromass_check/substrate_decay.csv
Writes necromass_decay.png next to the csv.
"""
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

plt.rcParams.update({"font.size": 9.5, "axes.titlesize": 10, "axes.labelsize": 9.5, "legend.fontsize": 8.5, "figure.dpi": 150})
d = pd.read_csv(sys.argv[1])
p1 = d[d.pulse == 1]
order = ["DeadMic", "DeadEnz", "OrgP2", "Chitin", "Protein1", "Protein2", "OrgP1", "Protein3", "Starch", "Hemicellulose", "Cellulose", "Lignin"]
per = {"days 1-100 (wet)": (1, 100), "days 150-250 (dry)": (150, 250), "days 300-365 (rewetted)": (300, 365)}
cols = ["#9ecae1", "#fd8d3c", "#4292c6"]
fig, axes = plt.subplots(1, 3, figsize=(17, 4.8))
ax = axes[0]
x = np.arange(len(order))
for k, ((lab, (lo, hi)), c) in enumerate(zip(per.items(), cols)):
    g = p1[(p1.doy >= lo) & (p1.doy <= hi)].groupby(["substrate", "seed"]).agg(m=("mass", "mean"), dc=("decay", "mean")).groupby("substrate").mean()
    ax.bar(x + (k - 1) * 0.27, (g.dc / g.m).loc[order], 0.27, color=c, label=lab)
ax.set_xticks(x)
ax.set_xticklabels(order, rotation=55, ha="right")
ax.set_ylabel("First-order decay rate (decay / mass, per day)")
ax.set_title("Necromass decays as fast as the labile substrates")
ax.legend(frameon=False)
ax = axes[1]
dm = p1[p1.substrate == "DeadMic"].groupby("doy")[["mass", "decay", "psi"]].mean()
ax.plot(dm.index, dm["mass"], color="#8c564b", label="DeadMic pool (C+N+P per cell)")
ax.set_xlabel("Day of pulse")
ax.set_ylabel("Pool per cell", color="#8c564b")
ax2 = ax.twinx()
ax2.plot(dm.index, dm["decay"].rolling(7, center=True, min_periods=1).mean(), color="k", lw=1.2, label="decay per day")
ax2.set_ylabel("Decay per cell per day")
ax3 = ax.twinx()
ax3.spines["right"].set_position(("axes", 1.18))
ax3.plot(dm.index, dm["psi"], color="#2166ac", lw=0.9, alpha=0.8)
ax3.set_ylabel("Ψ (MPa)", color="#2166ac")
ax.set_title("DeadMic fills in the dry season and is discarded at the reset")
ax = axes[2]
for s, c in (("DeadMic", "#8c564b"), ("DeadEnz", "#bcbd22"), ("Starch", "#17becf"), ("Hemicellulose", "#1f77b4")):
    m = p1[p1.substrate == s].groupby("doy")["monomer_after_uptake"].mean()
    ax.plot(m.index, m.rolling(7, center=True, min_periods=1).mean(), color=c, label=f"Mon_{s}")
ax.set_yscale("log")
ax.set_xlabel("Day of pulse")
ax.set_ylabel("Monomer pool left after uptake (per cell)")
ax.set_title("Monomers are taken up; the pools stay small")
ax.legend(frameon=False)
for a in (axes[0], axes[1], axes[2]):
    a.spines[["top"]].set_visible(False)
fig.tight_layout()
out = os.path.join(os.path.dirname(os.path.abspath(sys.argv[1])), "necromass_decay.png")
fig.savefig(out, dpi=150, bbox_inches="tight")
print("wrote", out)
