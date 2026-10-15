"""
Does raising the N supply remove the fungal advantage? Compares the baseline with runs in which the N of every initial substrate is
multiplied by k (fungal_tracking.py --sub-n-mult k), same communities, 20x20, 6 pulses.

Usage: python plot_n_supply.py ../analysis_output/fungal_tracking [--baseline-ts ../analysis_output/full_climate_km25_osmorange2/timeseries.csv]
Looks for baseline/ and subN_x<k>/ folders. Writes n_supply.png and prints the numbers.
"""
import argparse
import glob
import os
import re

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

plt.rcParams.update({"font.size": 9.5, "axes.titlesize": 10, "axes.labelsize": 9.5, "legend.fontsize": 8.5, "figure.dpi": 150})
ap = argparse.ArgumentParser()
ap.add_argument("root")
ap.add_argument("--baseline-ts", default="../analysis_output/full_climate_km25_osmorange2/timeseries.csv")
a = ap.parse_args()
pd.set_option("display.width", 220)
BAC, FUN = "#1f77b4", "#d6172a"
NP = 6

runs = {1.0: a.root + "/baseline"}
for d in glob.glob(os.path.join(a.root, "subN_x*")):
    runs[float(re.search(r"subN_x([0-9.]+)", d).group(1))] = d
ks = sorted(runs)
cm = dict(zip(ks, plt.cm.viridis(np.linspace(0.0, 0.9, len(ks)))))
base_ts = pd.read_csv(a.baseline_ts)
base_ts = base_ts[base_ts["seed"] < 5]


def shares(t):
    out = []
    for (s, p), x in t.groupby(["seed", "pulse"]):
        f = x[x["fungal"] == 1]
        out.append(dict(seed=s, pulse=p, rows1=f["cells_day1"].sum() / max(x["cells_day1"].sum(), 1),
                        rows365=f["cells_day365"].sum() / max(x["cells_day365"].sum(), 1),
                        C1=f["C_day1"].sum() / max(x["C_day1"].sum(), 1e-9), C365=f["C_day365"].sum() / max(x["C_day365"].sum(), 1e-9),
                        tot1=x["C_day1"].sum(), tot365=x["C_day365"].sum(),
                        fun_taxa365=int(((x.fungal == 1) & (x.cells_day365 > 0)).sum()),
                        bac_taxa365=int(((x.fungal == 0) & (x.cells_day365 > 0)).sum())))
    return pd.DataFrame(out)


res = {}
for k in ks:
    g = pd.read_csv(runs[k] + "/groups.csv")
    t = pd.read_csv(runs[k] + "/taxa_track.csv")
    g, t = g[g.pulse <= NP], t[t.pulse <= NP]
    sh = shares(t).groupby("pulse").mean()
    nc = g[g.doy >= 100].groupby(["pulse", "group"])["N_to_C"].mean().unstack("group")
    pc = g[g.doy >= 100].groupby(["pulse", "group"])["P_to_C"].mean().unstack("group")
    gr = g[g.pulse == 1].groupby(["group", "doy"])["biomass_C"].mean().unstack("group")
    growth = gr.loc[365] / gr.loc[1]
    if k == 1.0:
        d = base_ts[base_ts.pulse <= NP]
        ec = d[d.doy == 365].groupby("pulse")["energy_to_C"].mean()
        deg = 100 - d[d.doy == 365].groupby("pulse")["substrate_remaining_pct"].mean()
    else:
        dd = pd.read_csv(runs[k] + "/daily.csv")
        ec = dd[dd.doy == 365].groupby("pulse")["energy_to_C"].mean()
        deg = 100 - dd[dd.doy == 365].groupby("pulse")["remaining_pct"].mean()
    res[k] = dict(sh=sh, nc=nc, pc=pc, growth=growth, ec=ec, deg=deg)

fig, axes = plt.subplots(2, 3, figsize=(16, 8.6))
ax = axes[0][0]
x = np.arange(len(ks))
nb = [res[k]["nc"]["bacteria"].mean() for k in ks]
nf = [res[k]["nc"]["fungi"].mean() for k in ks]
ax.bar(x - 0.18, nb, 0.34, color=BAC, label="bacteria")
ax.bar(x + 0.18, nf, 0.34, color=FUN, label="fungi")
ax.axhline(0.16 / 0.825, color=BAC, ls=":", lw=1)
ax.axhline(0.09 / 0.9, color=FUN, ls=":", lw=1)
ax.set_xticks(x)
ax.set_xticklabels([f"x{k:g}" for k in ks])
ax.set_xlabel("Multiplier on the N of the initial substrates")
ax.set_ylabel("N : C of the biomass (day 100-365, pulses 1-6)")
ax.set_title("N:C reached (dotted = target of each group)")
ax.legend(frameon=False)
panels = [("C1", "Fungal share of biomass C on day 1 (%)", axes[0][1]), ("rows365", "Fungal share of cell-taxa on day 365 (%)", axes[0][2])]
for col, lab, ax in panels:
    for k in ks:
        m = res[k]["sh"][col] * 100
        ax.plot(m.index, m.values, "-o", ms=3.5, color=cm[k], label=f"N x{k:g}")
    ax.set_xlabel("Pulse")
    ax.set_ylabel(lab)
axes[0][1].legend(frameon=False)
ax = axes[1][0]
for k in ks:
    m = res[k]["sh"]["tot1"]
    ax.plot(m.index, m.values, "-o", ms=3.5, color=cm[k], label=f"N x{k:g}")
ax.set_xlabel("Pulse")
ax.set_ylabel("Total starting biomass C per cell (day 1)")
ax = axes[1][1]
for k in ks:
    m = res[k]["ec"]
    ax.plot(m.index, m.values, "-o", ms=3.5, color=cm[k], label=f"N x{k:g}")
ax.set_xlabel("Pulse")
ax.set_ylabel("Energy : C on day 365")
ax = axes[1][2]
gb = [res[k]["growth"]["bacteria"] for k in ks]
gf = [res[k]["growth"]["fungi"] for k in ks]
ax.bar(x - 0.18, gb, 0.34, color=BAC, label="bacteria")
ax.bar(x + 0.18, gf, 0.34, color=FUN, label="fungi")
ax.set_yscale("log")
ax.set_xticks(x)
ax.set_xticklabels([f"x{k:g}" for k in ks])
ax.set_xlabel("Multiplier on the N of the initial substrates")
ax.set_ylabel("Pulse-1 biomass growth (day 365 / day 1)")
ax.set_title("Biomass growth in pulse 1")
ax.legend(frameon=False)
for ax in axes.flat:
    ax.spines[["top", "right"]].set_visible(False)
fig.suptitle("Raising the N supply (20x20, 5 communities, 6 pulses, real climate)", fontsize=11, fontweight="bold")
fig.tight_layout(rect=(0, 0, 1, 0.97))
fig.savefig(os.path.join(a.root, "n_supply.png"), dpi=150, bbox_inches="tight")

print("N multiplier -> mean N:C (day 100-365, pulses 1-6): bacteria / fungi:  ", {f"x{k:g}": (round(res[k]["nc"]["bacteria"].mean(), 3), round(res[k]["nc"]["fungi"].mean(), 3)) for k in ks})
print("P:C: ", {f"x{k:g}": (round(res[k]["pc"]["bacteria"].mean(), 4), round(res[k]["pc"]["fungi"].mean(), 4)) for k in ks})
print("pulse-1 biomass growth (bacteria, fungi):", {f"x{k:g}": (round(res[k]["growth"]["bacteria"], 1), round(res[k]["growth"]["fungi"], 1)) for k in ks})
for col, lab in (("C1", "fungal share of biomass C, day 1 (%)"), ("C365", "fungal share of biomass C, day 365 (%)"),
                 ("rows1", "fungal share of cells, day 1 (%)"), ("rows365", "fungal share of cells, day 365 (%)"),
                 ("tot1", "total starting biomass C per cell"), ("tot365", "total biomass C per cell, day 365")):
    print(f"\n{lab}, pulses 1-{NP}:")
    print(pd.DataFrame({f"x{k:g}": res[k]["sh"][col] * (100 if col in ("C1", "C365", "rows1", "rows365") else 1) for k in ks}).round(1).T.to_string())
print("\nenergy:C on day 365:"); print(pd.DataFrame({f"x{k:g}": res[k]["ec"] for k in ks}).round(0).T.to_string())
print("\nsubstrate degraded % on day 365:"); print(pd.DataFrame({f"x{k:g}": res[k]["deg"] for k in ks}).round(0).T.to_string())
print("\ntaxa present on day 365 (fungal / bacterial), pulse 6:", {f"x{k:g}": (round(res[k]["sh"]["fun_taxa365"].loc[NP], 1), round(res[k]["sh"]["bac_taxa365"].loc[NP], 1)) for k in ks})
