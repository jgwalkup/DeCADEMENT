"""
Real climate, slopes 2.5/2.5: osmolyte allocation range (tolerance varies) vs identical allocation (tolerance 1 for all taxa).
Usage: python plot_full_climate_km25_range.py MAIN_DIR EQ_DIR [OTHER_DIR ...]
  e.g. python plot_full_climate_km25_range.py ../analysis_output/full_climate_km25_osmorange2 ../analysis_output/full_climate_km25_osmo0 ../analysis_output/full_climate_km25_osmorange
Writes <MAIN_DIR>/km25_range_community.png (columns: identical allocation | each OTHER_DIR | MAIN_DIR; line per pulse) and
<first dir>/km25_range_taxa.png (tolerance mapping, change in cell-taxa by allocation, selection across pulses) and prints tables.
"""
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from scenario_utils import scenario_color_map  # noqa: E402

plt.rcParams.update({"font.size": 10, "axes.titlesize": 10, "axes.labelsize": 10, "legend.fontsize": 9, "figure.dpi": 150})
rng_dir, eq_dir, others = sys.argv[1], sys.argv[2], sys.argv[3:]
rng = pd.read_csv(os.path.join(rng_dir, "timeseries.csv"))
tx = pd.read_csv(os.path.join(rng_dir, "taxa.csv"))
eq = pd.read_csv(os.path.join(eq_dir, "timeseries.csv"))
n = rng["seed"].nunique()


def label(d, ts):
    t = pd.read_csv(os.path.join(d, "taxa.csv"))
    g = t[(t["pulse"] == 1) & (t["seed"] == t["seed"].min())]
    grid = "50x50" if "50x50" in d else "20x20"
    return (f"Allocation {g['alloc'].min():.3f}-{g['alloc'].max():.3f}, tolerance {g['tol'].min():.2f}-1 (mean {g['tol'].mean():.2f}), "
            f"{ts['seed'].nunique()} communities, {grid}")


sets = [(f"Identical allocation 0.055 (tolerance 1), {eq['seed'].nunique()} communities", eq)]
sets += [(label(d, pd.read_csv(os.path.join(d, "timeseries.csv"))), pd.read_csv(os.path.join(d, "timeseries.csv"))) for d in others]
sets += [(label(rng_dir, rng), rng)]
pulses = sorted(rng["pulse"].unique())
if len(pulses) <= 6:
    col = scenario_color_map([str(p) for p in pulses], palette="discrete")
else:  # more pulses than the discrete palette has colours: light -> dark gradient
    col = {str(p): c for p, c in zip(pulses, plt.cm.viridis(np.linspace(0.95, 0.05, len(pulses))))}
rows = [("alive_rows", "Cell-taxa alive"), ("taxa_alive", "Taxa present"), ("biomass_C", "Biomass C per grid cell"),
        ("substrate_remaining_pct", "Substrate C remaining (% of initial)"), ("energy_to_C", "Energy : C")]
fig, axes = plt.subplots(len(rows), len(sets), figsize=(6 * len(sets), 13), sharex=True)
for j, (title, d) in enumerate(sets):
    m = d.groupby(["pulse", "doy"], as_index=False)[[r[0] for r in rows]].mean()
    for i, (cn, lab) in enumerate(rows):
        a = axes[i][j]
        for p in pulses:
            q = m[m["pulse"] == p].sort_values("doy")
            y = q[cn].rolling(7, center=True, min_periods=1).mean() if cn == "energy_to_C" else q[cn]
            a.plot(q["doy"], y, color=col[str(p)], lw=1.2, label=f"pulse {p}")
        if i == 0:
            a.set_title(title, fontsize=8)
        if j == 0:
            a.set_ylabel(lab)
        if i == len(rows) - 1:
            a.set_xlabel("Day of pulse")
        a.set_ylim(bottom=0)
        a.spines[["top", "right"]].set_visible(False)
h, l = axes[0][0].get_legend_handles_labels()
fig.legend(h, l, loc="lower center", ncols=len(pulses), framealpha=0.85, bbox_to_anchor=(0.5, -0.005))
fig.suptitle("Full model, real climate, slopes 2.5/2.5: identical vs varied osmolyte allocation (20x20)", fontsize=11, fontweight="bold")
fig.tight_layout(rect=(0, 0.02, 1, 0.98))
fig.savefig(os.path.join(rng_dir, "km25_range_community.png"), dpi=150, bbox_inches="tight")
plt.close(fig)

fig, axes = plt.subplots(1, 3, figsize=(15, 4.6))
g = tx[(tx["pulse"] == 1) & (tx["seed"] == tx["seed"].min())].sort_values("alloc")
axes[0].plot(g["alloc"], g["tol"], color="#1f77b4", lw=1.8)
axes[0].set_xlabel("Osmolyte allocation (consti. + induci.)")
axes[0].set_ylabel("Drought tolerance")
axes[0].set_title("Tolerance assigned to each taxon (lower bound 0)")
bins = np.linspace(tx["alloc"].min(), tx["alloc"].max() + 1e-9, 7)
mid = 0.5 * (bins[:-1] + bins[1:])
for p, c in ((1, "#9ecae1"), (3, "#4292c6"), (6, "#08306b")):
    q = tx[tx["pulse"] == p].copy()
    q["bin"] = pd.cut(q["alloc"], bins, labels=False, include_lowest=True)
    s = q.groupby("bin")[["cells_day1", "cells_day365"]].sum()
    axes[1].plot(mid[s.index.astype(int)], s["cells_day365"] / s["cells_day1"], "-o", ms=4, color=c, label=f"pulse {p}")
axes[1].set_xlabel("Osmolyte allocation (binned)")
axes[1].set_ylabel("Cell-taxa on day 365 / day 1")
axes[1].set_title("Change in cell-taxa through a pulse, by allocation")
axes[1].legend(frameon=False)
a = tx.groupby(["seed", "pulse"]).apply(lambda x: np.average(x["alloc"], weights=x["cells_day1"]) if x["cells_day1"].sum() else np.nan,
                                        include_groups=False).groupby("pulse").mean()
axes[2].plot(a.index, a.values, "-o", ms=4, color="#1f77b4")
axes[2].axhline(tx["alloc"].mean(), color="gray", lw=0.8, ls="--")
axes[2].set_xlabel("Pulse")
axes[2].set_ylabel("Abundance-weighted mean allocation on day 1")
axes[2].set_title("Selection on allocation (dashed = all-taxa mean)")
for ax in axes:
    ax.spines[["top", "right"]].set_visible(False)
fig.tight_layout()
fig.savefig(os.path.join(rng_dir, "km25_range_taxa.png"), dpi=150, bbox_inches="tight")
plt.close(fig)

pd.set_option("display.width", 220)
def table(d, label):
    t = d[d.doy == 365].groupby("pulse")[["alive_rows", "taxa_alive", "biomass_C", "energy_to_C", "CUE", "substrate_remaining_pct"]].mean()
    t["degraded"] = 100 - t.pop("substrate_remaining_pct")
    t["day1_taxa"] = d[d.doy == 1].groupby("pulse")["taxa_alive"].mean()
    t["day1_biomass"] = d[d.doy == 1].groupby("pulse")["biomass_C"].mean()
    print(f"\n{label}, day 365 (mean of seeds)"); print(t.round(2).to_string())
table(eq, "Identical allocation"); table(rng, "Allocation range")
print("\nabundance-weighted mean allocation at day 1 by pulse:", a.round(4).to_dict(), "| all-taxa mean", round(tx["alloc"].mean(), 4))
z = lambda d: d[d.seed == 0].groupby("pulse").apply(lambda x: int((x.energy_to_C < 5).sum()), include_groups=False).to_dict()
print("days with energy:C < 5 (seed 0): identical", z(eq), "| range", z(rng))
for lbl, d in (("identical", eq), ("range", rng)):
    dry = d[(d.doy >= 200) & (d.doy <= 330)]
    print(f"osmolyte share of energy cost days 200-330, {lbl}:", round(dry.osmo_E_cost.sum() / dry.E_cost_total.sum(), 3))
