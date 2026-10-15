"""
Figures on how fungi gain share and come to dominate across pulses, from fungal_tracking.py output.

Usage: python plot_fungal_dominance.py ../analysis_output/fungal_tracking
  expects baseline/ (12 pulses) and optionally stoich_f_eq_b/, size_f_eq_b/, stoich_and_size_eq_b/ (6 pulses) in that folder.
Writes fungal_within_pulse.png, fungal_ratchet.png, fungal_variants.png and prints the numbers used in the text.
"""
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

plt.rcParams.update({"font.size": 9.5, "axes.titlesize": 10, "axes.labelsize": 9.5, "legend.fontsize": 8.5, "figure.dpi": 150})
BAC, FUN = "#1f77b4", "#d6172a"
root = sys.argv[1]
g = pd.read_csv(os.path.join(root, "baseline", "groups.csv"))
t = pd.read_csv(os.path.join(root, "baseline", "taxa_track.csv"))
n = g["seed"].nunique()
PULSES = [1, 4, 12]
pd.set_option("display.width", 220)


def sm(x, w=7):
    return x.rolling(w, center=True, min_periods=1).mean()


# ---------------- figure 1: within-pulse dynamics by group
rows = [("cells", "Cell-taxa alive", None), ("biomass_C", "Biomass C per grid cell", None),
        ("mean_C_per_cell", "Mean C per cell", None), ("N_to_C", "N : C of the biomass", None),
        ("rate", "Births (solid), deaths (dashed)\nper 1,000 cells per day", None)]
fig, axes = plt.subplots(len(rows), len(PULSES), figsize=(4.6 * len(PULSES), 2.7 * len(rows)), sharex=True)
for j, p in enumerate(PULSES):
    d = g[g["pulse"] == p].groupby(["group", "doy"]).mean(numeric_only=True).reset_index()
    for i, (col, lab, _) in enumerate(rows):
        ax = axes[i][j]
        for grp, c in (("bacteria", BAC), ("fungi", FUN)):
            q = d[d["group"] == grp].sort_values("doy")
            if col == "rate":
                cells = q["cells"].replace(0, np.nan)
                ax.plot(q["doy"], sm(1000 * q["new_cells"] / cells), color=c, lw=1.3, label=f"{grp}: new cells")
                ax.plot(q["doy"], sm(1000 * (q["starve_deaths"] + q["stochastic_deaths"]) / cells), color=c, lw=1.3, ls="--",
                        label=f"{grp}: deaths")
            else:
                ax.plot(q["doy"], q[col] if col != "cells" else q[col], color=c, lw=1.4, label=grp)
        if col == "N_to_C":
            ax.axhline(0.16 / 0.825, color=BAC, lw=0.8, ls=":")
            ax.axhline(0.09 / 0.9, color=FUN, lw=0.8, ls=":")
        if col in ("cells", "biomass_C"):
            ax.set_yscale("log")
        if i == 0:
            ax.set_title(f"Pulse {p}")
        if j == 0:
            ax.set_ylabel(lab)
        if i == len(rows) - 1:
            ax.set_xlabel("Day of pulse")
        ax.spines[["top", "right"]].set_visible(False)
axes[0][0].legend(frameon=False)
axes[4][0].legend(frameon=False, fontsize=7.5, ncols=2)
fig.suptitle(f"Bacteria vs fungi within a pulse (20x20, mean of {n} communities; dotted = target N:C of each group)", fontsize=11,
             fontweight="bold")
fig.tight_layout(rect=(0, 0, 1, 0.985))
fig.savefig(os.path.join(root, "fungal_within_pulse.png"), dpi=150, bbox_inches="tight")
plt.close(fig)

# ---------------- figure 2: the ratchet across pulses
def shares(tt):
    out = []
    for (s, p), x in tt.groupby(["seed", "pulse"]):
        f = x[x["fungal"] == 1]
        out.append(dict(seed=s, pulse=p,
                        rows_day1=f["cells_day1"].sum() / max(x["cells_day1"].sum(), 1),
                        rows_day365=f["cells_day365"].sum() / max(x["cells_day365"].sum(), 1),
                        reseed=f["reseed_freq"].sum(),
                        C_day1=f["C_day1"].sum() / max(x["C_day1"].sum(), 1e-9),
                        C_day365=f["C_day365"].sum() / max(x["C_day365"].sum(), 1e-9),
                        taxa_fun_day365=int(((x["fungal"] == 1) & (x["cells_day365"] > 0)).sum()),
                        taxa_bac_day365=int(((x["fungal"] == 0) & (x["cells_day365"] > 0)).sum())))
    return pd.DataFrame(out)


sh = shares(t)
fig, axes = plt.subplots(2, 3, figsize=(15, 8.6))
ax = axes[0][0]
for col, lab, c, ls in (("rows_day1", "day 1 (reseeded)", "#9ecae1", "-"), ("rows_day365", "day 365", "#08519c", "-"),
                        ("reseed", "reseeding frequency for the next pulse", "#d6172a", "--")):
    m = sh.groupby("pulse")[col].agg(["mean", "min", "max"]) * 100
    ax.plot(m.index, m["mean"], ls, color=c, marker="o", ms=3.5, label=lab)
    ax.fill_between(m.index, m["min"], m["max"], color=c, alpha=0.12, lw=0)
ax.set_xlabel("Pulse")
ax.set_ylabel("Fungal share of cell-taxa (%)")
ax.set_title("Within each pulse fungi gain share; the next pulse is seeded from it")
ax.legend(frameon=False)
ax = axes[0][1]
for col, lab, c in (("C_day1", "day 1", "#9ecae1"), ("C_day365", "day 365", "#08519c")):
    m = sh.groupby("pulse")[col].mean() * 100
    ax.plot(m.index, m.values, "-o", ms=3.5, color=c, label=lab)
ax.set_xlabel("Pulse")
ax.set_ylabel("Fungal share of biomass C (%)")
ax.set_title("Fungal share of biomass")
ax.legend(frameon=False)
ax = axes[0][2]
tb = t.copy()
tb["growth"] = tb["cells_day365"] / tb["cells_day1"].replace(0, np.nan)
for grp, c, off in ((0, BAC, -0.18), (1, FUN, 0.18)):
    data = [tb[(tb["pulse"] == p) & (tb["fungal"] == grp)]["growth"].dropna().to_numpy() for p in range(1, 13)]
    ax.boxplot(data, positions=np.arange(1, 13) + off, widths=0.32, showfliers=False, patch_artist=True,
               boxprops=dict(facecolor=c, alpha=0.5, lw=0.6), medianprops=dict(color="k", lw=1), whiskerprops=dict(lw=0.6), capprops=dict(lw=0.6))
ax.axhline(1, color="gray", lw=0.8, ls="--")
ax.set_yscale("log")
ax.set_xticks(range(1, 13))
ax.set_xticklabels([str(k) for k in range(1, 13)])
ax.set_xlabel("Pulse")
ax.set_ylabel("Cells on day 365 / cells on day 1, per taxon")
ax.set_title("Growth in cell number of each taxon (blue bacteria, red fungi)")
# rank-abundance for one community in three pulses
seed = int(sh.groupby("seed")["rows_day1"].last().sort_values().index[len(sh["seed"].unique()) // 2])
for k, p in enumerate(PULSES):
    ax = axes[1][k]
    x = t[(t["seed"] == seed) & (t["pulse"] == p)].sort_values("reseed_freq", ascending=False).reset_index(drop=True)
    x = x[x["reseed_freq"] > 0]
    ax.bar(np.arange(len(x)), 100 * x["reseed_freq"], color=[FUN if f else BAC for f in x["fungal"]], width=0.9)
    ax.set_xlabel("Taxa ranked by their share of the cumulative cell counts")
    ax.set_ylabel("Share of cumulative cell counts (%)" if k == 0 else "")
    ax.set_title(f"Community {seed}, pulse {p}: {int(x['fungal'].sum())} fungal of {len(x)} taxa present, "
                 f"{100 * x.loc[x['fungal'] == 1, 'reseed_freq'].sum():.0f} % fungal")
    ax.spines[["top", "right"]].set_visible(False)
for ax in axes.flat[:3]:
    ax.spines[["top", "right"]].set_visible(False)
fig.suptitle(f"How fungal share builds up across pulses (20x20, {n} communities)", fontsize=11, fontweight="bold")
fig.tight_layout(rect=(0, 0, 1, 0.985))
fig.savefig(os.path.join(root, "fungal_ratchet.png"), dpi=150, bbox_inches="tight")
plt.close(fig)

print("fungal share by pulse (mean over communities):")
print((sh.groupby("pulse").mean().drop(columns="seed") * [100, 100, 100, 100, 100, 1, 1]).round(1).to_string())
print("\nmedian growth factor in cells (day 365 / day 1) by taxon, pulses 1-6:")
print(tb[tb.pulse <= 6].groupby(["pulse", "fungal"])["growth"].median().unstack().round(2).to_string())

# daily rates and totals per pulse by group
tot = g.groupby(["pulse", "seed", "group"])[["new_cells", "starve_deaths", "stochastic_deaths"]].sum()
c1 = g[g.doy == 1].set_index(["pulse", "seed", "group"])["cells"]
mc = g.groupby(["pulse", "seed", "group"])["cells"].mean()
rate = (tot.div(mc, axis=0) * 1).groupby(["pulse", "group"]).mean()
print("\nper cell present on average over the year (mean over communities): new cells / starvation deaths / drought+basal deaths")
print(rate.loc[[1, 4, 12]].round(3).to_string())
nc = g.groupby(["pulse", "group", "doy"]).mean(numeric_only=True).reset_index()
last = nc[nc.doy == 365].set_index(["pulse", "group"])[["cells", "biomass_C", "mean_C_per_cell", "N_to_C", "P_to_C"]]
first = nc[nc.doy == 1].set_index(["pulse", "group"])[["cells", "biomass_C", "mean_C_per_cell"]]
print("\nday 1 -> day 365 by group (mean over communities), pulses 1, 4, 12:")
for p in PULSES:
    for grp in ("bacteria", "fungi"):
        a, b = first.loc[(p, grp)], last.loc[(p, grp)]
        print(f"  pulse {p:2d} {grp:8s} cells {a.cells:8.0f} -> {b.cells:8.0f}  biomass {a.biomass_C:6.2f} -> {b.biomass_C:6.2f}  mean C/cell {a.mean_C_per_cell:5.2f} -> {b.mean_C_per_cell:5.2f}  N:C {b.N_to_C:.3f}  P:C {b.P_to_C:.4f}")

# ---------------- figure 3: variants
var = {"baseline": "baseline (fb 0.2)", "stoich_f_eq_b": "fungal C:N:P = bacterial", "size_f_eq_b": "fungal size = bacterial (max_size_f 2)",
       "stoich_and_size_eq_b": "both"}
cols = {"baseline": "#d6172a", "stoich_f_eq_b": "#ff7f0e", "size_f_eq_b": "#2ca02c", "stoich_and_size_eq_b": "#7f7f7f"}
have = {k: v for k, v in var.items() if os.path.exists(os.path.join(root, k, "taxa_track.csv"))}
if len(have) > 1:
    fig, axes = plt.subplots(1, 4, figsize=(18, 4.4))
    print("\nvariants (mean over communities):")
    for k, lab in have.items():
        tk = pd.read_csv(os.path.join(root, k, "taxa_track.csv"))
        gk = pd.read_csv(os.path.join(root, k, "groups.csv"))
        s = shares(tk).query("pulse <= 6")
        m = s.groupby("pulse").mean()
        e = gk[(gk.doy == 365)].groupby(["pulse", "seed"])["biomass_C"].sum().groupby("pulse").mean().loc[:6]
        b1 = gk[(gk.doy == 1)].groupby(["pulse", "seed"])["biomass_C"].sum().groupby("pulse").mean().loc[:6]
        ax_list = [("C_day1", "Fungal share of biomass C on day 1 (%)", 100), ("C_day365", "Fungal share of biomass C on day 365 (%)", 100)]
        for ax, (col, ylab, f) in zip(axes[:2], ax_list):
            ax.plot(m.index, f * m[col], "-o", ms=3.5, color=cols[k], label=lab)
            ax.set_ylabel(ylab)
        axes[2].plot(b1.index, b1.values, "-o", ms=3.5, color=cols[k], label=lab)
        axes[2].set_ylabel("Total biomass C per cell on day 1")
        axes[3].plot(e.index, e.values, "-o", ms=3.5, color=cols[k], label=lab)
        axes[3].set_ylabel("Total biomass C per cell on day 365")
        print(f"  {lab:40s} fungal share of cells day 1: {' '.join(f'{100*v:.0f}' for v in m.rows_day1)} | day 365: "
              f"{' '.join(f'{100*v:.0f}' for v in m.rows_day365)} || share of biomass C day 1: {' '.join(f'{100*v:.0f}' for v in m.C_day1)} | day 365: "
              f"{' '.join(f'{100*v:.0f}' for v in m.C_day365)} || total biomass day 1: {' '.join(f'{v:.1f}' for v in b1)} | day 365: {' '.join(f'{v:.0f}' for v in e)}")
    for ax in axes:
        ax.set_xlabel("Pulse")
        ax.spines[["top", "right"]].set_visible(False)
    axes[0].legend(frameon=False)
    fig.suptitle("Which fungal traits matter? Setting fungal parameters equal to the bacterial ones (20x20, 5 communities)",
                 fontsize=11, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(os.path.join(root, "fungal_variants.png"), dpi=150, bbox_inches="tight")
