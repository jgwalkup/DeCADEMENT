"""
Figures and summary tables for full_climate_osmo_range.py output.

Usage:
    python plot_full_climate_osmo_range.py ../analysis_output/full_climate_osmo_range
Writes <dir>/osmo_range_community.png (community time series), <dir>/osmo_range_taxa.png (tolerance mapping, survival by
allocation, selection across pulses) and prints summary tables.
"""
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

plt.rcParams.update({"font.size": 10, "axes.titlesize": 11, "axes.labelsize": 10, "legend.fontsize": 9, "figure.dpi": 150})
LABEL = {"rel_min": "relative, lower = community min (original)", "rel_lower0": "relative, lower = 0 (new default)",
         "abs0.055": "absolute, ref 0.055 (mean alloc.)", "abs0.095": "absolute, ref 0.095 (max alloc.)"}
COL = {"rel_min": "#d62728", "rel_lower0": "#ff7f0e", "abs0.055": "#1f77b4", "abs0.095": "#2ca02c"}
ORDER = list(LABEL)


def lab(r):
    return LABEL.get(r, str(r))


def add_cond(df):
    """Older output (before the Osmo_tol_lower option) has only `ref`: 0 = original relative method."""
    if "cond" not in df:
        df["cond"] = np.where(df["ref"] > 0, "abs" + df["ref"].map("{:g}".format), "rel_min")
    return df


def main(d):
    ts = add_cond(pd.read_csv(os.path.join(d, "timeseries.csv")))
    tx = add_cond(pd.read_csv(os.path.join(d, "taxa.csv")))
    refs = [c for c in ORDER if c in set(ts["cond"])]
    seeds = ts["seed"].nunique()

    # ---- community time series
    fig, axes = plt.subplots(2, 4, figsize=(16, 7), sharex=True)
    fig.suptitle(f"Full model, real climate, osmolyte allocation varied across taxa (0.015-0.095): relative vs absolute "
                 f"drought tolerance, mean of {seeds} communities", fontsize=11, fontweight="bold")
    cols = [("alive_rows", "Cell-taxa alive"), ("taxa_alive", "Taxa present"), ("biomass_C", "Biomass C per cell"),
            ("energy_to_C", "Energy : C")]
    for r, (name, sel) in enumerate((("Pulse 1", ts["pulse"] == 1), ("Pulses 2-6 (mean)", ts["pulse"] > 1))):
        for c, (col, title) in enumerate(cols):
            ax = axes[r][c]
            for ref in refs:
                m = ts[sel & (ts["cond"] == ref)].groupby("doy")[col].mean()
                ax.plot(m.index, m.values, color=COL.get(ref, "k"), lw=1.5, label=lab(ref))
            ax.set_title(f"{title}, {name}", fontsize=10)
            ax.spines[["top", "right"]].set_visible(False)
            if r == 1:
                ax.set_xlabel("Day of year")
    h, l = axes[0][0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncols=2, framealpha=0.85, bbox_to_anchor=(0.5, -0.02))
    fig.tight_layout()
    fig.savefig(os.path.join(d, "osmo_range_community.png"), dpi=150, bbox_inches="tight")
    plt.close(fig)

    # ---- taxon-level figure
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.8))
    ax = axes[0]
    p1 = tx[(tx["pulse"] == 1) & (tx["seed"] == tx["seed"].min())]
    for ref in refs:
        g = p1[p1["cond"] == ref].sort_values("alloc")
        ax.plot(g["alloc"], g["tol"], color=COL.get(ref, "k"), lw=1.8, label=lab(ref))
    ax.set_xlabel("Osmolyte allocation (consti. + induci.)")
    ax.set_ylabel("Drought tolerance")
    ax.set_title("Tolerance assigned to each taxon", fontsize=10)
    ax.legend(frameon=False)
    ax.spines[["top", "right"]].set_visible(False)

    ax = axes[1]
    bins = np.linspace(tx["alloc"].min(), tx["alloc"].max() + 1e-9, 7)
    mid = 0.5 * (bins[:-1] + bins[1:])
    for ref in refs:
        g = tx[(tx["pulse"] == 1) & (tx["cond"] == ref)].copy()
        g["bin"] = pd.cut(g["alloc"], bins, labels=False, include_lowest=True)
        s = g.groupby("bin")[["cells_day1", "cells_day365"]].sum()
        ax.plot(mid[s.index.astype(int)], s["cells_day365"] / s["cells_day1"], "-o", ms=4, color=COL.get(ref, "k"), label=lab(ref))
    ax.set_xlabel("Osmolyte allocation (binned)")
    ax.set_ylabel("Cell-taxa on day 365 / day 1 (pulse 1)")
    ax.set_title("Survival through the first year by allocation", fontsize=10)
    ax.spines[["top", "right"]].set_visible(False)

    ax = axes[2]
    for ref in refs:
        g = tx[tx["cond"] == ref].copy()
        g["w"] = g["cells_day1"]
        a = g.groupby(["seed", "pulse"]).apply(lambda x: np.average(x["alloc"], weights=x["w"]) if x["w"].sum() else np.nan,
                                               include_groups=False).groupby("pulse").mean()
        ax.plot(a.index, a.values, "-o", ms=4, color=COL.get(ref, "k"), label=lab(ref))
    ax.axhline(tx["alloc"].mean(), color="gray", lw=0.8, ls="--")
    ax.set_xlabel("Pulse")
    ax.set_ylabel("Abundance-weighted mean allocation on day 1")
    ax.set_title("Selection on allocation across pulses (dashed = all-taxa mean)", fontsize=10)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(os.path.join(d, "osmo_range_taxa.png"), dpi=150, bbox_inches="tight")
    plt.close(fig)

    # ---- tables
    end = ts[ts["doy"] == 365]
    t = end.assign(grp=np.where(end["pulse"] == 1, "pulse 1", "pulses 2-6")).groupby(["cond", "grp"])[
        ["alive_rows", "taxa_alive", "biomass_C", "energy_to_C", "substrate_remaining_pct"]].mean()
    t["degraded_pct"] = 100 - t.pop("substrate_remaining_pct")
    print("Day-365 values (mean of seeds; pulses 2-6 also averaged over pulses):")
    print(t.round(2).to_string())
    dry = ts[(ts["doy"] >= 200) & (ts["doy"] <= 330)]
    print("\nOsmolyte share of total energy cost, days 200-330:")
    print((dry.groupby("cond")["osmo_E_cost"].sum() / dry.groupby("cond")["E_cost_total"].sum()).round(3).to_string())
    print("\nAbundance-weighted mean allocation on day 1 of each pulse (all-taxa mean %.4f):" % tx["alloc"].mean())
    rows = {}
    for ref in refs:
        g = tx[tx["cond"] == ref]
        rows[ref] = g.groupby(["seed", "pulse"]).apply(
            lambda x: np.average(x["alloc"], weights=x["cells_day1"]) if x["cells_day1"].sum() else np.nan,
            include_groups=False).groupby("pulse").mean()
    print(pd.DataFrame(rows).round(4).to_string())
    print("\nTaxa present on day 1 of each pulse:")
    print(pd.DataFrame({ref: tx[tx["cond"] == ref].groupby(["seed", "pulse"]).apply(
        lambda x: int((x["cells_day1"] > 0).sum()), include_groups=False).groupby("pulse").mean() for ref in refs}).round(1).to_string())


if __name__ == "__main__":
    main(sys.argv[1])
