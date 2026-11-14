#!/usr/bin/env python3
"""
taxon_trait_relationships.py
============================
How the initialized traits of the taxa relate to what the taxa actually do,
from a per-taxon table written by taxon_trait_figures.py
(taxon_table_<treatment>_minus_<reference>.csv: one row per taxon and day,
values pooled over the communities a taxon occurs in) and taxon_traits.csv.
Reads CSV only (no pyarrow needed).

Figures  (next to the table unless --out-dir is given)
-------
  trait_correlations_<scenario>_<tag>.png
        rank (Spearman) correlation of every trait with every realized quantity,
        one panel per day; rows = traits, columns = realized quantities
  trait_scatter_<x>_vs_<y>_<scenario>_<tag>.png
        one trait against one realized quantity, one panel per day, one point
        per taxon, coloured by a second trait

Usage
-----
  python analysis/taxon_trait_relationships.py --table \\
      analysis_output/taxon_snapshots_summer_heatwaves_20261110/figures/taxa/taxon_table_climate_July_heatwave_minus_climate_ambient_temp.csv
      [--scenario treatment|reference] [--days 1845 2013 2055 2190]
      [--x uptake_trait] [--y hsp_share] [--color hsp_prod_consti]
      [--traits-file DIR/taxon_traits.csv] [--out-dir DIR]
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

TRAIT_LABELS = {
    "uptake_trait":     "Uptake (transporter per biomass)",
    "uptake_genes":     "Transporter types",
    "enz_genes":        "Enzyme types",
    "enz_prod_consti":  "Enzyme production, constitutive",
    "enz_prod_induci":  "Enzyme production, inducible",
    "hsp_prod_consti":  "HSP production, constitutive",
    "hsp_prod_induci":  "HSP production, inducible",
    "osmo_prod_consti": "Osmolyte production, constitutive",
    "osmo_prod_induci": "Osmolyte production, inducible",
}
# realized quantities: column stem in the table -> label
REALIZED = {
    "biomass_C":     "Biomass C",
    "alive_share":   "Communities\nwhere alive",
    "net_growth":    "Net specific\ngrowth rate",
    "growth_yield":  "Growth yield",
    "hsp_share":     "HSP pool share\nof biomass",
    "hsp_C_pb":      "HSP production\nper biomass",
    "enz_C_pb":      "Enzyme production\nper biomass",
    "osmo_C_pb":     "Osmolyte production\nper biomass",
    "osmo_share":    "Osmolyte pool share\nof biomass",
    "uptake_C_pb":   "C uptake\nper biomass",
    "energy_frac":   "Share of uptake C\nsent to energy",
    "prod_E_pb":     "Energy production\nper biomass",
    "transport_E_pb": "Transporter energy\ncost per biomass",
    "enz_E_pb":      "Enzyme energy\ncost per biomass",
    "hsp_E_pb":      "HSP energy\ncost per biomass",
    "osmo_E_pb":     "Osmolyte energy\ncost per biomass",
    "maint_E_pb":    "Maintenance cost\nper biomass",
    "maint_share_E": "Maintenance /\nenergy production",
    "endo_pb":       "Endogenous resp.\nper biomass",
}


def rank_corr(a: pd.Series, b: pd.Series) -> float:
    """Spearman correlation without scipy: Pearson correlation of the ranks."""
    ok = a.notna() & b.notna()
    if ok.sum() < 5 or a[ok].nunique() < 2 or b[ok].nunique() < 2:
        return np.nan
    return float(a[ok].rank().corr(b[ok].rank()))


def fig_correlations(table, traits, days, suffix, title, out):
    cols = [c for c in REALIZED if f"{c}_{suffix}" in table.columns]
    rows = [t for t in TRAIT_LABELS if t in traits.columns and traits[t].nunique() > 1]
    ncol = len(days) if len(cols) <= 10 else min(len(days), 2)      # wide grids: two panels per row
    nrow = -(-len(days) // ncol)
    fig, axes = plt.subplots(nrow, ncol, squeeze=False, layout="constrained", sharey=True,
                             figsize=(2.2 + (0.5 * len(cols) + 0.4) * ncol, (0.5 * len(rows) + 2.3) * nrow))
    for ax in axes.ravel()[len(days):]:
        ax.set_visible(False)
    for ax, day in zip(axes.ravel(), days):
        q = table[table.day == day].set_index("taxon")
        tr = traits.reindex(q.index)
        m = np.array([[rank_corr(tr[r], q[f"{c}_{suffix}"]) for c in cols] for r in rows])
        im = ax.imshow(m, cmap="RdBu", vmin=-1, vmax=1, aspect="auto")
        for i in range(len(rows)):
            for j in range(len(cols)):
                if np.isfinite(m[i, j]):
                    ax.text(j, i, f"{m[i, j]:+.2f}", ha="center", va="center", fontsize=7,
                            color="w" if abs(m[i, j]) > 0.6 else "k")
        ax.set_xticks(range(len(cols)))
        ax.set_xticklabels([REALIZED[c] for c in cols], rotation=60, ha="right", fontsize=8)
        ax.set_yticks(range(len(rows)))
        ax.set_yticklabels([TRAIT_LABELS[r] for r in rows], fontsize=8)
        ax.set_title(f"day {day}", fontsize=10)
    fig.colorbar(im, ax=list(axes.ravel()), label="rank correlation across taxa", shrink=0.8, pad=0.01)
    fig.suptitle(title)
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"Saved {out}")


def fig_scatter(table, traits, days, suffix, x, y, color, title, out):
    ycol = f"{y}_{suffix}"
    if ycol not in table.columns:
        sys.exit(f"ERROR: {ycol!r} is not a column of the table")
    norm = plt.Normalize(traits[color].min(), traits[color].max())
    fig, axes = plt.subplots(1, len(days), figsize=(4.4 * len(days) + 1, 4.3), squeeze=False,
                             layout="constrained")
    for ax, day in zip(axes.ravel(), days):
        q = table[table.day == day].set_index("taxon")
        tr = traits.reindex(q.index)
        sc = ax.scatter(tr[x], q[ycol], c=tr[color], cmap="viridis", norm=norm, s=16)
        ax.set_title(f"day {day}  (ρ = {rank_corr(tr[x], q[ycol]):+.2f})", fontsize=10)
        ax.set_xlabel(TRAIT_LABELS.get(x, x))
        ax.set_ylabel(REALIZED.get(y, y).replace("\n", " "))
    fig.colorbar(sc, ax=list(axes.ravel()), label=TRAIT_LABELS.get(color, color), shrink=0.85, pad=0.01)
    fig.suptitle(title)
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"Saved {out}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--table", required=True, help="taxon_table_<treatment>_minus_<reference>.csv")
    ap.add_argument("--traits-file", default=None,
                    help="taxon_traits.csv (default: in the snapshot folder two levels above the table)")
    ap.add_argument("--scenario", choices=("treatment", "reference"), default="treatment",
                    help="which scenario of the table's contrast to show")
    ap.add_argument("--days", type=int, nargs="+", default=[1845, 2013, 2055, 2190])
    ap.add_argument("--x", default="uptake_trait", help="trait on the x axis of the scatter")
    ap.add_argument("--y", default="hsp_share", help="realized quantity on the y axis of the scatter")
    ap.add_argument("--color", default="hsp_prod_consti", help="trait that colours the scatter points")
    ap.add_argument("--out-dir", default=None, help="default: the table's folder")
    a = ap.parse_args()

    path = Path(a.table)
    table = pd.read_csv(path)
    tfile = Path(a.traits_file) if a.traits_file else path.parents[2] / "taxon_traits.csv"
    traits = pd.read_csv(tfile, index_col=0)
    for t in (a.x, a.color):
        if t not in traits.columns:
            sys.exit(f"ERROR: trait {t!r} not in {tfile.name} ({', '.join(traits.columns)})")
    missing = [d for d in a.days if d not in set(table.day)]
    if missing:
        sys.exit(f"ERROR: day(s) {missing} not in the table")
    out_dir = Path(a.out_dir) if a.out_dir else path.parent
    out_dir.mkdir(parents=True, exist_ok=True)
    tag = path.stem.replace("taxon_table_", "")
    treat, ref = tag.split("_minus_")
    scen = treat if a.scenario == "treatment" else ref
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from scenario_utils import pretty_label
    name = pretty_label(scen)

    fig_correlations(table, traits, a.days, a.scenario, f"{name}: initialized traits against realized quantities "
                     f"(one value per taxon, pooled over its communities)",
                     out_dir / f"trait_correlations_{scen}_{tag}.png")
    fig_scatter(table, traits, a.days, a.scenario, a.x, a.y, a.color,
                f"{name}: {TRAIT_LABELS.get(a.x, a.x)} against {REALIZED.get(a.y, a.y).replace(chr(10), ' ')}, "
                f"one point per taxon",
                out_dir / f"trait_scatter_{a.x}_vs_{a.y}_{scen}_{tag}.png")


if __name__ == "__main__":
    main()
