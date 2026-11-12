#!/usr/bin/env python3
"""
taxon_snapshot_figures.py
=========================
Heatwave effects by taxon trait, from the per-taxon snapshot tables written by
export_taxon_snapshots.py.

Taxa are grouped into equal-sized classes of one trait (default: inducible HSP
production, hsp_prod_induci; class 1 = lowest).  For each contrast between two
scenarios that differ by one heatwave (treatment:reference), every quantity is
computed per community and trait class in both scenarios and the paired
difference (treatment − reference, same community) is averaged over the
communities; error bars are 95 % confidence intervals across communities.

Quantities per community and trait class (sums over the class's taxa):
  biomass C, living cells      % change relative to the reference scenario
  taxa alive                   difference in the number of taxa with living cells
  HSP pool share               HSP pool C / biomass C
  HSP, maintenance energy cost and energy production     per unit biomass C
  endogenous respiration       per unit biomass C
  CUE                          retained C / C taken up over the trailing window
                               of the snapshot export (7 days by default)
and per community (all taxa): the biomass-weighted mean of the trait, which
shows whether a heatwave shifts the community towards high- or low-trait taxa.

Output  (<out-dir>, default <snap-dir>/figures)
------
  by_day_<trait>_<treatment>_minus_<reference>.png   one per contrast: each
                             quantity against day, one line per trait class
  by_class_<trait>_day<D>.png   each quantity against trait class on day D
                             (--at-day), one line per contrast
  taxon_contrasts_<trait>.csv   everything plotted (long format)

Usage
-----
  python analysis/taxon_snapshot_figures.py \\
      --snap-dir analysis_output/taxon_snapshots_summer_heatwaves_20261110
      [--contrasts TREATMENT:REFERENCE ...] [--trait hsp_prod_induci] [--classes 5]
      [--at-day 2055] [--days START END] [--heatwave june:2020-2026 ...] [--out-dir DIR]

Reads .parquet tables (needs pyarrow) or .csv.gz ones, whichever the snapshot
folder holds.
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

_HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(_HERE))
from heatwave_delta_figures import parse_heatwaves      # noqa: E402
from scenario_utils import pretty_label                 # noqa: E402

# the July wave of the summer heatwave-timing experiment, by heatwave history
DEFAULT_CONTRASTS = ["climate_July_heatwave:climate_ambient_temp",
                     "climate_June_July_heatwave:climate_June_heatwave",
                     "climate_July_August_heatwave:climate_August_heatwave",
                     "climate_all_three_heatwave:climate_June_August_heatwave"]
CONTRAST_STYLES = [dict(color="#1f77b4", ls="-"), dict(color="#d62728", ls="-"),
                   dict(color="#1f77b4", ls="--"), dict(color="#d62728", ls="--")]

SUMMED = ["cells", "biomass_C", "hsp_pool_C", "hsp_E", "maintenance_E", "production_E", "endo_resp_C"]

# (key, panel title, y label of the paired difference, relative change in %)
QUANTITIES = [
    ("biomass_C",     "Biomass C",                "% of reference",      True),
    ("cells",         "Living cells",             "% of reference",      True),
    ("taxa_alive",    "Taxa alive",               "taxa",                False),
    ("hsp_share",     "HSP pool share of biomass", "C C⁻¹",              False),
    ("hsp_E_pb",      "HSP energy cost",          "E d⁻¹ per biomass C", False),
    ("maint_E_pb",    "Maintenance energy cost",  "E d⁻¹ per biomass C", False),
    ("prod_E_pb",     "Energy production",        "E d⁻¹ per biomass C", False),
    ("endo_pb",       "Endogenous respiration",   "C d⁻¹ per biomass C", False),
    ("cue",           "CUE ({w}-day)",            "CUE",                 False),
    ("weighted_trait", "Biomass-weighted mean trait (all taxa)", "trait units", False),
]


# ── Loading ───────────────────────────────────────────────────────────────────

def read_snapshots(snap_dir: Path, scenario: str, columns: list[str], days) -> pd.DataFrame:
    pq, gz = (snap_dir / f"taxon_snapshots_{scenario}.{e}" for e in ("parquet", "csv.gz"))
    if pq.exists():
        df = pd.read_parquet(pq, columns=columns)
    elif gz.exists():
        df = pd.read_csv(gz, usecols=columns)
    else:
        sys.exit(f"ERROR: no snapshot table for scenario {scenario!r} in {snap_dir}")
    if days:
        df = df[df.day.between(days[0], days[1])]
    return df


def class_table(snap_dir: Path, scenario: str, klass: pd.Series, trait: pd.Series,
                window: int, days) -> pd.DataFrame:
    """Quantities per community x trait class x day for one scenario, plus the
    community-level biomass-weighted trait mean (class 0)."""
    w = [f"uptake_C_{window}d", f"retained_C_{window}d"]
    df = read_snapshots(snap_dir, scenario, ["community", "taxon", "day"] + SUMMED + w, days)
    df["klass"] = df.taxon.map(klass)
    df["taxa_alive"] = (df.cells > 0).astype(np.int32)
    df["trait_x_biomass"] = df.taxon.map(trait) * df.biomass_C
    g = df.groupby(["community", "klass", "day"], sort=False)[SUMMED + w + ["taxa_alive"]].sum()
    with np.errstate(divide="ignore", invalid="ignore"):
        out = g[["biomass_C", "cells", "taxa_alive"]].astype(float)
        out["hsp_share"] = g.hsp_pool_C / g.biomass_C
        out["hsp_E_pb"] = g.hsp_E / g.biomass_C
        out["maint_E_pb"] = g.maintenance_E / g.biomass_C
        out["prod_E_pb"] = g.production_E / g.biomass_C
        out["endo_pb"] = g.endo_resp_C / g.biomass_C
        out["cue"] = g[w[1]] / g[w[0]]
        c = df.groupby(["community", "day"], sort=False)[["trait_x_biomass", "biomass_C"]].sum()
        comm = (c.trait_x_biomass / c.biomass_C).rename("weighted_trait").to_frame()
    comm["klass"] = 0
    out = pd.concat([out, comm.set_index("klass", append=True).reorder_levels(["community", "klass", "day"])])
    return out.replace([np.inf, -np.inf], np.nan)


def contrast_stats(treat: pd.DataFrame, ref: pd.DataFrame) -> pd.DataFrame:
    """Mean and 95 % CI across communities of the paired difference, per class,
    day and quantity (long format)."""
    ref = ref.reindex(treat.index)
    rows = []
    for key, _, _, relative in QUANTITIES:
        with np.errstate(divide="ignore", invalid="ignore"):
            d = (treat[key] - ref[key]) / ref[key] * 100 if relative else treat[key] - ref[key]
        d = d.replace([np.inf, -np.inf], np.nan).dropna()
        s = d.groupby(level=["klass", "day"]).agg(["mean", "std", "count"])
        s["ci95"] = 1.96 * s["std"] / np.sqrt(s["count"])
        rows.append(s[["mean", "ci95", "count"]].assign(quantity=key).reset_index())
    return pd.concat(rows, ignore_index=True)


# ── Figures ───────────────────────────────────────────────────────────────────

def _shade(ax, heatwaves, lo, hi):
    for _, a, b in heatwaves:                 # climate rows a..b = output days a+1..b+1
        if b + 1 >= lo and a <= hi:
            ax.axvspan(a, b + 1, color="0.88", zorder=0)


def _grid(n, ncol=5):
    nrow = -(-n // ncol)
    fig, axes = plt.subplots(nrow, ncol, figsize=(4.3 * ncol, 3.2 * nrow), squeeze=False)
    for ax in axes.ravel()[n:]:
        ax.set_visible(False)
    return fig, axes.ravel()[:n]


def fig_by_day(stats, label, class_labels, heatwaves, window, trait, out):
    days = np.sort(stats.day.unique())
    lo, hi = days[0], days[-1]
    colors = plt.cm.viridis(np.linspace(0, 0.92, len(class_labels)))
    fig, axes = _grid(len(QUANTITIES))
    for ax, (key, title, ylab, _) in zip(axes, QUANTITIES):
        q = stats[stats.quantity == key]
        for k in sorted(q.klass.unique()):
            s = q[q.klass == k].sort_values("day")
            kw = (dict(color="k", label="all taxa") if k == 0
                  else dict(color=colors[k - 1], label=class_labels[k - 1]))
            ax.errorbar(s.day, s["mean"], yerr=s.ci95, lw=1.2, marker="o", ms=3, capsize=0, **kw)
        _shade(ax, heatwaves, lo, hi)
        ax.axhline(0, color="k", lw=0.7)
        ax.set_title(f"Δ {title.format(w=window)}", fontsize=10)
        ax.set_ylabel(ylab, fontsize=9)
        ax.set_xlabel("day")
        ax.tick_params(labelsize=8)
    h, l = axes[0].get_legend_handles_labels()
    fig.suptitle(f"{label}, by class of {trait} (paired difference per community, mean ± 95 % CI)", y=1.0)
    fig.tight_layout()
    fig.legend(h, l, title=f"{trait} class (range)", loc="upper center", bbox_to_anchor=(0.5, 0.0),
               ncols=len(class_labels), frameon=False, fontsize=8)
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {out}")


def fig_by_class(all_stats, day, class_labels, window, trait, out):
    quantities = [q for q in QUANTITIES if q[0] != "weighted_trait"]
    fig, axes = _grid(len(quantities))
    x = np.arange(1, len(class_labels) + 1)
    for ax, (key, title, ylab, _) in zip(axes, quantities):
        for (label, stats), style in zip(all_stats.items(), CONTRAST_STYLES):
            s = stats[(stats.quantity == key) & (stats.day == day)].set_index("klass").reindex(x)
            ax.errorbar(x, s["mean"], yerr=s.ci95, lw=1.3, marker="o", ms=4, label=label, **style)
        ax.axhline(0, color="k", lw=0.7)
        ax.set_title(f"Δ {title.format(w=window)}", fontsize=10)
        ax.set_ylabel(ylab, fontsize=9)
        ax.set_xlabel(f"{trait} class (1 = lowest)")
        ax.set_xticks(x)
        ax.tick_params(labelsize=8)
    h, l = axes[0].get_legend_handles_labels()
    fig.suptitle(f"Heatwave effect on day {day} by class of {trait} "
                 f"(paired difference per community, mean ± 95 % CI)", y=1.0)
    fig.tight_layout()
    fig.legend(h, l, loc="upper center", bbox_to_anchor=(0.5, 0.0), ncols=2, frameon=False)
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {out}")


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--snap-dir", required=True, help="export_taxon_snapshots.py output folder")
    ap.add_argument("--contrasts", nargs="+", default=DEFAULT_CONTRASTS, metavar="TREATMENT:REFERENCE",
                    help="scenario pairs (up to 4); default: the July wave by heatwave history")
    ap.add_argument("--trait", default="hsp_prod_induci", help="column of taxon_traits.csv")
    ap.add_argument("--classes", type=int, default=5, help="number of equal-sized trait classes")
    ap.add_argument("--at-day", type=int, default=2055,
                    help="day of the by-class figure (default 2055, last day of the July wave)")
    ap.add_argument("--days", type=int, nargs=2, metavar=("START", "END"), default=None,
                    help="restrict the by-day figures to this window (default: all snapshot days)")
    ap.add_argument("--heatwave", nargs="+", default=["june:2020-2026", "july:2048-2054", "august:2076-2082"],
                    metavar="NAME:START-END", help="heatwave windows to shade; 'none' for none")
    ap.add_argument("--out-dir", default=None, help="default <snap-dir>/figures")
    a = ap.parse_args()

    snap = Path(a.snap_dir)
    out_dir = Path(a.out_dir) if a.out_dir else snap / "figures"
    out_dir.mkdir(parents=True, exist_ok=True)
    meta = pd.read_csv(snap / "taxon_snapshots_meta.csv", index_col=0).iloc[:, 0]
    window = int(meta["window"])
    traits = pd.read_csv(snap / "taxon_traits.csv", index_col=0)
    if a.trait not in traits.columns:
        sys.exit(f"ERROR: trait {a.trait!r} not in taxon_traits.csv ({', '.join(traits.columns)})")
    trait = traits[a.trait]
    klass = (pd.qcut(trait, a.classes, labels=False, duplicates="drop") + 1).astype(int)
    edges = trait.groupby(klass).agg(["min", "max"])
    class_labels = [f"{k}: {r['min']:.3g}–{r['max']:.3g}" for k, r in edges.iterrows()]
    heatwaves = parse_heatwaves(a.heatwave)
    contrasts = [tuple(c.split(":")) for c in a.contrasts]
    print(f"{a.trait}: {len(trait)} taxa in {len(class_labels)} classes; "
          f"{len(contrasts)} contrasts; window {window} d")

    tables, all_stats, rows = {}, {}, []
    for treat, ref in contrasts:
        for s in (treat, ref):
            if s not in tables:
                print(f"  reading {s}", flush=True)
                tables[s] = class_table(snap, s, klass, trait, window, a.days)
        label = f"{pretty_label(treat)} − {pretty_label(ref)}"
        stats = contrast_stats(tables[treat], tables[ref])
        all_stats[label] = stats
        rows.append(stats.assign(contrast=label))
        fig_by_day(stats, label, class_labels, heatwaves, window, a.trait,
                   out_dir / f"by_day_{a.trait}_{treat}_minus_{ref}.png")
    if a.at_day in set(rows[0].day):
        fig_by_class(all_stats, a.at_day, class_labels, window, a.trait,
                     out_dir / f"by_class_{a.trait}_day{a.at_day}.png")
    else:
        print(f"  day {a.at_day} is not a snapshot day: by-class figure skipped")
    table = pd.concat(rows, ignore_index=True)
    table["class_range"] = table.klass.map(dict(enumerate(class_labels, start=1))).fillna("all taxa")
    table.to_csv(out_dir / f"taxon_contrasts_{a.trait}.csv", index=False)
    print(f"Wrote {out_dir / f'taxon_contrasts_{a.trait}.csv'}")


if __name__ == "__main__":
    main()
