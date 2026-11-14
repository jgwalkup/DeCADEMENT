#!/usr/bin/env python3
"""
taxon_trait_figures.py
======================
Heatwave effects taxon by taxon, against a trait, from the per-taxon snapshot
tables written by export_taxon_snapshots.py.  One contrast (treatment:reference,
two scenarios that differ by one heatwave) per run.

Each taxon is one point or line.  Its values pool the communities it occurs in:
sums over those communities for pools and fluxes (so rates per biomass C are
ratios of sums, not means of noisy per-community ratios), and the share of
those communities in which it still has living cells.

Figures  (<out-dir>, default <snap-dir>/figures/taxa)
-------
  taxon_deltas_<contrast>.png      treatment − reference against day, one line
                                   per taxon, coloured by the trait
  taxon_deltas_biomass_cells_<contrast>.png   the biomass C and living cells
                                   panels of that figure on their own
  taxon_at_day<D>_<contrast>.png   treatment − reference on day D against the
                                   trait, with the rank correlation; the y axes
                                   scale to the spread among taxa, so small
                                   trait effects on a large common response show
  potential_vs_realized_<contrast>.png   the trait (potential) against realized
                                   HSP production and HSP pool share, in both
                                   scenarios, on the first and the last day of
                                   the heatwave
  maintenance_<contrast>.png       on day D: maintenance cost, its share of
                                   energy production, growth yield (share of
                                   the day's C uptake that becomes biomass) and
                                   net specific growth rate (C retained per
                                   biomass C and day, over the export's window)
                                   against the realized HSP pool share
  taxon_table_<contrast>.csv       per taxon and day: both scenarios and the
                                   difference, for every quantity

Usage
-----
  python analysis/taxon_trait_figures.py \\
      --snap-dir analysis_output/taxon_snapshots_summer_heatwaves_20261110
      [--contrast climate_July_heatwave:climate_ambient_temp]
      [--first-day 2049] [--at-day 2055] [--trait hsp_prod_induci]
      [--heatwave june:2020-2026 july:2048-2054 august:2076-2082] [--out-dir DIR]
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
from heatwave_delta_figures import parse_heatwaves              # noqa: E402
from scenario_utils import pretty_label                         # noqa: E402
from taxon_snapshot_figures import read_snapshots               # noqa: E402

SUMMED = ["cells", "biomass_C", "hsp_pool_C", "hsp_C", "hsp_E", "maintenance_E", "production_E",
          "endo_resp_C", "uptake_C", "biosynthesis_C", "enzyme_C", "osmolyte_C", "osmo_pool_C",
          "transport_E", "enzyme_E", "osmolyte_E"]
CMAP = "viridis"

# (key, title, unit, difference shown as % of the reference)
QUANTITIES = [
    ("biomass_C",   "Biomass C",                 "% change from {ref}",  True),
    ("cells",       "Living cells",              "% change from {ref}",  True),
    ("alive_share", "Communities where alive",   "share of communities", False),
    ("hsp_share",   "HSP pool share of biomass", "C C⁻¹",                False),
    ("hsp_C_pb",    "HSP production",            "C d⁻¹ per biomass C",  False),
    ("maint_E_pb",  "Maintenance energy cost",   "E d⁻¹ per biomass C",  False),
    ("endo_pb",     "Endogenous respiration",    "C d⁻¹ per biomass C",  False),
    ("cue",         "CUE ({w}-day)",             "CUE",                  False),
]


def taxon_table(snap: Path, scenario: str, window: int) -> pd.DataFrame:
    """Quantities per taxon x day for one scenario, pooled over communities."""
    w = [f"uptake_C_{window}d", f"retained_C_{window}d"]
    df = read_snapshots(snap, scenario, ["community", "taxon", "day"] + SUMMED + w, None)
    df["alive"] = df.cells > 0
    g = df.groupby(["taxon", "day"], sort=False)
    s = g[SUMMED + w].sum()
    with np.errstate(divide="ignore", invalid="ignore"):
        out = s[["biomass_C", "cells"]].astype(float)
        out["alive_share"] = g.alive.mean()
        out["n_communities"] = g.size()
        out["hsp_share"] = s.hsp_pool_C / s.biomass_C
        out["hsp_C_pb"] = s.hsp_C / s.biomass_C
        out["hsp_E_pb"] = s.hsp_E / s.biomass_C
        out["maint_E_pb"] = s.maintenance_E / s.biomass_C
        out["maint_share_E"] = s.maintenance_E / s.production_E
        out["endo_pb"] = s.endo_resp_C / s.biomass_C
        # realized production and energy budget terms, for the per-taxon table
        out["uptake_C_pb"] = s.uptake_C / s.biomass_C
        out["enz_C_pb"] = s.enzyme_C / s.biomass_C
        out["osmo_C_pb"] = s.osmolyte_C / s.biomass_C
        out["osmo_share"] = s.osmo_pool_C / s.biomass_C
        out["prod_E_pb"] = s.production_E / s.biomass_C
        out["transport_E_pb"] = s.transport_E / s.biomass_C
        out["enz_E_pb"] = s.enzyme_E / s.biomass_C
        out["osmo_E_pb"] = s.osmolyte_E / s.biomass_C
        out["cue"] = s[w[1]] / s[w[0]]
        # share of the day's uptake that becomes biomass (0-1), and net change of
        # biomass C per unit biomass and day over the window (negative = shrinking)
        out["growth_yield"] = (s.biosynthesis_C - s.enzyme_C) / s.uptake_C
        out["energy_frac"] = (s.uptake_C - s.biosynthesis_C) / s.uptake_C   # share of uptake C sent to energy
        out["net_growth"] = s[w[1]] / window / s.biomass_C
    return out.replace([np.inf, -np.inf], np.nan)


def _rho(x: pd.Series, y: pd.Series) -> str:
    ok = x.notna() & y.notna()
    if ok.sum() < 5:
        return ""
    return f"  (Spearman ρ = {x[ok].corr(y[ok], method='spearman'):+.2f})"


def _colorbar(fig, axes, norm, trait):
    sm = plt.cm.ScalarMappable(norm=norm, cmap=CMAP)
    fig.colorbar(sm, ax=list(axes), label=trait, shrink=0.8, pad=0.015)


def fig_deltas(delta, trait, norm, heatwaves, window, label, ref, out, keys=None):
    """Treatment − reference against day, one line per taxon; `keys` picks the
    quantities (default all)."""
    quantities = [q for q in QUANTITIES if keys is None or q[0] in keys]
    ncol = min(len(quantities), 4)
    nrow = -(-len(quantities) // ncol)
    days = np.sort(delta.index.get_level_values("day").unique())
    cmap = plt.get_cmap(CMAP)
    order = trait.sort_values().index                    # draw low-trait taxa first
    fig, axes = plt.subplots(nrow, ncol, figsize=(4.75 * ncol, 3.6 * nrow), squeeze=False,
                             layout="constrained")
    for ax, (key, title, unit, _) in zip(axes.ravel(), quantities):
        wide = delta[key].unstack("day").reindex(columns=days)
        for t in order:
            if t in wide.index:
                ax.plot(days, wide.loc[t].to_numpy(), color=cmap(norm(trait[t])), lw=0.6, alpha=0.7)
        for _, a, b in heatwaves:                        # climate rows a..b = output days a+1..b+1
            ax.axvspan(a, b + 1, color="0.88", zorder=0)
        ax.axhline(0, color="k", lw=0.7)
        ax.set_title(f"Δ {title.format(w=window)}", fontsize=10)
        ax.set_ylabel(unit.format(ref=ref), fontsize=9)
        ax.set_xlabel("day")
        ax.set_xlim(days[0], days[-1])
    _colorbar(fig, axes.ravel(), norm, trait.name)
    fig.suptitle(f"{label}: one line per taxon (pooled over its communities)")
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"Saved {out}")


def fig_at_day(delta, trait, norm, day, window, label, ref, out):
    d = delta.xs(day, level="day")
    x = trait.reindex(d.index)
    fig, axes = plt.subplots(2, 4, figsize=(19, 7.2), squeeze=False, layout="constrained")
    for ax, (key, title, unit, _) in zip(axes.ravel(), QUANTITIES):
        ax.scatter(x, d[key], c=x, cmap=CMAP, norm=norm, s=14)
        ok = x.notna() & d[key].notna()
        if ok.sum() > 2:
            b, a = np.polyfit(x[ok], d[key][ok], 1)
            xs = np.array([x.min(), x.max()])
            ax.plot(xs, a + b * xs, color="k", lw=1)
        ax.set_title(f"Δ {title.format(w=window)}{_rho(x, d[key])}", fontsize=10)
        ax.set_ylabel(unit.format(ref=ref), fontsize=9)
        ax.set_xlabel(trait.name)
    fig.suptitle(f"{label} on day {day}: one point per taxon, against {trait.name}")
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"Saved {out}")


def fig_potential(tables, names, trait, norm, days, label, out):
    """Trait against realized HSP production and pool share, per scenario and day."""
    rows = [("hsp_C_pb", "Realized HSP production (C d⁻¹ per biomass C)"),
            ("hsp_share", "Realized HSP pool share of biomass")]
    cols = [(s, d) for d in days for s in names]
    fig, axes = plt.subplots(len(rows), len(cols), figsize=(4.6 * len(cols), 7.2), squeeze=False,
                             layout="constrained", sharey="row")
    for i, (key, ylab) in enumerate(rows):
        for j, (s, d) in enumerate(cols):
            ax = axes[i, j]
            y = tables[s][key].xs(d, level="day")
            x = trait.reindex(y.index)
            ax.scatter(x, y, c=x, cmap=CMAP, norm=norm, s=14)
            ax.set_title(f"{pretty_label(s)}, day {d}{_rho(x, y)}", fontsize=9)
            ax.set_xlabel(f"{trait.name} (potential)")
            if j == 0:
                ax.set_ylabel(ylab, fontsize=9)
    fig.suptitle(f"Potential against realized HSP production — {label}: one point per taxon")
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"Saved {out}")


def fig_maintenance(treat, ref, trait, norm, day, window, label, out):
    t, r = treat.xs(day, level="day"), ref.xs(day, level="day").reindex(treat.xs(day, level="day").index)
    c = trait.reindex(t.index)
    panels = [("maint_E_pb", "Maintenance energy cost (E d⁻¹ per biomass C)"),
              ("maint_share_E", "Maintenance cost / energy production"),
              ("growth_yield", "Growth yield: (biosynthesis C − enzyme C) / C uptake, that day"),
              ("net_growth", f"Net specific growth rate (d⁻¹, {window}-day mean)")]
    fig, axes = plt.subplots(1, 4, figsize=(22, 4.8), squeeze=False, layout="constrained")
    for ax, (key, ylab) in zip(axes.ravel(), panels):
        ax.scatter(r.hsp_share, r[key], color="0.6", s=12, marker="x", label="reference scenario")
        ax.scatter(t.hsp_share, t[key], c=c, cmap=CMAP, norm=norm, s=16, label="treatment scenario")
        ax.set_xlabel("Realized HSP pool share of biomass")
        ax.set_ylabel(ylab, fontsize=9)
        ax.set_title(f"treatment{_rho(t.hsp_share, t[key])}", fontsize=10)
        if key == "net_growth":
            ax.axhline(0, color="k", lw=0.7)
    axes[0, 0].legend(frameon=False, fontsize=8)
    _colorbar(fig, axes.ravel(), norm, trait.name)
    fig.suptitle(f"{label} on day {day}: maintenance cost and efficiency against realized HSP, one point per taxon")
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"Saved {out}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--snap-dir", required=True, help="export_taxon_snapshots.py output folder")
    ap.add_argument("--contrast", default="climate_July_heatwave:climate_ambient_temp",
                    metavar="TREATMENT:REFERENCE")
    ap.add_argument("--trait", default="hsp_prod_induci", help="column of taxon_traits.csv")
    ap.add_argument("--first-day", type=int, default=2049, help="first day of the heatwave (default 2049, July)")
    ap.add_argument("--at-day", type=int, default=2055, help="last day of the heatwave (default 2055, July)")
    ap.add_argument("--heatwave", nargs="+", default=["june:2020-2026", "july:2048-2054", "august:2076-2082"],
                    metavar="NAME:START-END", help="heatwave windows to shade; 'none' for none")
    ap.add_argument("--out-dir", default=None, help="default <snap-dir>/figures/taxa")
    a = ap.parse_args()

    snap = Path(a.snap_dir)
    out_dir = Path(a.out_dir) if a.out_dir else snap / "figures" / "taxa"
    out_dir.mkdir(parents=True, exist_ok=True)
    window = int(pd.read_csv(snap / "taxon_snapshots_meta.csv", index_col=0).iloc[:, 0]["window"])
    traits = pd.read_csv(snap / "taxon_traits.csv", index_col=0)
    if a.trait not in traits.columns:
        sys.exit(f"ERROR: trait {a.trait!r} not in taxon_traits.csv ({', '.join(traits.columns)})")
    trait = traits[a.trait]
    norm = plt.Normalize(trait.min(), trait.max())
    treat_name, ref_name = a.contrast.split(":")
    ref_label = pretty_label(ref_name)
    label = f"{pretty_label(treat_name)} − {ref_label}"
    tag = f"{treat_name}_minus_{ref_name}"

    tables = {}
    for s in (treat_name, ref_name):
        print(f"  reading {s}", flush=True)
        tables[s] = taxon_table(snap, s, window)
    treat, ref = tables[treat_name], tables[ref_name].reindex(tables[treat_name].index)
    keys = [q[0] for q in QUANTITIES]
    with np.errstate(divide="ignore", invalid="ignore"):
        delta = pd.DataFrame({k: ((treat[k] - ref[k]) / ref[k] * 100 if rel else treat[k] - ref[k])
                              for k, _, _, rel in QUANTITIES}).replace([np.inf, -np.inf], np.nan)
    for d in (a.first_day, a.at_day):
        if d not in delta.index.get_level_values("day"):
            sys.exit(f"ERROR: day {d} is not a snapshot day")

    heatwaves = parse_heatwaves(a.heatwave)
    fig_deltas(delta, trait, norm, heatwaves, window, label, ref_label, out_dir / f"taxon_deltas_{tag}.png")
    fig_deltas(delta, trait, norm, heatwaves, window, label, ref_label,
               out_dir / f"taxon_deltas_biomass_cells_{tag}.png",
               keys=("biomass_C", "cells"))
    fig_at_day(delta, trait, norm, a.at_day, window, label, ref_label,
               out_dir / f"taxon_at_day{a.at_day}_{tag}.png")
    fig_potential(tables, [ref_name, treat_name], trait, norm, [a.first_day, a.at_day], label,
                  out_dir / f"potential_vs_realized_{tag}.png")
    fig_maintenance(treat, ref, trait, norm, a.at_day, window, label, out_dir / f"maintenance_{tag}.png")

    extra = ["hsp_E_pb", "maint_share_E", "growth_yield", "net_growth", "energy_frac", "uptake_C_pb", "enz_C_pb",
             "osmo_C_pb", "osmo_share", "prod_E_pb", "transport_E_pb", "enz_E_pb", "osmo_E_pb",
             "n_communities"]
    table = pd.concat([treat[keys + extra].add_suffix("_treatment"), ref[keys + extra].add_suffix("_reference"),
                       delta.add_suffix("_delta")], axis=1).reset_index()
    table.insert(1, a.trait, table.taxon.map(trait))
    table.to_csv(out_dir / f"taxon_table_{tag}.csv", index=False)
    print(f"Wrote {out_dir / f'taxon_table_{tag}.csv'}")


if __name__ == "__main__":
    main()
