#!/usr/bin/env python3
"""
heatwave_delta_figures.py
=========================
Treatment − baseline differences over time, one panel per tracked quantity,
from export_summary_data.py output (community-level files only).  Every panel
is scenario − baseline of the community means; because the mean is linear this
equals the mean of the per-community paired differences.

Figures (one per group; each panel one quantity; heatwave weeks shaded):
  delta1_pools.png        plant-litter substrate C (excl. dead microbes / enzymes),
                          total monomer C, microbial biomass C, living cells
  delta2_respiration.png  daily respiration, mass-specific respiration, and the
                          two accumulated over the plotted window
  delta3_efficiency.png   daily CUE and CAE (model values, floored)
  delta4_energy_<basis>.png     energy production, total energy cost, net energy gain
  delta5_carbon_<basis>.png     C uptake, C to energy, C to biosynthesis
  delta6_e_pathways_<basis>.png energy cost of maintenance, enzymes, HSP, osmolytes, transporters
  delta7_c_pathways_<basis>.png C to enzymes, HSP, osmolytes, endogenous respiration
  <basis> = total (grid sum, default), per_biomass (total / biomass C) or
  per_microbe (mean across taxa of the per-alive-microbe rate); --basis takes
  one or more.  These need an export made with the three-basis comm_* series.
  heatwave_deltas.csv     every Δ series (day × quantity × scenario)

Usage
-----
  python analysis/heatwave_delta_figures.py --summary-dir analysis_output/summary_20261026
      [--baseline climate_ambient_temp] [--year 6 | --days 2030 2110]
      [--heatwave winter:1860-1866 summer:2048-2054] [--smooth 7] [--fig-dir DIR]
      [--basis total per_biomass per_microbe]
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
from scenario_utils import scenario_selected, set_scenario_filter, pretty_label, scenario_color_map  # noqa: E402

plt.rcParams.update({"font.size": 9, "axes.titlesize": 10, "axes.labelsize": 9,
                     "legend.fontsize": 8, "figure.dpi": 150})


# ── Quantities ────────────────────────────────────────────────────────────────
# Each entry: (panel title, y label, function(series_getter) -> Series or None)
# where series_getter(stem) returns that stem's daily mean series for a scenario.

def _get(g, *stems):
    out = [g(s) for s in stems]
    return None if any(o is None for o in out) else out


def _sum(g, *stems):
    out = _get(g, *stems)
    return None if out is None else sum(out[1:], out[0])


GROUPS = [
    ("delta1_pools", "C pools", [
        ("Plant-litter substrate C", "C",
         lambda g: (lambda o: None if o is None else o[0] - o[1] - o[2])(
             _get(g, "system_substrate_C", "system_deadmic_substrate_C", "system_deadenz_substrate_C"))),
        ("Total monomer C", "C", lambda g: g("system_monomer_C")),
        ("Microbial biomass C", "C", lambda g: g("system_microbial_biomass_C")),
        ("Living cells", "cells", lambda g: g("system_taxon_count")),
    ]),
    ("delta2_respiration", "Respiration", [
        ("Daily respiration", "C d⁻¹", lambda g: g("system_resp")),
        ("Mass-specific respiration", "C C⁻¹ d⁻¹", lambda g: g("system_resp_per_biomass")),
        ("Cumulative respiration (from window start)", "C", lambda g: _cum(g("system_resp"))),
        ("Cumulative mass-specific respiration (from window start)", "C C⁻¹",
         lambda g: _cum(g("system_resp_per_biomass"))),
    ]),
    ("delta3_efficiency", "Efficiency", [
        ("Daily CUE (model, floored)", "CUE", lambda g: g("system_CUE")),
        ("Daily CAE", "CAE", lambda g: g("system_CAE")),
    ]),
    ("delta4_energy", "Energy budget", [
        ("Energy production", "E d⁻¹", lambda g: g("comm_production_E")),
        ("Total energy cost", "E d⁻¹", lambda g: _sum(g, *[f"comm_{p}_E" for p in
                                                           ("transport", "maintenance", "hsp", "enzyme", "osmolyte")])),
        ("Net energy gain (production − cost)", "E d⁻¹",
         lambda g: (lambda p, c: None if (p is None or c is None) else p - c)(
             g("comm_production_E"),
             _sum(g, *[f"comm_{p}_E" for p in ("transport", "maintenance", "hsp", "enzyme", "osmolyte")]))),
    ]),
    ("delta5_carbon", "C allocation", [
        ("Total C uptake", "C d⁻¹", lambda g: g("comm_uptake_C")),
        ("C to energy (uptake − biosynthesis)", "C d⁻¹",
         lambda g: (lambda o: None if o is None else o[0] - o[1])(
             _get(g, "comm_uptake_C", "comm_biosynthesis_C"))),
        ("C to biosynthesis", "C d⁻¹", lambda g: g("comm_biosynthesis_C")),
    ]),
    ("delta6_e_pathways", "Energy cost by pathway", [
        ("Maintenance", "E d⁻¹", lambda g: g("comm_maintenance_E")),
        ("Enzymes", "E d⁻¹", lambda g: g("comm_enzyme_E")),
        ("HSP", "E d⁻¹", lambda g: g("comm_hsp_E")),
        ("Osmolytes", "E d⁻¹", lambda g: g("comm_osmolyte_E")),
        ("Transporters", "E d⁻¹", lambda g: g("comm_transport_E")),
    ]),
    ("delta7_c_pathways", "C cost by pathway", [
        ("Enzymes", "C d⁻¹", lambda g: g("comm_enzyme_C")),
        ("HSP", "C d⁻¹", lambda g: g("comm_hsp_C")),
        ("Osmolytes", "C d⁻¹", lambda g: g("comm_osmolyte_C")),
        ("Endogenous respiration", "C d⁻¹", lambda g: g("comm_endo_resp_C")),
    ]),
]


def _cum(s):
    return None if s is None else s.cumsum()


BASIS_SUFFIX = {"per_microbe": "per_microbe", "total": "total", "per_biomass": "per_biomass"}
BASIS_UNIT = {"per_microbe": " per alive microbe", "total": " (grid total)", "per_biomass": " per biomass C"}
BASIS_TITLE = {"per_microbe": "per alive microbe", "total": "grid total", "per_biomass": "per unit biomass C"}


# ── I/O ───────────────────────────────────────────────────────────────────────

def read_mean(scen_dir: Path, stem: str, lo: int, hi: int):
    f = scen_dir / f"{stem}_mean.csv"
    if not f.exists():
        return None
    s = pd.read_csv(f, index_col=0).iloc[:, 0]
    s.index = s.index.astype(int)
    return s.loc[lo:hi]


def read_meta(scen_dir: Path) -> dict:
    m = pd.read_csv(scen_dir / "summary_meta.csv", index_col=0).iloc[:, 0]
    return {k: int(float(m[k])) for k in ("pulse", "end_time", "gridsize")}


def parse_heatwaves(items):
    if len(items) == 1 and items[0].lower() == "none":    # e.g. year-round warming: nothing to shade
        return []
    out = []
    for it in items:
        name, span = it.split(":")
        a, b = span.split("-")
        out.append((name, int(a), int(b)))
    return out


# ── Plotting ──────────────────────────────────────────────────────────────────

def make_group(root, fname, gtitle, quantities, scen, baseline, lo, hi, heatwaves, colors, smooth, fig_dir, rows,
               basis="total"):
    def _stem(stem):                       # community series come in three bases
        return f"{stem}_{BASIS_SUFFIX[basis]}" if stem.startswith("comm_") else stem
    getters = {s: (lambda stem, s=s: read_mean(root / s, _stem(stem), lo, hi)) for s in scen}
    is_comm = fname in ("delta4_energy", "delta5_carbon", "delta6_e_pathways", "delta7_c_pathways")
    n = len(quantities)
    ncol = 2 if n <= 4 else 3
    nrow = int(np.ceil(n / ncol))
    # Skip the figure outright when the export lacks the series for this basis
    # (comm_<x>_total / _per_biomass need an export made after 2026-10-29).
    available = [q for q in quantities if q[2](getters[baseline]) is not None]
    if not available:
        print(f"Skipped {fname} ({basis}): none of its series are in the export"
              + (" — re-run export_summary_data.py for the total / per_biomass series" if is_comm else ""))
        return
    missing = [q[0] for q in quantities if q not in available]
    if missing:
        print(f"  {fname} ({basis}): not in export, panels left out: {', '.join(missing)}")
    fig, axes = plt.subplots(nrow, ncol, figsize=(5.6 * ncol, 3.3 * nrow), squeeze=False)
    treatments = [s for s in scen if s != baseline]
    for ax, (title, ylab, fn) in zip(axes.ravel(), quantities):
        base = fn(getters[baseline])
        if base is None:
            ax.set_visible(False)
            continue
        for s in treatments:
            other = fn(getters[s])
            if other is None:
                continue
            d = (other - base.reindex(other.index))
            if smooth and smooth > 1:
                d = d.rolling(smooth, center=True, min_periods=1).mean()
            ax.plot(d.index, d.values, color=colors[s], lw=1.2, label=pretty_label(s))
            rows.append(pd.DataFrame({"day": d.index, "quantity": title, "group": gtitle,
                                      "basis": basis if is_comm else "total",
                                      "scenario": s, "delta": d.values}))
        for _, a, b in heatwaves:
            if b >= lo and a <= hi:
                ax.axvspan(a, b, color="0.85", zorder=0)
        ax.axhline(0, color="k", lw=0.7, ls="--")
        ax.set_title(f"Δ {title}")
        ax.set_ylabel(ylab + (BASIS_UNIT[basis] if is_comm else ""))
        ax.set_xlabel("day")
        ax.set_xlim(lo, hi)
    for ax in axes.ravel()[n:]:
        ax.set_visible(False)
    axes.ravel()[0].legend(frameon=False)
    fig.suptitle(f"{gtitle}{', ' + BASIS_TITLE[basis] if is_comm else ''}: scenario − {pretty_label(baseline)} "
                 f"(difference of community means), days {lo}–{hi}"
                 + (f", {smooth}-day running mean" if smooth and smooth > 1 else ""), y=1.01)
    fig.tight_layout()
    out = fig_dir / (f"{fname}_{basis}.png" if is_comm else f"{fname}.png")
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {out}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--summary-dir", required=True)
    ap.add_argument("--scenarios", nargs="+", default=None, metavar="NAME",
                    help="plot only these scenario folders (stems or display labels; "
                         "include the baseline explicitly)")
    ap.add_argument("--baseline", default="climate_ambient_temp")
    ap.add_argument("--year", type=int, default=None, help="year to plot (default: last pulse)")
    ap.add_argument("--days", type=int, nargs=2, metavar=("START", "END"), default=None,
                    help="absolute day window (overrides --year)")
    ap.add_argument("--heatwave", nargs="+", default=["winter:1860-1866", "summer:2048-2054"],
                    # summer heatwave-timing runs: --heatwave june:2020-2026 july:2048-2054 august:2076-2082
                    metavar="NAME:START-END",
                    help="heatwave windows to shade (absolute days); 'none' for no shading")
    ap.add_argument("--smooth", type=int, default=0, help="running-mean window in days (0 = none)")
    ap.add_argument("--basis", choices=["total", "per_biomass", "per_microbe"], nargs="+", default=["total"],
                    help="basis for the per-community energy / C / pathway series (several allowed; "
                         "pools, respiration and efficiency panels are unaffected)")
    ap.add_argument("--fig-dir", default=None, help="default <summary-dir>/figures/deltas")
    a = ap.parse_args()
    set_scenario_filter(a.scenarios)

    root = Path(a.summary_dir)
    scen = sorted(d.name for d in root.iterdir() if d.is_dir() and (d / "summary_meta.csv").exists()
                  and scenario_selected(d.name))
    if a.baseline not in scen:
        sys.exit(f"ERROR: baseline {a.baseline!r} not among scenarios {scen}")
    scen = [a.baseline] + [s for s in scen if s != a.baseline]
    meta = read_meta(root / a.baseline)
    if a.days:
        lo, hi = a.days
    else:
        year = a.year or meta["pulse"]
        lo, hi = (year - 1) * meta["end_time"] + 1, year * meta["end_time"]
    heatwaves = parse_heatwaves(a.heatwave)
    fig_dir = Path(a.fig_dir) if a.fig_dir else root / "figures" / "deltas"
    fig_dir.mkdir(parents=True, exist_ok=True)
    labels = [pretty_label(s) for s in scen]
    cmap = scenario_color_map(labels)
    colors = {s: cmap[l] for s, l in zip(scen, labels)}
    print(f"Scenarios: {labels}; days {lo}–{hi}; gridsize {meta['gridsize']}")

    rows = []
    for fname, gtitle, quantities in GROUPS:
        comm = fname in ("delta4_energy", "delta5_carbon", "delta6_e_pathways", "delta7_c_pathways")
        for basis in (a.basis if comm else ["total"]):
            make_group(root, fname, gtitle, quantities, scen, a.baseline, lo, hi, heatwaves, colors, a.smooth,
                       fig_dir, rows, basis=basis)
    if rows:
        pd.concat(rows, ignore_index=True).to_csv(fig_dir / "heatwave_deltas.csv", index=False)
        print(f"Wrote {fig_dir / 'heatwave_deltas.csv'}")


if __name__ == "__main__":
    main()
