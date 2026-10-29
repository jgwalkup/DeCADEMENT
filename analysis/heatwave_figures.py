#!/usr/bin/env python3
"""
heatwave_figures.py
===================
Figures for the heatwave experiment from export_summary_data.py output
(community-level files only) and heatwave_compare.py output.

  fig0_all_years.png        daily series over the whole run, mean ± SD across
                            communities, by scenario (years before the branch
                            are shared, so the scenarios coincide there)
  fig1_year_series.png      branch-year daily series by scenario, mean ± SD:
                            daily CUE (floored), respiration, microbial biomass
                            C, endogenous respiration per biomass C, unpaid
                            energy, living cells; heatwave windows shaded
  fig2_year_differences.png the same series as scenario − baseline (difference
                            of community means = mean of paired differences)
  fig3_net_cue_pulses.png   net CUE per pulse by scenario (mean ± SE across
                            communities): spin-up years, then the branch year
  fig4_paired_differences.png
                            per-community paired differences for the branch
                            year (heatwave_compare_pairs.csv): net CUE,
                            retained C, endogenous respiration, biomass C at
                            year end; mean and 95 % CI from heatwave_compare.csv

Input
-----
  --summary-dir DIR  export root: <scenario>/system_*, comm_*, net_cue_pulse_*,
                     summary_meta.csv, plus heatwave_compare*.csv at the root
                     (run heatwave_compare.py --summary-dir DIR first; fig4 is
                     skipped without it).  Per-taxon tables are not needed.
  --baseline NAME    baseline scenario folder (default climate_ambient_temp)
  --heatwave NAME:START-END ...
                     heatwave windows to shade, as absolute days (default: the
                     Oct 2026 design, winter 1860-1866 and summer 2048-2054)
  --fig-dir DIR      default <summary-dir>/figures

Usage
-----
  python analysis/heatwave_figures.py --summary-dir analysis_output/summary_20261026
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

plt.rcParams.update({"font.size": 9, "axes.titlesize": 10, "axes.labelsize": 9,
                     "legend.fontsize": 8, "figure.dpi": 150})

COLORS = {"climate_ambient_temp": "#000000", "climate_summer_warm": "#d62728",
          "climate_winter_warm": "#1f77b4"}   # baseline black, dashed, on top (see scenario_utils)
_FALLBACK = ["#2ca02c", "#9467bd", "#ff7f0e", "#8c564b"]

SERIES = [  # (stem, title, ylabel)
    ("system_CUE",                  "Daily CUE (model, floored)",            "CUE"),
    ("system_resp",                 "Respiration",                           "C d⁻¹"),
    ("system_microbial_biomass_C",  "Microbial biomass C",                   "C"),
    ("comm_endo_resp_C_per_biomass", "Endogenous respiration per biomass C",  "C C⁻¹ d⁻¹"),
    ("comm_unpaid_E_per_microbe",   "Unpaid energy demand per microbe",      "E d⁻¹"),
    ("system_taxon_count",          "Living cells",                          "cells"),
]
ALL_YEARS_SERIES = [
    ("system_microbial_biomass_C",  "Microbial biomass C",                   "C"),
    ("system_resp",                 "Respiration",                           "C d⁻¹"),
    ("system_substrate_C",          "Substrate C",                           "C"),
    ("system_CUE",                  "Daily CUE (model, floored)",            "CUE"),
    ("system_microbial_energy",     "Microbial energy",                      "E"),
    ("comm_production_E_per_microbe", "Energy production per microbe",       "E d⁻¹"),
    ("comm_endo_resp_C_per_biomass", "Endogenous respiration per biomass C",  "C C⁻¹ d⁻¹"),
    ("comm_retained_C_per_biomass", "Retained C per biomass C",              "C C⁻¹ d⁻¹"),
    ("system_taxon_count",          "Living cells",                          "cells"),
]
PAIR_METRICS = [("net_cue", "Net CUE, branch year"), ("c_retained", "C retained"),
                ("c_endo_resp", "Endogenous respiration C"), ("biomass_C_end", "Biomass C, year end")]


sys.path.insert(0, str(Path(__file__).resolve().parent))
from scenario_utils import pretty_label as label, line_style, plot_order, legend_baseline_first  # noqa: E402


def color(s: str, i: int) -> str:
    return COLORS.get(s, _FALLBACK[i % len(_FALLBACK)])


def read_series(scen_dir: Path, stem: str):
    f = scen_dir / f"{stem}_mean.csv"
    if not f.exists():
        return None, None
    def _idx(x):
        try:
            x.index = x.index.astype(int)          # day labels; pulse labels stay as they are
        except (ValueError, TypeError):
            pass
        return x
    m = _idx(pd.read_csv(f, index_col=0).iloc[:, 0])
    sf = scen_dir / f"{stem}_sd.csv"
    s = _idx(pd.read_csv(sf, index_col=0).iloc[:, 0]) if sf.exists() else None
    return m, s


def read_meta(scen_dir: Path) -> dict:
    m = pd.read_csv(scen_dir / "summary_meta.csv", index_col=0).iloc[:, 0]
    return {k: int(float(m[k])) for k in ("pulse", "end_time", "gridsize")}


def parse_heatwaves(items) -> list[tuple[str, int, int]]:
    out = []
    for it in items:
        name, span = it.split(":")
        a, b = span.split("-")
        out.append((name, int(a), int(b)))
    return out


def shade(ax, heatwaves, lo, hi):
    for name, a, b in heatwaves:
        if b >= lo and a <= hi:
            ax.axvspan(a, b, color="0.85", zorder=0)


def fig_year_series(scen, root, year, dpp, heatwaves, out, diff_against=None,
                    series=SERIES, lo=None, hi=None, ribbon=True, n_pulses=None, title=None):
    if lo is None:
        lo, hi = (year - 1) * dpp + 1, year * dpp
    ncol = 3
    nrow = int(np.ceil(len(series) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(12, 3.1 * nrow))
    axes = np.atleast_2d(axes)
    for ax, (stem, ttl, ylab) in zip(axes.ravel(), series):
        base_m = read_series(root / diff_against, stem)[0] if diff_against else None
        plotted = False
        for s in plot_order(scen):                     # baseline last, on top
            i = scen.index(s)
            if diff_against and s == diff_against:
                continue
            m, sd = read_series(root / s, stem)
            if m is None:
                continue
            m = m.loc[lo:hi]
            if diff_against:
                m = m - base_m.reindex(m.index)
            elif ribbon and sd is not None:
                sd = sd.reindex(m.index)
                ax.fill_between(m.index, (m - sd).values, (m + sd).values, color=color(s, i), alpha=0.15, lw=0)
            ax.plot(m.index, m.values, color=color(s, i), lw=1.1, label=label(s), **line_style(s, scen))
            plotted = True
        if not plotted:
            ax.set_visible(False)
            continue
        shade(ax, heatwaves, lo, hi)
        if diff_against:
            ax.axhline(0, color="k", lw=0.6, ls="--")
            ax.set_title(f"Δ {ttl}")
        else:
            ax.set_title(ttl)
        if n_pulses:                                   # year boundaries over a multi-year span
            for yb in range(1, n_pulses):
                ax.axvline(yb * dpp + 0.5, color="0.75", lw=0.6, ls=":")
        ax.set_ylabel(ylab)
        ax.set_xlabel("day" if n_pulses else f"day (year {year})")
        ax.set_xlim(lo, hi)
    for ax in axes.ravel()[len(series):]:
        ax.set_visible(False)
    legend_baseline_first(axes[0, 0], [label(s) for s in scen], frameon=False)
    if diff_against:
        fig.suptitle(f"Scenario − {label(diff_against)} (difference of community means), year {year}", y=1.01)
    elif title:
        fig.suptitle(title, y=1.01)
    fig.tight_layout()
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {out}")


def fig_net_cue_pulses(scen, root, n_pulses, out):
    fig, ax = plt.subplots(figsize=(5.2, 3.6))
    x = np.arange(1, n_pulses + 1)
    for i, s in enumerate(scen):
        m, sd = read_series(root / s, "net_cue_pulse")
        if m is None:
            continue
        m.index = [int(str(p).split()[-1]) for p in m.index]
        n = pd.read_csv(root / "communities_included.csv").shape[0] if (root / "communities_included.csv").exists() else np.nan
        se = (sd.values / np.sqrt(n)) if (sd is not None and np.isfinite(n)) else None
        off = (i - (len(scen) - 1) / 2) * 0.08
        ax.errorbar(m.index + off, m.values, yerr=se, fmt="o", ms=3.5, lw=1, capsize=2,
                    color=color(s, i), label=label(s), zorder=5 if s == scen[0] else 3)
        ax.plot(m.index + off, m.values, color=color(s, i), lw=1, **line_style(s, scen))
    ax.axvline(n_pulses - 0.5, color="0.6", lw=0.8, ls=":")
    ax.set_xticks(x)
    ax.set_xlabel("pulse (year)")
    ax.set_ylabel("net CUE  (Σ retained / Σ uptake)")
    ax.set_title("Net CUE per year, mean ± SE across communities")
    legend_baseline_first(ax, [label(s) for s in scen], frameon=False)
    fig.tight_layout()
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {out}")


def fig_paired(root, baseline, out):
    pf, tf = root / "heatwave_compare_pairs.csv", root / "heatwave_compare.csv"
    if not (pf.exists() and tf.exists()):
        print("heatwave_compare*.csv not found — run heatwave_compare.py first; fig4 skipped")
        return
    pairs, table = pd.read_csv(pf), pd.read_csv(tf)
    scen = list(pairs.scenario.unique())
    fig, axes = plt.subplots(1, len(PAIR_METRICS), figsize=(3.4 * len(PAIR_METRICS), 3.6))
    for ax, (m, title) in zip(axes, PAIR_METRICS):
        sub = pairs[pairs.metric == m]
        vals = sub["diff"].to_numpy(float)
        vals = vals[np.isfinite(vals)]
        if vals.size == 0:
            ax.set_visible(False)
            continue
        q = np.nanpercentile(vals, [0.5, 99.5])
        bins = np.linspace(q[0], q[1], 51)
        for i, s in enumerate(scen):
            d = sub[sub.scenario == s]["diff"].to_numpy(float)
            d = d[np.isfinite(d)]
            c = color(s, i)
            ax.hist(np.clip(d, q[0], q[1]), bins=bins, alpha=0.5, color=c, label=label(s))
            r = table[(table.scenario == s) & (table.metric == m)]
            if len(r):
                r = r.iloc[0]
                ax.axvline(r.diff_mean, color=c, lw=1.6)
                ax.axvspan(r.ci95_lo, r.ci95_hi, color=c, alpha=0.25)
        ax.axvline(0, color="k", lw=0.7, ls="--")
        ax.set_title(title)
        ax.set_xlabel(f"scenario − {label(baseline)}, per community")
    axes[0].set_ylabel("communities")
    axes[0].legend(frameon=False)
    fig.suptitle("Paired branch-year differences (mean and 95 % CI shaded)", y=1.03)
    fig.tight_layout()
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {out}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--summary-dir", required=True)
    ap.add_argument("--baseline", default="climate_ambient_temp")
    ap.add_argument("--heatwave", nargs="+", default=["winter:1860-1866", "summer:2048-2054"],
                    metavar="NAME:START-END")
    ap.add_argument("--year", type=int, default=None, help="year to plot (default: last pulse)")
    ap.add_argument("--fig-dir", default=None)
    a = ap.parse_args()

    root = Path(a.summary_dir)
    scen = sorted(d.name for d in root.iterdir() if d.is_dir() and (d / "summary_meta.csv").exists())
    if a.baseline not in scen:
        sys.exit(f"ERROR: baseline {a.baseline!r} not among scenarios {scen}")
    scen = [a.baseline] + [s for s in scen if s != a.baseline]
    meta = read_meta(root / a.baseline)
    year = a.year or meta["pulse"]
    heatwaves = parse_heatwaves(a.heatwave)
    fig_dir = Path(a.fig_dir) if a.fig_dir else root / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)
    print(f"Scenarios: {scen}; year {year} of {meta['pulse']}, {meta['end_time']} days/year, "
          f"gridsize {meta['gridsize']}")

    fig_year_series(scen, root, year, meta["end_time"], heatwaves, fig_dir / "fig0_all_years.png",
                    series=ALL_YEARS_SERIES, lo=1, hi=meta["pulse"] * meta["end_time"],
                    n_pulses=meta["pulse"], title="Daily series, mean ± SD across communities, all years")
    fig_year_series(scen, root, year, meta["end_time"], heatwaves, fig_dir / "fig1_year_series.png",
                    title=f"Year {year}, mean ± SD across communities")
    fig_year_series(scen, root, year, meta["end_time"], heatwaves, fig_dir / "fig2_year_differences.png",
                    diff_against=a.baseline)
    fig_net_cue_pulses(scen, root, meta["pulse"], fig_dir / "fig3_net_cue_pulses.png")
    fig_paired(root, a.baseline, fig_dir / "fig4_paired_differences.png")


if __name__ == "__main__":
    main()
