#!/usr/bin/env python3
"""
heatwave_carryover_figure.py
============================
Does an earlier heatwave change the response to a later one?

The effect of one heatwave (the "focal" wave) is isolated as a contrast between
two scenarios that differ only by that wave, and the same contrast is drawn for
communities with a different heatwave history.  Defaults: the July wave of the
summer heatwave-timing experiment,

    July − Ambient                      no earlier wave
    June + July − June                  after a June wave
    July + August − August              no earlier wave, August wave follows
    All three − June + August           after a June wave, August wave follows

Lines that coincide before the focal wave and separate during or after it show
a carry-over of the earlier wave; the pairs with and without a later wave are
identical until that wave starts.

Pools (biomass C, living cells, plant-litter substrate C) are shown as a
percentage of the reference scenario on the same day, because a community that
has been through an earlier wave is smaller; rates are per unit biomass C.
All values are differences of community means (= means of the paired
differences, the communities being paired across scenarios).

Usage
-----
  python analysis/heatwave_carryover_figure.py --summary-dir analysis_output/summary_summer_heatwaves_20261110
      [--contrasts TREATMENT:REFERENCE ...] [--days START END]
      [--focal july:2048-2054] [--heatwave june:2020-2026 august:2076-2082]
      [--smooth 0] [--out FILE.png]

Output: <summary-dir>/figures/carryover/carryover_<focal>.png and .csv
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

_HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(_HERE))
from heatwave_delta_figures import read_mean, read_meta, parse_heatwaves   # noqa: E402
from scenario_utils import pretty_label                                    # noqa: E402

DEFAULT_CONTRASTS = ["climate_July_heatwave:climate_ambient_temp",
                     "climate_June_July_heatwave:climate_June_heatwave",
                     "climate_July_August_heatwave:climate_August_heatwave",
                     "climate_all_three_heatwave:climate_June_August_heatwave"]
# one colour per heatwave history before the focal wave, dashed when another wave follows
STYLES = [dict(color="#1f77b4", ls="-"), dict(color="#d62728", ls="-"),
          dict(color="#1f77b4", ls="--"), dict(color="#d62728", ls="--")]


def _litter(g):
    parts = [g(s) for s in ("system_substrate_C", "system_deadmic_substrate_C", "system_deadenz_substrate_C")]
    return None if any(p is None for p in parts) else parts[0] - parts[1] - parts[2]


# (title, y label, series getter, relative to the reference scenario in %)
PANELS = [
    ("Microbial biomass C", "% of reference", lambda g: g("system_microbial_biomass_C"), True),
    ("Living cells", "% of reference", lambda g: g("system_taxon_count"), True),
    ("Plant-litter substrate C", "% of reference", _litter, True),
    ("Daily CUE (model, floored)", "CUE", lambda g: g("system_CUE"), False),
    ("C uptake", "C d⁻¹ per biomass C", lambda g: g("comm_uptake_C_per_biomass"), False),
    ("C retained", "C d⁻¹ per biomass C", lambda g: g("comm_retained_C_per_biomass"), False),
    ("Endogenous respiration", "C d⁻¹ per biomass C", lambda g: g("comm_endo_resp_C_per_biomass"), False),
    ("Energy production", "E d⁻¹ per biomass C", lambda g: g("comm_production_E_per_biomass"), False),
    ("HSP energy cost", "E d⁻¹ per biomass C", lambda g: g("comm_hsp_E_per_biomass"), False),
    ("Maintenance energy cost", "E d⁻¹ per biomass C", lambda g: g("comm_maintenance_E_per_biomass"), False),
    ("Enzyme energy cost", "E d⁻¹ per biomass C", lambda g: g("comm_enzyme_E_per_biomass"), False),
]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--summary-dir", required=True)
    ap.add_argument("--contrasts", nargs="+", default=DEFAULT_CONTRASTS, metavar="TREATMENT:REFERENCE",
                    help="scenario folder pairs that differ only by the focal wave (up to 4)")
    ap.add_argument("--focal", default="july:2048-2054", metavar="NAME:START-END",
                    help="the heatwave the contrasts isolate (shaded dark)")
    ap.add_argument("--heatwave", nargs="+", default=["june:2020-2026", "august:2076-2082"],
                    metavar="NAME:START-END", help="other heatwaves to shade (light); 'none' for none")
    ap.add_argument("--days", type=int, nargs=2, metavar=("START", "END"), default=None,
                    help="day window (default: 3 weeks before the focal wave to the end of its year)")
    ap.add_argument("--smooth", type=int, default=0, help="running-mean window in days (0 = none)")
    ap.add_argument("--out", default=None, help="default <summary-dir>/figures/carryover/carryover_<focal>.png")
    a = ap.parse_args()

    root = Path(a.summary_dir)
    contrasts = [tuple(c.split(":")) for c in a.contrasts]
    for s in {s for c in contrasts for s in c}:
        if not (root / s / "summary_meta.csv").exists():
            sys.exit(f"ERROR: scenario {s!r} not found under {root}")
    (fname, f0, f1), = parse_heatwaves([a.focal])
    others = parse_heatwaves(a.heatwave)
    meta = read_meta(root / contrasts[0][1])
    if a.days:
        lo, hi = a.days
    else:
        lo, hi = f0 - 21, ((f0 - 1) // meta["end_time"] + 1) * meta["end_time"]
    out = Path(a.out) if a.out else root / "figures" / "carryover" / f"carryover_{fname}.png"
    out.parent.mkdir(parents=True, exist_ok=True)

    ncol = 4
    nrow = -(-len(PANELS) // ncol)
    fig, axes = plt.subplots(nrow, ncol, figsize=(4.6 * ncol, 3.1 * nrow), squeeze=False)
    rows = []
    for ax, (title, ylab, fn, relative) in zip(axes.ravel(), PANELS):
        for (treat, ref), style in zip(contrasts, STYLES):
            x = fn(lambda stem: read_mean(root / treat, stem, lo, hi))
            r = fn(lambda stem: read_mean(root / ref, stem, lo, hi))
            if x is None or r is None:
                continue
            d = (x - r) / r * 100 if relative else x - r
            if a.smooth > 1:
                d = d.rolling(a.smooth, center=True, min_periods=1).mean()
            label = f"{pretty_label(treat)} − {pretty_label(ref)}"
            ax.plot(d.index, d.values, lw=1.4, label=label, **style)
            rows.append(pd.DataFrame({"day": d.index, "quantity": title, "unit": ylab,
                                      "contrast": label, "delta": d.values}))
        for _, s0, s1 in others:
            if s1 >= lo and s0 <= hi:
                ax.axvspan(s0, s1, color="0.92", zorder=0)
        ax.axvspan(f0, f1, color="0.78", zorder=0)
        ax.axhline(0, color="k", lw=0.7)
        ax.set_title(f"Δ {title}", fontsize=10)
        ax.set_ylabel(ylab, fontsize=9)
        ax.set_xlim(lo, hi)
        ax.tick_params(labelsize=8)
    for ax in axes[-1]:
        ax.set_xlabel("day")
    for ax in axes.ravel()[len(PANELS):]:
        ax.set_visible(False)
    fig.suptitle(f"Effect of the {fname.capitalize()} heatwave (dark band) by heatwave history: "
                 f"scenario with the wave − the same scenario without it", y=1.0)
    fig.tight_layout()
    h, l = axes[0, 0].get_legend_handles_labels()
    fig.legend(h, l, loc="upper center", bbox_to_anchor=(0.5, 0.0), ncols=2, frameon=False)
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    pd.concat(rows, ignore_index=True).to_csv(out.with_suffix(".csv"), index=False)
    print(f"Saved {out} and {out.with_suffix('.csv').name}")


if __name__ == "__main__":
    main()
