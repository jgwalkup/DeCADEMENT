#!/usr/bin/env python3
"""
substrate_quartile_figures.py

Temperature-gradient figures from substrate_quartile_metrics.py output.

Reads the tidy metrics CSV (one row per scenario × pulse × quartile) and plots
each metric as a response curve across temperature, with one line per plant-litter loss
loss level (25 / 50 / 75 % lost).  Mirrors the style of gradient_response_figures.py.

Error bands
-----------
The metrics CSV holds one value per (temperature, pulse, quartile) — computed on
the ensemble-mean trajectory — so there is no across-community spread in it.
By default this script shows the **across-pulse** spread (mean ± 1 SD over the
pulses at each temperature), which captures pulse-to-pulse variability.

If the CSV also contains matching '<metric>_sd' columns (e.g. produced by a
future per-community version of the metrics script), pass --error community to
use those as the ensemble spread instead.  --error none draws lines only.

Usage
-----
  python analysis/substrate_quartile_figures.py \\
      --metrics-csv analysis_out/substrate_quartile_metrics.csv \\
      --fig-dir     figures/substrate_quartiles \\
      [--error pulse|community|none]
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Quartile (substrate REMAINING) → colour; displayed as % LOST (little lost = cool, much lost = warm)
_Q_COLORS = {75: "#2166ac", 50: "#f4a836", 25: "#d6172a"}
_Q_ORDER = [75, 50, 25]

# Level/rate variables → the four metric suffixes they each carry
_LEVEL_RATE_VARS = [
    ("substrate",        "Substrate"),
    ("monomer",          "Monomer"),
    ("mbc",              "Microbial biomass C"),
    ("enzymes",          "Enzymes"),
    ("respiration",      "Respiration (cumulative)"),
    ("resp_per_biomass", "Mass-specific respiration (cumulative)"),
]
_RATE_SUFFIXES = [
    ("total",         "Total / level"),
    ("inst_rate",     "Instantaneous rate"),
    ("cum_rate",      "Cumulative rate"),
    ("quartile_rate", "Quartile rate"),
]

# Instant-only variable groups → (column, label)
_INSTANT_GROUPS = {
    "efficiency": [("cue", "CUE"), ("cae", "CAE")],
    "c_pathway":  [("c_uptake", "C uptake"), ("c_biosynthesis", "C biosynthesis"),
                   ("c_enzyme", "C enzyme"), ("c_hsp", "C HSP"),
                   ("c_osmolyte", "C osmolyte")],
    "e_balance":  [("e_production", "E production"), ("e_gain", "E gain"),
                   ("e_cost", "E cost")],
    "e_cost":     [("e_transport", "E transport"), ("e_maintenance", "E maintenance"),
                   ("e_hsp", "E HSP"), ("e_enzyme", "E enzyme"),
                   ("e_osmolyte", "E osmolyte")],
}


# ── Aggregation ─────────────────────────────────────────────────────────────────

def _aggregate(df: pd.DataFrame, value_col: str, error_mode: str):
    """Return {quartile: (temps, mean, lo, hi, n)} for one metric column.

    error_mode:
      'pulse'     → mean ± SD across pulses at each temperature
      'community' → use df[value_col] as the mean and df[value_col+'_sd'] as SD
                    (averaged across pulses); falls back to 'pulse' if absent
      'none'      → mean across pulses, no band
    """
    out = {}
    sd_col = f"{value_col}_sd"
    use_community = error_mode == "community" and sd_col in df.columns

    for q in _Q_ORDER:
        sub = df[df["quartile_pct"] == q]
        if sub.empty or value_col not in sub.columns:
            continue
        grp = sub.groupby("temperature")[value_col]
        temps = np.array(sorted(sub["temperature"].dropna().unique()))
        mean = grp.mean().reindex(temps).values
        n    = grp.count().reindex(temps).values

        if error_mode == "none":
            sd = np.full_like(mean, np.nan, dtype=float)
        elif use_community:
            sd = sub.groupby("temperature")[sd_col].mean().reindex(temps).values
        else:  # across-pulse SD
            sd = grp.std(ddof=1).reindex(temps).values

        lo = mean - sd
        hi = mean + sd
        out[q] = (temps, mean, lo, hi, n)
    return out


def _plot_metric(ax, df, value_col, title, error_mode):
    """Plot one metric vs temperature, one line per quartile."""
    series = _aggregate(df, value_col, error_mode)
    any_data = False
    for q in _Q_ORDER:
        if q not in series:
            continue
        temps, mean, lo, hi, n = series[q]
        mask = ~np.isnan(mean)
        if not mask.any():
            continue
        any_data = True
        color = _Q_COLORS[q]
        ax.plot(temps[mask], mean[mask], "-o", color=color, ms=4,
                lw=1.8, label=f"{100 - q}%")   # display as % plant litter lost
        if error_mode != "none":
            band = ~(np.isnan(lo) | np.isnan(hi))
            if band.any():
                ax.fill_between(temps[band], lo[band], hi[band],
                                color=color, alpha=0.18, linewidth=0)
    ax.set_title(title, fontsize=10)
    ax.set_xlabel("Temperature (°C)")
    ax.grid(alpha=0.25)
    if not any_data:
        ax.text(0.5, 0.5, "No data", transform=ax.transAxes,
                ha="center", va="center")
    return any_data


# ── Figure builders ─────────────────────────────────────────────────────────────

def fig_crossing_days(df, fig_dir, error_mode):
    fig, ax = plt.subplots(figsize=(7, 5))
    fig.suptitle("Plant Litter Loss Crossing Day vs Temperature",
                 fontsize=12, fontweight="bold")
    _plot_metric(ax, df, "crossing_pulse_day",
                 "Day within pulse at threshold", error_mode)
    ax.set_ylabel("Crossing day (within pulse)")
    ax.legend(title="Plant litter lost", fontsize=8)
    fig.tight_layout()
    out = fig_dir / "fig_litter_loss_crossing_days.png"
    fig.savefig(out, dpi=150, bbox_inches="tight"); plt.close(fig)
    print(f"Saved {out.name}")


def fig_level_rate(df, var, var_label, fig_dir, error_mode):
    fig, axes = plt.subplots(2, 2, figsize=(12, 9))
    fig.suptitle(f"{var_label} at Plant Litter Loss Quartiles vs Temperature\n"
                 f"(lines = 25/50/75% plant litter lost; band = ±1 SD across pulses)",
                 fontsize=12, fontweight="bold")
    for ax, (suffix, slabel) in zip(axes.ravel(), _RATE_SUFFIXES):
        col = f"{var}_{suffix}"
        ok = _plot_metric(ax, df, col, slabel, error_mode)
        ax.set_ylabel(slabel)
    axes.ravel()[0].legend(title="Plant litter lost", fontsize=8)
    fig.tight_layout()
    out = fig_dir / f"fig_{var}.png"
    fig.savefig(out, dpi=150, bbox_inches="tight"); plt.close(fig)
    print(f"Saved {out.name}")


def fig_instant_group(df, group_name, members, fig_dir, error_mode):
    n = len(members)
    ncol = min(3, n)
    nrow = int(np.ceil(n / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(5 * ncol, 4.2 * nrow),
                             squeeze=False)
    fig.suptitle(f"{group_name.replace('_', ' ').title()} at Plant Litter Loss Quartiles "
                 f"vs Temperature",
                 fontsize=12, fontweight="bold")
    axflat = axes.ravel()
    for ax, (col, label) in zip(axflat, members):
        _plot_metric(ax, df, f"{col}_value", label, error_mode)
        ax.set_ylabel(label)
        if group_name == "efficiency":   # CUE / CAE are bounded [0, 1]
            ax.set_ylim(0, 1)
    for ax in axflat[len(members):]:
        ax.axis("off")
    axflat[0].legend(title="Plant litter lost", fontsize=8)
    fig.tight_layout()
    out = fig_dir / f"fig_instant_{group_name}.png"
    fig.savefig(out, dpi=150, bbox_inches="tight"); plt.close(fig)
    print(f"Saved {out.name}")


# ── CLI / main ───────────────────────────────────────────────────────────────

def fig_taxon_survival(df, fig_dir, error_mode):
    """Surviving-taxon count (biomass > 0 at the loss crossing day) vs temperature."""
    if "survival_value" not in df.columns or df["survival_value"].notna().sum() == 0:
        print("Skipped fig_taxon_survival — no survival_value data "
              "(per-community/--base mode only)")
        return

    # Community taxon-pool size (M) for "N of M" context.
    pool_txt = ""
    if "pool_size" in df.columns:
        pool = df["pool_size"].dropna().unique()
        if len(pool) == 1:
            pool_txt = f" (of {int(round(pool[0]))} taxa per community)"
        elif len(pool) > 1:
            pool_txt = (f" (of {int(round(np.nanmin(pool)))}–"
                        f"{int(round(np.nanmax(pool)))} taxa per community)")

    fig, ax = plt.subplots(figsize=(7, 5))
    fig.suptitle("Taxon Survival at Plant Litter Loss Quartiles vs Temperature",
                 fontsize=12, fontweight="bold")
    _plot_metric(ax, df, "survival_value",
                 f"Surviving taxa (biomass > 0){pool_txt}", error_mode)
    ax.set_ylabel("Surviving taxa (count)")
    ax.legend(title="Plant litter lost", fontsize=8)
    fig.tight_layout()
    out = fig_dir / "fig_taxon_survival.png"
    fig.savefig(out, dpi=150, bbox_inches="tight"); plt.close(fig)
    print(f"Saved {out.name}")


def _parse_args():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--metrics-csv", required=True, metavar="PATH",
                   help="substrate_quartile_metrics.py output CSV.")
    p.add_argument("--fig-dir", default=None, metavar="DIR",
                   help="Directory for figures (default: alongside the CSV).")
    p.add_argument("--error", choices=["pulse", "community", "none"],
                   default="pulse",
                   help="Error band source: 'pulse' = SD across pulses (default); "
                        "'community' = use <metric>_sd columns if present; "
                        "'none' = lines only.")
    return p.parse_args()


def main():
    args = _parse_args()
    csv = Path(args.metrics_csv)
    if not csv.exists():
        sys.exit(f"CSV not found: {csv}")
    df = pd.read_csv(csv)

    fig_dir = Path(args.fig_dir) if args.fig_dir else csv.parent / "substrate_quartile_figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    n_temps = df["temperature"].nunique()
    print(f"Metrics CSV : {csv}")
    print(f"Temperatures: {sorted(df['temperature'].dropna().unique())}")
    print(f"Error mode  : {args.error}")
    print(f"Figure dir  : {fig_dir}\n")

    if args.error == "community" and not any(c.endswith("_sd") for c in df.columns):
        print("  NOTE: --error community requested but no '<metric>_sd' columns "
              "found; falling back to across-pulse SD.\n")
        args.error = "pulse"

    fig_crossing_days(df, fig_dir, args.error)
    fig_taxon_survival(df, fig_dir, args.error)
    for var, label in _LEVEL_RATE_VARS:
        if f"{var}_total" in df.columns:
            fig_level_rate(df, var, label, fig_dir, args.error)
    for group, members in _INSTANT_GROUPS.items():
        if any(f"{c}_value" in df.columns for c, _ in members):
            fig_instant_group(df, group, members, fig_dir, args.error)

    print(f"\nAll figures saved to: {fig_dir}")


if __name__ == "__main__":
    main()
