#!/usr/bin/env python3
"""
realized_vs_potential_figures.py

Regression figures and summary table comparing realized trait scores
(from realized_trait_scores.py) against genomic potential trait scores
(from community_trait_table.py traits_combined.csv).

For each trait × normalization combination, produces one figure with three
panels (one per substrate quartile). Each panel shows realized (y) vs
genomic potential (x) as a colored cloud of points (one per taxon per
temperature), colored by the blue→orange→red temperature gradient.

Genomic potential scalars
--------------------------
Derived from traits_combined.csv by summing across gene columns:
  enzyme potential    : sum(EnzConsti_*) + sum(EnzInduci_*)
  transporter potential: sum(Upt_*)
  HSP potential       : HSPConsti_HSP1 + HSPInduci_HSP1

Gene presence/absence columns (EnzGene_*, HSPGene_*, and the binary
UptakeGenes columns) are excluded from the aggregation.

Taxon ID alignment
------------------
The super_community traits_combined.csv should contain every taxon in
the realized_trait_scores.csv. Mismatches are flagged with a warning
and the unmatched taxa are listed before processing continues.

Regression options
------------------
  --regression none        : scatter only (default)
  --regression all         : one fitted line across all temperatures per panel
  --regression per-temp    : one fitted line per temperature per panel,
                             each drawn in its gradient color

  --show-equation          : annotate each line with slope, intercept, and R²

Output
------
  fig_regression_<trait>_<norm>.png  (9 figures: 3 traits × 3 norms)
  regression_summary.csv             (one row per fit)

Regression summary columns
--------------------------
  trait, normalization, quartile_pct, temperature (NaN for all-points fit),
  n_points, slope, intercept, r_squared, p_value, std_err

Usage
-----
    python analysis/realized_vs_potential_figures.py \\
        --traits-csv    analysis_out/realized_trait_scores.csv \\
        --potential-csv community_pool/super_community/traits/traits_combined.csv \\
        --fig-dir       figures/realized_vs_potential \\
        --regression    per-temp \\
        --show-equation

    # Subset to specific pulses or normalizations
    python analysis/realized_vs_potential_figures.py \\
        --traits-csv    analysis_out/realized_trait_scores.csv \\
        --potential-csv community_pool/super_community/traits/traits_combined.csv \\
        --fig-dir       figures/realized_vs_potential \\
        --pulses        5 \\
        --norms         per_biomass \\
        --regression    per-temp \\
        --show-equation
"""

import argparse
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from scipy import stats

_HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(_HERE))
try:
    from scenario_utils import _TEMP_CMAP_OBJ
except ImportError:
    _TEMP_CMAP_OBJ = None

warnings.filterwarnings("ignore")

# ── Constants ─────────────────────────────────────────────────────────────────

_Q_ORDER  = [75, 50, 25]
_Q_LABELS = {75: "25% lost", 50: "50% lost", 25: "75% lost"}

_TRAITS = [
    ("enzyme_C",    "Enzyme C investment",      "enzyme_potential"),
    ("transport_C", "Transporter C investment", "transport_potential"),
    ("hsp_C",       "HSP C investment",         "hsp_potential"),
]

_NORMS = [
    ("raw",         "C cost (model units)"),
    ("per_biomass", "C cost per unit biomass C"),
    ("per_uptake",  "C cost per unit uptake C"),
]

plt.rcParams.update({
    "font.size": 10, "axes.titlesize": 11, "axes.labelsize": 10,
    "legend.fontsize": 8, "figure.dpi": 150,
})


# ── Gradient colormap ─────────────────────────────────────────────────────────

def _temp_colors(temps: list) -> dict:
    temps = sorted(temps)
    n = len(temps)
    if n == 1:
        return {temps[0]: "#f4a836"}
    if _TEMP_CMAP_OBJ is not None:
        cmap = _TEMP_CMAP_OBJ
    else:
        try:
            cmap = plt.get_cmap("blue_orange_red")
        except Exception:
            cmap = plt.get_cmap("RdYlBu_r")
    positions = np.linspace(0, 1, n)
    return {t: mcolors.to_hex(cmap(p)) for t, p in zip(temps, positions)}


# ── Genomic potential ─────────────────────────────────────────────────────────

def load_genomic_potential(potential_csv: Path) -> pd.DataFrame:
    """
    Load traits_combined.csv and derive three scalar potential scores per taxon.
    Gene presence/absence columns (EnzGene_*, HSPGene_*) are excluded.
    Returns a DataFrame indexed by taxon with columns:
      enzyme_potential, transport_potential, hsp_potential
    """
    df = pd.read_csv(potential_csv, index_col=0)

    enz_consti = [c for c in df.columns if c.startswith("EnzConsti_")]
    enz_induci = [c for c in df.columns if c.startswith("EnzInduci_")]
    upt_cols   = [c for c in df.columns if c.startswith("Upt_")]
    hsp_consti = [c for c in df.columns if c.startswith("HSPConsti_")]
    hsp_induci = [c for c in df.columns if c.startswith("HSPInduci_")]

    pot = pd.DataFrame(index=df.index)
    pot["enzyme_potential"]    = df[enz_consti].sum(axis=1) + df[enz_induci].sum(axis=1)
    pot["transport_potential"] = df[upt_cols].sum(axis=1)
    pot["hsp_potential"]       = (df[hsp_consti].sum(axis=1) +
                                   df[hsp_induci].sum(axis=1))
    return pot


# ── Realized trait aggregation ────────────────────────────────────────────────

def aggregate_realized(df: pd.DataFrame, trait_col: str) -> pd.DataFrame:
    """
    Aggregate realized traits to one value per taxon × temperature × quartile
    by averaging across communities and pulses (alive taxa only).
    Returns columns: taxon, temperature, quartile_pct, <trait_col>
    """
    alive = df[df["alive"] == True].copy()
    if alive.empty or trait_col not in alive.columns:
        return pd.DataFrame()
    grp = (alive
           .groupby(["taxon", "temperature", "quartile_pct"])[trait_col]
           .mean()
           .reset_index(name=trait_col))
    return grp


# ── Alignment and mismatch check ──────────────────────────────────────────────

def check_alignment(realized_taxa: set, potential_taxa: set) -> bool:
    """
    Flag any taxon ID mismatches between the two datasets.
    Returns True if alignment is perfect, False otherwise.
    """
    in_real_not_pot = realized_taxa - potential_taxa
    in_pot_not_real = potential_taxa - realized_taxa

    ok = True
    if in_real_not_pot:
        ok = False
        print(f"  WARNING: {len(in_real_not_pot)} taxa in realized scores "
              f"but NOT in genomic potential:")
        for t in sorted(in_real_not_pot)[:20]:
            print(f"    {t}")
        if len(in_real_not_pot) > 20:
            print(f"    ... and {len(in_real_not_pot) - 20} more")

    if in_pot_not_real:
        # Taxa in potential but not in realized — expected if some taxa
        # never appear in any community simulation
        print(f"  NOTE: {len(in_pot_not_real)} taxa in genomic potential "
              f"but not observed in any realized score "
              f"(may not have been sampled into any community).")

    if ok and not in_real_not_pot:
        print(f"  Taxon alignment: OK ({len(realized_taxa)} taxa matched)")
    return ok


# ── Regression ────────────────────────────────────────────────────────────────

def _fit(x: np.ndarray, y: np.ndarray) -> dict | None:
    """OLS regression of y on x. Returns stats dict or None if insufficient data."""
    mask = np.isfinite(x) & np.isfinite(y)
    if mask.sum() < 3:
        return None
    result = stats.linregress(x[mask], y[mask])
    return {
        "n_points":   int(mask.sum()),
        "slope":      result.slope,
        "intercept":  result.rvalue,
        "r_squared":  result.rvalue ** 2,
        "p_value":    result.pvalue,
        "std_err":    result.stderr,
    }


def _draw_regression(ax, x: np.ndarray, y: np.ndarray, color: str,
                     show_equation: bool, label: str | None = None,
                     eq_y_offset: float = 0.0):
    """Fit and draw one regression line on ax."""
    mask = np.isfinite(x) & np.isfinite(y)
    if mask.sum() < 3:
        return None
    result = stats.linregress(x[mask], y[mask])
    x_range = np.linspace(x[mask].min(), x[mask].max(), 200)
    y_hat   = result.slope * x_range + result.intercept
    ax.plot(x_range, y_hat, color=color, lw=1.4, alpha=0.85, label=label)

    if show_equation:
        eq = (f"y = {result.slope:.3f}x + {result.intercept:.3f}\n"
              f"R² = {result.rvalue**2:.3f}  p = {result.pvalue:.2e}")
        ax.text(0.03, 0.97 - eq_y_offset, eq,
                transform=ax.transAxes,
                fontsize=7, va="top", ha="left",
                color=color, alpha=0.9,
                bbox=dict(boxstyle="round,pad=0.2", fc="white",
                          ec=color, alpha=0.5, lw=0.6))
    return {
        "n_points":  int(mask.sum()),
        "slope":     result.slope,
        "intercept": result.intercept,
        "r_squared": result.rvalue ** 2,
        "p_value":   result.pvalue,
        "std_err":   result.stderr,
    }


# ── Figure builder ────────────────────────────────────────────────────────────

def make_figure(merged: pd.DataFrame,
                trait_prefix: str, trait_label: str,
                potential_col: str,
                norm_suffix: str, norm_label: str,
                temps: list, color_map: dict,
                regression: str, show_equation: bool,
                fig_dir: Path) -> list[dict]:
    """
    Build one 1×3 figure (three quartile panels) and return regression rows.
    """
    realized_col = f"{trait_prefix}_{norm_suffix}"
    if realized_col not in merged.columns:
        print(f"  Skipped {trait_prefix}_{norm_suffix} — column not in merged data")
        return []

    fig, axes = plt.subplots(1, 3, figsize=(5.5 * 3, 5), sharey=False)
    fig.suptitle(
        f"Realized vs Genomic Potential: {trait_label}\n"
        f"({norm_suffix.replace('_', ' ')} normalisation)",
        fontsize=11, fontweight="bold",
    )

    reg_rows = []

    for ax, q in zip(axes, _Q_ORDER):
        sub = merged[merged["quartile_pct"] == q].copy()

        # ── Scatter: one point per taxon per temperature ───────────────────
        for t in temps:
            t_sub = sub[sub["temperature"] == t]
            x = t_sub[potential_col].values
            y = t_sub[realized_col].values
            mask = np.isfinite(x) & np.isfinite(y)
            if mask.any():
                ax.scatter(x[mask], y[mask],
                           color=color_map[t], s=14, alpha=0.45,
                           linewidths=0, zorder=2)

        # ── Regression lines ───────────────────────────────────────────────
        if regression == "per-temp":
            for i, t in enumerate(temps):
                t_sub = sub[sub["temperature"] == t]
                x = t_sub[potential_col].values
                y = t_sub[realized_col].values
                eq_offset = i * 0.14 if show_equation else 0.0
                stats_d = _draw_regression(
                    ax, x, y,
                    color=color_map[t],
                    show_equation=show_equation,
                    label=f"{int(t)}°C",
                    eq_y_offset=eq_offset,
                )
                if stats_d:
                    reg_rows.append({
                        "trait":         trait_prefix,
                        "normalization": norm_suffix,
                        "quartile_pct":  q,
                        "temperature":   t,
                        **stats_d,
                    })

        elif regression == "all":
            x = sub[potential_col].values
            y = sub[realized_col].values
            stats_d = _draw_regression(
                ax, x, y,
                color="#444444",
                show_equation=show_equation,
                label="All temperatures",
            )
            if stats_d:
                reg_rows.append({
                    "trait":         trait_prefix,
                    "normalization": norm_suffix,
                    "quartile_pct":  q,
                    "temperature":   np.nan,
                    **stats_d,
                })

        # All-points fit always goes into the summary table regardless of mode
        if regression != "all":
            x = sub[potential_col].values
            y = sub[realized_col].values
            fit = _fit(x, y)
            if fit:
                reg_rows.append({
                    "trait":         trait_prefix,
                    "normalization": norm_suffix,
                    "quartile_pct":  q,
                    "temperature":   np.nan,
                    **fit,
                })

        # ── Axes formatting ────────────────────────────────────────────────
        ax.set_title(_Q_LABELS[q], fontsize=10)
        ax.set_xlabel(f"Genomic potential\n({trait_label})", fontsize=9)
        ax.set_ylabel(norm_label if ax is axes[0] else "", fontsize=9)
        ax.spines[["top", "right"]].set_visible(False)
        ax.grid(alpha=0.2)

        # Temperature legend on first panel only
        if ax is axes[0]:
            from matplotlib.lines import Line2D
            handles = [
                Line2D([0], [0], marker="o", color="w",
                       markerfacecolor=color_map[t],
                       markersize=7, label=f"{int(t)}°C")
                for t in temps
            ]
            ax.legend(handles=handles, title="Temperature",
                      fontsize=7, loc="upper left",
                      framealpha=0.7)

    fig.tight_layout()
    out = fig_dir / f"fig_regression_{trait_prefix}_{norm_suffix}.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {out.name}")
    return reg_rows


# ── CLI ───────────────────────────────────────────────────────────────────────

def _parse_args():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--traits-csv",    required=True, metavar="PATH",
                   help="realized_trait_scores.py output CSV.")
    p.add_argument("--potential-csv", required=True, metavar="PATH",
                   help="traits_combined.csv from community_trait_table.py.")
    p.add_argument("--fig-dir",       default=None,  metavar="DIR",
                   help="Output directory for figures and summary table "
                        "(default: alongside --traits-csv).")
    p.add_argument("--regression",
                   choices=["none", "all", "per-temp"], default="none",
                   help="Regression line mode: 'none' = scatter only (default); "
                        "'all' = one line across all temperatures; "
                        "'per-temp' = one line per temperature.")
    p.add_argument("--show-equation", action="store_true",
                   help="Annotate regression lines with slope, intercept, R².")
    p.add_argument("--pulses", nargs="+", type=int, default=None, metavar="N",
                   help="Pulse numbers to include (default: all).")
    p.add_argument("--norms", nargs="+",
                   choices=["raw", "per_biomass", "per_uptake"],
                   default=None,
                   help="Normalizations to plot (default: all three).")
    p.add_argument("--traits", nargs="+",
                   choices=["enzyme_C", "transport_C", "hsp_C"],
                   default=None,
                   help="Traits to plot (default: all three).")
    return p.parse_args()


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    args = _parse_args()

    traits_csv    = Path(args.traits_csv)
    potential_csv = Path(args.potential_csv)

    for path, name in [(traits_csv, "--traits-csv"),
                       (potential_csv, "--potential-csv")]:
        if not path.exists():
            sys.exit(f"ERROR: {name} file not found: {path}")

    fig_dir = (Path(args.fig_dir) if args.fig_dir
               else traits_csv.parent / "realized_vs_potential_figures")
    fig_dir.mkdir(parents=True, exist_ok=True)

    # ── Load data ─────────────────────────────────────────────────────────────
    print("Loading realized trait scores ...", flush=True)
    realized = pd.read_csv(traits_csv)

    if args.pulses:
        all_pulses = sorted(realized["pulse"].unique())
        realized   = realized[realized["pulse"].isin(args.pulses)].copy()
        if realized.empty:
            sys.exit(f"ERROR: no data after filtering to pulses {args.pulses}. "
                     f"Available: {all_pulses}")

    print("Loading genomic potential ...", flush=True)
    potential = load_genomic_potential(potential_csv)

    # ── Alignment check ───────────────────────────────────────────────────────
    print("\nChecking taxon ID alignment ...")
    realized_taxa = set(realized["taxon"].dropna().unique())
    potential_taxa = set(potential.index)
    check_alignment(realized_taxa, potential_taxa)

    # ── Temperature and color setup ───────────────────────────────────────────
    temps     = sorted(realized["temperature"].dropna().unique())
    color_map = _temp_colors(temps)

    print(f"\nRealized CSV  : {traits_csv}")
    print(f"Potential CSV : {potential_csv}")
    print(f"Rows (realized): {len(realized):,}")
    print(f"Temperatures  : {temps}")
    print(f"Pulses        : {args.pulses if args.pulses else 'all'}")
    print(f"Regression    : {args.regression}")
    print(f"Show equation : {args.show_equation}")
    print(f"Figure dir    : {fig_dir}\n")

    # ── Apply subsets ─────────────────────────────────────────────────────────
    traits_to_plot = [(p, l, pc) for p, l, pc in _TRAITS
                      if args.traits is None or p in args.traits]
    norms_to_plot  = [(s, l) for s, l in _NORMS
                      if args.norms is None or s in args.norms]

    # ── Build merged dataset once per trait (potential col shared across norms)
    all_reg_rows = []

    for trait_prefix, trait_label, potential_col in traits_to_plot:
        print(f"── {trait_label} ──")

        # Aggregate realized: one value per taxon × temperature × quartile
        realized_cols = ([f"{trait_prefix}_{s}" for s, _ in norms_to_plot]
                         + ["taxon", "temperature", "quartile_pct", "alive"])
        realized_sub  = realized[[c for c in realized_cols
                                   if c in realized.columns]].copy()

        agg_dict = {}
        for norm_suffix, _ in norms_to_plot:
            col = f"{trait_prefix}_{norm_suffix}"
            if col not in realized_sub.columns:
                continue
            alive_sub = realized_sub[realized_sub["alive"] == True]
            agg = (alive_sub
                   .groupby(["taxon", "temperature", "quartile_pct"])[col]
                   .mean()
                   .reset_index(name=col))
            if agg_dict:
                agg_dict[col] = agg.set_index(
                    ["taxon", "temperature", "quartile_pct"]
                )[col]
            else:
                base = agg.copy()
                agg_dict[col] = agg.set_index(
                    ["taxon", "temperature", "quartile_pct"]
                )[col]

        if not agg_dict:
            print(f"  No data for {trait_prefix} — skipping")
            continue

        # Build merged DataFrame: join all norm columns onto base index
        base = base.set_index(["taxon", "temperature", "quartile_pct"])
        for col, s in agg_dict.items():
            base[col] = s
        base = base.reset_index()

        # Join genomic potential
        base = base.join(
            potential[[potential_col]],
            on="taxon", how="inner",
        )

        n_unmatched = len(agg_dict[list(agg_dict.keys())[0]]) - len(base)
        if n_unmatched > 0:
            print(f"  WARNING: {n_unmatched} taxon-temperature rows dropped "
                  f"due to no genomic potential match.")

        # ── Generate figures per normalization ────────────────────────────────
        for norm_suffix, norm_label in norms_to_plot:
            rows = make_figure(
                merged       = base,
                trait_prefix = trait_prefix,
                trait_label  = trait_label,
                potential_col= potential_col,
                norm_suffix  = norm_suffix,
                norm_label   = norm_label,
                temps        = temps,
                color_map    = color_map,
                regression   = args.regression,
                show_equation= args.show_equation,
                fig_dir      = fig_dir,
            )
            all_reg_rows.extend(rows)

    # ── Regression summary table ──────────────────────────────────────────────
    # Regression summary always written — all-points fits are computed
    # regardless of --regression mode; per-temperature fits are added on top
    # when --regression per-temp or --regression all is specified.
    if all_reg_rows:
        reg_df = pd.DataFrame(all_reg_rows)
        reg_df = reg_df.sort_values(
            ["trait", "normalization", "quartile_pct", "temperature"],
            ascending=[True, True, False, True],
            na_position="last",
        ).reset_index(drop=True)
        out_csv = fig_dir / "regression_summary.csv"
        reg_df.to_csv(out_csv, index=False)
        print(f"\nRegression summary: {len(reg_df)} rows → {out_csv}")
        n_all   = (reg_df["temperature"].isna()).sum()
        n_temp  = (~reg_df["temperature"].isna()).sum()
        print(f"  All-points fits  : {n_all}")
        print(f"  Per-temp fits    : {n_temp}")
    else:
        print("\nNo regression rows produced — check that data loaded correctly.")

    print(f"\nAll figures saved to: {fig_dir}")


if __name__ == "__main__":
    main()
