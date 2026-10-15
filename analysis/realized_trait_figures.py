#!/usr/bin/env python3
"""
realized_trait_figures.py

Figures from realized_trait_scores.py output.

Produces two figure types for each trait (enzyme, transporter, HSP) and each
normalization (raw, per_biomass, per_uptake):

Type 1 — Temperature response curves
--------------------------------------
  X axis   : temperature
  Y axis   : community-mean trait value (averaged across taxa within each
              community, then mean ± SD across communities)
  Lines    : one per substrate quartile (75 / 50 / 25 % remaining)
  Coloring : quartile color scheme matching substrate_quartile_figures.py
  Output   : fig_response_<trait>_<norm>.png  (9 figures total)

Type 2 — Per-taxon violin plots
---------------------------------
  X axis   : temperature treatment (one violin per temperature)
  Y axis   : per-taxon mean trait value (averaged across communities,
              alive taxa only)
  Panels   : one per substrate quartile (75 / 50 / 25 %), shared y-axis
  Coloring : blue→orange→red gradient (same as other figure scripts)
  Output   : fig_violin_<trait>_<norm>.png  (9 figures total)

Usage
-----
    python analysis/realized_trait_figures.py \\
        --traits-csv  analysis_out/realized_trait_scores.csv \\
        --fig-dir     figures/realized_traits \\
        [--error      community|none]
        [--palette    temperature|discrete]
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

_HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(_HERE))
try:
    from scenario_utils import _TEMP_CMAP_OBJ
except ImportError:
    _TEMP_CMAP_OBJ = None

# ── Constants ─────────────────────────────────────────────────────────────────

# Substrate quartile colors — % remaining → color; displayed as % lost
_Q_COLORS = {75: "#2166ac", 50: "#f4a836", 25: "#d6172a"}
_Q_ORDER   = [75, 50, 25]    # quartile_pct values (substrate remaining %)
_Q_LABELS  = {75: "25% lost", 50: "50% lost", 25: "75% lost"}

# Traits — (column_prefix, display_label)
_TRAITS = [
    ("enzyme_C",    "Enzyme C investment"),
    ("transport_C", "Transporter C investment"),
    ("hsp_C",       "HSP C investment"),
]

# Normalizations — (suffix, y-axis label)
_NORMS = [
    ("raw",          "C cost (model units)"),
    ("per_biomass",  "C cost per unit biomass C"),
    ("per_uptake",   "C cost per unit uptake C"),
]

plt.rcParams.update({
    "font.size": 10, "axes.titlesize": 11, "axes.labelsize": 10,
    "legend.fontsize": 9, "figure.dpi": 150,
})


# ── Gradient colormap helpers ─────────────────────────────────────────────────

def _temp_colors(temps: list) -> dict:
    """Map sorted temperature values to colors from the blue→orange→red gradient."""
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


# ── Aggregation ───────────────────────────────────────────────────────────────

def _community_means(df: pd.DataFrame, col: str) -> pd.DataFrame:
    """
    Step 1: average across alive taxa within each community × temperature ×
    quartile × pulse. Returns one row per community × temperature × quartile
    × pulse.
    """
    alive = df[df["alive"] == True].copy()
    if alive.empty or col not in alive.columns:
        return pd.DataFrame()
    grp = alive.groupby(
        ["scenario", "temperature", "community", "pulse", "quartile_pct"]
    )[col].mean().reset_index(name=col)
    return grp


def _response_agg(comm_means: pd.DataFrame, col: str, error_mode: str):
    """
    Step 2 (Type 1): aggregate community means across communities and pulses
    to get mean ± SD per temperature × quartile.

    Returns {quartile_pct: (temps, mean, lo, hi)}.
    """
    out = {}
    sd_col = f"{col}_sd"
    use_community_sd = (error_mode == "community" and
                        sd_col in comm_means.columns)

    for q in _Q_ORDER:
        sub = comm_means[comm_means["quartile_pct"] == q]
        if sub.empty:
            continue
        grp   = sub.groupby("temperature")[col]
        temps = np.array(sorted(sub["temperature"].dropna().unique()))
        mean  = grp.mean().reindex(temps).values

        if error_mode == "none":
            sd = np.full_like(mean, np.nan, dtype=float)
        elif use_community_sd:
            sd = sub.groupby("temperature")[sd_col].mean().reindex(temps).values
        else:
            sd = grp.std(ddof=1).reindex(temps).values

        out[q] = (temps, mean, mean - sd, mean + sd)
    return out


def _taxon_means(df: pd.DataFrame, col: str) -> pd.DataFrame:
    """
    Compute per-taxon mean trait values averaged across communities,
    alive taxa only. Returns one row per taxon × temperature × quartile × pulse.
    """
    alive = df[df["alive"] == True].copy()
    if alive.empty or col not in alive.columns:
        return pd.DataFrame()
    grp = alive.groupby(
        ["taxon", "temperature", "quartile_pct", "pulse"]
    )[col].mean().reset_index(name=col)
    return grp


# ── Figure Type 1: Temperature response curves ────────────────────────────────

def fig_response(df: pd.DataFrame, trait_prefix: str, trait_label: str,
                 norm_suffix: str, norm_label: str,
                 fig_dir: Path, error_mode: str):
    col = f"{trait_prefix}_{norm_suffix}"
    comm_means = _community_means(df, col)
    if comm_means.empty:
        print(f"  Skipped fig_response_{trait_prefix}_{norm_suffix} — no data")
        return

    series = _response_agg(comm_means, col, error_mode)
    if not series:
        print(f"  Skipped fig_response_{trait_prefix}_{norm_suffix} — no data")
        return

    fig, ax = plt.subplots(figsize=(7, 5))
    fig.suptitle(
        f"{trait_label} vs Temperature\n"
        f"(community mean ± 1 SD, {norm_suffix.replace('_', ' ')} normalisation)",
        fontsize=11, fontweight="bold",
    )

    for q in _Q_ORDER:
        if q not in series:
            continue
        temps, mean, lo, hi = series[q]
        mask = ~np.isnan(mean)
        if not mask.any():
            continue
        color = _Q_COLORS[q]
        ax.plot(temps[mask], mean[mask], "-o", color=color, ms=4,
                lw=1.8, label=_Q_LABELS[q])
        if error_mode != "none":
            band = ~(np.isnan(lo) | np.isnan(hi))
            if band.any():
                ax.fill_between(temps[band], lo[band], hi[band],
                                color=color, alpha=0.18, linewidth=0)

    ax.set_xlabel("Temperature (°C)")
    ax.set_ylabel(norm_label)
    ax.legend(title="Plant litter lost", fontsize=9)
    ax.grid(alpha=0.25)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()

    out = fig_dir / f"fig_response_{trait_prefix}_{norm_suffix}.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {out.name}")


# ── Figure Type 2: Per-taxon violin plots ─────────────────────────────────────

def fig_violin(df: pd.DataFrame, trait_prefix: str, trait_label: str,
               norm_suffix: str, norm_label: str,
               fig_dir: Path, temp_color_map: dict):
    col = f"{trait_prefix}_{norm_suffix}"
    tax_means = _taxon_means(df, col)
    if tax_means.empty:
        print(f"  Skipped fig_violin_{trait_prefix}_{norm_suffix} — no data")
        return

    temps   = sorted(tax_means["temperature"].dropna().unique())
    n_temps = len(temps)
    colors  = temp_color_map

    fig, axes = plt.subplots(
        1, 3, figsize=(5 * 3, 5),
        sharey=True,
    )
    fig.suptitle(
        f"{trait_label} — Per-Taxon Distribution vs Temperature\n"
        f"({norm_suffix.replace('_', ' ')} normalisation; "
        f"mean across communities, alive taxa only)",
        fontsize=11, fontweight="bold",
    )

    for ax, q in zip(axes, _Q_ORDER):
        sub = tax_means[tax_means["quartile_pct"] == q]

        # Collect data per temperature — one list of taxon-mean values per violin
        data_by_temp = []
        positions    = []
        valid_temps  = []
        for i, t in enumerate(temps):
            vals = sub[sub["temperature"] == t][col].dropna().values
            if len(vals) > 0:
                data_by_temp.append(vals)
                positions.append(i)
                valid_temps.append(t)

        if not data_by_temp:
            ax.set_title(_Q_LABELS[q])
            ax.text(0.5, 0.5, "No data", transform=ax.transAxes,
                    ha="center", va="center")
            continue

        # Draw violins
        vp = ax.violinplot(
            data_by_temp,
            positions=positions,
            showmedians=True,
            showextrema=False,
            widths=0.7,
        )

        # Color each violin body by temperature gradient
        for body, t in zip(vp["bodies"], valid_temps):
            body.set_facecolor(colors[t])
            body.set_edgecolor("none")
            body.set_alpha(0.75)

        vp["cmedians"].set_color("black")
        vp["cmedians"].set_linewidth(1.2)

        ax.set_title(_Q_LABELS[q], fontsize=10)
        ax.set_xticks(range(len(temps)))
        ax.set_xticklabels(
            [f"{int(t)}°C" for t in temps],
            rotation=45, ha="right", fontsize=8,
        )
        ax.set_xlabel("Temperature")
        ax.grid(axis="y", alpha=0.25)
        ax.spines[["top", "right"]].set_visible(False)

    axes[0].set_ylabel(norm_label)

    fig.tight_layout()
    out = fig_dir / f"fig_violin_{trait_prefix}_{norm_suffix}.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {out.name}")


# ── CLI / main ────────────────────────────────────────────────────────────────

def _parse_args():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--traits-csv", required=True, metavar="PATH",
                   help="realized_trait_scores.py output CSV.")
    p.add_argument("--fig-dir", default=None, metavar="DIR",
                   help="Output directory for figures "
                        "(default: alongside the CSV).")
    p.add_argument("--error", choices=["community", "none"],
                   default="community",
                   help="Error band source for response curves: "
                        "'community' = SD across communities (default); "
                        "'none' = lines only.")
    p.add_argument("--traits", nargs="+",
                   choices=["enzyme_C", "transport_C", "hsp_C"],
                   default=None,
                   help="Subset of traits to plot (default: all three).")
    p.add_argument("--norms", nargs="+",
                   choices=["raw", "per_biomass", "per_uptake"],
                   default=None,
                   help="Subset of normalizations to plot (default: all three).")
    p.add_argument("--pulses", nargs="+", type=int, default=None, metavar="N",
                   help="Pulse numbers to include (1-indexed). "
                        "Default: all pulses. E.g. --pulses 5 for the final pulse.")
    return p.parse_args()


def main():
    args   = _parse_args()
    csv    = Path(args.traits_csv)
    if not csv.exists():
        sys.exit(f"CSV not found: {csv}")

    fig_dir = (Path(args.fig_dir) if args.fig_dir
               else csv.parent / "realized_trait_figures")
    fig_dir.mkdir(parents=True, exist_ok=True)

    print(f"Loading {csv} ...", flush=True)
    df = pd.read_csv(csv)

    if args.pulses:
        all_pulses = sorted(df["pulse"].unique())
        df = df[df["pulse"].isin(args.pulses)].copy()
        if df.empty:
            sys.exit(f"ERROR: no data after filtering to pulses {args.pulses}. "
                     f"Available: {all_pulses}")

    temps = sorted(df["temperature"].dropna().unique())
    print(f"Traits CSV  : {csv}")
    print(f"Rows        : {len(df):,}")
    print(f"Temperatures: {temps}")
    print(f"Pulses      : {args.pulses if args.pulses else 'all'}")
    print(f"Communities : {df['community'].nunique()}")
    print(f"Unique taxa : {df['taxon'].nunique()}")
    print(f"Error mode  : {args.error}")
    print(f"Figure dir  : {fig_dir}\n")

    # Build temperature → color map once for violin plots
    color_map = _temp_colors(temps)

    # Apply trait/norm subsets if specified
    traits = [(p, l) for p, l in _TRAITS
              if args.traits is None or p in args.traits]
    norms  = [(s, l) for s, l in _NORMS
              if args.norms is None or s in args.norms]

    for trait_prefix, trait_label in traits:
        print(f"── {trait_label} ──")
        for norm_suffix, norm_label in norms:
            fig_response(df, trait_prefix, trait_label,
                         norm_suffix, norm_label,
                         fig_dir, args.error)
            fig_violin(df, trait_prefix, trait_label,
                       norm_suffix, norm_label,
                       fig_dir, color_map)

    print(f"\nAll figures saved to: {fig_dir}")


if __name__ == "__main__":
    main()
