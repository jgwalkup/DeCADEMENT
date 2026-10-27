#!/usr/bin/env python3
"""
hsp_annotated_e_heatmap.py

Per-taxon energy-allocation heatmaps with taxa annotated and sorted by HSP production.

Layout (GridSpec 2 rows × 3 cols):

  ┌──────────────────┬──────────────────────────┬──────────────────────────┐
  │                  │  Control                 │  Warm                    │
  │  HSP production  │  E cost / E production   │  E cost / E production   │
  │  sidebar bar     │  heatmap  (taxa × time)  │  heatmap  (taxa × time)  │
  │  (spans rows)    ├──────────────────────────┼──────────────────────────┤
  │                  │  HSP fraction of E cost  │  HSP fraction of E cost  │
  │                  │  heatmap  (taxa × time)  │  heatmap  (taxa × time)  │
  └──────────────────┴──────────────────────────┴──────────────────────────┘

Taxa are sorted descending by mean HSP C production (highest HSP producer at top).

HSP annotation source (in order of preference):
  1. HSP_C_Taxon.csv   – daily HSP C production cost per alive microbe (new tracking)
  2. Taxon_HSP_Series.csv – cumulative HSP pool per taxon (proxy for old runs)
"""

import argparse
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.colors import Normalize
from pathlib import Path

warnings.filterwarnings("ignore")
import sys
sys.path.insert(0, str(__import__('pathlib').Path(__file__).parent))
from scenario_utils import (
    load_series, pick_runs, scenario_color_map,
    rep_mean, rep_std, read_runtime, year_range, filter_cols,
    add_year_labels, add_year_labels_heatmap,
    add_scenario_args, resolve_scenario_list,
    normalise_by_gridsize, get_scenario_gridsize,
    resolve_runs,
)

# Raw per-run tables are read with scenario_utils.load_series, which joins split
# branched output (branch_output.json) into the full series, keeps the
# (taxon, type) index of Taxon_HSP_Series / Taxon_Osmo_Series aligned to the
# right day (type level dropped), and renumbers day columns 1..n.
load = load_series


# ── Paths ─────────────────────────────────────────────────────────────────────
BASE = Path(
    "/Users/jswalkup/Library/CloudStorage/"
    "GoogleDrive-jwalkup@uci.edu/My Drive/"
    "Decoupled-Energy-DEMENTpy-20250902/"
    "Decoupled-Energy-DEMENTpy/output"
)

# ── Run selection ──────────────────────────────────────────────────────────────
# Three ways to choose which run folders to use, in order of precedence:
#
#   1. Command-line flags (highest priority):
#        python hsp_annotated_e_heatmap.py \
#            --control control_001_20260305093045 control_002_20260305094512 \
#            --warm    warm_001_20260305093045    warm_002_20260305094512
#      Folder names resolved relative to BASE/control and BASE/warm.
#      Pass absolute paths to use folders anywhere on disk.
#      Override scenario parent dirs with --control-dir / --warm-dir.
#      Override BASE itself with --base.
#
#   2. Hard-coded lists below (used when no CLI flags given):
#
#   3. Auto-pick (default when both lists are empty).
#
CTRL_RUNS: list[str] = []
WARM_RUNS: list[str] = []
MIN_COLS  = 100

# ── Style ─────────────────────────────────────────────────────────────────────
CMAP_EHEAT  = "RdYlBu_r"   # E cost fraction heatmap
CMAP_HSPHEAT = "YlOrRd"    # HSP E fraction heatmap
CMAP_HSPBAR  = "Greens"    # HSP production sidebar bars
SMOOTH = 14

plt.rcParams.update({
    "font.size": 9, "axes.titlesize": 10, "axes.labelsize": 9,
    "legend.fontsize": 8, "figure.dpi": 150,
})

# ── Helpers ───────────────────────────────────────────────────────────────────
def year_range(years_arg, n_years: int, days_per_year: int) -> tuple[int, int]:
    if years_arg is None:
        return 1, n_years * days_per_year
    if len(years_arg) == 1:
        yr = years_arg[0]
        return (yr - 1) * days_per_year + 1, yr * days_per_year
    yr_start, yr_end = years_arg[0], years_arg[-1]
    return (yr_start - 1) * days_per_year + 1, yr_end * days_per_year

def filter_cols(df, col_start: int, col_end: int):
    if df is None:
        return None
    if isinstance(df, pd.Series):
        return df.loc[col_start:col_end]
    cols = [c for c in df.columns if col_start <= c <= col_end]
    return df[cols] if cols else df

def add_year_labels_heatmap(ax, mat_ncols: int, col_start: int, col_end: int,
                             n_years: int, days_per_year: int) -> None:
    total = col_end - col_start + 1
    for yr in range(1, n_years + 1):
        boundary = yr * days_per_year
        if col_start < boundary < col_end:
            frac  = (boundary - col_start) / total
            x_pos = frac * (mat_ncols - 1)
            ax.axvline(x_pos, color="white", lw=0.8, ls=":", alpha=0.6, zorder=5)

def safe_div(a: pd.DataFrame, b: pd.DataFrame) -> pd.DataFrame:
    with np.errstate(divide="ignore", invalid="ignore"):
        r = a / b
    r.replace([np.inf, -np.inf], np.nan, inplace=True)
    return r

def taxon_mean(df: pd.DataFrame | None) -> pd.Series | None:
    if df is None:
        return None
    return df.mean(axis=1, skipna=True)

def subsample(df: pd.DataFrame | None, n: int = 200) -> pd.DataFrame | None:
    if df is None or df.shape[1] <= n:
        return df
    idx = np.linspace(0, df.shape[1] - 1, n, dtype=int)
    return df.iloc[:, idx]

# ── Per-run data loading ───────────────────────────────────────────────────────
def compute_run(run_dir: Path) -> dict | None:
    EP  = load(run_dir, "E_Prod_Taxon")
    ET  = load(run_dir, "Transport_E_Taxon")
    EM  = load(run_dir, "Maint_E_Taxon")
    EH  = load(run_dir, "HSP_E_Taxon")
    EE  = load(run_dir, "Enzyme_E_Taxon")
    EO  = load(run_dir, "Osmo_E_Taxon")

    # HSP C production: prefer new per-day tracking, fall back to pool series
    HSP_C   = load(run_dir, "HSP_C_Taxon")
    HSP_pool = load(run_dir, "Taxon_HSP_Series")

    if EP is None and HSP_C is None and HSP_pool is None:
        return None

    # Total pathway E cost
    total_cost = None
    if ET is not None and EM is not None:
        total_cost = ET.copy()
        for df in [EM, EH, EE, EO]:
            if df is not None:
                cc = sorted(set(total_cost.columns) & set(df.columns))
                total_cost[cc] = total_cost[cc].add(df[cc], fill_value=0)

    cost_frac = safe_div(total_cost, EP).clip(0, 1) if (total_cost is not None and EP is not None) else None

    # HSP fraction of total E cost (per taxon × day)
    hsp_pf = None
    if EH is not None and total_cost is not None:
        cc = sorted(set(EH.columns) & set(total_cost.columns))
        hsp_pf = safe_div(EH[cc], total_cost[cc]).clip(0, 1)

    # Scalar HSP annotation per taxon:
    #   Option 1 – mean daily HSP C production cost (HSP_C_Taxon.csv)
    #   Option 2 – mean HSP pool size             (Taxon_HSP_Series.csv)
    hsp_scalar = None
    hsp_label  = None
    if HSP_C is not None:
        hsp_scalar = HSP_C.mean(axis=1, skipna=True)
        hsp_label  = "Mean HSP C production\n(per alive microbe per day)"
    elif HSP_pool is not None:
        hsp_scalar = HSP_pool.mean(axis=1, skipna=True)
        hsp_label  = "Mean HSP pool C\n(community total; proxy)"

    return dict(
        cost_frac=cost_frac,
        hsp_pf=hsp_pf,
        hsp_scalar=hsp_scalar,
        hsp_label=hsp_label,
    )

# ── Parse CLI ──────────────────────────────────────────────────────────────────
import argparse as _ap
_p = _ap.ArgumentParser(description="hsp_annotated_e_heatmap.py (multi-scenario)",
                        formatter_class=_ap.RawDescriptionHelpFormatter)
add_scenario_args(_p)
_args = _p.parse_args()

scenario_list = resolve_scenario_list(_args, BASE, min_cols=MIN_COLS,
                                      sentinel="E_Prod_Taxon")
if not scenario_list:
    print("No scenarios found. Use --base or --scenario-dirs.")
    import sys; sys.exit(1)

scen_labels = [s[0] for s in scenario_list]
COLOR       = scenario_color_map(scen_labels)
N_SCEN      = len(scenario_list)

from pathlib import Path as _Path
base_path = _Path(_args.base) if _args.base else BASE
FIG_DIR   = _Path(_args.fig_dir) if _args.fig_dir else base_path / "hsp_annotated_e_heatmap"
FIG_DIR.mkdir(parents=True, exist_ok=True)

KEYS = ["cost_frac", "hsp_pf", "hsp_scalar", "hsp_label"]

print(f"Scenarios ({N_SCEN}):")
scen_data = {}
for _label, _run_dirs in scenario_list:
    print(f"  {_label}: {len(_run_dirs)} run(s)")
    _runs = [compute_run(d) for d in _run_dirs]
    scen_data[_label]         = {k: rep_mean(_runs, k) for k in KEYS}
    scen_data[_label + "_sd"] = {k: rep_std(_runs,  k) for k in KEYS}

_first_run = scenario_list[0][1][0] if scenario_list[0][1] else None
_rt = read_runtime(_first_run)
N_YEARS       = _rt["pulse"]
DAYS_PER_YEAR = int(_rt["end_time"])
COL_START, COL_END = year_range(_args.years, N_YEARS, DAYS_PER_YEAR)
print(f"Simulation: {N_YEARS} year(s), {DAYS_PER_YEAR} days/year  |  "
      f"Displaying: day {COL_START}-{COL_END}")

# Day filtering applies to the taxa × day frames only; hsp_scalar is indexed by
# taxon and hsp_label is a string.
for _label in scen_labels:
    scen_data[_label]         = {k: (filter_cols(v, COL_START, COL_END) if isinstance(v, pd.DataFrame) else v)
                                   for k, v in scen_data[_label].items()}
    scen_data[_label + "_sd"] = {k: (filter_cols(v, COL_START, COL_END) if isinstance(v, pd.DataFrame) else v)
                                   for k, v in scen_data[_label + "_sd"].items()}

# ── Backward-compat aliases (ctrl/warm map to first two scenarios) ─────────────
_l0 = scen_labels[0]
_l1 = scen_labels[1] if N_SCEN > 1 else scen_labels[0]
ctrl    = scen_data[_l0]
warm    = scen_data[_l1]
ctrl_sd = scen_data[_l0 + "_sd"]
warm_sd = scen_data[_l1 + "_sd"]
scenarios = [(lab, scen_data[lab], scen_data[lab + "_sd"]) for lab in scen_labels]


# ── Normalise absolute series by gridsize ────────────────────────────────────
# Fractions and per-alive-microbe scalars are not grid totals; labels are strings.
_NORM_SKIP = {"cost_frac", "hsp_pf", "hsp_scalar", "hsp_label"}
GRIDSIZE = get_scenario_gridsize(scenario_list)
print(f"Grid size     : {GRIDSIZE} cells  (normalising absolute series by {GRIDSIZE})")
for _label in scen_labels:
    scen_data[_label] = {
        k: (v if k in _NORM_SKIP else normalise_by_gridsize(v, GRIDSIZE))
        for k, v in scen_data[_label].items()
    }

# ── Per-scenario hsp variables (backward compat) ──────────────────────────────
ctrl_cf  = scen_data[_l0].get("cost_frac")
ctrl_hpf = scen_data[_l0].get("hsp_pf")
warm_cf  = scen_data[_l1].get("cost_frac")
warm_hpf = scen_data[_l1].get("hsp_pf")


# HSP scalar: prefer control, fall back to warm
def _first_hsp(runs: list) -> tuple[pd.Series | None, str | None]:
    for r in runs:
        if r and r.get("hsp_scalar") is not None:
            return r["hsp_scalar"], r["hsp_label"]
    return None, None

ctrl_hsp_scalar, hsp_label = _first_hsp(ctrl_runs)
if ctrl_hsp_scalar is None:
    ctrl_hsp_scalar, hsp_label = _first_hsp(warm_runs)

# Average HSP scalar across replicates that have it
hsp_frames = [r["hsp_scalar"] for r in ctrl_runs if r and r.get("hsp_scalar") is not None]
if len(hsp_frames) > 1:
    hsp_idx = hsp_frames[0].index
    ctrl_hsp_scalar = pd.Series(
        np.nanmean(np.stack([s.reindex(hsp_idx).values for s in hsp_frames], axis=1), axis=1),
        index=hsp_idx,
    )

if ctrl_hsp_scalar is None:
    print("WARNING: No HSP data found in any run — HSP sidebar will be blank.")
    # Still attempt to build heatmaps sorted by E cost fraction
    hsp_order = taxon_mean(ctrl_cf)
    hsp_order = hsp_order.sort_values(ascending=False).index if hsp_order is not None else None
    hsp_label = "No HSP data (sorted by E cost fraction)"
else:
    hsp_order = ctrl_hsp_scalar.sort_values(ascending=False).index

if hsp_order is None:
    print("ERROR: No data available for any taxa — exiting.")
    raise SystemExit(1)

n_taxa = len(hsp_order)
day_labels = np.linspace(COL_START, COL_END, 6, dtype=int)

# ── Build HSP sidebar color array ─────────────────────────────────────────────
if ctrl_hsp_scalar is not None:
    hsp_vals = ctrl_hsp_scalar.reindex(hsp_order).fillna(0).values
    vmin_h, vmax_h = hsp_vals.min(), hsp_vals.max()
    norm_h = Normalize(vmin=vmin_h, vmax=max(vmax_h, 1e-12))
    bar_colors = plt.get_cmap(CMAP_HSPBAR)(norm_h(hsp_vals))
else:
    hsp_vals   = np.zeros(n_taxa)
    bar_colors = ["#cccccc"] * n_taxa

# ── Figure ────────────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(18, 12))
gs  = GridSpec(
    2, 3, figure=fig,
    width_ratios=[0.38, 1, 1],
    wspace=0.06, hspace=0.35,
    left=0.04, right=0.97, top=0.92, bottom=0.07,
)

ax_hsp  = fig.add_subplot(gs[:, 0])            # HSP sidebar — spans both rows
ax_e_c  = fig.add_subplot(gs[0, 1], sharey=ax_hsp)   # Control E cost heatmap
ax_e_w  = fig.add_subplot(gs[0, 2], sharey=ax_hsp)   # Warm    E cost heatmap
ax_hp_c = fig.add_subplot(gs[1, 1], sharey=ax_hsp)   # Control HSP E frac heatmap
ax_hp_w = fig.add_subplot(gs[1, 2], sharey=ax_hsp)   # Warm    HSP E frac heatmap

ytick_pos = np.arange(n_taxa)

fig.suptitle(
    "Per-Taxon Energy Allocation Heatmaps — Taxa Sorted and Annotated by HSP Production\n"
    "(rows = taxa, highest HSP producer at top; cols = time subsampled to 200 days)",
    fontsize=11, fontweight="bold",
)

# ── HSP sidebar ───────────────────────────────────────────────────────────────
ax_hsp.barh(
    ytick_pos, hsp_vals,
    color=bar_colors, edgecolor="white", linewidth=0.3, height=0.85,
)
ax_hsp.set_yticks(ytick_pos)
ax_hsp.set_yticklabels(hsp_order, fontsize=5.5)
ax_hsp.invert_yaxis()          # highest HSP at top, matching imshow row 0
ax_hsp.set_xlabel(hsp_label if hsp_label else "HSP production", fontsize=8)
ax_hsp.set_title(
    "HSP production\n(control mean)" if ctrl_hsp_scalar is not None else "No HSP data",
    fontsize=9,
)
ax_hsp.tick_params(axis="y", which="both", length=0)
ax_hsp.spines[["top", "right"]].set_visible(False)

# Color each y-tick label to match its bar color
for lbl, col in zip(ax_hsp.get_yticklabels(), bar_colors):
    lbl.set_color(col[:3])   # RGB, ignore alpha for text

# Colorbar for HSP sidebar
sm_hsp = plt.cm.ScalarMappable(cmap=CMAP_HSPBAR,
                                norm=Normalize(vmin=hsp_vals.min(), vmax=hsp_vals.max()))
sm_hsp.set_array([])
cbar_hsp = fig.colorbar(sm_hsp, ax=ax_hsp, orientation="horizontal", pad=0.12, shrink=0.9)
cbar_hsp.ax.tick_params(labelsize=7)
cbar_hsp.set_label("HSP C production (model units)", fontsize=7)

# ── Helper: draw one heatmap panel ────────────────────────────────────────────
def draw_heatmap(ax, df, title, cmap, vmin, vmax, cbar_label, xtick_pos=None):
    """Draw a (taxa × time) heatmap with taxa in hsp_order."""
    ax.set_title(title, fontsize=9)
    if df is None:
        ax.text(0.5, 0.5,
                "Data not available\n(re-run model with\npathway tracking)",
                transform=ax.transAxes, ha="center", va="center", fontsize=8,
                bbox=dict(facecolor="lightyellow", edgecolor="gray", alpha=0.85))
        ax.set_yticks([]); ax.set_xticks([])
        return None

    sub = subsample(df.reindex(hsp_order))
    mat = sub.values
    im  = ax.imshow(mat, aspect="auto", cmap=cmap, vmin=vmin, vmax=vmax,
                    interpolation="nearest", origin="upper")

    ax.set_yticks(ytick_pos)
    ax.set_yticklabels([])          # shared y-axis labels on HSP sidebar
    ax.tick_params(axis="y", length=0)

    n_xticks = 6
    xt_pos = np.linspace(0, mat.shape[1] - 1, n_xticks, dtype=int)
    ax.set_xticks(xt_pos)
    ax.set_xticklabels(day_labels[:n_xticks], fontsize=7)
    ax.set_xlabel("Day (subsampled to 200)", fontsize=8)
    add_year_labels_heatmap(ax, mat.shape[1], COL_START, COL_END, N_YEARS, DAYS_PER_YEAR)

    cb = fig.colorbar(im, ax=ax, shrink=0.75, pad=0.02)
    cb.set_label(cbar_label, fontsize=7)
    cb.ax.tick_params(labelsize=6)
    return im

# ── Row 0: E cost fraction heatmaps ──────────────────────────────────────────
draw_heatmap(ax_e_c, ctrl_cf,
             "Control — E cost / E production",
             CMAP_EHEAT, 0, 1, "E cost fraction")
draw_heatmap(ax_e_w, warm_cf,
             "Warm — E cost / E production",
             CMAP_EHEAT, 0, 1, "E cost fraction")

# ── Row 1: HSP fraction of E cost heatmaps ───────────────────────────────────
# Determine shared color scale across control and warm
hpf_vals = [v for df in [ctrl_hpf, warm_hpf] if df is not None
            for v in df.values.ravel() if np.isfinite(v)]
hpf_vmax = np.percentile(hpf_vals, 95) if hpf_vals else 1.0

draw_heatmap(ax_hp_c, ctrl_hpf,
             "Control — HSP E cost / total E cost",
             CMAP_HSPHEAT, 0, hpf_vmax, "HSP E fraction")
draw_heatmap(ax_hp_w, warm_hpf,
             "Warm — HSP E cost / total E cost",
             CMAP_HSPHEAT, 0, hpf_vmax, "HSP E fraction")

# ── Save ──────────────────────────────────────────────────────────────────────
out = FIG_DIR / "hsp_annotated_e_heatmap.png"
fig.savefig(out, dpi=150, bbox_inches="tight")
plt.close(fig)
print(f"Saved: {out}")
