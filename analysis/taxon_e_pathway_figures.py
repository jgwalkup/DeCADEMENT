#!/usr/bin/env python3
"""
taxon_e_pathway_figures.py

Per-taxon absolute energy expenditure with HSP:Biomass annotation.

Fig 1 – Per-taxon HSP:Biomass ratio time series (ctrl vs warm)
         All taxa plotted as lines, coloured by mean control HSP:Biomass
         ratio (YlOrRd colourmap).  The top/bottom TOP_ANNOT taxa by mean
         ratio are annotated at the right margin.

Fig 2 – Per-taxon absolute E cost per pathway
         Grid layout: available pathways (rows) × 2 scenarios (ctrl/warm cols).
         Each taxon is one line, coloured by its mean ctrl HSP:Biomass ratio.
         Top/bottom TOP_ANNOT taxa are annotated "taxon  r=X.XXX" at the
         right margin of each panel.

HSP:Biomass ratio source (in order of preference):
  1. Taxon_HSP_Series.csv + MicrobesSeries.csv
       Taxon_HSP_Series has a (taxon, HSP_type) MultiIndex saved to CSV; after
       loading with index_col=0 the HSP_type column is dropped and taxon names
       repeat once per HSP type.  All HSP types are summed per taxon via
       groupby before computing:
           ratio = total_HSP_pool_C / biomass_C   [unbounded; > 0]
  2. HSP_C_Taxon.csv + MicrobesSeries.csv
       ratio = HSP_C_production / biomass_C       [production rate proxy; unbounded]
"""

import argparse
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
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
    get_palette_from_args,
)


# ── Paths ─────────────────────────────────────────────────────────────────────
BASE = Path(
    "/Users/jswalkup/Library/CloudStorage/"
    "GoogleDrive-jwalkup@uci.edu/My Drive/"
    "Decoupled-Energy-DEMENTpy-20250902/"
    "Decoupled-Energy-DEMENTpy/output"
)

# ── Run selection ──────────────────────────────────────────────────────────────
CTRL_RUNS: list[str] = []
WARM_RUNS: list[str] = []
MIN_COLS  = 100
TOP_ANNOT = 5          # number of highest- and lowest-ratio taxa to annotate

# ── Pathway style ─────────────────────────────────────────────────────────────
C_CTRL   = "#1f77b4"
C_WARM   = "#d62728"
PATHWAY_ORDER  = ["Transport", "Maintenance", "Enzymes", "Osmolytes", "HSP"]
PATHWAY_COLORS = ["#ff7f0e",  "#9467bd",      "#e377c2",  "#bcbd22",   "#8c564b"]
CMAP_RATIO = "plasma"    # colourmap for HSP:Biomass ratio (dark at low, bright at high)

plt.rcParams.update({
    "font.size": 9, "axes.titlesize": 10, "axes.labelsize": 9,
    "legend.fontsize": 7, "figure.dpi": 150,
})

# ── Helpers ───────────────────────────────────────────────────────────────────
def load(run_dir: Path, name: str) -> pd.DataFrame | None:
    f = run_dir / f"{name}.csv"
    if not f.exists():
        return None
    df = pd.read_csv(f, index_col=0)
    # Positional renaming: handles both correct (int) and broken ('C','C.1',...) column names
    df.columns = range(len(df.columns))        # 0 = init state, 1..n = days
    if 0 in df.columns:
        df = df.drop(columns=[0])              # drop init state
    df.columns = range(1, len(df.columns) + 1) # renumber 1..n_days
    return df

def pick_runs(scenario_dir: Path, min_cols: int = MIN_COLS) -> list[Path]:
    if not scenario_dir.exists():
        return []
    by_rep: dict[str, Path] = {}
    for d in sorted(scenario_dir.iterdir()):
        if not d.is_dir():
            continue
        f = d / "C_Uptake_Taxon.csv"
        if not f.exists():
            continue
        try:
            nc = pd.read_csv(f, index_col=0, nrows=1).shape[1]
        except Exception:
            continue
        if nc < min_cols:
            continue
        parts = d.name.split("_")
        rep = parts[1] if len(parts) > 1 else "000"
        if rep not in by_rep or d.name > by_rep[rep].name:
            by_rep[rep] = d
    return list(by_rep.values())

def resolve_runs(names: list[str], base_dir: Path) -> list[Path]:
    dirs = []
    for n in names:
        p = Path(n)
        if p.is_absolute():
            if p.is_dir():
                dirs.append(p)
            else:
                print(f"  WARNING: absolute path not found: {p}")
        elif (base_dir / n).is_dir():
            dirs.append(base_dir / n)
        else:
            print(f"  WARNING: run folder not found: {n!r} (looked in {base_dir})")
    return dirs

def safe_div(a: pd.DataFrame, b: pd.DataFrame) -> pd.DataFrame:
    with np.errstate(divide="ignore", invalid="ignore"):
        r = a / b
    r.replace([np.inf, -np.inf], np.nan, inplace=True)
    return r

def rep_mean_df(pool_list: list, key: str) -> pd.DataFrame | None:
    """Average a (taxa × days) DataFrame across replicate runs."""
    frames = [a[key] for a in pool_list if a and a.get(key) is not None]
    if not frames:
        return None
    idx, cols = frames[0].index, frames[0].columns
    try:
        arr = np.stack([f.reindex(index=idx, columns=cols).values
                        for f in frames], axis=2)
        return pd.DataFrame(np.nanmean(arr, axis=2), index=idx, columns=cols)
    except Exception:
        return frames[0]

# ── Year / pulse helpers ───────────────────────────────────────────────────────
def read_runtime(run_dir: Path | None) -> dict:
    if run_dir is None:
        return {"pulse": 1, "end_time": 365}
    f = run_dir / "Runtime.csv"
    if not f.exists():
        return {"pulse": 1, "end_time": 365}
    try:
        rt  = pd.read_csv(f, index_col=0, header=0)
        col = rt.columns[0]
        return {"pulse":    int(float(rt.loc["pulse",    col])),
                "end_time": int(float(rt.loc["end_time", col]))}
    except Exception:
        return {"pulse": 1, "end_time": 365}

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

def add_year_labels(ax, col_start: int, col_end: int,
                    n_years: int, days_per_year: int) -> None:
    trans = ax.get_xaxis_transform()
    for yr in range(1, n_years + 1):
        boundary = yr * days_per_year
        yr_start = (yr - 1) * days_per_year + 1
        yr_end   = yr * days_per_year
        if col_start < boundary < col_end:
            ax.axvline(boundary, color="gray", lw=0.6, ls=":", alpha=0.45, zorder=0)
        label_day = (max(yr_start, col_start) + min(yr_end, col_end)) / 2
        if col_start <= label_day <= col_end:
            ax.text(label_day, 0.98, f"Year {yr}",
                    transform=trans, ha="center", va="top",
                    fontsize=6.5, color="dimgray", alpha=0.75)

# ── Per-run data loading ───────────────────────────────────────────────────────
def compute_run(run_dir: Path) -> dict | None:
    ET  = load(run_dir, "Transport_E_Taxon")   # taxa × days DataFrames
    EM  = load(run_dir, "Maint_E_Taxon")
    EH  = load(run_dir, "HSP_E_Taxon")
    EE  = load(run_dir, "Enzyme_E_Taxon")
    EO  = load(run_dir, "Osmo_E_Taxon")
    MIC = load(run_dir, "MicrobesSeries")

    # HSP pool data — prefer cumulative pool, fall back to daily production
    HSP_pool = load(run_dir, "Taxon_HSP_Series")
    HSP_C    = load(run_dir, "HSP_C_Taxon")

    if all(x is None for x in [ET, EM, EH, EE, EO]):
        return None

    # ── HSP:Biomass ratio per taxon per day ────────────────────────────────────
    # Taxon_HSP_Series is saved from a (taxon, HSP_type) MultiIndex DataFrame.
    # After loading with index_col=0 the HSP_type becomes a dropped data column,
    # leaving repeated taxon-name rows — one per HSP type.  We sum across them.
    hsp_biomass_ratio = None
    ratio_label = None
    if MIC is not None:
        if HSP_pool is not None:
            # Sum all HSP types per taxon to get total HSP pool C per taxon
            hsp_total = HSP_pool.groupby(HSP_pool.index).sum()
            common_taxa = sorted(set(hsp_total.index) & set(MIC.index))
            common_days = sorted(set(hsp_total.columns) & set(MIC.columns))
            if common_taxa and common_days:
                h = hsp_total.loc[common_taxa, common_days].astype(float)
                m = MIC.loc[common_taxa, common_days].astype(float)
                hsp_biomass_ratio = safe_div(h, m).clip(lower=0)
                ratio_label = "Total HSP pool C / biomass C"
        elif HSP_C is not None:
            common_taxa = sorted(set(HSP_C.index) & set(MIC.index))
            common_days = sorted(set(HSP_C.columns) & set(MIC.columns))
            if common_taxa and common_days:
                h = HSP_C.loc[common_taxa, common_days].astype(float)
                m = MIC.loc[common_taxa, common_days].astype(float)
                hsp_biomass_ratio = safe_div(h, m).clip(lower=0)
                ratio_label = "HSP C production / biomass C (proxy)"

    return dict(
        Transport=ET, Maintenance=EM, HSP=EH, Enzymes=EE, Osmolytes=EO,
        hsp_biomass_ratio=hsp_biomass_ratio,
        ratio_label=ratio_label,
    )

# ── Parse CLI ──────────────────────────────────────────────────────────────────
import argparse as _ap
_p = _ap.ArgumentParser(description="taxon_e_pathway_figures.py (multi-scenario)",
                        formatter_class=_ap.RawDescriptionHelpFormatter)
add_scenario_args(_p)
_args = _p.parse_args()

scenario_list = resolve_scenario_list(_args, BASE, min_cols=MIN_COLS,
                                      sentinel="E_Prod_Taxon")
if not scenario_list:
    print("No scenarios found. Use --base or --scenario-dirs.")
    import sys; sys.exit(1)

scen_labels = [s[0] for s in scenario_list]
COLOR       = scenario_color_map(scen_labels, palette=get_palette_from_args(_args))
N_SCEN      = len(scenario_list)

from pathlib import Path as _Path
base_path = _Path(_args.base) if _args.base else BASE
FIG_DIR   = _Path(_args.fig_dir) if _args.fig_dir else base_path / "taxon_e_pathway_figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

KEYS = ["EP", "Transport", "Maintenance", "HSP", "Enzymes", "Osmolytes", "MIC", "hsp_pool_frac"]

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

for _label in scen_labels:
    scen_data[_label]         = {k: filter_cols(v, COL_START, COL_END)
                                   for k, v in scen_data[_label].items()}
    scen_data[_label + "_sd"] = {k: filter_cols(v, COL_START, COL_END)
                                   for k, v in scen_data[_label + "_sd"].items()}

# ── Backward-compat aliases (ctrl/warm/ctrl_sd/warm_sd/scenarios) ─────────────
# These allow existing two-scenario figure code to work unchanged.
# For 3+ scenarios, ctrl = first, warm = second; extra scenarios available via scen_data.
_l0 = scen_labels[0]
_l1 = scen_labels[1] if N_SCEN > 1 else scen_labels[0]
ctrl    = scen_data[_l0]
warm    = scen_data[_l1]
ctrl_sd = scen_data[_l0 + "_sd"]
warm_sd = scen_data[_l1 + "_sd"]
scenarios = [(lab, scen_data[lab], scen_data[lab + "_sd"]) for lab in scen_labels]

# hsp_pool_frac ratio DataFrames for fig1 (ctrl = first scenario, warm = second)
ctrl_ratio = scen_data[_l0].get("hsp_pool_frac")
warm_ratio = scen_data[_l1].get("hsp_pool_frac")

# ratio_label: description of the HSP:biomass ratio used for colouring
# Stored per-run in compute_run(); extract from first scenario's runs
_first_runs = scenario_list[0][1] if scenario_list else []
_ratio_labels = [compute_run(d).get('ratio_label')
                 for d in _first_runs
                 if compute_run(d) is not None]
ratio_label = next((r for r in _ratio_labels if r is not None),
                   'HSP C / biomass C')

# Build the full taxon list from whichever pathway DataFrame is available first
_ref_df = next(
    (ctrl[k] for k in PATHWAY_ORDER if ctrl.get(k) is not None),
    next((warm[k] for k in PATHWAY_ORDER if warm.get(k) is not None), None)
)
all_taxa = list(_ref_df.index) if _ref_df is not None else []

# mean_ratio: mean biomass fraction per taxon in first scenario (used for colouring)
_mic_ref = scen_data[scen_labels[0]].get('MIC')
mean_ratio = _mic_ref.mean(axis=1) / _mic_ref.mean(axis=1).sum() if _mic_ref is not None else None

if mean_ratio is not None:
    mean_ratio  = mean_ratio.reindex(all_taxa).fillna(0.0)
    taxa_order  = mean_ratio.sort_values(ascending=False).index.tolist()
    vmin_r      = float(mean_ratio.min())
    vmax_r      = float(max(mean_ratio.max(), 1e-9))
else:
    taxa_order  = all_taxa
    vmin_r, vmax_r = 0.0, 1.0

norm_ratio = Normalize(vmin=vmin_r, vmax=vmax_r)
cmap_fn    = cm.get_cmap(CMAP_RATIO)

def ratio_color(taxon):
    if mean_ratio is None or taxon not in mean_ratio.index:
        return "grey"
    return cmap_fn(norm_ratio(float(mean_ratio[taxon])))

# Annotation set: top-N and bottom-N taxa by mean ratio
if mean_ratio is not None and len(mean_ratio) >= 2:
    _sorted = mean_ratio.sort_values(ascending=False)
    annot_taxa = set(_sorted.index[:TOP_ANNOT]) | set(_sorted.index[-TOP_ANNOT:])
else:
    annot_taxa = set()

def annot_label(taxon: str) -> str:
    if mean_ratio is not None and taxon in mean_ratio.index:
        return f"{taxon}  r={mean_ratio[taxon]:.3f}"
    return taxon

# Shared ScalarMappable for colourbars
sm = plt.cm.ScalarMappable(cmap=CMAP_RATIO, norm=norm_ratio)
sm.set_array([])

# ── Helper: draw all taxon lines on an axes ────────────────────────────────────
def plot_taxa(ax, df: pd.DataFrame | None, annotate: bool = True) -> bool:
    """Plot one line per taxon, coloured by mean ctrl HSP:Biomass ratio.
    Returns True if any lines were drawn."""
    if df is None or not taxa_order:
        ax.text(0.5, 0.5, "No data", transform=ax.transAxes,
                ha="center", va="center", fontsize=9)
        return False
    days = df.columns
    for taxon in taxa_order:
        if taxon not in df.index:
            continue
        in_annot = taxon in annot_taxa
        ax.plot(days, df.loc[taxon].values,
                color=ratio_color(taxon),
                lw=1.5 if in_annot else 0.6,
                alpha=0.90 if in_annot else 0.35)
    ax.set_xlim(days.min(), days.max())
    ax.set_ylim(bottom=0)
    ax.spines[["top", "right"]].set_visible(False)
    add_year_labels(ax, COL_START, COL_END, N_YEARS, DAYS_PER_YEAR)

    if annotate:
        for taxon in annot_taxa:
            if taxon not in df.index:
                continue
            series = df.loc[taxon].dropna()
            if series.empty:
                continue
            ax.annotate(annot_label(taxon),
                        xy=(days[-1], series.iloc[-1]),
                        xytext=(4, 0), textcoords="offset points",
                        fontsize=5.5, color=ratio_color(taxon),
                        va="center", clip_on=False)
    return True

# ── Normalise absolute series by gridsize ────────────────────────────────────
# Divides absolute pool/pathway quantities by grid cells so results are
# expressed per grid cell. Ratio/fraction/efficiency series are excluded
# (they are already normalised by biomass, uptake, or total energy).
_NORM_SKIP = {"hsp_pool_frac"}
GRIDSIZE = get_scenario_gridsize(scenario_list)
print(f"Grid size     : {GRIDSIZE} cells  (normalising absolute series by {GRIDSIZE})")
for _label in scen_labels:
    scen_data[_label] = {
        k: (v if k in _NORM_SKIP else normalise_by_gridsize(v, GRIDSIZE))
        for k, v in scen_data[_label].items()
    }
    if _label + '_sd' in scen_data:
        scen_data[_label + '_sd'] = {
            k: (v if k in _NORM_SKIP else normalise_by_gridsize(v, GRIDSIZE))
            for k, v in scen_data[_label + '_sd'].items()
        }

# ────────────────────────────────────────────────────────────────────────────
# Fig 1: Per-taxon HSP:Biomass ratio time series
# ────────────────────────────────────────────────────────────────────────────
fig1, axes1 = plt.subplots(1, N_SCEN, figsize=(7 * N_SCEN, 5), sharey=True)
if N_SCEN == 1: axes1 = [axes1]
fig1.suptitle(
    "Per-Taxon HSP:Biomass Ratio Over Time\n"
    f"({ratio_label}; lines coloured by mean ctrl ratio; "
    f"top/bottom {TOP_ANNOT} taxa annotated)",
    fontsize=11, fontweight="bold",
)

for ax, label in zip(axes1, scen_labels):
    ratio_df = scen_data[label].get("hsp_pool_frac")
    ax.set_title(label)
    plot_taxa(ax, ratio_df)
    ax.set_xlabel("Day")
    ax.set_ylabel(ratio_label)

cbar1 = fig1.colorbar(sm, ax=list(axes1) if hasattr(axes1, "__iter__") else [axes1], shrink=0.85, pad=0.02)
cbar1.set_label(f"Mean ctrl {ratio_label}", fontsize=8)

fig1.tight_layout()
fig1.savefig(FIG_DIR / "fig1_taxon_hsp_biomass_ratio.png", dpi=150, bbox_inches="tight")
plt.close(fig1); print("Saved fig1")

# ────────────────────────────────────────────────────────────────────────────
# Fig 2: Per-taxon absolute E cost per pathway
# ────────────────────────────────────────────────────────────────────────────
available_pathways = [p for p in PATHWAY_ORDER
                      if ctrl.get(p) is not None or warm.get(p) is not None]
n_paths = len(available_pathways)

if n_paths == 0:
    print("Skipped fig2 — no pathway CSVs found in any run.")
else:
    path_color_map = dict(zip(PATHWAY_ORDER, PATHWAY_COLORS))
    fig2, axes2 = plt.subplots(n_paths, 2,
                                figsize=(13, 4 * n_paths),
                                squeeze=False)
    fig2.suptitle(
        "Per-Taxon Absolute E Cost by Pathway — Control vs Warm\n"
        f"(lines coloured by mean ctrl HSP:Biomass ratio; "
        f"top/bottom {TOP_ANNOT} taxa annotated)",
        fontsize=11, fontweight="bold",
    )

    for row, name in enumerate(available_pathways):
        path_color = path_color_map[name]
        for col, (title, dat) in enumerate([("Control", ctrl), ("Warm", warm)]):
            ax = axes2[row, col]
            ax.set_title(f"{name} — {title}", color=path_color, fontweight="bold")
            plot_taxa(ax, dat.get(name))
            ax.set_xlabel("Day")
            ax.set_ylabel("E cost (model units per alive microbe)")

    # Reserve space on the right for a single shared colorbar, then pin it there
    fig2.subplots_adjust(right=0.87, hspace=0.40, wspace=0.30)
    cbar_ax2 = fig2.add_axes([0.89, 0.15, 0.015, 0.70])  # [left, bottom, width, height]
    cbar2 = fig2.colorbar(sm, cax=cbar_ax2)
    cbar2.set_label(f"Mean ctrl {ratio_label}", fontsize=8)

    fig2.savefig(FIG_DIR / "fig2_taxon_e_pathway_lines.png", dpi=150, bbox_inches="tight")
    plt.close(fig2); print("Saved fig2")

print(f"\nAll figures saved to: {FIG_DIR}")
