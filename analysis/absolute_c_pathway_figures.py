#!/usr/bin/env python3
"""
absolute_c_pathway_figures.py

Supports raw runs (--base / --scenario-dirs) or pre-aggregated summary CSVs
(--summary-dir): pathways come from comm_*_C_per_microbe stems (community mean
per alive microbe), whose cross-community mean/SD match raw mode exactly.

Absolute carbon allocation pathway breakdown for DEMENTpy.
Analogous to absolute_e_pathway_figures.py but for C use rather than energy use.

C pathways tracked per taxon (per alive microbe, community mean):
  • Enzymes   – Enzyme_C_Taxon.csv   (Taxon_Enzyme_Cost_C in grid.py)
  • HSP       – HSP_C_Taxon.csv      (Taxon_HSP_Cost_C in grid.py)
  • Osmolytes – Osmo_C_Taxon.csv     (Taxon_Osmo_Cost_C in grid.py)

Overlay series:
  • C uptake        – C_Uptake_Taxon.csv   (total C taken up)
  • Biosynthesis C  – BC_Uptake_Taxon.csv  (C allocated to biomass / growth)

Note: Transport and Maintenance are energy-only costs in this model formulation;
      they do not appear as explicit C-cost pathways.

Fig 1  – Stacked absolute C cost pathways + C uptake & biosynthesis overlays (ctrl vs warm)
Fig 2  – Absolute C cost per pathway as individual lines (ctrl vs warm)
Fig 3  – Warm − Control Δ absolute C cost per pathway
"""

import argparse
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
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
    add_summary_arg, resolve_summary_scenarios, load_summary, read_summary_meta,
    _parse_treatment_map,
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
SMOOTH    = 14

# ── Pathway colours ───────────────────────────────────────────────────────────
C_HSP    = "#8c564b"
C_ENZ    = "#e377c2"
C_OSMO   = "#bcbd22"
C_UPTAKE = "#2ca02c"    # total C uptake overlay
C_BIOSYN = "#1f77b4"    # biosynthesis C overlay
C_CTRL   = "#1f77b4"
C_WARM   = "#d62728"

PATHWAY_ORDER  = ["Enzymes", "Osmolytes", "HSP"]
PATHWAY_COLORS = [C_ENZ,     C_OSMO,      C_HSP]

plt.rcParams.update({
    "font.size": 10, "axes.titlesize": 11, "axes.labelsize": 10,
    "legend.fontsize": 9, "figure.dpi": 150,
})

# ── Helpers ───────────────────────────────────────────────────────────────────
def load(run_dir: Path, name: str) -> pd.DataFrame | None:
    f = run_dir / f"{name}.csv"
    if not f.exists():
        return None
    df = pd.read_csv(f, index_col=0)
    # Positional renaming: handles both correct (int) and broken ('C','C.1',...) column names
    df.columns = range(len(df.columns))        # 0 = init state, 1..n = days
    df = df.drop(columns=[0], errors="ignore") # drop init state column
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
            dirs.append(p) if p.is_dir() else print(f"  WARNING: not found: {p}")
        elif (base_dir / n).is_dir():
            dirs.append(base_dir / n)
        else:
            print(f"  WARNING: run folder not found: {n!r} (looked in {base_dir})")
    return dirs

def rep_mean(pool_list: list, key: str) -> pd.Series | None:
    """Average a per-day community Series across replicate runs."""
    frames = [a[key] for a in pool_list if a and a.get(key) is not None]
    if not frames:
        return None
    if len(frames) == 1:
        return frames[0]
    idx = frames[0].index
    arr = np.stack([s.reindex(idx).values for s in frames], axis=1)
    return pd.Series(np.nanmean(arr, axis=1), index=idx)

def rep_std(pool_list: list, key: str) -> pd.Series | None:
    frames = [a[key] for a in pool_list if a and a.get(key) is not None]
    if len(frames) < 2:
        return None
    idx = frames[0].index
    arr = np.stack([s.reindex(idx).values for s in frames], axis=1)
    return pd.Series(np.nanstd(arr, axis=1, ddof=1), index=idx)

def smooth(s: pd.Series) -> pd.Series:
    return s

# ── Year / pulse helpers ───────────────────────────────────────────────────────
def read_runtime(run_dir: Path | None) -> dict:
    """Read Runtime.csv → {'pulse': n_years, 'end_time': days_per_year}."""
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
    """Convert --years argument to (col_start, col_end) inclusive day range."""
    if years_arg is None:
        return 1, n_years * days_per_year
    if len(years_arg) == 1:
        yr = years_arg[0]
        return (yr - 1) * days_per_year + 1, yr * days_per_year
    yr_start, yr_end = years_arg[0], years_arg[-1]
    return (yr_start - 1) * days_per_year + 1, yr_end * days_per_year

def filter_series(s, col_start: int, col_end: int):
    """Slice a per-day Series to [col_start, col_end]."""
    if s is None:
        return None
    return s.loc[col_start:col_end]

def add_year_labels(ax, col_start: int, col_end: int,
                    n_years: int, days_per_year: int) -> None:
    """Draw dashed year-boundary lines and 'Year N' labels on a time-series axes."""
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

# ── Per-run loading ────────────────────────────────────────────────────────────
def compute_run(run_dir: Path) -> dict | None:
    """
    Load per-taxon C pathway CSVs and reduce to community-mean per-alive-microbe
    time series.

    Returns a dict with keys:
      Enzymes, HSP, Osmolytes   – explicit C cost pathways
      C_uptake                  – total C taken up (per alive microbe)
      Biosynthesis              – C allocated to biomass (per alive microbe)
      total_cost                – sum of tracked pathway costs
    """
    EE = load(run_dir, "Enzyme_C_Taxon")
    EH = load(run_dir, "HSP_C_Taxon")
    EO = load(run_dir, "Osmo_C_Taxon")
    CU = load(run_dir, "C_Uptake_Taxon")
    BC = load(run_dir, "BC_Uptake_Taxon")

    if all(x is None for x in [EE, EH, EO]):
        return None

    def _comm(df: pd.DataFrame | None) -> pd.Series | None:
        """Community mean across taxa (rows already normalised per alive microbe)."""
        return df.mean(axis=0, skipna=True) if df is not None else None

    ee = _comm(EE)
    eh = _comm(EH)
    eo = _comm(EO)
    cu = _comm(CU)
    bc = _comm(BC)

    pathway_series = [s for s in [ee, eh, eo] if s is not None]
    if pathway_series:
        idx = pathway_series[0].index
        total_cost = pd.Series(
            np.nansum(np.stack([s.reindex(idx).values for s in pathway_series], axis=1), axis=1),
            index=idx,
        )
    else:
        total_cost = None

    return dict(
        Enzymes=ee, HSP=eh, Osmolytes=eo,
        C_uptake=cu, Biosynthesis=bc,
        total_cost=total_cost,
    )

# ── Parse CLI ──────────────────────────────────────────────────────────────────
import argparse as _ap
_p = _ap.ArgumentParser(description="absolute_c_pathway_figures.py (multi-scenario)",
                        formatter_class=_ap.RawDescriptionHelpFormatter)
add_scenario_args(_p)
add_summary_arg(_p)
_p.add_argument("--per-biomass", action="store_true",
                help="Normalise all pathway series by total microbial biomass C "
                     "(system_microbial_biomass_C) instead of alive microbe count. "
                     "Changes y-axis to model units per unit biomass C.")
_args = _p.parse_args()

from pathlib import Path as _Path

KEYS = ["Enzymes", "HSP", "Osmolytes", "C_uptake", "Biosynthesis", "total_cost"]

# key → export stem.  Community-mean-per-microbe series whose cross-community
# mean/SD match raw-mode compute_run exactly.  total_cost is derived-but-unplotted
# (vestigial in KEYS) → None.
_SUMMARY_STEMS = {
    "Enzymes":      "comm_enzyme_C_per_microbe",
    "HSP":          "comm_hsp_C_per_microbe",
    "Osmolytes":    "comm_osmolyte_C_per_microbe",
    "C_uptake":     "comm_uptake_C_per_microbe",
    "Biosynthesis": "comm_biosynthesis_C_per_microbe",
    "total_cost":   None,
}
_MBC_STEM = "system_microbial_biomass_C"

if _args.summary_dir:
    # ── Summary mode (all discovered scenarios → full gradient / N-scenario) ────
    summary_root = _Path(_args.summary_dir)
    _tmap = _parse_treatment_map(getattr(_args, "treatment_map", None))
    scenario_list = resolve_summary_scenarios(summary_root, treatment_map=_tmap)
    if not scenario_list:
        print(f"No summary scenarios found under {summary_root}.")
        import sys; sys.exit(1)

    scen_labels = [s[0] for s in scenario_list]
    COLOR       = scenario_color_map(scen_labels, palette=get_palette_from_args(_args))
    N_SCEN      = len(scenario_list)

    base_path = summary_root
    FIG_DIR   = _Path(_args.fig_dir) if _args.fig_dir else base_path / "absolute_c_pathway_figures"
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    print(f"Summary mode — scenarios ({N_SCEN}):")
    scen_data = {}
    mbc_data  = {}   # label → MBC mean Series (for --per-biomass)
    for _label, _sdir in scenario_list:
        print(f"  {_label}: {_sdir}")
        _pairs = {k: (load_summary(_sdir, _SUMMARY_STEMS[k]) if _SUMMARY_STEMS[k]
                      else (None, None)) for k in KEYS}
        scen_data[_label]         = {k: _pairs[k][0] for k in KEYS}
        scen_data[_label + "_sd"] = {k: _pairs[k][1] for k in KEYS}
        _mbc_mean, _ = load_summary(_sdir, _MBC_STEM)
        mbc_data[_label] = _mbc_mean

    _meta = read_summary_meta(scenario_list[0][1])
    N_YEARS       = _meta["pulse"]
    DAYS_PER_YEAR = int(_meta["end_time"])
    GRIDSIZE      = _meta["gridsize"]
    if _args.years:
        COL_START, COL_END = year_range(_args.years, N_YEARS, DAYS_PER_YEAR)
    elif _meta["col_start"] is not None:
        COL_START, COL_END = _meta["col_start"], _meta["col_end"]
    else:
        _cols = [c for v in scen_data[scen_labels[0]].values()
                 if v is not None
                 for c in getattr(v, "index", [])]
        COL_START, COL_END = (min(_cols), max(_cols)) if _cols else (1, DAYS_PER_YEAR)
    print(f"Simulation: {N_YEARS} year(s), {DAYS_PER_YEAR} days/year  |  "
          f"Displaying: day {COL_START}-{COL_END}  |  gridsize {GRIDSIZE}")

else:
    # ── Raw-run mode (original behaviour) ───────────────────────────────────────
    scenario_list = resolve_scenario_list(_args, BASE, min_cols=MIN_COLS,
                                          sentinel="C_Uptake_Taxon")
    if not scenario_list:
        print("No scenarios found. Use --base or --scenario-dirs.")
        import sys; sys.exit(1)

    scen_labels = [s[0] for s in scenario_list]
    COLOR       = scenario_color_map(scen_labels, palette=get_palette_from_args(_args))
    N_SCEN      = len(scenario_list)

    base_path = _Path(_args.base) if _args.base else BASE
    FIG_DIR   = _Path(_args.fig_dir) if _args.fig_dir else base_path / "absolute_c_pathway_figures"
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    print(f"Scenarios ({N_SCEN}):")
    scen_data = {}
    mbc_data  = {}   # label → MBC mean Series (for --per-biomass)
    for _label, _run_dirs in scenario_list:
        print(f"  {_label}: {len(_run_dirs)} run(s)")
        _runs = [compute_run(d) for d in _run_dirs]
        scen_data[_label]         = {k: rep_mean(_runs, k) for k in KEYS}
        scen_data[_label + "_sd"] = {k: rep_std(_runs,  k) for k in KEYS}
        # MBC not available in raw mode — per-biomass requires summary mode
        mbc_data[_label] = None

    _first_run = scenario_list[0][1][0] if scenario_list[0][1] else None
    _rt = read_runtime(_first_run)
    N_YEARS       = _rt["pulse"]
    DAYS_PER_YEAR = int(_rt["end_time"])
    COL_START, COL_END = year_range(_args.years, N_YEARS, DAYS_PER_YEAR)
    GRIDSIZE = get_scenario_gridsize(scenario_list)
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


# ── Normalise all series by gridsize ─────────────────────────────────────────
print(f"Grid size     : {GRIDSIZE} cells  (normalising all series by {GRIDSIZE})")
for _label in scen_labels:
    scen_data[_label] = {
        k: normalise_by_gridsize(v, GRIDSIZE)
        for k, v in scen_data[_label].items()
    }
    if _label + '_sd' in scen_data:
        scen_data[_label + '_sd'] = {
            k: normalise_by_gridsize(v, GRIDSIZE)
            for k, v in scen_data[_label + '_sd'].items()
        }

# ── Optional per-biomass normalisation ───────────────────────────────────────
PER_BIOMASS = getattr(_args, "per_biomass", False)
if PER_BIOMASS:
    _missing_mbc = [l for l in scen_labels if mbc_data.get(l) is None]
    if _missing_mbc:
        print(f"WARNING: --per-biomass requested but MBC series not available for: "
              f"{_missing_mbc}. Per-biomass normalisation requires --summary-dir.")
        PER_BIOMASS = False
    else:
        print("Applying per-biomass normalisation (dividing by system_microbial_biomass_C)")
        for _label in scen_labels:
            _mbc = mbc_data[_label]
            def _div_by_mbc(s, mbc=_mbc):
                if s is None:
                    return None
                return s.div(mbc.reindex(s.index)).replace([float("inf"),
                                                            float("-inf")], float("nan"))
            scen_data[_label] = {k: _div_by_mbc(v)
                                 for k, v in scen_data[_label].items()}
            # SD: propagate as s_sd / mbc (ignore mbc uncertainty — mbc_sd not loaded)
            scen_data[_label + "_sd"] = {k: _div_by_mbc(v)
                                         for k, v in scen_data[_label + "_sd"].items()}

NORM_LABEL = "per unit biomass C" if PER_BIOMASS else "per alive microbe"
NORM_SUFFIX = "_per_biomass" if PER_BIOMASS else ""

# ────────────────────────────────────────────────────────────────────────────
# Fig 1: Stacked absolute C cost pathways — one panel per scenario
# ────────────────────────────────────────────────────────────────────────────
fig1, axes = plt.subplots(1, N_SCEN, figsize=(7 * N_SCEN, 5), sharey=True)
if N_SCEN == 1: axes = [axes]
fig1.suptitle(
    "Absolute C Cost by Pathway\n"
    f"(community mean {NORM_LABEL}; stacked area = tracked pathway costs; "
    "dashed = C uptake, dotted = biosynthesis C)",
    fontsize=12, fontweight="bold",
)

for ax, label in zip(axes, scen_labels):
    dat    = scen_data[label]
    dat_sd = scen_data[label + "_sd"]
    ax.set_title(label)

    stack_vals   = []
    stack_labels = []
    stack_colors = []
    for name, color in zip(PATHWAY_ORDER, PATHWAY_COLORS):
        s = dat.get(name)
        if s is not None:
            stack_vals.append(smooth(s).values)
            stack_labels.append(name)
            stack_colors.append(color)

    if not stack_vals:
        ax.text(0.5, 0.5, "No C pathway data", transform=ax.transAxes,
                ha="center", va="center")
        continue

    days = next(dat[k] for k in KEYS if dat.get(k) is not None).index
    ax.stackplot(days, stack_vals,
                 labels=stack_labels, colors=stack_colors, alpha=0.85)

    if dat.get("C_uptake") is not None:
        cu_s = smooth(dat["C_uptake"])
        ax.plot(days, cu_s.values, color="black", lw=2, ls="--",
                label="C uptake", zorder=5)
        if dat_sd.get("C_uptake") is not None:
            cu_sd = smooth(dat_sd["C_uptake"])
            ax.fill_between(days, (cu_s - cu_sd).clip(0), cu_s + cu_sd,
                            color="black", alpha=0.10, linewidth=0)

    if dat.get("Biosynthesis") is not None:
        bc_s = smooth(dat["Biosynthesis"])
        ax.plot(days, bc_s.values, color=C_BIOSYN, lw=1.8, ls=":",
                label="Biosynthesis C", zorder=5)
        if dat_sd.get("Biosynthesis") is not None:
            bc_sd = smooth(dat_sd["Biosynthesis"])
            ax.fill_between(days, (bc_s - bc_sd).clip(0), bc_s + bc_sd,
                            color=C_BIOSYN, alpha=0.10, linewidth=0)

    ax.set_xlabel("Day")
    ax.set_ylabel(f"C (model units {NORM_LABEL})")
    ax.set_xlim(days.min(), days.max())
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.22),
              ncols=3, fontsize=8, framealpha=0.85)
    ax.spines[["top", "right"]].set_visible(False)
    add_year_labels(ax, COL_START, COL_END, N_YEARS, DAYS_PER_YEAR)

fig1.tight_layout()
fig1.savefig(FIG_DIR / f"fig1_absolute_c_pathway_stacked{NORM_SUFFIX}.png", dpi=150, bbox_inches="tight")
plt.close(fig1); print("Saved fig1")

# ────────────────────────────────────────────────────────────────────────────
# Fig 2: Individual pathway lines — all scenarios on shared axes per pathway
# ────────────────────────────────────────────────────────────────────────────
available_pathways = [name for name in PATHWAY_ORDER
                      if any(scen_data[l].get(name) is not None for l in scen_labels)]
n_paths = len(available_pathways)

if n_paths == 0:
    print("Skipped fig2 — no C pathway data")
else:
    ncols = min(n_paths, 3)
    nrows = (n_paths + ncols - 1) // ncols
    fig2, axes2 = plt.subplots(nrows, ncols,
                               figsize=(5 * ncols, 4 * nrows),
                               squeeze=False)
    fig2.suptitle(
        "Absolute C Cost Per Pathway — All Scenarios\n"
        f"(community mean {NORM_LABEL})",
        fontsize=12, fontweight="bold",
    )
    path_color_map = dict(zip(PATHWAY_ORDER, PATHWAY_COLORS))

    for idx, name in enumerate(available_pathways):
        ax    = axes2[idx // ncols][idx % ncols]
        color = path_color_map[name]
        days  = None

        for label in scen_labels:
            dat    = scen_data[label]
            dat_sd = scen_data[label + "_sd"]
            s = dat.get(name)
            if s is None:
                continue
            days = s.index
            s_sm = smooth(s)
            ax.plot(days, s_sm.values, color=COLOR[label], lw=1.8, label=label)
            if dat_sd.get(name) is not None:
                sd_sm = smooth(dat_sd[name])
                ax.fill_between(days,
                                (s_sm - sd_sm).clip(0), s_sm + sd_sm,
                                color=COLOR[label], alpha=0.18, linewidth=0)

        ax.set_title(name, color=color, fontweight="bold")
        ax.set_xlabel("Day")
        ax.set_ylabel(f"C cost (model units {NORM_LABEL})")
        if days is not None:
            ax.set_xlim(days.min(), days.max())
        ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.22),
                  ncols=min(N_SCEN, 3), fontsize=8, framealpha=0.85)
        ax.spines[["top", "right"]].set_visible(False)
        if days is not None:
            add_year_labels(ax, COL_START, COL_END, N_YEARS, DAYS_PER_YEAR)

    for idx in range(n_paths, nrows * ncols):
        axes2[idx // ncols][idx % ncols].set_visible(False)

    fig2.tight_layout()
    fig2.savefig(FIG_DIR / f"fig2_absolute_c_pathway_lines{NORM_SUFFIX}.png", dpi=150, bbox_inches="tight")
    plt.close(fig2); print("Saved fig2")

# ────────────────────────────────────────────────────────────────────────────
# Fig 3: Pairwise Δ absolute C cost per pathway (adjacent scenario pairs)
# ────────────────────────────────────────────────────────────────────────────
pairs = [(scen_labels[i], scen_labels[i + 1]) for i in range(N_SCEN - 1)]
path_color_map = dict(zip(PATHWAY_ORDER, PATHWAY_COLORS))

for a_label, b_label in pairs:
    a = scen_data[a_label]
    b = scen_data[b_label]
    delta_paths = [name for name in PATHWAY_ORDER
                   if a.get(name) is not None and b.get(name) is not None]
    if not delta_paths:
        print(f"Skipped fig3 pair {a_label}/{b_label} — no matching pathway data")
        continue

    fig3, ax3 = plt.subplots(figsize=(10, 4))
    fig3.suptitle(
        f"{b_label} − {a_label}: Δ Absolute C Cost Per Pathway\n"
        f"(community mean {NORM_LABEL}; positive = higher in warmer scenario)",
        fontsize=12, fontweight="bold",
    )
    days = None
    for name in delta_paths:
        delta = smooth(b[name] - a[name].reindex(b[name].index))
        days  = delta.index
        ax3.plot(days, delta.values, color=path_color_map[name], lw=1.8, label=name)

    ax3.axhline(0, color="black", lw=0.8, ls="--")
    ax3.set_xlabel("Day")
    ax3.set_ylabel(f"Δ C cost ({b_label} − {a_label}, model units {NORM_LABEL})")
    if days is not None:
        ax3.set_xlim(days.min(), days.max())
    ax3.legend(loc="upper center", bbox_to_anchor=(0.5, -0.22),
               ncols=3, fontsize=9, framealpha=0.85)
    ax3.spines[["top", "right"]].set_visible(False)
    if days is not None:
        add_year_labels(ax3, COL_START, COL_END, N_YEARS, DAYS_PER_YEAR)

    safe_a = a_label.replace("/", "-"); safe_b = b_label.replace("/", "-")
    fig3.tight_layout()
    fig3.savefig(FIG_DIR / f"fig3_delta_c_pathway_{safe_b}_minus_{safe_a}{NORM_SUFFIX}.png",
                 dpi=150, bbox_inches="tight")
    plt.close(fig3)
    print(f"Saved fig3  ({b_label} − {a_label})")

print(f"\nAll figures saved to: {FIG_DIR}")
