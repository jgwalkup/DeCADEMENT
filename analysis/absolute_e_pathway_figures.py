#!/usr/bin/env python3
"""
absolute_e_pathway_figures.py

Supports raw runs (--base / --scenario-dirs / --treatment-map) or pre-aggregated
summary CSVs (--summary-dir): pathways come from comm_*_E_per_microbe stems
(community mean per alive microbe), whose cross-community mean/SD match raw mode.

Absolute energy cost pathway breakdown for DEMENTpy.
Supports an arbitrary number of climate scenarios (ambient, moderate warming,
extreme warming, ...), each plotted as a separate line.

Fig 1  - Stacked absolute E cost pathways + E production overlay (one panel per scenario)
Fig 2  - Absolute E cost per pathway -- one line per scenario
Fig 3  - Delta absolute E cost per pathway relative to the first (ambient) scenario

Usage
-----
  python absolute_e_pathway_figures.py \
      --base /path/to/output \
      --treatment-map ambient_dir=Ambient \
                      mod_dir="Moderate Warming" \
                      ext_dir="Extreme Warming" \
      --fig-dir /path/to/figures
"""

import argparse
import warnings
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

warnings.filterwarnings("ignore")

sys.path.insert(0, str(Path(__file__).parent))
from scenario_utils import (
    resolve_summary_scenarios, load_summary, read_summary_meta, add_summary_arg,
    scenario_color_map, get_palette_from_args, _parse_treatment_map,
    load_series, pick_runs, resolve_runs, read_runtime,
    line_style, plot_order, legend_baseline_first,
)

# Raw per-run tables are read with scenario_utils.load_series, which joins split
# branched output (branch_output.json) into the full series, keeps the
# (taxon, type) index of Taxon_HSP_Series / Taxon_Osmo_Series aligned to the
# right day (type level dropped), and renumbers day columns 1..n.
load = load_series

# -- Paths ---------------------------------------------------------------------
BASE = Path(
    "/Users/jswalkup/Library/CloudStorage/"
    "GoogleDrive-jwalkup@uci.edu/My Drive/"
    "Decoupled-Energy-DEMENTpy-20250902/"
    "Decoupled-Energy-DEMENTpy/output"
)

# -- Run selection -------------------------------------------------------------
MIN_COLS = 100
SMOOTH   = 14

# -- Pathway colours -----------------------------------------------------------
C_TRANS = "#ff7f0e"
C_MAINT = "#9467bd"
C_HSP   = "#8c564b"
C_ENZ   = "#e377c2"
C_OSMO  = "#bcbd22"
C_PROD  = "#2ca02c"

PATHWAY_ORDER  = ["Transport", "Maintenance", "Enzymes", "Osmolytes", "HSP"]
PATHWAY_COLORS = [C_TRANS, C_MAINT, C_ENZ, C_OSMO, C_HSP]

plt.rcParams.update({
    "font.size": 10, "axes.titlesize": 11, "axes.labelsize": 10,
    "legend.fontsize": 9, "figure.dpi": 150,
})

# -- Helpers -------------------------------------------------------------------
def rep_mean(pool_list, key):
    frames = [a[key] for a in pool_list if a and a.get(key) is not None]
    if not frames:
        return None
    if len(frames) == 1:
        return frames[0]
    idx = frames[0].index
    arr = np.stack([s.reindex(idx).values for s in frames], axis=1)
    return pd.Series(np.nanmean(arr, axis=1), index=idx)

def rep_std(pool_list, key):
    frames = [a[key] for a in pool_list if a and a.get(key) is not None]
    if len(frames) < 2:
        return None
    idx = frames[0].index
    arr = np.stack([s.reindex(idx).values for s in frames], axis=1)
    return pd.Series(np.nanstd(arr, axis=1, ddof=1), index=idx)

def smooth(s):
    return s

def year_range(years_arg, n_years, days_per_year):
    if years_arg is None:
        return 1, n_years * days_per_year
    if len(years_arg) == 1:
        yr = years_arg[0]
        return (yr - 1) * days_per_year + 1, yr * days_per_year
    yr_start, yr_end = years_arg[0], years_arg[-1]
    return (yr_start - 1) * days_per_year + 1, yr_end * days_per_year

def filter_series(s, col_start, col_end):
    if s is None:
        return None
    return s.loc[col_start:col_end]

def add_year_labels(ax, col_start, col_end, n_years, days_per_year):
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

# -- Per-run loading -----------------------------------------------------------
def compute_run(run_dir):
    EP = load(run_dir, "E_Prod_Taxon")
    ET = load(run_dir, "Transport_E_Taxon")
    EM = load(run_dir, "Maint_E_Taxon")
    EH = load(run_dir, "HSP_E_Taxon")
    EE = load(run_dir, "Enzyme_E_Taxon")
    EO = load(run_dir, "Osmo_E_Taxon")

    if all(x is None for x in [ET, EM, EH, EE, EO]):
        return None

    def _comm(df):
        return df.mean(axis=0, skipna=True) if df is not None else None

    ep = _comm(EP)
    et = _comm(ET)
    em = _comm(EM)
    eh = _comm(EH)
    ee = _comm(EE)
    eo = _comm(EO)

    pathway_series = [s for s in [et, em, eh, ee, eo] if s is not None]
    if pathway_series:
        idx = pathway_series[0].index
        total_cost = pd.Series(
            np.nansum(np.stack([s.reindex(idx).values for s in pathway_series], axis=1), axis=1),
            index=idx,
        )
    else:
        total_cost = None

    net_gain = (ep - total_cost) if (ep is not None and total_cost is not None) else None

    return dict(
        EP=ep, Transport=et, Maintenance=em, HSP=eh, Enzymes=ee, Osmolytes=eo,
        total_cost=total_cost, net_gain=net_gain,
    )

# -- Parse CLI -----------------------------------------------------------------
_p = argparse.ArgumentParser(
    description="Absolute E pathway figures for DeCA-DEMENT",
    formatter_class=argparse.RawDescriptionHelpFormatter)
_p.add_argument("--base",            default=None,
                help="Base directory containing scenario subdirectories")
_p.add_argument("--scenario-dirs",   nargs="+", default=None, metavar="DIR",
                help="Explicit scenario directories (ordered: ambient first)")
_p.add_argument("--scenario-labels", nargs="+", default=None, metavar="LABEL",
                help="Display labels matching --scenario-dirs order")
_p.add_argument("--treatment-map",   nargs="+", default=None, metavar="DIR=LABEL",
                help="Map directory stems to labels (DIR_STEM=LABEL). "
                     "Order determines scenario order.")
_p.add_argument("--palette",         default="temperature", metavar="PALETTE",
                help="Colour palette: 'temperature' (default, blue->orange->red) "
                     "or 'discrete' (fixed categorical colours).")
_p.add_argument("--per-biomass",     action="store_true",
                help="Normalise all pathway series by total microbial biomass C "
                     "(system_microbial_biomass_C) instead of alive microbe count. "
                     "Requires --summary-dir.")
_p.add_argument("--fig-dir",         default=None,
                help="Output directory for figures")
_p.add_argument("--years",           nargs="+", type=int, default=None, metavar="YEAR",
                help="Year(s) to display. Single = that year; two = inclusive range.")
add_summary_arg(_p)
_args = _p.parse_args()

# -- Resolve scenario directories and labels -----------------------------------
def _parse_tmap(entries):
    if not entries:
        return {}
    result = {}
    for entry in entries:
        if "=" not in entry:
            print(f"  WARNING: --treatment-map entry has no '=': {entry!r}  (ignored)")
            continue
        stem, label = entry.split("=", 1)
        result[stem.strip()] = label.strip()
    return result

KEYS = ["EP", "Transport", "Maintenance", "HSP", "Enzymes", "Osmolytes",
        "total_cost", "net_gain"]

# key → export stem.  Community-mean-per-microbe E series whose cross-community
# mean/SD match raw-mode compute_run exactly.  total_cost / net_gain are
# derived-but-unplotted (vestigial in KEYS) → None.
_SUMMARY_STEMS = {
    "EP":          "comm_production_E_per_microbe",
    "Transport":   "comm_transport_E_per_microbe",
    "Maintenance": "comm_maintenance_E_per_microbe",
    "HSP":         "comm_hsp_E_per_microbe",
    "Enzymes":     "comm_enzyme_E_per_microbe",
    "Osmolytes":   "comm_osmolyte_E_per_microbe",
    "total_cost":  None,
    "net_gain":    None,
}
_MBC_STEM = "system_microbial_biomass_C"

if _args.summary_dir:
    # ── Summary mode (all discovered scenarios → full gradient / N-scenario) ────
    summary_root = Path(_args.summary_dir)
    _tmap = _parse_treatment_map(getattr(_args, "treatment_map", None))
    scenario_dirs_labels = resolve_summary_scenarios(summary_root, treatment_map=_tmap)
    if not scenario_dirs_labels:
        print(f"No summary scenarios found under {summary_root}.")
        raise SystemExit(1)

    FIG_DIR = (Path(_args.fig_dir) if _args.fig_dir
               else summary_root / "absolute_e_pathway_figures")
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    print(f"Summary mode — scenarios ({len(scenario_dirs_labels)}):")
    scenarios = []
    mbc_data  = {}
    for i, (label, sdir) in enumerate(scenario_dirs_labels):
        print(f"  {label}: {sdir}")
        pairs = {k: (load_summary(sdir, _SUMMARY_STEMS[k]) if _SUMMARY_STEMS[k]
                     else (None, None)) for k in KEYS}
        dat   = {k: pairs[k][0] for k in KEYS}
        sd    = {k: pairs[k][1] for k in KEYS}
        scenarios.append((label, None, dat, sd))
        _mbc_mean, _ = load_summary(sdir, _MBC_STEM)
        mbc_data[label] = _mbc_mean

    scen_labels = [label for label, _, _, _ in scenarios]
    COLOR = scenario_color_map(scen_labels, palette=get_palette_from_args(_args))
    scenarios = [(label, COLOR[label], dat, sd) for label, _, dat, sd in scenarios]

    _meta = read_summary_meta(scenario_dirs_labels[0][1])
    N_YEARS       = _meta["pulse"]
    DAYS_PER_YEAR = int(_meta["end_time"])
    if _args.years:
        COL_START, COL_END = year_range(_args.years, N_YEARS, DAYS_PER_YEAR)
    elif _meta["col_start"] is not None:
        COL_START, COL_END = _meta["col_start"], _meta["col_end"]
    else:
        _idx = [c for _, _, d, _ in scenarios for v in d.values()
                if v is not None for c in getattr(v, "index", [])]
        COL_START, COL_END = (min(_idx), max(_idx)) if _idx else (1, DAYS_PER_YEAR)
    print(f"Simulation: {N_YEARS} year(s), {DAYS_PER_YEAR} days/year  |  "
          f"Displaying: day {COL_START}-{COL_END}")

else:
    # ── Raw-run mode (original behaviour) ───────────────────────────────────────
    tmap = _parse_tmap(_args.treatment_map)

    if _args.scenario_dirs:
        _paths = [Path(p) for p in _args.scenario_dirs]
        _raw_labels = (_args.scenario_labels
                       if _args.scenario_labels and len(_args.scenario_labels) == len(_paths)
                       else [p.stem for p in _paths])
        scenario_dirs_labels = [(tmap.get(p.stem, lab), p)
                                for p, lab in zip(_paths, _raw_labels)]
    elif tmap:
        base = Path(_args.base) if _args.base else BASE
        scenario_dirs_labels = []
        for stem, label in tmap.items():
            sdir = base / stem
            if not sdir.exists():
                print(f"  WARNING: scenario directory not found: {sdir}")
                continue
            scenario_dirs_labels.append((label, sdir))
    else:
        base = Path(_args.base) if _args.base else BASE
        found = sorted([d for d in base.iterdir() if d.is_dir()])
        pos_labels = _args.scenario_labels or []
        scenario_dirs_labels = [(pos_labels[i] if i < len(pos_labels) else d.stem, d)
                                for i, d in enumerate(found)]

    FIG_DIR = (Path(_args.fig_dir) if _args.fig_dir else
               (Path(_args.base) if _args.base else BASE) / "absolute_e_pathway_figures")
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    scenarios = []   # list of (label, color, dat_dict, sd_dict)
    all_run_dirs = []
    mbc_data = {}
    for i, (label, sdir) in enumerate(scenario_dirs_labels):
        run_dirs = pick_runs(sdir)
        all_run_dirs.extend(run_dirs)
        print(f"  {label}: {len(run_dirs)} run(s) in {sdir.name}")
        runs = [compute_run(d) for d in run_dirs]
        dat  = {k: rep_mean(runs, k) for k in KEYS}
        sd   = {k: rep_std(runs,  k) for k in KEYS}
        scenarios.append((label, None, dat, sd))
        mbc_data[label] = None   # not available in raw mode

    scen_labels_raw = [label for label, _, _, _ in scenarios]
    COLOR_RAW = scenario_color_map(scen_labels_raw, palette=get_palette_from_args(_args))
    scenarios = [(label, COLOR_RAW[label], dat, sd) for label, _, dat, sd in scenarios]

    _rt = read_runtime(all_run_dirs[0] if all_run_dirs else None)
    N_YEARS       = _rt["pulse"]
    DAYS_PER_YEAR = int(_rt["end_time"])
    COL_START, COL_END = year_range(_args.years, N_YEARS, DAYS_PER_YEAR)
    print(f"Simulation: {N_YEARS} year(s), {DAYS_PER_YEAR} days/year  |  "
          f"Displaying: day {COL_START}-{COL_END}")

if not scenarios:
    print("ERROR: no scenarios found.")
    raise SystemExit(1)

if all(dat["Transport"] is None for _, _, dat, _ in scenarios):
    print("ERROR: No pathway CSVs found in any run.")
    raise SystemExit(1)

scenarios = [
    (label, color,
     {k: filter_series(v, COL_START, COL_END) for k, v in dat.items()},
     {k: filter_series(v, COL_START, COL_END) for k, v in sd.items()})
    for label, color, dat, sd in scenarios
]

# ── Optional per-biomass normalisation ───────────────────────────────────────
PER_BIOMASS = getattr(_args, "per_biomass", False)
if PER_BIOMASS:
    _missing_mbc = [l for l, _, _, _ in scenarios if mbc_data.get(l) is None]
    if _missing_mbc:
        print(f"WARNING: --per-biomass requested but MBC series not available for: "
              f"{_missing_mbc}. Per-biomass normalisation requires --summary-dir.")
        PER_BIOMASS = False
    else:
        print("Applying per-biomass normalisation (dividing by system_microbial_biomass_C)")
        def _div_by_mbc(s, mbc):
            if s is None:
                return None
            _mbc_r = mbc.reindex(s.index)
            return s.div(_mbc_r).replace([float("inf"), float("-inf")], float("nan"))
        scenarios = [
            (label, color,
             {k: _div_by_mbc(v, filter_series(mbc_data[label], COL_START, COL_END))
              for k, v in dat.items()},
             {k: _div_by_mbc(v, filter_series(mbc_data[label], COL_START, COL_END))
              for k, v in sd.items()})
            for label, color, dat, sd in scenarios
        ]

NORM_LABEL  = "per unit biomass C" if PER_BIOMASS else "per alive microbe"
NORM_SUFFIX = "_per_biomass" if PER_BIOMASS else ""

n_scen = len(scenarios)
# =============================================================================
# Fig 1: Stacked absolute E cost pathways + E production overlay
# =============================================================================
fig1, axes = plt.subplots(1, n_scen, figsize=(7 * n_scen, 5), sharey=True)
if n_scen == 1:
    axes = [axes]
fig1.suptitle(
    "Absolute Energy Cost by Pathway\n"
    f"(community mean {NORM_LABEL}; stacked area = total cost; line = E production)",
    fontsize=12, fontweight="bold",
)

for ax, (label, _scen_color, dat, dat_sd) in zip(axes, scenarios):
    ax.set_title(label)

    stack_vals   = []
    stack_labels = []
    stack_colors = []
    for name, color in zip(PATHWAY_ORDER, PATHWAY_COLORS):
        s = dat[name]
        if s is not None:
            stack_vals.append(smooth(s).values)
            stack_labels.append(name)
            stack_colors.append(color)

    if not stack_vals:
        ax.text(0.5, 0.5, "No pathway data", transform=ax.transAxes,
                ha="center", va="center")
        continue

    days = dat[[k for k in KEYS if dat[k] is not None][0]].index
    ax.stackplot(days, stack_vals,
                 labels=stack_labels, colors=stack_colors, alpha=0.85)

    if dat["EP"] is not None:
        ep_s = smooth(dat["EP"])
        ax.plot(days, ep_s.values, color="black", lw=2, ls="--",
                label="E production", zorder=5)
        if dat_sd["EP"] is not None:
            ep_sd = smooth(dat_sd["EP"])
            ax.fill_between(days,
                            (ep_s - ep_sd).clip(0), ep_s + ep_sd,
                            color="black", alpha=0.12, linewidth=0)

    ax.set_xlabel("Day")
    ax.set_ylabel(f"Energy (model units {NORM_LABEL})")
    ax.set_xlim(days.min(), days.max())
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.20),
              ncols=3, fontsize=8, framealpha=0.85)
    add_year_labels(ax, COL_START, COL_END, N_YEARS, DAYS_PER_YEAR)

fig1.tight_layout()
fig1.savefig(FIG_DIR / f"fig1_absolute_e_pathway_stacked{NORM_SUFFIX}.png", dpi=150, bbox_inches="tight")
plt.close(fig1)
print("Saved fig1")

# =============================================================================
# Fig 2: Individual pathway lines -- one line per scenario
# =============================================================================
available_pathways = [name for name in PATHWAY_ORDER
                      if any(dat[name] is not None for _, _, dat, _ in scenarios)]
n_paths = len(available_pathways)

if n_paths == 0:
    print("Skipped fig2 -- no pathway data")
else:
    ncols = min(n_paths, 3)
    nrows = (n_paths + ncols - 1) // ncols
    fig2, axes2 = plt.subplots(nrows, ncols,
                               figsize=(5 * ncols, 4 * nrows),
                               squeeze=False)
    scen_labels = [label for label, _, _, _ in scenarios]
    fig2.suptitle(
        "Absolute E Cost Per Pathway -- " + " vs ".join(scen_labels) + "\n"
        f"(community mean {NORM_LABEL})",
        fontsize=12, fontweight="bold",
    )

    path_color_map = dict(zip(PATHWAY_ORDER, PATHWAY_COLORS))

    for idx, name in enumerate(available_pathways):
        ax = axes2[idx // ncols][idx % ncols]
        pcolor = path_color_map[name]
        days   = None

        _labels = [lab for lab, _, _, _ in scenarios]
        for label, scen_color, dat, dat_sd in sorted(scenarios, key=lambda t: plot_order(_labels).index(t[0])):
            s = dat[name]                              # baseline last, on top
            if s is None:
                continue
            days = s.index
            s_sm = smooth(s)
            ax.plot(days, s_sm.values, color=scen_color, lw=1.8,
                    label=label, **line_style(label, _labels))
            if dat_sd[name] is not None:
                sd_sm = smooth(dat_sd[name])
                ax.fill_between(days,
                                (s_sm - sd_sm).clip(0), s_sm + sd_sm,
                                color=scen_color, alpha=0.18, linewidth=0)

        ax.set_title(name, color=pcolor, fontweight="bold")
        ax.set_xlabel("Day")
        ax.set_ylabel(f"E cost (model units {NORM_LABEL})")
        if days is not None:
            ax.set_xlim(days.min(), days.max())
        legend_baseline_first(ax, _labels, loc="upper center", bbox_to_anchor=(0.5, -0.20),
                              ncols=min(n_scen, 3), fontsize=8, framealpha=0.85)
        ax.spines[["top", "right"]].set_visible(False)
        if days is not None:
            add_year_labels(ax, COL_START, COL_END, N_YEARS, DAYS_PER_YEAR)

    for idx in range(n_paths, nrows * ncols):
        axes2[idx // ncols][idx % ncols].set_visible(False)

    fig2.tight_layout()
    fig2.savefig(FIG_DIR / f"fig2_absolute_e_pathway_lines{NORM_SUFFIX}.png", dpi=150, bbox_inches="tight")
    plt.close(fig2)
    print("Saved fig2")

# =============================================================================
# Fig 3: Delta absolute E cost per pathway -- each scenario relative to ambient
# =============================================================================
ambient_label, _ac, ambient_dat, _as_dict = scenarios[0]
comparison_scenarios = scenarios[1:]

if not comparison_scenarios:
    print("Skipped fig3 -- need at least two scenarios for delta plot")
else:
    delta_paths = [name for name in PATHWAY_ORDER
                   if ambient_dat[name] is not None and
                   any(dat[name] is not None for _, _, dat, _ in comparison_scenarios)]

    if not delta_paths:
        print("Skipped fig3 -- no matching pathways between ambient and other scenarios")
    else:
        path_color_map = dict(zip(PATHWAY_ORDER, PATHWAY_COLORS))

        fig3, ax3 = plt.subplots(figsize=(10, 4))
        fig3.suptitle(
            f"Delta Absolute E Cost Per Pathway Relative to {ambient_label}\n"
            f"(community mean {NORM_LABEL}; positive = higher than ambient)",
            fontsize=12, fontweight="bold",
        )

        days = None
        for name in delta_paths:
            pcolor = path_color_map[name]
            for si, (label, scen_color, dat, _sd) in enumerate(comparison_scenarios):
                if dat[name] is None:
                    continue
                ref   = ambient_dat[name].reindex(dat[name].index)
                delta = smooth(dat[name] - ref)
                days  = delta.index
                ax3.plot(days, delta.values, color=pcolor, lw=1.8,
                         ls=["-", "--", ":", "-."][si % 4],      # line style tells the scenarios apart
                         label=f"{name} ({label})")

        ax3.axhline(0, color="black", lw=0.8, ls="--")
        ax3.set_xlabel("Day")
        ax3.set_ylabel(
            f"Delta E cost (scenario - {ambient_label}, model units {NORM_LABEL})")
        if days is not None:
            ax3.set_xlim(days.min(), days.max())
        ax3.legend(loc="upper center", bbox_to_anchor=(0.5, -0.22),
                   ncols=min(len(delta_paths), 3), fontsize=8, framealpha=0.85)
        ax3.spines[["top", "right"]].set_visible(False)
        if days is not None:
            add_year_labels(ax3, COL_START, COL_END, N_YEARS, DAYS_PER_YEAR)

        fig3.tight_layout()
        fig3.savefig(FIG_DIR / f"fig3_absolute_e_pathway_delta{NORM_SUFFIX}.png",
                     dpi=150, bbox_inches="tight")
        plt.close(fig3)
        print("Saved fig3")

print(f"\nAll figures saved to: {FIG_DIR}")
