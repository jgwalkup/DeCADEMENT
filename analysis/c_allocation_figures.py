#!/usr/bin/env python3
"""
c_allocation_figures.py

Proportional and absolute C-allocation figures for DEMENTpy.

Fig 1  – Proportional allocation time series  (energy vs. biosynthesis)
Fig 2  – Absolute allocation time series      (per alive microbe, community mean)
Fig 3  – Per-taxon mean allocation bar charts  (energy vs. biosynthesis, sorted)
Fig 4  – Within-biosynthesis breakdown         (HSP / enzyme / osmolyte / net growth)
Fig 5  – Per-taxon allocation heatmaps         (energy fraction + HSP pool fraction)
Fig 6  – Warm − Control Δ allocation           (per taxon)
"""

import argparse
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
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
CTRL_RUNS: list[str] = []
WARM_RUNS: list[str] = []
MIN_COLS  = 100

# ── Style ─────────────────────────────────────────────────────────────────────
C_ENERGY = "#d62728"
C_BIOSYN = "#1f77b4"
C_HSP    = "#2ca02c"
C_ENZ    = "#ff7f0e"
C_OSMO   = "#9467bd"
C_GROW   = "#aec7e8"
C_TRANS  = "#17becf"   # cyan/teal — transporter production C
C_CTRL   = "#1f77b4"
C_WARM   = "#d62728"

plt.rcParams.update({
    "font.size": 10, "axes.titlesize": 11, "axes.labelsize": 10,
    "legend.fontsize": 9, "figure.dpi": 150,
})

# ── Helpers ───────────────────────────────────────────────────────────────────
def safe_div(a, b):
    with np.errstate(divide="ignore", invalid="ignore"):
        r = a / b
    r.replace([np.inf, -np.inf], np.nan, inplace=True)
    return r

def compute_alloc(run_dir: Path):
    C   = load(run_dir, "C_Uptake_Taxon")
    BC  = load(run_dir, "BC_Uptake_Taxon")
    MIC = load(run_dir, "MicrobesSeries")
    HSP = load(run_dir, "Taxon_HSP_Series")
    ENZ = load(run_dir, "EnzymesProduction")
    H_C = load(run_dir, "HSP_C_Taxon")
    E_C = load(run_dir, "Enzyme_C_Taxon")
    O_C = load(run_dir, "Osmo_C_Taxon")

    if C is None or BC is None:
        return None

    energy_C = (C - BC).clip(lower=0)
    biosyn_C  = BC.clip(lower=0)
    energy_frac = safe_div(energy_C, C).clip(0, 1)
    biosyn_frac = safe_div(biosyn_C, C).clip(0, 1)

    hsp_pool_frac = None
    if HSP is not None and MIC is not None:
        day_cols = [c for c in HSP.columns if isinstance(c, int)]
        hv = HSP[day_cols].astype(float) if day_cols else HSP.astype(float)
        common = sorted(set(hv.columns) & set(MIC.columns))
        if common:
            hsp_pool_frac = safe_div(hv[common], hv[common] + MIC[common]).clip(0, 1)

    hsp_bf = safe_div(H_C, BC).clip(0, 1) if H_C is not None else None
    enz_bf = safe_div(E_C, BC).clip(0, 1) if E_C is not None else None
    osm_bf = safe_div(O_C, BC).clip(0, 1) if O_C is not None else None
    T_C    = load(run_dir, "Transport_C_Taxon")
    trans_bf = safe_div(T_C, BC).clip(0, 1) if T_C is not None else None

    enz_comm = ENZ.sum(axis=0) if ENZ is not None else None

    return dict(
        C=C, BC=BC, MIC=MIC,
        energy_C=energy_C, biosyn_C=biosyn_C,
        energy_frac=energy_frac, biosyn_frac=biosyn_frac,
        hsp_pool_frac=hsp_pool_frac,
        hsp_bf=hsp_bf, enz_bf=enz_bf, osm_bf=osm_bf, trans_bf=trans_bf,
        enz_comm=enz_comm,
    )

def rep_mean(allocs, key):
    frames = [a[key] for a in allocs if a and a.get(key) is not None]
    if not frames:
        return None
    # Handle Series (e.g. enz_comm) vs DataFrame
    if isinstance(frames[0], pd.Series):
        try:
            arr = np.stack([f.reindex(frames[0].index).values for f in frames], axis=1)
            return pd.Series(np.nanmean(arr, axis=1), index=frames[0].index)
        except Exception:
            return frames[0]
    idx, cols = frames[0].index, frames[0].columns
    try:
        arr = np.stack([f.reindex(index=idx, columns=cols).values for f in frames], axis=2)
        return pd.DataFrame(np.nanmean(arr, axis=2), index=idx, columns=cols)
    except Exception:
        return frames[0]

def comm_ts(df):
    if df is None:
        return None, None
    mu  = df.mean(axis=0, skipna=True)
    sig = df.std(axis=0,  skipna=True)
    return mu, sig

def taxon_mean(df):
    if df is None:
        return None
    return df.mean(axis=1, skipna=True)

def subsample(df, n=200):
    if df is None or df.shape[1] <= n:
        return df
    idx = np.linspace(0, df.shape[1] - 1, n, dtype=int)
    return df.iloc[:, idx]

# ── Year / pulse helpers ───────────────────────────────────────────────────────
def year_range(years_arg, n_years: int, days_per_year: int) -> tuple[int, int]:
    """Convert --years argument to (col_start, col_end) inclusive day range."""
    if years_arg is None:
        return 1, n_years * days_per_year
    if len(years_arg) == 1:
        yr = years_arg[0]
        return (yr - 1) * days_per_year + 1, yr * days_per_year
    yr_start, yr_end = years_arg[0], years_arg[-1]
    return (yr_start - 1) * days_per_year + 1, yr_end * days_per_year

def filter_cols(df, col_start: int, col_end: int):
    """Slice a DataFrame (rows × days) or Series (days,) to [col_start, col_end]."""
    if df is None:
        return None
    if isinstance(df, pd.Series):
        return df.loc[col_start:col_end]
    cols = [c for c in df.columns if col_start <= c <= col_end]
    return df[cols] if cols else df

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

def add_year_labels_heatmap(ax, mat_ncols: int, col_start: int, col_end: int,
                             n_years: int, days_per_year: int) -> None:
    """Draw year-boundary markers on a subsampled imshow axes."""
    total = col_end - col_start + 1
    for yr in range(1, n_years + 1):
        boundary = yr * days_per_year
        if col_start < boundary < col_end:
            frac  = (boundary - col_start) / total
            x_pos = frac * (mat_ncols - 1)
            ax.axvline(x_pos, color="white", lw=0.8, ls=":", alpha=0.6, zorder=5)

# ── Parse CLI ──────────────────────────────────────────────────────────────────
import argparse as _ap
_p = _ap.ArgumentParser(description="c_allocation_figures.py (multi-scenario)",
                        formatter_class=_ap.RawDescriptionHelpFormatter)
add_scenario_args(_p)
add_summary_arg(_p)
_args = _p.parse_args()

from pathlib import Path as _Path

KEYS = ["C", "BC", "MIC", "energy_C", "biosyn_C", "energy_frac", "biosyn_frac",
        "hsp_pool_frac", "hsp_bf", "enz_bf", "osm_bf", "trans_bf", "enz_comm"]

# key → export_summary_data.py stem.  enz_comm is vestigial (never plotted) → None.
_SUMMARY_STEMS = {
    "C":            "c_uptake_taxon",
    "BC":           "c_biosynthesis_taxon",
    "MIC":          "microbes_by_taxon",
    "energy_C":     "c_energy_taxon",
    "biosyn_C":     "c_biosynthesis_taxon",
    "energy_frac":  "c_energy_fraction_taxon",
    "biosyn_frac":  "c_biosyn_fraction_taxon",
    "hsp_pool_frac":"hsp_pool_fraction",
    "hsp_bf":       "c_hsp_biosyn_frac",
    "enz_bf":       "c_enzyme_biosyn_frac",
    "osm_bf":       "c_osmolyte_biosyn_frac",
    "trans_bf":     "c_transport_biosyn_frac",
    "enz_comm":     None,
}

if _args.summary_dir:
    # ── Summary mode (all discovered scenarios → full gradient / N-scenario) ────
    summary_root = _Path(_args.summary_dir)
    _tmap = _parse_treatment_map(getattr(_args, "treatment_map", None))
    scenario_list = resolve_summary_scenarios(summary_root, treatment_map=_tmap)
    if not scenario_list:
        print(f"No summary scenarios found under {summary_root}.")
        import sys; sys.exit(1)
    scen_labels = [s[0] for s in scenario_list]
    N_SCEN      = len(scenario_list)
    COLOR       = scenario_color_map(scen_labels, palette=get_palette_from_args(_args))

    base_path = summary_root
    FIG_DIR   = _Path(_args.fig_dir) if _args.fig_dir else base_path / "c_allocation_figures"
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    print(f"Summary mode — scenarios ({N_SCEN}):")
    scen_data = {}
    for _label, _sdir in scenario_list:
        print(f"  {_label}: {_sdir}")
        scen_data[_label]         = {k: (load_summary(_sdir, _SUMMARY_STEMS[k])[0]
                                         if _SUMMARY_STEMS[k] else None) for k in KEYS}
        scen_data[_label + "_sd"] = {k: (load_summary(_sdir, _SUMMARY_STEMS[k])[1]
                                         if _SUMMARY_STEMS[k] else None) for k in KEYS}

    _meta = read_summary_meta(scenario_list[0][1])
    N_YEARS       = _meta["pulse"]
    DAYS_PER_YEAR = int(_meta["end_time"])
    GRIDSIZE      = _meta["gridsize"]
    # The summary CSVs hold the full exported window; honor --years by slicing to
    # the requested year(s) at plot time.  Without --years, show the full range
    # the export recorded (or the data range if meta lacks it).
    if _args.years:
        COL_START, COL_END = year_range(_args.years, N_YEARS, DAYS_PER_YEAR)
    elif _meta["col_start"] is not None:
        COL_START, COL_END = _meta["col_start"], _meta["col_end"]
    else:
        _cols = [c for v in scen_data[scen_labels[0]].values()
                 if v is not None and getattr(v, "ndim", 0) == 2
                 for c in v.columns if isinstance(c, int)]
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
    FIG_DIR   = _Path(_args.fig_dir) if _args.fig_dir else base_path / "c_allocation_figures"
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    print(f"Scenarios ({N_SCEN}):")
    scen_data = {}
    for _label, _run_dirs in scenario_list:
        print(f"  {_label}: {len(_run_dirs)} run(s)")
        _runs = [compute_alloc(d) for d in _run_dirs]
        scen_data[_label]         = {k: rep_mean(_runs, k) for k in KEYS}
        scen_data[_label + "_sd"] = {k: rep_std(_runs,  k) for k in KEYS}

    _first_run = scenario_list[0][1][0] if scenario_list[0][1] else None
    _rt = read_runtime(_first_run)
    N_YEARS       = _rt["pulse"]
    DAYS_PER_YEAR = int(_rt["end_time"])
    COL_START, COL_END = year_range(_args.years, N_YEARS, DAYS_PER_YEAR)
    print(f"Simulation: {N_YEARS} year(s), {DAYS_PER_YEAR} days/year  |  "
          f"Displaying: day {COL_START}-{COL_END}")
    GRIDSIZE = get_scenario_gridsize(scenario_list)

for _label in scen_labels:
    scen_data[_label]         = {k: filter_cols(v, COL_START, COL_END)
                                   for k, v in scen_data[_label].items()}
    scen_data[_label + "_sd"] = {k: filter_cols(v, COL_START, COL_END)
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
# Divides absolute pool/pathway quantities by grid cells so results are
# expressed per grid cell. Ratio/fraction/efficiency series are excluded
# (they are already normalised by biomass, uptake, or total energy).
_NORM_SKIP = {"biosyn_frac", "energy_frac", "enz_bf", "hsp_bf", "hsp_pool_frac", "osm_bf", "trans_bf"}
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
# Fig 1: Proportional C allocation time series — one panel per scenario
# ────────────────────────────────────────────────────────────────────────────
fig1, axes = plt.subplots(1, N_SCEN, figsize=(7 * N_SCEN, 5), sharey=True)
if N_SCEN == 1: axes = [axes]
fig1.suptitle("Proportional C Allocation — Energy vs. Biosynthesis\n"
              "(community mean across alive taxa)",
              fontsize=12, fontweight="bold")

for ax, label in zip(axes, scen_labels):
    dat = scen_data[label]
    ef_mu, ef_sd = comm_ts(dat["energy_frac"])
    bf_mu, bf_sd = comm_ts(dat["biosyn_frac"])
    if ef_mu is None:
        ax.set_title(f"{label} — no data"); continue
    days = ef_mu.index
    ax.stackplot(days, ef_mu, bf_mu, colors=[C_ENERGY, C_BIOSYN], alpha=0.82,
                 labels=["Energy substrate C", "Biosynthesis C"])
    if ef_sd is not None:
        ax.fill_between(days, (ef_mu - ef_sd).clip(0), (ef_mu + ef_sd).clip(0, 1),
                        color=C_ENERGY, alpha=0.25, linewidth=0)
    ax.set_xlabel("Day"); ax.set_ylabel("Fraction of total C uptake")
    ax.set_title(label); ax.set_xlim(days.min(), days.max()); ax.set_ylim(0, 1)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.20),
              ncols=2, framealpha=0.85, fontsize=8)
    ax.axhline(0.5, ls="--", lw=0.6, color="gray", alpha=0.45)
    add_year_labels(ax, COL_START, COL_END, N_YEARS, DAYS_PER_YEAR)

fig1.tight_layout()
fig1.savefig(FIG_DIR / "fig1_proportional_alloc_time_series.png", dpi=150, bbox_inches="tight")
plt.close(fig1); print("Saved fig1")

# ────────────────────────────────────────────────────────────────────────────
# Fig 2: Absolute C allocation time series — one panel per scenario
# ────────────────────────────────────────────────────────────────────────────
fig2, axes = plt.subplots(1, N_SCEN, figsize=(7 * N_SCEN, 5), sharey=True)
if N_SCEN == 1: axes = [axes]
fig2.suptitle("Absolute C Allocation Per Alive Microbe\n(community mean ± 1 SD)",
              fontsize=12, fontweight="bold")

for ax, label in zip(axes, scen_labels):
    dat = scen_data[label]
    ec_mu, ec_sd = comm_ts(dat["energy_C"])
    bc_mu, bc_sd = comm_ts(dat["biosyn_C"])
    if ec_mu is None:
        ax.set_title(f"{label} — no data"); continue
    days = ec_mu.index
    ax.stackplot(days, ec_mu, bc_mu, colors=[C_ENERGY, C_BIOSYN], alpha=0.82,
                 labels=["Energy substrate C", "Biosynthesis C"])
    if ec_sd is not None:
        ax.fill_between(days, (ec_mu - ec_sd).clip(0), ec_mu + ec_sd,
                        color=C_ENERGY, alpha=0.2, linewidth=0)
    ax.set_xlabel("Day"); ax.set_ylabel("C per alive microbe (model units)")
    ax.set_title(label); ax.set_xlim(days.min(), days.max())
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.20),
              ncols=2, framealpha=0.85, fontsize=8)
    add_year_labels(ax, COL_START, COL_END, N_YEARS, DAYS_PER_YEAR)

fig2.tight_layout()
fig2.savefig(FIG_DIR / "fig2_absolute_alloc_time_series.png", dpi=150, bbox_inches="tight")
plt.close(fig2); print("Saved fig2")

# ────────────────────────────────────────────────────────────────────────────
# Fig 3: Per-taxon mean allocation bar chart — one panel per scenario
# (taxa sorted by first-scenario energy fraction)
# ────────────────────────────────────────────────────────────────────────────
_ref_ef = taxon_mean(scen_data[scen_labels[0]]["energy_frac"])
order = _ref_ef.sort_values(ascending=True).index if _ref_ef is not None else None

fig3, axes = plt.subplots(1, N_SCEN, figsize=(7 * N_SCEN, 9), sharey=True)
if N_SCEN == 1: axes = [axes]
fig3.suptitle("Mean Proportional C Allocation Per Taxon\n"
              f"(averaged over all alive days; taxa sorted by {scen_labels[0]} energy fraction)",
              fontsize=12, fontweight="bold")

for ax, label in zip(axes, scen_labels):
    dat = scen_data[label]
    if dat["energy_frac"] is None or order is None:
        ax.set_title(f"{label} — no data"); continue
    ef = taxon_mean(dat["energy_frac"]).reindex(order)
    bf = taxon_mean(dat["biosyn_frac"]).reindex(order)
    yticks = np.arange(len(order))
    ax.barh(yticks, ef,           color=C_ENERGY, alpha=0.85, label="Energy substrate C")
    ax.barh(yticks, bf, left=ef,  color=C_BIOSYN, alpha=0.85, label="Biosynthesis C")
    ax.set_yticks(yticks); ax.set_yticklabels(order, fontsize=7)
    ax.set_xlabel("Fraction of total C uptake"); ax.set_xlim(0, 1)
    ax.set_title(label); ax.axvline(0.5, ls="--", lw=0.8, color="gray", alpha=0.5)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.20),
              ncols=2, framealpha=0.85, fontsize=8)

fig3.tight_layout()
fig3.savefig(FIG_DIR / "fig3_taxon_allocation_bars.png", dpi=150, bbox_inches="tight")
plt.close(fig3); print("Saved fig3")

# ────────────────────────────────────────────────────────────────────────────
# Fig 4: Within-biosynthesis breakdown — one panel per scenario
# ────────────────────────────────────────────────────────────────────────────
has_new = any(scen_data[scen_labels[0]].get(k) is not None
              for k in ["hsp_bf", "enz_bf", "osm_bf", "trans_bf"])

fig4, axes = plt.subplots(1, N_SCEN, figsize=(7 * N_SCEN, 5), sharey=True)
if N_SCEN == 1: axes = [axes]

if has_new:
    fig4.suptitle("Within-Biosynthesis C Allocation Per Alive Microbe\n"
                  "(HSP + enzyme + osmolyte + transporter + net growth = biosynthesis C)",
                  fontsize=12, fontweight="bold")
    for ax, label in zip(axes, scen_labels):
        dat = scen_data[label]
        hsp_m, _ = comm_ts(dat["hsp_bf"])
        enz_m, _ = comm_ts(dat["enz_bf"])
        osm_m, _ = comm_ts(dat["osm_bf"])
        trans_m, _ = comm_ts(dat["trans_bf"])
        if hsp_m is None:
            ax.set_title(f"{label} — no data"); continue
        enz_m   = enz_m   if enz_m   is not None else pd.Series(0, index=hsp_m.index)
        osm_m   = osm_m   if osm_m   is not None else pd.Series(0, index=hsp_m.index)
        trans_m = trans_m if trans_m is not None else pd.Series(0, index=hsp_m.index)
        grow_m  = (1 - hsp_m - enz_m - osm_m - trans_m).clip(0)
        days = hsp_m.index
        ax.stackplot(days, hsp_m, enz_m, osm_m, trans_m, grow_m,
                     colors=[C_HSP, C_ENZ, C_OSMO, C_TRANS, C_GROW], alpha=0.85,
                     labels=["HSP", "Enzymes (exuded)", "Osmolytes",
                             "Transporter production", "Net growth biomass"])
        ax.set_xlabel("Day"); ax.set_ylabel("Fraction of biosynthesis C")
        ax.set_title(label); ax.set_xlim(days.min(), days.max()); ax.set_ylim(0, 1)
        ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.20),
                  ncols=3, fontsize=8, framealpha=0.85)
        add_year_labels(ax, COL_START, COL_END, N_YEARS, DAYS_PER_YEAR)
else:
    fig4.suptitle("HSP Pool Fraction of Total Cellular C\n"
                  "(proxy = HSP pool / (HSP pool + biomass C))",
                  fontsize=11, fontweight="bold")
    for ax, label in zip(axes, scen_labels):
        dat = scen_data[label]
        hp_m, hp_s = comm_ts(dat["hsp_pool_frac"])
        if hp_m is None:
            ax.set_title(f"{label} — no data"); continue
        days = hp_m.index
        ax.stackplot(days, hp_m, (1 - hp_m).clip(0),
                     colors=[C_HSP, C_BIOSYN], alpha=0.82,
                     labels=["HSP pool", "Non-HSP cellular C (biomass)"])
        if hp_s is not None:
            ax.fill_between(days, (hp_m - hp_s).clip(0), (hp_m + hp_s).clip(0, 1),
                            color=C_HSP, alpha=0.25, linewidth=0)
        ax.set_xlabel("Day"); ax.set_ylabel("Fraction of total cellular C")
        ax.set_title(label); ax.set_xlim(days.min(), days.max()); ax.set_ylim(0, 1)
        ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.20),
                  ncols=2, fontsize=8, framealpha=0.85)
        add_year_labels(ax, COL_START, COL_END, N_YEARS, DAYS_PER_YEAR)

fig4.tight_layout()
fig4.savefig(FIG_DIR / "fig4_within_biosynthesis_breakdown.png", dpi=150, bbox_inches="tight")
plt.close(fig4); print("Saved fig4")

# ────────────────────────────────────────────────────────────────────────────
# Fig 5: Per-taxon allocation heatmaps — N_SCEN columns × 2 rows
# ────────────────────────────────────────────────────────────────────────────
order_e = taxon_mean(scen_data[scen_labels[0]]["energy_frac"])
order_e = order_e.sort_values(ascending=False).index if order_e is not None else None

fig5, axes5 = plt.subplots(2, N_SCEN, figsize=(8 * N_SCEN, 11), squeeze=False)
fig5.suptitle(f"Per-Taxon C Allocation Heatmaps\n"
              f"(rows = taxa sorted by {scen_labels[0]} energy fraction; cols = days)",
              fontsize=12, fontweight="bold")
day_labels = np.linspace(COL_START, COL_END, 6, dtype=int)

if order_e is not None:
    for col, label in enumerate(scen_labels):
        dat = scen_data[label]
        xtick_pos = None

        for row, (key, cmap_name, cbar_label) in enumerate([
                ("energy_frac",   "RdYlBu_r", "Energy fraction"),
                ("hsp_pool_frac", "Greens",   "HSP pool fraction"),
        ]):
            ax  = axes5[row, col]
            df  = dat.get(key)
            ax.set_title(f"{label} — {cbar_label}", fontsize=10)
            if df is None:
                ax.axis("off"); continue
            df_sub = subsample(df)
            mat    = df_sub.reindex(order_e).fillna(np.nan).values
            vmax   = np.nanmax(mat) if np.nanmax(mat) > 0 else 1.0
            vmin   = 0
            im = ax.imshow(mat, aspect="auto", cmap=cmap_name,
                           vmin=vmin, vmax=vmax, interpolation="nearest")
            ax.set_yticks(np.arange(len(order_e)))
            ax.set_yticklabels(order_e, fontsize=5)
            ax.set_xlabel("Day (subsampled to 200)")
            ax.set_ylabel("Taxon")
            xtick_pos = np.linspace(0, mat.shape[1] - 1, 6, dtype=int)
            ax.set_xticks(xtick_pos); ax.set_xticklabels(day_labels)
            plt.colorbar(im, ax=ax, label=cbar_label, shrink=0.8)
            add_year_labels_heatmap(ax, mat.shape[1], COL_START, COL_END, N_YEARS, DAYS_PER_YEAR)

fig5.tight_layout()
fig5.savefig(FIG_DIR / "fig5_taxon_allocation_heatmaps.png", dpi=150, bbox_inches="tight")
plt.close(fig5); print("Saved fig5")

# ────────────────────────────────────────────────────────────────────────────
# Fig 6: Pairwise Δ allocation per taxon (adjacent scenario pairs)
# ────────────────────────────────────────────────────────────────────────────
pairs = [(scen_labels[i], scen_labels[i + 1]) for i in range(N_SCEN - 1)]

for a_label, b_label in pairs:
    a = scen_data[a_label]
    b = scen_data[b_label]
    if a["energy_frac"] is None or b["energy_frac"] is None:
        print(f"Skipped fig6 pair {a_label}/{b_label} — no energy_frac data")
        continue

    delta_ef = (taxon_mean(b["energy_frac"]) - taxon_mean(a["energy_frac"])).dropna()
    delta_hp = None
    if a.get("hsp_pool_frac") is not None and b.get("hsp_pool_frac") is not None:
        delta_hp = (taxon_mean(b["hsp_pool_frac"]) - taxon_mean(a["hsp_pool_frac"])).dropna()

    order6  = delta_ef.sort_values(ascending=True).index
    yticks6 = np.arange(len(order6))
    n_panels = 2 if delta_hp is not None else 1

    fig6, axs = plt.subplots(1, n_panels, figsize=(7 * n_panels, 9), sharey=True)
    if n_panels == 1: axs = [axs]
    fig6.suptitle(f"{b_label} − {a_label}: Difference in C Allocation Per Taxon\n"
                  "(positive = higher in warmer scenario)",
                  fontsize=12, fontweight="bold")

    c_pos = COLOR[b_label]
    c_neg = COLOR[a_label]

    vals6  = delta_ef.reindex(order6).fillna(0).values
    colors6 = [c_pos if v > 0 else c_neg for v in vals6]
    axs[0].barh(yticks6, vals6, color=colors6, alpha=0.85, edgecolor="white", linewidth=0.3)
    axs[0].axvline(0, color="black", lw=0.9)
    axs[0].set_yticks(yticks6); axs[0].set_yticklabels(order6, fontsize=7)
    axs[0].set_xlabel(f"Δ Energy fraction ({b_label} − {a_label})")
    axs[0].set_title("Δ Energy substrate fraction")

    if delta_hp is not None:
        vals_hp  = delta_hp.reindex(order6).fillna(0).values
        colors_hp = [c_pos if v > 0 else c_neg for v in vals_hp]
        axs[1].barh(yticks6, vals_hp, color=colors_hp, alpha=0.85, edgecolor="white", linewidth=0.3)
        axs[1].axvline(0, color="black", lw=0.9)
        axs[1].set_xlabel(f"Δ HSP pool fraction ({b_label} − {a_label})")
        axs[1].set_title("Δ HSP pool fraction of cellular C")

    legend_handles = [Patch(facecolor=c_pos, label=f"Higher in {b_label}"),
                      Patch(facecolor=c_neg, label=f"Higher in {a_label}")]
    axs[-1].legend(handles=legend_handles, loc="upper center",
                   bbox_to_anchor=(0.5, -0.20), ncols=2, framealpha=0.85, fontsize=8)
    safe_a = a_label.replace("/", "-"); safe_b = b_label.replace("/", "-")
    fig6.tight_layout()
    fig6.savefig(FIG_DIR / f"fig6_delta_alloc_{safe_b}_minus_{safe_a}.png",
                 dpi=150, bbox_inches="tight")
    plt.close(fig6)
    print(f"Saved fig6  ({b_label} − {a_label})")

print(f"\nAll figures saved to: {FIG_DIR}")
