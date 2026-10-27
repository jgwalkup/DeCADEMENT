#!/usr/bin/env python3
"""
energy_allocation_figures.py

Energy production and allocation figures for DEMENTpy.

Supports raw scenario runs (--base / --scenario-dirs) or pre-aggregated summary
CSVs from export_summary_data.py (--summary-dir, full gradient / N-scenario).
All keys map to e_*_taxon / e_*_frac stems whose definitions mirror compute_e.

Fig 1  – E production, total cost, net gain  (community time series)
Fig 2  – Proportional E allocation            (cost fraction vs. net-gain fraction)
Fig 3  – Per-taxon mean E allocation bar chart
Fig 4  – Within-cost pathway breakdown        (transport / maintenance / HSP / enzyme / osmolyte)
Fig 5  – Per-taxon E allocation heatmaps
Fig 6  – Warm − Control Δ E allocation        (per taxon)
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
    add_summary_arg, resolve_summary_scenarios, load_summary, read_summary_meta,
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
C_PROD   = "#2ca02c"
C_COST   = "#d62728"
C_NET    = "#1f77b4"
C_TRANS  = "#ff7f0e"
C_MAINT  = "#9467bd"
C_HSP    = "#8c564b"
C_ENZ    = "#e377c2"
C_OSMO   = "#bcbd22"
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

def compute_e(run_dir: Path):
    EP  = load(run_dir, "E_Prod_Taxon")
    EC  = load(run_dir, "Microbe_E_Cost")
    EG  = load(run_dir, "Microbe_E_Gain")
    ET  = load(run_dir, "Transport_E_Taxon")
    EM  = load(run_dir, "Maint_E_Taxon")
    EH  = load(run_dir, "HSP_E_Taxon")
    EE  = load(run_dir, "Enzyme_E_Taxon")
    EO  = load(run_dir, "Osmo_E_Taxon")

    if EP is None:
        return None

    total_cost = None
    if ET is not None and EM is not None:
        total_cost = ET.copy()
        for df in [EM, EH, EE, EO]:
            if df is not None:
                cc = sorted(set(total_cost.columns) & set(df.columns))
                total_cost[cc] = total_cost[cc].add(df[cc], fill_value=0)

    cost_frac = safe_div(total_cost, EP).clip(0, 1) if total_cost is not None else None
    net_frac  = (1 - cost_frac).clip(0, 1)          if cost_frac  is not None else None

    def _pf(x):
        if x is None or total_cost is None:
            return None
        cc = sorted(set(x.columns) & set(total_cost.columns))
        return safe_div(x[cc], total_cost[cc]).clip(0, 1)

    return dict(
        EP=EP, EC=EC, EG=EG,
        ET=ET, EM=EM, EH=EH, EE=EE, EO=EO,
        total_cost=total_cost,
        cost_frac=cost_frac, net_frac=net_frac,
        trans_pf=_pf(ET), maint_pf=_pf(EM),
        hsp_pf=_pf(EH),   enz_pf=_pf(EE), osm_pf=_pf(EO),
    )

def rep_mean(allocs, key):
    frames = [a[key] for a in allocs if a and a.get(key) is not None]
    if not frames:
        return None
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

def add_year_labels_heatmap(ax, mat_ncols: int, col_start: int, col_end: int,
                             n_years: int, days_per_year: int) -> None:
    total = col_end - col_start + 1
    for yr in range(1, n_years + 1):
        boundary = yr * days_per_year
        if col_start < boundary < col_end:
            frac  = (boundary - col_start) / total
            x_pos = frac * (mat_ncols - 1)
            ax.axvline(x_pos, color="white", lw=0.8, ls=":", alpha=0.6, zorder=5)

# ── Parse CLI ──────────────────────────────────────────────────────────────────
import argparse as _ap
_p = _ap.ArgumentParser(description="energy_allocation_figures.py (multi-scenario)",
                        formatter_class=_ap.RawDescriptionHelpFormatter)
add_scenario_args(_p)
add_summary_arg(_p)
_args = _p.parse_args()

from pathlib import Path as _Path

# All keys returned by compute_e() — used to build scen_data dicts
KEYS = ["EP", "EC", "EG", "ET", "EM", "EH", "EE", "EO",
        "total_cost", "cost_frac", "net_frac",
        "trans_pf", "maint_pf", "hsp_pf", "enz_pf", "osm_pf"]

# key → export_summary_data.py stem.  compute_run mirrors compute_e exactly
# (same total_cost sum, cost_frac = total_cost/EP, *_pf = part/total_cost).
_SUMMARY_STEMS = {
    "EP":         "e_production_taxon",
    "EC":         "e_cost_taxon",
    "EG":         "e_gain_taxon",
    "ET":         "e_transport_taxon",
    "EM":         "e_maintenance_taxon",
    "EH":         "e_hsp_taxon",
    "EE":         "e_enzyme_taxon",
    "EO":         "e_osmolyte_taxon",
    "total_cost": "e_total_cost_taxon",
    "cost_frac":  "e_cost_fraction_taxon",
    "net_frac":   "e_net_fraction",
    "trans_pf":   "e_transport_frac",
    "maint_pf":   "e_maintenance_frac",
    "hsp_pf":     "e_hsp_frac",
    "enz_pf":     "e_enzyme_frac",
    "osm_pf":     "e_osmolyte_frac",
}

if _args.summary_dir:
    # ── Summary mode (all discovered scenarios → full gradient / N-scenario) ────
    summary_root = _Path(_args.summary_dir)
    scenario_list = resolve_summary_scenarios(summary_root)
    if not scenario_list:
        print(f"No summary scenarios found under {summary_root}.")
        import sys; sys.exit(1)

    scen_labels = [s[0] for s in scenario_list]
    COLOR       = scenario_color_map(scen_labels)
    N_SCEN      = len(scenario_list)

    base_path = summary_root
    FIG_DIR   = _Path(_args.fig_dir) if _args.fig_dir else base_path / "energy_allocation_figures"
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    print(f"Summary mode — scenarios ({N_SCEN}):")
    scen_data = {}
    for _label, _sdir in scenario_list:
        print(f"  {_label}: {_sdir}")
        _pairs = {k: load_summary(_sdir, _SUMMARY_STEMS[k]) for k in KEYS}
        scen_data[_label]         = {k: _pairs[k][0] for k in KEYS}
        scen_data[_label + "_sd"] = {k: _pairs[k][1] for k in KEYS}

    _meta = read_summary_meta(scenario_list[0][1])
    N_YEARS       = _meta["pulse"]
    DAYS_PER_YEAR = int(_meta["end_time"])
    GRIDSIZE      = _meta["gridsize"]
    # Summary CSVs hold the full exported window; honor --years by slicing at
    # plot time, else fall back to the recorded window or the data range.
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
                                          sentinel="E_Prod_Taxon")
    if not scenario_list:
        print("No scenarios found. Use --base or --scenario-dirs.")
        import sys; sys.exit(1)

    scen_labels = [s[0] for s in scenario_list]
    COLOR       = scenario_color_map(scen_labels)
    N_SCEN      = len(scenario_list)

    base_path = _Path(_args.base) if _args.base else BASE
    FIG_DIR   = _Path(_args.fig_dir) if _args.fig_dir else base_path / "energy_allocation_figures"
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    print(f"Scenarios ({N_SCEN}):")
    scen_data = {}
    for _label, _run_dirs in scenario_list:
        print(f"  {_label}: {len(_run_dirs)} run(s)")
        _runs = [compute_e(d) for d in _run_dirs]
        scen_data[_label]         = {k: rep_mean(_runs, k) for k in KEYS}
        scen_data[_label + "_sd"] = {k: rep_std(_runs,  k) for k in KEYS}

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
_NORM_SKIP = {"cost_frac", "enz_pf", "hsp_pf", "maint_pf", "net_frac", "osm_pf", "trans_pf"}
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
# Fig 1: E production, cost, net gain — community time series
# ────────────────────────────────────────────────────────────────────────────
fig1, axes = plt.subplots(1, 2, figsize=(13, 5), sharey=True)
fig1.suptitle("Energy Production, Total Cost, and Net Gain Per Alive Microbe\n"
              "(community mean)",
              fontsize=12, fontweight="bold")

for ax, dat, title in zip(axes, [scen_data[l] for l in scen_labels[:2]], scen_labels[:2]):
    ep_mu, ep_sd = comm_ts(dat["EP"])
    tc_mu, tc_sd = comm_ts(dat["total_cost"])
    if ep_mu is None:
        ax.set_title(f"{title} — no data"); continue
    days = ep_mu.index
    ax.plot(days, ep_mu, color=C_PROD, lw=1.8, label="E production")
    ax.fill_between(days, ep_mu - ep_sd, ep_mu + ep_sd, color=C_PROD, alpha=0.2, linewidth=0)
    if tc_mu is not None:
        ax.plot(days, tc_mu, color=C_COST, lw=1.8, label="Total E cost")
        ax.fill_between(days, tc_mu, ep_mu, where=(ep_mu >= tc_mu),
                        color=C_NET, alpha=0.30, label="Net E gain", linewidth=0)
    ax.set_xlabel("Day"); ax.set_ylabel("Energy (model units per alive microbe)")
    ax.set_title(title); ax.set_xlim(days.min(), days.max())
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.20),
              ncols=3, framealpha=0.85, fontsize=8)
    add_year_labels(ax, COL_START, COL_END, N_YEARS, DAYS_PER_YEAR)

fig1.tight_layout()
fig1.savefig(FIG_DIR / "fig1_e_production_cost_gain.png", dpi=150, bbox_inches="tight")
plt.close(fig1); print("Saved fig1")

# ────────────────────────────────────────────────────────────────────────────
# Fig 2: Proportional E allocation
# ────────────────────────────────────────────────────────────────────────────
fig2, axes = plt.subplots(1, 2, figsize=(13, 5), sharey=True)
fig2.suptitle("Proportional Energy Allocation — Cost vs. Net Gain\n"
              "(community mean across alive taxa)",
              fontsize=12, fontweight="bold")

for ax, dat, title in zip(axes, [scen_data[l] for l in scen_labels[:2]], scen_labels[:2]):
    cf_mu, cf_sd = comm_ts(dat["cost_frac"])
    nf_mu, nf_sd = comm_ts(dat["net_frac"])
    if cf_mu is None:
        ax.set_title(f"{title} — no per-pathway data\n(re-run model to generate pathway CSVs)")
        continue
    days = cf_mu.index
    ax.stackplot(days, cf_mu, nf_mu.clip(0),
                 colors=[C_COST, C_NET], alpha=0.82,
                 labels=["Total E cost", "Net E gain"])
    ax.fill_between(days, (cf_mu - cf_sd).clip(0), (cf_mu + cf_sd).clip(0, 1),
                    color=C_COST, alpha=0.20, linewidth=0)
    ax.set_xlabel("Day"); ax.set_ylabel("Fraction of E production")
    ax.set_title(title); ax.set_xlim(days.min(), days.max()); ax.set_ylim(0, 1)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.20),
              ncols=2, framealpha=0.85, fontsize=8)
    ax.axhline(0.5, ls="--", lw=0.6, color="gray", alpha=0.45)
    add_year_labels(ax, COL_START, COL_END, N_YEARS, DAYS_PER_YEAR)

fig2.tight_layout()
fig2.savefig(FIG_DIR / "fig2_proportional_e_alloc.png", dpi=150, bbox_inches="tight")
plt.close(fig2); print("Saved fig2")

# ────────────────────────────────────────────────────────────────────────────
# Fig 3: Per-taxon mean E allocation bar chart
# ────────────────────────────────────────────────────────────────────────────
ctrl_cf_m = taxon_mean(ctrl["cost_frac"])
if ctrl_cf_m is not None:
    order = ctrl_cf_m.sort_values(ascending=True).index
else:
    order = None

fig3, axes = plt.subplots(1, 2, figsize=(14, 9), sharey=True)
fig3.suptitle("Mean Proportional E Allocation Per Taxon\n"
              "(averaged over all alive days; taxa sorted by Control cost fraction)",
              fontsize=12, fontweight="bold")

for ax, dat, title in zip(axes, [scen_data[l] for l in scen_labels[:2]], scen_labels[:2]):
    cf_df = dat["cost_frac"]
    nf_df = dat["net_frac"]
    if cf_df is None or order is None:
        ax.set_title(f"{title} — no per-pathway data"); continue
    cf = taxon_mean(cf_df).reindex(order)
    nf = taxon_mean(nf_df).reindex(order)
    yticks = np.arange(len(order))
    ax.barh(yticks, cf,           color=C_COST, alpha=0.85, label="Total E cost")
    ax.barh(yticks, nf.clip(0), left=cf.clip(0, 1), color=C_NET,  alpha=0.85, label="Net E gain")
    ax.set_yticks(yticks); ax.set_yticklabels(order, fontsize=7)
    ax.set_xlabel("Fraction of E production"); ax.set_xlim(0, 1)
    ax.set_title(title); ax.axvline(0.5, ls="--", lw=0.8, color="gray", alpha=0.5)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.20),
              ncols=2, framealpha=0.85, fontsize=8)

fig3.tight_layout()
fig3.savefig(FIG_DIR / "fig3_taxon_e_alloc_bars.png", dpi=150, bbox_inches="tight")
plt.close(fig3); print("Saved fig3")

# ────────────────────────────────────────────────────────────────────────────
# Fig 4: Within-cost pathway breakdown
# ────────────────────────────────────────────────────────────────────────────
fig4, axes = plt.subplots(1, 2, figsize=(13, 5), sharey=True)

has_pathways = any(scen_data[scen_labels[0]].get(k) is not None
                   for k in ['trans_pf', 'maint_pf', 'hsp_pf', 'enz_pf', 'osm_pf'])
if has_pathways:
    fig4.suptitle("Within-Cost E Pathway Breakdown\n"
                  "(fraction of total E cost; transport / maintenance / HSP / enzyme / osmolyte)",
                  fontsize=12, fontweight="bold")
    for ax, dat, title in zip(axes, [ctrl, warm], ["Control", "Warm"]):
        tr_m, _ = comm_ts(dat["trans_pf"])
        ma_m, _ = comm_ts(dat["maint_pf"])
        hs_m, _ = comm_ts(dat["hsp_pf"])
        en_m, _ = comm_ts(dat["enz_pf"])
        os_m, _ = comm_ts(dat["osm_pf"])
        if tr_m is None:
            ax.set_title(f"{title} — no data"); continue
        zero = pd.Series(0, index=tr_m.index)
        ma_m = ma_m if ma_m is not None else zero
        hs_m = hs_m if hs_m is not None else zero
        en_m = en_m if en_m is not None else zero
        os_m = os_m if os_m is not None else zero
        res_m = (1 - tr_m - ma_m - hs_m - en_m - os_m).clip(0)
        days = tr_m.index
        ax.stackplot(days, tr_m, ma_m, hs_m, en_m, os_m, res_m,
                     colors=[C_TRANS, C_MAINT, C_HSP, C_ENZ, C_OSMO, "#cccccc"],
                     alpha=0.85,
                     labels=["Transport", "Maintenance", "HSP",
                             "Enzymes", "Osmolytes", "Other/residual"])
        ax.set_xlabel("Day"); ax.set_ylabel("Fraction of total E cost")
        ax.set_title(title); ax.set_xlim(days.min(), days.max()); ax.set_ylim(0, 1)
        ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.20),
                  ncols=3, fontsize=8, framealpha=0.85)
        add_year_labels(ax, COL_START, COL_END, N_YEARS, DAYS_PER_YEAR)
else:
    fig4.suptitle("Within-Cost E Breakdown — pathway CSVs not available\n"
                  "Re-run the model after adding per-pathway tracking to output.py",
                  fontsize=11, fontweight="bold")
    for ax, dat, title in zip(axes, [ctrl, warm], ["Control", "Warm"]):
        ep_mu, _ = comm_ts(dat["EP"])
        if ep_mu is None:
            ax.set_title(f"{title} — no data"); continue
        ax.set_title(f"{title} — pathway breakdown not available")
        ax.text(0.5, 0.5,
                "Transport_E_Taxon, Maint_E_Taxon,\n"
                "HSP_E_Taxon, Enzyme_E_Taxon,\n"
                "Osmo_E_Taxon CSVs not found.\n\n"
                "These are produced by the per-pathway\n"
                "tracking added in the latest grid.py / output.py update.",
                transform=ax.transAxes, ha="center", va="center", fontsize=9,
                bbox=dict(facecolor="lightyellow", edgecolor="gray", alpha=0.8))

fig4.tight_layout()
fig4.savefig(FIG_DIR / "fig4_within_cost_pathway_breakdown.png", dpi=150, bbox_inches="tight")
plt.close(fig4); print("Saved fig4")

# ────────────────────────────────────────────────────────────────────────────
# Fig 5: Per-taxon E allocation heatmaps
# ────────────────────────────────────────────────────────────────────────────
cf_order = (taxon_mean(ctrl["cost_frac"]).sort_values(ascending=False).index
            if ctrl["cost_frac"] is not None else None)

fig5, axes = plt.subplots(2, 2, figsize=(16, 11))
fig5.suptitle("Per-Taxon E Allocation Heatmaps  "
              "(rows = taxa sorted by ctrl cost fraction; cols = days)",
              fontsize=12, fontweight="bold")

datasets = [
    (ctrl["cost_frac"], ctrl["trans_pf"], "Control"),
    (warm["cost_frac"], warm["trans_pf"], "Warm"),
]

day_labels = np.linspace(COL_START, COL_END, 6, dtype=int)

if cf_order is not None:
    for col, (cf, tp, scen) in enumerate(datasets):
        ax = axes[0, col]
        cf_sub = subsample(cf)
        if cf_sub is not None:
            mat = cf_sub.reindex(cf_order).fillna(np.nan).values
            im = ax.imshow(mat, aspect="auto", cmap="RdYlBu_r", vmin=0, vmax=1,
                           interpolation="nearest")
            ax.set_title(f"{scen} — E cost fraction of production", fontsize=10)
            ax.set_xlabel("Day (subsampled to 200)")
            ax.set_ylabel("Taxon")
            ax.set_yticks(np.arange(len(cf_order)))
            ax.set_yticklabels(cf_order, fontsize=5)
            xtick_pos = np.linspace(0, mat.shape[1] - 1, 6, dtype=int)
            ax.set_xticks(xtick_pos); ax.set_xticklabels(day_labels)
            plt.colorbar(im, ax=ax, label="E cost fraction", shrink=0.8)
            add_year_labels_heatmap(ax, mat.shape[1], COL_START, COL_END, N_YEARS, DAYS_PER_YEAR)

        ax2 = axes[1, col]
        tp_sub = subsample(tp)
        if tp_sub is not None:
            mat2 = tp_sub.reindex(cf_order).fillna(np.nan).values
            vmax2 = max(np.nanmax(mat2), 0.01)
            im2 = ax2.imshow(mat2, aspect="auto", cmap="Oranges", vmin=0, vmax=vmax2,
                             interpolation="nearest")
            ax2.set_title(f"{scen} — Transport fraction of E cost", fontsize=10)
            ax2.set_xlabel("Day (subsampled to 200)")
            ax2.set_ylabel("Taxon")
            ax2.set_yticks(np.arange(len(cf_order)))
            ax2.set_yticklabels(cf_order, fontsize=5)
            ax2.set_xticks(xtick_pos); ax2.set_xticklabels(day_labels)
            plt.colorbar(im2, ax=ax2, label="Transport E fraction", shrink=0.8)
            add_year_labels_heatmap(ax2, mat2.shape[1], COL_START, COL_END, N_YEARS, DAYS_PER_YEAR)
        else:
            axes[1, col].set_title(f"{scen} — Transport E data not available")
            axes[1, col].text(0.5, 0.5, "Re-run model with pathway tracking",
                              transform=axes[1, col].transAxes, ha="center", va="center")
else:
    for ax in axes.flat:
        ax.set_visible(False)
    fig5.text(0.5, 0.5, "No E cost fraction data available.",
              ha="center", va="center", fontsize=14)

fig5.tight_layout()
fig5.savefig(FIG_DIR / "fig5_taxon_e_alloc_heatmaps.png", dpi=150, bbox_inches="tight")
plt.close(fig5); print("Saved fig5")

# ────────────────────────────────────────────────────────────────────────────
# Fig 6: Warm − Control Δ E allocation per taxon
# ────────────────────────────────────────────────────────────────────────────
if ctrl["cost_frac"] is not None and warm["cost_frac"] is not None:
    ctrl_cf_m = taxon_mean(ctrl["cost_frac"])
    warm_cf_m = taxon_mean(warm["cost_frac"])
    delta_cf  = (warm_cf_m - ctrl_cf_m).dropna()

    ctrl_tr_m = taxon_mean(ctrl["trans_pf"]) if ctrl["trans_pf"] is not None else None
    warm_tr_m = taxon_mean(warm["trans_pf"]) if warm["trans_pf"] is not None else None
    delta_tr  = ((warm_tr_m - ctrl_tr_m).dropna()
                 if (ctrl_tr_m is not None and warm_tr_m is not None) else None)

    order6   = delta_cf.sort_values(ascending=True).index
    yticks6  = np.arange(len(order6))
    n_panels = 2 if delta_tr is not None else 1

    fig6, axs = plt.subplots(1, n_panels, figsize=(7 * n_panels, 9), sharey=True)
    if n_panels == 1:
        axs = [axs]
    fig6.suptitle("Warm − Control Difference in E Allocation Per Taxon\n"
                  "(positive = higher under warming)", fontsize=12, fontweight="bold")

    vals6   = delta_cf.reindex(order6).fillna(0).values
    colors6 = [C_WARM if v > 0 else C_CTRL for v in vals6]
    axs[0].barh(yticks6, vals6, color=colors6, alpha=0.85, edgecolor="white", linewidth=0.3)
    axs[0].axvline(0, color="black", lw=0.9)
    axs[0].set_yticks(yticks6); axs[0].set_yticklabels(order6, fontsize=7)
    axs[0].set_xlabel("Δ E cost fraction (Warm − Control)")
    axs[0].set_title("Δ E cost / production")

    if delta_tr is not None and n_panels > 1:
        vals_tr   = delta_tr.reindex(order6).fillna(0).values
        colors_tr = [C_WARM if v > 0 else C_CTRL for v in vals_tr]
        axs[1].barh(yticks6, vals_tr, color=colors_tr, alpha=0.85,
                    edgecolor="white", linewidth=0.3)
        axs[1].axvline(0, color="black", lw=0.9)
        axs[1].set_xlabel("Δ Transport E fraction (Warm − Control)")
        axs[1].set_title("Δ Transport fraction of E cost")

    legend_handles = [Patch(facecolor=C_WARM, label="Higher under warming"),
                      Patch(facecolor=C_CTRL, label="Lower under warming")]
    axs[-1].legend(handles=legend_handles, loc="upper center",
                   bbox_to_anchor=(0.5, -0.20), ncols=2, framealpha=0.85, fontsize=8)
    fig6.tight_layout()
    fig6.savefig(FIG_DIR / "fig6_warm_ctrl_delta_e_alloc.png", dpi=150, bbox_inches="tight")
    plt.close(fig6); print("Saved fig6")
else:
    print("Skipped fig6 — missing cost_frac data for control or warm")

print(f"\nAll figures saved to: {FIG_DIR}")
