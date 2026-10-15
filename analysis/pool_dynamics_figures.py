#!/usr/bin/env python3
"""
pool_dynamics_figures.py

Substrate, monomer, and microbial biomass pool dynamics for DEMENTpy.

Supports any number of climate scenarios produced by run_multi_climate_simulations.py.
Scenario directories are auto-discovered from --base or supplied explicitly via
--scenario-dirs.  Alternatively, pass --summary-dir to read pre-aggregated summary
CSVs from export_summary_data.py (full gradient / N-scenario), in which case pools
come from substrates_by_type / monomers_by_type / microbes_by_taxon / system_resp.
Cumulative-respiration-per-pulse (Fig 8) has no summary stem and is skipped there.

Fig 1  – Total C pool partitioning over time    (one panel per scenario)
Fig 2  – Substrate pools by type                (stacked area, one panel per scenario)
Fig 3  – Monomer pools by type                  (organic stacked; mineral NH4/PO4 line)
Fig 4  – Microbial biomass                      (total + top-N taxa per scenario)
Fig 4b – Per-taxon microbial biomass            (taxa only, no community total)
Fig 5  – Per-taxon biomass heatmap              (taxa × time, one panel per scenario)
Fig 6  – Pairwise Δ C pools                     (scenario B − scenario A, all pairs)
Fig 7  – Daily respiration                      (per grid cell and per unit biomass)
Fig 8  – Cumulative respiration per pulse        (per grid cell and per unit biomass; raw only)
Fig 9  – Cumulative respiration over full run    (per grid cell and per unit biomass)
"""

import argparse
import warnings
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.patches import Patch

warnings.filterwarnings("ignore")

# ── Add parent dir to path so scenario_utils is importable ───────────────────
sys.path.insert(0, str(Path(__file__).parent))
from scenario_utils import (
    load_series, pick_runs, resolve_runs, discover_scenarios,
    scenario_color_map, rep_mean, rep_std,
    read_runtime, year_range, filter_cols,
    add_year_labels, add_year_labels_heatmap,
    add_scenario_args, resolve_scenario_list,
    normalise_by_gridsize, get_scenario_gridsize,
    get_palette_from_args,
    add_summary_arg, resolve_summary_scenarios, load_summary, read_summary_meta,
)

# ── Hardcoded base path (override with --base) ────────────────────────────────
BASE = Path(
    "/Users/jswalkup/Library/CloudStorage/"
    "GoogleDrive-jwalkup@uci.edu/My Drive/"
    "Decoupled-Energy-DEMENTpy-20250902/"
    "Decoupled-Energy-DEMENTpy/output"
)

MIN_COLS = 100
TOP_TAXA = 10

# ── Display labels ─────────────────────────────────────────────────────────────
SUB_LABELS = {
    "DeadMic":      "Dead microbes",
    "DeadEnz":      "Dead enzymes",
    "Cellulose":    "Cellulose",
    "Hemicellulose":"Hemicellulose",
    "Starch":       "Starch",
    "Chitin":       "Chitin",
    "Lignin":       "Lignin",
    "Protein1":     "Protein 1",
    "Protein2":     "Protein 2",
    "Protein3":     "Protein 3",
    "OrgP1":        "Organic P 1",
    "OrgP2":        "Organic P 2",
}
MINERAL_MONOMERS = {"NH4", "PO4"}

# ── Style ─────────────────────────────────────────────────────────────────────
C_SUB  = "#8c564b"
C_MON  = "#ff7f0e"
C_MIC  = "#1f77b4"
SMOOTH = 14

plt.rcParams.update({
    "font.size": 10, "axes.titlesize": 11, "axes.labelsize": 10,
    "legend.fontsize": 8, "figure.dpi": 150,
})

# ── Helpers ───────────────────────────────────────────────────────────────────
def smooth(s: pd.Series, window: int = SMOOTH) -> pd.Series:
    return s.rolling(window, center=True, min_periods=1).mean()

def smooth_df(df: pd.DataFrame, window: int = SMOOTH) -> pd.DataFrame:
    return df.T.rolling(window, center=True, min_periods=1).mean().T

def total_ts(df: pd.DataFrame) -> pd.Series:
    return df.sum(axis=0)

def resp_per_biomass(resp, MIC):
    """Mass-specific respiration: respiration / total microbial biomass.

    Both inputs are already gridsize-normalised, so the per-cell factor cancels
    and the result is respiration per unit biomass.  Returns a day-indexed
    Series aligned on the days the two share, or None if either is missing.
    """
    if resp is None or MIC is None:
        return None
    tot = total_ts(MIC)
    common = resp.index.intersection(tot.index)
    if len(common) == 0:
        return None
    denom = tot.reindex(common).replace(0, np.nan)
    return resp.reindex(common) / denom

def pulse_mean_biomass(MIC, n_years, days_per_year):
    """Mean total biomass within each pulse window, indexed by pulse number."""
    if MIC is None:
        return None
    tot = total_ts(MIC)
    out = {}
    for p in range(1, n_years + 1):
        lo, hi = (p - 1) * days_per_year + 1, p * days_per_year
        days_in = [d for d in tot.index if lo <= d <= hi]
        if days_in:
            out[p] = tot.loc[days_in].mean()
    return pd.Series(out) if out else None

def monomer_label(name: str) -> str:
    if name.startswith("Mon_"):
        return "Mon. " + name[4:]
    return {"NH4": "NH₄⁺", "PO4": "PO₄³⁻", "FermProd": "Ferm. products"}.get(name, name)

def sub_label(name: str) -> str:
    return SUB_LABELS.get(name, name)

def tab_colors(n: int, cmap: str = "tab20") -> list:
    return [cm.get_cmap(cmap)(i / max(n - 1, 1)) for i in range(n)]

def subsample(df: pd.DataFrame, n: int = 200) -> pd.DataFrame:
    if df is None or df.shape[1] <= n:
        return df
    idx = np.linspace(0, df.shape[1] - 1, n, dtype=int)
    return df.iloc[:, idx]

# ── Pool loading ──────────────────────────────────────────────────────────────
def compute_pools(run_dir: Path) -> dict | None:
    SUB  = load_series(run_dir, "SubstratesSeries")
    MON  = load_series(run_dir, "MonomersSeries")
    MIC  = load_series(run_dir, "MicrobesSeries")
    RESP = load_series(run_dir, "RespSeries")       # 1 row × days
    # CumRespByPulse is a pulse-indexed Series exported as (n_pulses × 1) CSV.
    # load_series drops the single column thinking it's an init-state column,
    # so we read it directly with pd.read_csv instead.
    _crbp_path = run_dir / "CumRespByPulse.csv"
    if _crbp_path.exists():
        try:
            CRBP = pd.read_csv(_crbp_path, index_col=0).iloc[:, 0]
        except Exception:
            CRBP = None
    else:
        CRBP = None
    if SUB is None and MON is None and MIC is None:
        return None
    # Flatten 1-row DataFrames → Series for respiration
    if RESP is not None and hasattr(RESP, 'iloc'):
        resp_s = RESP.iloc[0] if RESP.ndim == 2 else RESP
    else:
        resp_s = RESP
    # CumRespByPulse is stored as a Series (index = pulse numbers).
    # load_series may return it as a DataFrame with various shapes — normalise.
    crbp_s = None
    if CRBP is not None:
        if isinstance(CRBP, pd.Series):
            crbp_s = CRBP
        elif hasattr(CRBP, 'iloc'):
            # Could be (n_pulses, 1) or (1, n_pulses) — find the long axis
            if CRBP.shape[1] == 1:
                crbp_s = CRBP.iloc[:, 0]   # (n_pulses, 1) → Series
            elif CRBP.shape[0] == 1:
                crbp_s = CRBP.iloc[0, :]   # (1, n_pulses) → Series
            elif CRBP.shape[1] > 1:
                crbp_s = CRBP.iloc[:, 0]   # take first column
            else:
                crbp_s = None              # empty DataFrame — skip
    return dict(SUB=SUB, MON=MON, MIC=MIC, RESP=resp_s, CRBP=crbp_s)

# ── Parse CLI ─────────────────────────────────────────────────────────────────
_p = argparse.ArgumentParser(
    description="Pool dynamics figures for DEMENTpy (multi-scenario)",
    formatter_class=argparse.RawDescriptionHelpFormatter,
)
add_scenario_args(_p)
add_summary_arg(_p)
_args = _p.parse_args()

# key → export_summary_data.py stem.  CRBP (cumulative resp by pulse) has no
# summary stem, so it is None in summary mode and Fig 8 auto-skips.
_SUMMARY_STEMS = {
    "SUB":  "substrates_by_type",
    "MON":  "monomers_by_type",
    "MIC":  "microbes_by_taxon",
    "RESP": "system_resp",
}

if _args.summary_dir:
    # ── Summary mode (all discovered scenarios → full gradient / N-scenario) ────
    summary_root = Path(_args.summary_dir)
    scenario_list = resolve_summary_scenarios(summary_root)
    if not scenario_list:
        print(f"No summary scenarios found under {summary_root}.")
        sys.exit(1)

    scen_labels = [s[0] for s in scenario_list]
    COLOR       = scenario_color_map(scen_labels, palette=get_palette_from_args(_args))
    N_SCEN      = len(scenario_list)

    base_path = summary_root
    FIG_DIR   = Path(_args.fig_dir) if _args.fig_dir else base_path / "pool_dynamics_figures"
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    print(f"Summary mode — scenarios ({N_SCEN}):")
    scen_data = {}   # label -> {SUB, MON, MIC, RESP, CRBP, *_sd}
    for label, sdir in scenario_list:
        print(f"  {label}: {sdir}")
        _sub  = load_summary(sdir, _SUMMARY_STEMS["SUB"])
        _mon  = load_summary(sdir, _SUMMARY_STEMS["MON"])
        _mic  = load_summary(sdir, _SUMMARY_STEMS["MIC"])
        _resp = load_summary(sdir, _SUMMARY_STEMS["RESP"])
        scen_data[label] = dict(
            SUB=_sub[0], MON=_mon[0], MIC=_mic[0], RESP=_resp[0], CRBP=None,
            SUB_sd=_sub[1], MON_sd=_mon[1], MIC_sd=_mic[1], RESP_sd=_resp[1],
        )
        if _mic[0] is None:
            print(f"    note: no microbes_by_taxon for {label}; "
                  f"per-taxon biomass figures (4, 5) will be empty for it.")

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
                                          sentinel="SubstratesSeries")
    if not scenario_list:
        print("No scenarios found. Use --base or --scenario-dirs.")
        sys.exit(1)

    scen_labels = [s[0] for s in scenario_list]
    COLOR       = scenario_color_map(scen_labels, palette=get_palette_from_args(_args))
    N_SCEN      = len(scenario_list)

    base_path = Path(_args.base) if _args.base else BASE
    FIG_DIR   = Path(_args.fig_dir) if _args.fig_dir else base_path / "pool_dynamics_figures"
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    # ── Load data per scenario ────────────────────────────────────────────────
    print(f"Scenarios ({N_SCEN}):")
    scen_data = {}   # label -> {SUB, MON, MIC, RESP, CRBP, *_sd}
    for label, run_dirs in scenario_list:
        print(f"  {label}: {len(run_dirs)} run(s)")
        pools = [compute_pools(d) for d in run_dirs]
        scen_data[label] = dict(
            SUB    = rep_mean(pools, "SUB"),
            MON    = rep_mean(pools, "MON"),
            MIC    = rep_mean(pools, "MIC"),
            RESP   = rep_mean(pools, "RESP"),
            CRBP   = rep_mean(pools, "CRBP"),
            SUB_sd = rep_std(pools,  "SUB"),
            MON_sd = rep_std(pools,  "MON"),
            MIC_sd = rep_std(pools,  "MIC"),
            RESP_sd= rep_std(pools,  "RESP"),
        )

    # ── Runtime / year filtering ──────────────────────────────────────────────
    _first_run = scenario_list[0][1][0] if scenario_list[0][1] else None
    _rt = read_runtime(_first_run)
    N_YEARS       = _rt["pulse"]
    DAYS_PER_YEAR = int(_rt["end_time"])
    COL_START, COL_END = year_range(_args.years, N_YEARS, DAYS_PER_YEAR)
    GRIDSIZE = get_scenario_gridsize(scenario_list)
    print(f"Simulation: {N_YEARS} year(s), {DAYS_PER_YEAR} days/year  |  "
          f"Displaying: day {COL_START}–{COL_END}")

# ── Year filtering (shared) ───────────────────────────────────────────────────
for label in scen_labels:
    d = scen_data[label]
    for k in list(d.keys()):
        # CRBP is indexed by pulse number, not day — don't filter by day range
        if k == "CRBP":
            continue
        d[k] = filter_cols(d[k], COL_START, COL_END)


# ── Normalise all series by gridsize ─────────────────────────────────────────
# Divides pool/pathway values by the number of grid cells so results
# are expressed per grid cell rather than as community totals.
# This makes runs on different grid sizes directly comparable.
print(f"Grid size     : {GRIDSIZE} cells  (normalising all series by {GRIDSIZE})")
for _label in scen_labels:
    scen_data[_label] = {
        k: normalise_by_gridsize(v, GRIDSIZE)
        for k, v in scen_data[_label].items()
    }

# ────────────────────────────────────────────────────────────────────────────
# Fig 1: Total C pool partitioning — one panel per scenario
# ────────────────────────────────────────────────────────────────────────────
fig1, axes = plt.subplots(1, N_SCEN, figsize=(7 * N_SCEN, 6), sharey=True)
if N_SCEN == 1:
    axes = [axes]
fig1.suptitle("Total C Pool Partitioning Over Time\n"
              "(substrates / monomers / microbial biomass)",
              fontsize=12, fontweight="bold")

for ax, label in zip(axes, scen_labels):
    d = scen_data[label]
    ax.set_title(label)
    plotted = False
    SUB, MON, MIC = d["SUB"], d["MON"], d["MIC"]
    SUBsd, MONsd, MICsd = d["SUB_sd"], d["MON_sd"], d["MIC_sd"]
    if SUB is not None:
        s_tot = total_ts(SUB)
        ax.plot(s_tot.index, s_tot.values, color=C_SUB, lw=2, label="Substrates")
        if SUBsd is not None:
            s_sd = total_ts(SUBsd)
            ax.fill_between(s_tot.index, (s_tot - s_sd).clip(0), s_tot + s_sd,
                            color=C_SUB, alpha=0.20, linewidth=0)
        plotted = True
    if MON is not None:
        m_tot = total_ts(MON)
        ax.plot(m_tot.index, m_tot.values, color=C_MON, lw=2, label="Monomers")
        if MONsd is not None:
            m_sd = total_ts(MONsd)
            ax.fill_between(m_tot.index, (m_tot - m_sd).clip(0), m_tot + m_sd,
                            color=C_MON, alpha=0.20, linewidth=0)
        plotted = True
    if MIC is not None:
        b_tot = total_ts(MIC)
        ax.plot(b_tot.index, b_tot.values, color=C_MIC, lw=2, label="Microbial biomass")
        if MICsd is not None:
            b_sd = total_ts(MICsd)
            ax.fill_between(b_tot.index, (b_tot - b_sd).clip(0), b_tot + b_sd,
                            color=C_MIC, alpha=0.20, linewidth=0)
        plotted = True
    if plotted:
        ax.set_xlabel("Day"); ax.set_ylabel("Total C (model units)")
        ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.20),
                  ncols=3, framealpha=0.85, fontsize=8)
        add_year_labels(ax, COL_START, COL_END, N_YEARS, DAYS_PER_YEAR)
    else:
        ax.text(0.5, 0.5, "No data", transform=ax.transAxes, ha="center", va="center")

fig1.tight_layout()
fig1.savefig(FIG_DIR / "fig1_total_pool_partitioning.png", dpi=150, bbox_inches="tight")
plt.close(fig1); print("Saved fig1")

# ────────────────────────────────────────────────────────────────────────────
# Fig 2: Substrate pools by type — stacked area, one panel per scenario
# ────────────────────────────────────────────────────────────────────────────
ref_SUB = next((scen_data[l]["SUB"] for l in scen_labels
                if scen_data[l]["SUB"] is not None), None)
sub_types  = list(ref_SUB.index) if ref_SUB is not None else []
sub_colors = tab_colors(len(sub_types), "tab20")
sub_lbls   = [sub_label(t) for t in sub_types]

fig2, axes = plt.subplots(1, N_SCEN, figsize=(7 * N_SCEN, 6), sharey=True)
if N_SCEN == 1:
    axes = [axes]
fig2.suptitle("Substrate C Mass by Type\n(lines per substrate type)",
              fontsize=12, fontweight="bold")

for ax, label in zip(axes, scen_labels):
    SUB = scen_data[label]["SUB"]
    ax.set_title(label)
    if SUB is None or not sub_types:
        ax.text(0.5, 0.5, "No data", transform=ax.transAxes, ha="center", va="center")
        continue
    df   = SUB.reindex(sub_types).fillna(0)
    days = df.columns
    for i, stype in enumerate(sub_types):
        ax.plot(days, df.loc[stype].values, color=sub_colors[i],
                lw=1.5, label=sub_lbls[i], alpha=0.85)
    ax.set_xlabel("Day"); ax.set_ylabel("C mass per cell (model units)")
    ax.set_xlim(days.min(), days.max())
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.20),
              ncols=2, fontsize=7, framealpha=0.85)
    add_year_labels(ax, COL_START, COL_END, N_YEARS, DAYS_PER_YEAR)
    ax.spines[["top", "right"]].set_visible(False)

fig2.tight_layout()
fig2.savefig(FIG_DIR / "fig2_substrate_pools_by_type.png", dpi=150, bbox_inches="tight")
plt.close(fig2); print("Saved fig2")

# ────────────────────────────────────────────────────────────────────────────
# Fig 3: Organic monomer pools by type — stacked area, one panel per scenario
# ────────────────────────────────────────────────────────────────────────────
ref_MON = next((scen_data[l]["MON"] for l in scen_labels
                if scen_data[l]["MON"] is not None), None)
if ref_MON is not None:
    organic_idx = [t for t in ref_MON.index if t not in MINERAL_MONOMERS]
    org_colors  = tab_colors(len(organic_idx), "tab20b")
else:
    organic_idx = []
    org_colors  = []

fig3, axes3 = plt.subplots(1, N_SCEN, figsize=(7 * N_SCEN, 6), sharey=True)
if N_SCEN == 1:
    axes3 = [axes3]
fig3.suptitle("Organic Monomer C Mass by Type\n"
              "(lines per monomer type; inorganic NH\u2084\u207a and PO\u2084\u00b3\u207b excluded)",
              fontsize=12, fontweight="bold")

for ax, label in zip(axes3, scen_labels):
    MON = scen_data[label]["MON"]
    ax.set_title(label)
    if MON is None or not organic_idx:
        ax.text(0.5, 0.5, "No data", transform=ax.transAxes, ha="center", va="center")
        continue
    org_df = MON.reindex(organic_idx).fillna(0)
    days   = org_df.columns
    for i, mtype in enumerate(organic_idx):
        ax.plot(days, org_df.loc[mtype].values, color=org_colors[i],
                lw=1.5, label=monomer_label(mtype), alpha=0.85)
    ax.set_xlabel("Day"); ax.set_ylabel("C mass per cell (model units)")
    ax.set_xlim(days.min(), days.max())
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.20),
              ncols=2, fontsize=7, framealpha=0.85)
    add_year_labels(ax, COL_START, COL_END, N_YEARS, DAYS_PER_YEAR)
    ax.spines[["top", "right"]].set_visible(False)

fig3.tight_layout()
fig3.savefig(FIG_DIR / "fig3_organic_monomer_pools_by_type.png", dpi=150, bbox_inches="tight")
plt.close(fig3); print("Saved fig3")

# ────────────────────────────────────────────────────────────────────────────
# Fig 4: Microbial biomass — total + top-N taxa, one panel per scenario
# ────────────────────────────────────────────────────────────────────────────
ref_MIC = next((scen_data[l]["MIC"] for l in scen_labels
                if scen_data[l]["MIC"] is not None), None)
top_taxa   = []
if ref_MIC is not None:
    top_taxa = list(ref_MIC.mean(axis=1).sort_values(ascending=False).index[:TOP_TAXA])
top_colors = tab_colors(TOP_TAXA, "tab10")

fig4, axes = plt.subplots(1, N_SCEN, figsize=(7 * N_SCEN, 6), sharey=True)
if N_SCEN == 1:
    axes = [axes]
fig4.suptitle(f"Microbial Biomass Over Time\n"
              f"(bold = community total; coloured = top {TOP_TAXA} taxa "
              f"by mean biomass in first scenario; grey = remaining)",
              fontsize=11, fontweight="bold")

for ax, label in zip(axes, scen_labels):
    MIC = scen_data[label]["MIC"]
    ax.set_title(label)
    if MIC is None:
        ax.text(0.5, 0.5, "No data", transform=ax.transAxes, ha="center", va="center")
        continue
    days = MIC.columns
    for taxon in [t for t in MIC.index if t not in top_taxa]:
        s = MIC.loc[taxon]
        ax.plot(s.index, s.values, color="grey", lw=0.5, alpha=0.25)
    for i, taxon in enumerate(top_taxa):
        if taxon not in MIC.index:
            continue
        s = MIC.loc[taxon]
        ax.plot(s.index, s.values, color=top_colors[i], lw=1.2, alpha=0.85, label=taxon)
    tot = total_ts(MIC)
    ax.plot(tot.index, tot.values, color="black", lw=2.5, label="Community total", zorder=5)
    ax.set_xlabel("Day"); ax.set_ylabel("Biomass C (model units)")
    ax.set_xlim(days.min(), days.max())
    handles, labels_leg = ax.get_legend_handles_labels()
    tidx = [j for j, lb in enumerate(labels_leg) if lb == "Community total"]
    if tidx:
        handles = [handles[tidx[0]]] + [h for j, h in enumerate(handles) if j not in tidx]
        labels_leg = ["Community total"] + [lb for j, lb in enumerate(labels_leg) if j not in tidx]
    ax.legend(handles, labels_leg, loc="upper center", bbox_to_anchor=(0.5, -0.20),
              ncols=4, fontsize=7, framealpha=0.85)
    add_year_labels(ax, COL_START, COL_END, N_YEARS, DAYS_PER_YEAR)

fig4.tight_layout()
fig4.savefig(FIG_DIR / "fig4_microbial_biomass.png", dpi=150, bbox_inches="tight")
plt.close(fig4); print("Saved fig4")

# ────────────────────────────────────────────────────────────────────────────
# Fig 4b: Per-taxon microbial biomass — taxa only (no community total)
# ────────────────────────────────────────────────────────────────────────────
fig4b, axes = plt.subplots(1, N_SCEN, figsize=(7 * N_SCEN, 6), sharey=True)
if N_SCEN == 1:
    axes = [axes]
fig4b.suptitle(f"Per-Taxon Microbial Biomass Over Time\n"
               f"(coloured = top {TOP_TAXA} taxa by mean biomass in first scenario; "
               f"grey = remaining; no community total)",
               fontsize=11, fontweight="bold")

for ax, label in zip(axes, scen_labels):
    MIC = scen_data[label]["MIC"]
    ax.set_title(label)
    if MIC is None:
        ax.text(0.5, 0.5, "No data", transform=ax.transAxes, ha="center", va="center")
        continue
    days = MIC.columns
    for taxon in [t for t in MIC.index if t not in top_taxa]:
        s = MIC.loc[taxon]
        ax.plot(s.index, s.values, color="grey", lw=0.5, alpha=0.25)
    for i, taxon in enumerate(top_taxa):
        if taxon not in MIC.index:
            continue
        s = MIC.loc[taxon]
        ax.plot(s.index, s.values, color=top_colors[i], lw=1.2, alpha=0.85, label=taxon)
    ax.set_xlabel("Day"); ax.set_ylabel("Biomass C (model units)")
    ax.set_xlim(days.min(), days.max())
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.20),
              ncols=5, fontsize=7, framealpha=0.85)
    add_year_labels(ax, COL_START, COL_END, N_YEARS, DAYS_PER_YEAR)

fig4b.tight_layout()
fig4b.savefig(FIG_DIR / "fig4b_taxon_biomass_only.png", dpi=150, bbox_inches="tight")
plt.close(fig4b); print("Saved fig4b")

# ────────────────────────────────────────────────────────────────────────────
# Fig 5: Per-taxon biomass heatmap — one panel per scenario
# ────────────────────────────────────────────────────────────────────────────
taxa_order = (ref_MIC.mean(axis=1).sort_values(ascending=False).index
              if ref_MIC is not None else None)

fig5, axes = plt.subplots(1, N_SCEN, figsize=(8 * N_SCEN, 10))
if N_SCEN == 1:
    axes = [axes]
fig5.suptitle("Per-Taxon Microbial Biomass Heatmap\n"
              "(rows = taxa sorted by mean biomass in first scenario; subsampled to 200 days)",
              fontsize=12, fontweight="bold")

day_labels = np.linspace(COL_START, COL_END, 6, dtype=int)
for ax, label in zip(axes, scen_labels):
    MIC = scen_data[label]["MIC"]
    ax.set_title(label)
    if MIC is None or taxa_order is None:
        ax.text(0.5, 0.5, "No data", transform=ax.transAxes, ha="center", va="center")
        continue
    MIC_sub = subsample(MIC.reindex(taxa_order))
    mat  = MIC_sub.values
    vmax = np.nanpercentile(mat[mat > 0], 95) if np.any(mat > 0) else 1.0
    im   = ax.imshow(mat, aspect="auto", cmap="YlOrBr", vmin=0, vmax=vmax,
                     interpolation="nearest")
    ax.set_yticks(np.arange(len(taxa_order)))
    ax.set_yticklabels(taxa_order, fontsize=5)
    ax.set_ylabel("Taxon"); ax.set_xlabel("Day (subsampled to 200)")
    xtick_pos = np.linspace(0, mat.shape[1] - 1, 6, dtype=int)
    ax.set_xticks(xtick_pos); ax.set_xticklabels(day_labels)
    plt.colorbar(im, ax=ax, label="Biomass C (model units)", shrink=0.75)
    add_year_labels_heatmap(ax, mat.shape[1], COL_START, COL_END, N_YEARS, DAYS_PER_YEAR)

fig5.tight_layout()
fig5.savefig(FIG_DIR / "fig5_taxon_biomass_heatmap.png", dpi=150, bbox_inches="tight")
plt.close(fig5); print("Saved fig5")

# ────────────────────────────────────────────────────────────────────────────
# Fig 6: Pairwise Δ C pools (scenario B − scenario A) for all adjacent pairs
# ────────────────────────────────────────────────────────────────────────────
pairs = [(scen_labels[i], scen_labels[i + 1]) for i in range(N_SCEN - 1)]

for a_label, b_label in pairs:
    a = scen_data[a_label]
    b = scen_data[b_label]
    pool_pairs = [("SUB", "Substrates",     sub_label,    "Δ C mass"),
                  ("MON", "Monomers",        monomer_label, "Δ C mass"),
                  ("MIC", "Microbial taxa",  lambda x: x,  "Δ Biomass C")]
    available  = [(k, t, lf, xl) for k, t, lf, xl in pool_pairs
                  if a[k] is not None and b[k] is not None]
    if not available:
        print(f"Skipped fig6 pair {a_label}/{b_label} — no overlapping data")
        continue

    n_panels = len(available)
    fig6, axs = plt.subplots(1, n_panels, figsize=(6 * n_panels, 8))
    if n_panels == 1:
        axs = [axs]
    fig6.suptitle(f"{b_label} − {a_label}  Mean C Mass Difference\n"
                  f"(positive = higher in {b_label}; averaged over all days)",
                  fontsize=12, fontweight="bold")

    c_pos = COLOR[b_label]
    c_neg = COLOR[a_label]

    for ax, (key, title, label_fn, xlabel) in zip(axs, available):
        a_m    = a[key].mean(axis=1)
        b_m    = b[key].mean(axis=1)
        common = a_m.index.intersection(b_m.index)
        delta  = (b_m - a_m).reindex(common).sort_values()
        colors = [c_pos if v > 0 else c_neg for v in delta.values]
        yticks = np.arange(len(delta))
        ax.barh(yticks, delta.values, color=colors, alpha=0.85,
                edgecolor="white", linewidth=0.4)
        ax.set_yticks(yticks)
        ax.set_yticklabels([label_fn(t) for t in delta.index], fontsize=9)
        ax.axvline(0, color="black", lw=0.9)
        ax.set_xlabel(f"{xlabel} ({b_label} − {a_label})")
        ax.set_title(title)

    legend_handles = [Patch(facecolor=c_pos, label=f"Higher in {b_label}"),
                      Patch(facecolor=c_neg, label=f"Higher in {a_label}")]
    axs[-1].legend(handles=legend_handles, loc="upper center",
                   bbox_to_anchor=(0.5, -0.20), ncols=2, framealpha=0.85, fontsize=8)
    fig6.tight_layout()
    safe_a = a_label.replace("/", "-")
    safe_b = b_label.replace("/", "-")
    fname  = FIG_DIR / f"fig6_delta_pools_{safe_b}_minus_{safe_a}.png"
    fig6.savefig(fname, dpi=150, bbox_inches="tight")
    plt.close(fig6)
    print(f"Saved fig6  ({b_label} − {a_label})")


# ────────────────────────────────────────────────────────────────────────────
# Figs 7–9: Respiration, drawn in two denominators each —
#   per grid cell   (respiration as-is; the gridsize-normalised series)
#   per unit biomass (respiration ÷ total microbial biomass; gridsize cancels)
# Each helper is called once per denominator.  `_daily_resp_series` returns the
# per-day series and its ±SD band amount for a scenario under either denominator.
# ────────────────────────────────────────────────────────────────────────────
def _daily_resp_series(label, per_biomass):
    """(series, band) for one scenario, or (None, None) if unavailable.

    per_biomass=False → respiration per grid cell (series=RESP, band=RESP_sd).
    per_biomass=True  → respiration ÷ total biomass (band = RESP_sd / biomass).
    """
    resp    = scen_data[label].get("RESP")
    resp_sd = scen_data[label].get("RESP_sd")
    MIC     = scen_data[label].get("MIC")
    if not per_biomass:
        if resp is None:
            return None, None
        return resp, (resp_sd if resp_sd is not None else None)
    series = resp_per_biomass(resp, MIC)
    if series is None:
        return None, None
    if resp_sd is None:
        return series, None
    tot = total_ts(MIC).reindex(series.index).replace(0, np.nan)
    return series, resp_sd.reindex(series.index) / tot


def draw_daily_respiration(per_biomass):
    tag    = "per biomass" if per_biomass else "per cell"
    title  = "Daily Mass-Specific Respiration" if per_biomass else "Daily Community Respiration"
    sub    = "respiration per unit microbial biomass" if per_biomass else "per grid cell"
    ylabel = ("Respiration per unit biomass (model units day⁻¹)" if per_biomass
              else "Respiration per grid cell (model units day⁻¹)")
    fname  = ("fig7_daily_respiration_per_biomass.png" if per_biomass
              else "fig7_daily_respiration_per_cell.png")

    fig, ax = plt.subplots(figsize=(10, 5))
    fig.suptitle(f"{title}\n({sub}; mean ± 1 SD across replicates; all scenarios overlaid)",
                 fontsize=12, fontweight="bold")
    any_resp = False
    for label in scen_labels:
        series, band = _daily_resp_series(label, per_biomass)
        if series is None:
            continue
        days = series.index
        ax.plot(days, series.values, color=COLOR[label], lw=1.8, label=label)
        if band is not None:
            ax.fill_between(days,
                            (series - band).clip(0).values,
                            (series + band).values,
                            color=COLOR[label], alpha=0.18, linewidth=0)
        any_resp = True

    if any_resp:
        ax.set_xlabel("Day"); ax.set_ylabel(ylabel)
        ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.15),
                  ncols=min(N_SCEN, 4), framealpha=0.85, fontsize=9)
        ax.spines[["top", "right"]].set_visible(False)
        add_year_labels(ax, COL_START, COL_END, N_YEARS, DAYS_PER_YEAR)
    else:
        ax.text(0.5, 0.5, "No respiration / biomass data found",
                transform=ax.transAxes, ha="center", va="center")
    fig.tight_layout()
    fig.savefig(FIG_DIR / fname, dpi=150, bbox_inches="tight")
    plt.close(fig); print(f"Saved fig7 ({tag})")


def draw_cumulative_respiration(per_biomass):
    tag    = "per biomass" if per_biomass else "per cell"
    title  = ("Cumulative Mass-Specific Respiration Over Full Simulation" if per_biomass
              else "Cumulative Respiration Over Full Simulation")
    sub    = ("running total of daily respiration per unit biomass" if per_biomass
              else "running daily total, per grid cell")
    ylabel = ("Cumulative respiration per unit biomass (model units)" if per_biomass
              else "Cumulative respiration per grid cell (model units)")
    fname  = ("fig9_cumresp_full_run_per_biomass.png" if per_biomass
              else "fig9_cumresp_full_run_per_cell.png")

    fig, ax = plt.subplots(figsize=(10, 5))
    fig.suptitle(f"{title}\n({sub}; mean ± 1 SD)", fontsize=12, fontweight="bold")
    any_cum = False
    for label in scen_labels:
        series, band = _daily_resp_series(label, per_biomass)
        if series is None:
            continue
        cum  = series.cumsum()
        days = cum.index
        ax.plot(days, cum.values, color=COLOR[label], lw=1.8, label=label)
        if band is not None:
            # Propagate SD through cumsum: Var(sum) = sum(Var) for independent days
            cum_sd = band.fillna(0).pow(2).cumsum().pow(0.5)
            ax.fill_between(days,
                            (cum - cum_sd).clip(0).values,
                            (cum + cum_sd).values,
                            color=COLOR[label], alpha=0.18, linewidth=0)
        any_cum = True

    if any_cum:
        ax.set_xlabel("Day"); ax.set_ylabel(ylabel)
        ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.15),
                  ncols=min(N_SCEN, 4), framealpha=0.85, fontsize=9)
        ax.spines[["top", "right"]].set_visible(False)
        add_year_labels(ax, COL_START, COL_END, N_YEARS, DAYS_PER_YEAR)
    else:
        ax.text(0.5, 0.5, "No respiration / biomass data found",
                transform=ax.transAxes, ha="center", va="center")
    fig.tight_layout()
    fig.savefig(FIG_DIR / fname, dpi=150, bbox_inches="tight")
    plt.close(fig); print(f"Saved fig9 ({tag})")


def draw_cumresp_per_pulse(per_biomass):
    # CRBP exists only in raw mode; the per-biomass variant also needs biomass.
    tag = "per biomass" if per_biomass else "per cell"
    if per_biomass:
        have = any(scen_data[l].get("CRBP") is not None
                   and scen_data[l].get("MIC") is not None for l in scen_labels)
    else:
        have = any(scen_data[l].get("CRBP") is not None for l in scen_labels)
    if not have:
        need = "CumRespByPulse + biomass" if per_biomass else "CumRespByPulse"
        print(f"Skipped fig8 ({tag}) — no {need} data found")
        return

    _ref_crbp = next(scen_data[l]["CRBP"] for l in scen_labels
                     if scen_data[l].get("CRBP") is not None)
    ref_pulses   = list(_ref_crbp.index)
    pulse_labels = [f"Pulse {int(p)}" for p in ref_pulses]
    n_pulses     = len(pulse_labels)
    bar_width    = 0.8 / N_SCEN
    x            = np.arange(n_pulses)

    title  = ("Cumulative Respiration per Pulse, per Unit Biomass" if per_biomass
              else "Cumulative Respiration per Pulse")
    sub    = ("pulse cumulative respiration ÷ mean biomass over the pulse; mean across replicates"
              if per_biomass else "per grid cell; mean across replicates")
    ylabel = ("Cumulative respiration per unit biomass (model units)" if per_biomass
              else "Cumulative respiration per grid cell (model units)")
    fname  = ("fig8_cumresp_per_pulse_per_biomass.png" if per_biomass
              else "fig8_cumresp_per_pulse_per_cell.png")

    fig, ax = plt.subplots(figsize=(max(7, n_pulses * 2), 5))
    fig.suptitle(f"{title}\n({sub})", fontsize=12, fontweight="bold")

    for i, label in enumerate(scen_labels):
        crbp = scen_data[label].get("CRBP")
        if crbp is None:
            continue
        if per_biomass:
            MIC = scen_data[label].get("MIC")
            pbm = pulse_mean_biomass(MIC, N_YEARS, DAYS_PER_YEAR) if MIC is not None else None
            if pbm is None:
                continue
            common = crbp.index.intersection(pbm.index)
            vals   = crbp.reindex(common) / pbm.reindex(common).replace(0, np.nan)
            idx    = common
        else:
            vals, idx = crbp, crbp.index
        pos    = [ref_pulses.index(p) for p in idx if p in ref_pulses]
        offset = (i - N_SCEN / 2 + 0.5) * bar_width
        ax.bar(np.array(pos) + offset, vals.values, width=bar_width * 0.9,
               color=COLOR[label], label=label, alpha=0.85)

    ax.set_xticks(x); ax.set_xticklabels(pulse_labels)
    ax.set_xlabel("Pulse"); ax.set_ylabel(ylabel)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.15),
              ncols=min(N_SCEN, 4), framealpha=0.85, fontsize=9)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(FIG_DIR / fname, dpi=150, bbox_inches="tight")
    plt.close(fig); print(f"Saved fig8 ({tag})")


for _pb in (False, True):
    draw_daily_respiration(_pb)
    draw_cumresp_per_pulse(_pb)
    draw_cumulative_respiration(_pb)

print(f"\nAll figures saved to: {FIG_DIR}")

