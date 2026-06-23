#!/usr/bin/env python3
"""
cue_figures.py

Carbon Use Efficiency (CUE) and Carbon Assimilation Efficiency (CAE) figures
for DEMENTpy — N-scenario gradient mode.

Definitions
-----------
CUE  Carbon Use Efficiency
     = growth C / assimilated C
     Fraction of assimilated C incorporated into new biomass (0–1).
CAE  Carbon Assimilation Efficiency
     = assimilated C / uptaken C
     Fraction of uptaken C that is assimilated (not immediately respired for
     transport costs) (0–1).

Data sources
------------
Raw-run mode  : CUE_system.csv, CAE_system.csv (1 row × days; system-level)
                CUE_TaxonSeries.csv, CAE_TaxonSeries.csv (taxa × days)
Summary mode  : system_CUE, system_CAE, cue_taxon, cae_taxon stems from
                export_summary_data.py

Figures
-------
Fig 1 – System-level CUE and CAE — all scenarios overlaid (1 × 2 panels)
Fig 2 – Community-mean taxon CUE and CAE — all scenarios with ±1 SD (1 × 2 panels)

Usage
-----
    # Summary mode — auto-discovers all scenarios under --summary-dir
    python analysis/cue_figures.py \\
        --summary-dir /path/to/summary \\
        --fig-dir     /path/to/figures \\
        --treatment-map static_T3_Pn0p5="3°C" static_T9_Pn0p5="9°C"

    # Raw-run mode — specify scenario directories explicitly
    python analysis/cue_figures.py \\
        --scenario-dirs /path/to/control /path/to/warm \\
        --treatment-map control=Control warm=Warm \\
        --fig-dir /path/to/figures
"""

import argparse
import sys
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

warnings.filterwarnings("ignore")

# Make scenario_utils importable
_HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(_HERE))
try:
    import scenario_utils as su
    from scenario_utils import (
        pick_runs, rep_mean, rep_std, read_runtime, year_range, filter_cols,
        add_year_labels,
        add_scenario_args, resolve_scenario_list,
        add_summary_arg, resolve_summary_scenarios, load_summary, read_summary_meta,
        scenario_color_map, get_palette_from_args, _parse_treatment_map,
    )
    _SU_ERR = None
except ImportError as _e:
    su = None
    _SU_ERR = _e

# ── Paths ─────────────────────────────────────────────────────────────────────
BASE = Path(".")

# ── Style ─────────────────────────────────────────────────────────────────────
plt.rcParams.update({
    "font.size": 10, "axes.titlesize": 11, "axes.labelsize": 10,
    "legend.fontsize": 9, "figure.dpi": 150,
})

# ── Per-run helpers ────────────────────────────────────────────────────────────
def _load(run_dir: Path, name: str):
    f = run_dir / f"{name}.csv"
    if not f.exists():
        return None
    df = pd.read_csv(f, index_col=0)
    df.columns = range(len(df.columns))
    df = df.drop(columns=[0], errors="ignore")
    df.columns = range(1, len(df.columns) + 1)
    df = df.replace([np.inf, -np.inf], np.nan)
    return df


def compute_run(run_dir: Path) -> dict | None:
    CUE_sys = _load(run_dir, "CUE_system")
    CAE_sys = _load(run_dir, "CAE_system")
    CUE_tax = _load(run_dir, "CUE_TaxonSeries")
    CAE_tax = _load(run_dir, "CAE_TaxonSeries")
    if all(x is None for x in [CUE_sys, CAE_sys, CUE_tax, CAE_tax]):
        return None
    def clip01(df):
        return df.clip(lower=0, upper=1) if df is not None else None
    return dict(
        CUE_sys=clip01(CUE_sys), CAE_sys=clip01(CAE_sys),
        CUE_tax=clip01(CUE_tax), CAE_tax=clip01(CAE_tax),
    )


def comm_ts(df):
    """Community mean ± SD across taxa (axis=0) at each day."""
    if df is None:
        return None, None
    mu  = df.mean(axis=0, skipna=True)
    sig = df.std(axis=0,  skipna=True)
    return mu, sig


# ── Summary helpers ────────────────────────────────────────────────────────────
_SUMMARY_STEMS = {
    "CUE_sys": "system_CUE",
    "CAE_sys": "system_CAE",
    "CUE_tax": "cue_taxon",
    "CAE_tax": "cae_taxon",
}
_MBC_STEM = "system_microbial_biomass_C"


def _build_from_summary(sdir: Path) -> dict:
    out = {}
    for key, stem in _SUMMARY_STEMS.items():
        mean, _sd = load_summary(sdir, stem)
        if mean is None:
            out[key] = None
        elif isinstance(mean, pd.Series):
            out[key] = pd.DataFrame([mean.values], columns=list(mean.index))
        else:
            out[key] = mean
    for key in out:
        if out[key] is not None:
            out[key] = out[key].clip(lower=0, upper=1)
    return out


def _sd_from_summary(sdir: Path) -> dict:
    """Load the SD arrays for taxon series (needed for ±SD ribbons in fig 2)."""
    out = {}
    for key, stem in _SUMMARY_STEMS.items():
        _mean, sd = load_summary(sdir, stem)
        if sd is None:
            out[key] = None
        elif isinstance(sd, pd.Series):
            out[key] = pd.DataFrame([sd.values], columns=list(sd.index))
        else:
            out[key] = sd.clip(lower=0, upper=1) if sd is not None else None
    return out


# ── CLI ────────────────────────────────────────────────────────────────────────
_p = argparse.ArgumentParser(
    description="CUE / CAE figures for DEMENTpy (N-scenario gradient mode)",
    formatter_class=argparse.RawDescriptionHelpFormatter,
)
add_scenario_args(_p)   # --base, --scenario-dirs, --treatment-map, --palette,
                        # --fig-dir, --years, --scenario-labels
add_summary_arg(_p)     # --summary-dir
_p.add_argument("--days-per-year", type=int, default=365, metavar="N",
                help="Days per simulation year (default: 365).")
_args = _p.parse_args()

if _args.base:
    BASE = Path(_args.base)

FIG_DIR = Path(_args.fig_dir) if _args.fig_dir else BASE / "cue_figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

KEYS = ["CUE_sys", "CAE_sys", "CUE_tax", "CAE_tax"]

# ── Data loading ───────────────────────────────────────────────────────────────
# scenarios: list of (label, data_dict, sd_dict)
# data_dict: {key: DataFrame | None}
# sd_dict:   {key: DataFrame | None}  — SD across replicates / taxa

scenarios = []   # [(label, data, sd), ...]

if _args.summary_dir:
    # ── Summary mode ────────────────────────────────────────────────────────────
    if su is None:
        sys.exit(f"ERROR: --summary-dir requires scenario_utils: {_SU_ERR}")

    summary_root = Path(_args.summary_dir)
    _tmap = _parse_treatment_map(getattr(_args, "treatment_map", None))
    scenario_dirs_labels = resolve_summary_scenarios(summary_root,
                                                     treatment_map=_tmap)
    if not scenario_dirs_labels:
        sys.exit(f"ERROR: no scenario subdirectories found in {summary_root}")

    print(f"Summary mode — {len(scenario_dirs_labels)} scenario(s):")
    for label, sdir in scenario_dirs_labels:
        print(f"  {label}: {sdir}")
        data = _build_from_summary(sdir)
        sd   = _sd_from_summary(sdir)
        scenarios.append((label, data, sd))

    _meta = read_summary_meta(scenario_dirs_labels[0][1])
    DAYS_PER_YEAR = int(_args.days_per_year)
    N_YEARS       = _meta.get("n_years", 1)

    # Infer day range from loaded data
    _spans = []
    for _, data, _ in scenarios:
        for v in data.values():
            if v is not None and v.shape[1] > 0:
                cols = [c for c in v.columns if isinstance(c, (int, np.integer))]
                if cols:
                    _spans.append((min(cols), max(cols)))
    if _spans:
        _data_start = min(s for s, _ in _spans)
        _data_end   = max(e for _, e in _spans)
    else:
        _data_start, _data_end = 1, DAYS_PER_YEAR

    # Apply --years filtering (same logic as raw-run mode)
    COL_START, COL_END = year_range(_args.years, N_YEARS, DAYS_PER_YEAR)
    # Clamp to actual data range
    COL_START = max(COL_START, _data_start)
    COL_END   = min(COL_END,   _data_end)

    # Apply column filtering
    scenarios = [
        (label,
         {k: filter_cols(v, COL_START, COL_END) for k, v in data.items()},
         {k: filter_cols(v, COL_START, COL_END) for k, v in sd.items()})
        for label, data, sd in scenarios
    ]

    print(f"Simulation: {N_YEARS} year(s), {DAYS_PER_YEAR} days/year  |  "
          f"Displaying: day {COL_START}–{COL_END}")

else:
    # ── Raw-run mode ─────────────────────────────────────────────────────────────
    if su is None:
        sys.exit(f"ERROR: scenario_utils required: {_SU_ERR}")

    scenario_list = resolve_scenario_list(_args, BASE, min_cols=100)
    if not scenario_list:
        sys.exit("ERROR: no scenarios found — use --base or --scenario-dirs.")

    print(f"Raw-run mode — {len(scenario_list)} scenario(s):")
    all_run_dirs = []
    for label, run_dirs in scenario_list:
        print(f"  {label}: {len(run_dirs)} run(s)")
        all_run_dirs.extend(run_dirs)
        runs = [compute_run(d) for d in run_dirs]
        runs = [r for r in runs if r is not None]
        data = {k: rep_mean(runs, k) for k in KEYS}
        # SD across replicates for taxon series
        sd   = {k: rep_std(runs,  k) for k in KEYS}
        scenarios.append((label, data, sd))

    _rt_dir = all_run_dirs[0] if all_run_dirs else None
    _rt = read_runtime(_rt_dir)
    N_YEARS       = _rt["pulse"]
    DAYS_PER_YEAR = int(_rt["end_time"])
    COL_START, COL_END = year_range(_args.years, N_YEARS, DAYS_PER_YEAR)

    scenarios = [
        (label,
         {k: filter_cols(v, COL_START, COL_END) for k, v in data.items()},
         {k: filter_cols(v, COL_START, COL_END) for k, v in sd.items()})
        for label, data, sd in scenarios
    ]

    print(f"Simulation: {N_YEARS} year(s), {DAYS_PER_YEAR} days/year  |  "
          f"Displaying: day {COL_START}–{COL_END}")

if not scenarios:
    sys.exit("ERROR: no scenario data loaded.")

# ── Colours ────────────────────────────────────────────────────────────────────
scen_labels = [label for label, _, _ in scenarios]
COLOR = scenario_color_map(scen_labels, palette=get_palette_from_args(_args))

# ── Legend layout helper ───────────────────────────────────────────────────────
def _legend_ncols(n: int) -> int:
    if n <= 4:  return n
    if n <= 8:  return 4
    return 5


# ────────────────────────────────────────────────────────────────────────────────
# Fig 1: System-level CUE and CAE — all scenarios overlaid
# ────────────────────────────────────────────────────────────────────────────────
fig1, axes1 = plt.subplots(1, 2, figsize=(13, 5))
fig1.suptitle(
    "System-Level CUE and CAE\n"
    "(community mean ± 1 SD across replicates)",
    fontsize=11, fontweight="bold",
)

for ax, key, ylabel, title in [
        (axes1[0], "CUE_sys", "Carbon Use Efficiency",
         "CUE  (growth C / assimilated C)"),
        (axes1[1], "CAE_sys", "Carbon Assimilation Efficiency",
         "CAE  (assimilated C / uptaken C)"),
]:
    has_data = False
    for label, data, sd in scenarios:
        df    = data.get(key)
        df_sd = sd.get(key)
        if df is None:
            continue

        s    = df.iloc[0]
        c    = COLOR[label]
        ax.plot(s.index, s.values, color=c, lw=1.6, label=label)

        # ±1 SD ribbon across replicates (from summary SD or rep_std)
        if df_sd is not None:
            s_sd = df_sd.iloc[0] if hasattr(df_sd, "iloc") else df_sd
            ax.fill_between(s.index,
                            (s - s_sd).clip(lower=0).values,
                            (s + s_sd).clip(upper=1).values,
                            color=c, alpha=0.15, linewidth=0)
        has_data = True

    if not has_data:
        ax.set_title(f"{title} — no data")
        continue

    ax.axhline(0.5, ls="--", lw=0.7, color="gray", alpha=0.5)
    ax.set_xlabel("Day")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.set_ylim(0, 1)
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(loc="upper center",
              bbox_to_anchor=(0.5, -0.20),
              ncols=_legend_ncols(len(scen_labels)),
              framealpha=0.85, fontsize=8)
    add_year_labels(ax, COL_START, COL_END, N_YEARS, DAYS_PER_YEAR)

fig1.tight_layout()
fig1.savefig(FIG_DIR / "fig1_system_cue_cae.png", dpi=150, bbox_inches="tight")
plt.close(fig1); print("Saved fig1")

# ────────────────────────────────────────────────────────────────────────────────
# Fig 2: Community-mean taxon CUE and CAE — all scenarios with ±1 SD
# ────────────────────────────────────────────────────────────────────────────────
fig2, axes2 = plt.subplots(1, 2, figsize=(13, 5))
fig2.suptitle(
    "Community-Mean Taxon CUE and CAE\n"
    "(mean ± 1 SD across alive taxa)",
    fontsize=11, fontweight="bold",
)

for ax, key, ylabel, title in [
        (axes2[0], "CUE_tax", "CUE (per taxon)",  "Per-Taxon CUE"),
        (axes2[1], "CAE_tax", "CAE (per taxon)",  "Per-Taxon CAE"),
]:
    has_data = False
    for label, data, sd in scenarios:
        df   = data.get(key)
        df_sd = sd.get(key)
        if df is None:
            continue

        mu, sig = comm_ts(df)
        if mu is None:
            continue

        # If per-scenario SD was loaded from summary, use it for the ribbon;
        # otherwise fall back to across-taxon SD from comm_ts.
        if df_sd is not None:
            _, sig = comm_ts(df_sd)

        c = COLOR[label]
        ax.plot(mu.index, mu.values, color=c, lw=1.6, label=label)
        if sig is not None:
            ax.fill_between(mu.index,
                            (mu - sig).clip(lower=0).values,
                            (mu + sig).clip(upper=1).values,
                            color=c, alpha=0.15, linewidth=0)
        has_data = True

    if not has_data:
        ax.set_title(f"{title} — no data")
        continue

    ax.axhline(0.5, ls="--", lw=0.7, color="gray", alpha=0.5)
    ax.set_xlabel("Day")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.set_ylim(0, 1)
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(loc="upper center",
              bbox_to_anchor=(0.5, -0.20),
              ncols=_legend_ncols(len(scen_labels)),
              framealpha=0.85, fontsize=8)
    add_year_labels(ax, COL_START, COL_END, N_YEARS, DAYS_PER_YEAR)

fig2.tight_layout()
fig2.savefig(FIG_DIR / "fig2_comm_mean_cue_cae.png", dpi=150, bbox_inches="tight")
plt.close(fig2); print("Saved fig2")

print(f"\nAll figures saved to: {FIG_DIR}")
