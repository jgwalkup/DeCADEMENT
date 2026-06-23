#!/usr/bin/env python3
"""
cue_figures.py

Carbon Use Efficiency (CUE) and Carbon Assimilation Efficiency (CAE) figures
for DEMENTpy.

Definitions
-----------
CUE  Carbon Use Efficiency
     = growth C / assimilated C
     Fraction of assimilated C that is incorporated into new biomass (0–1).
CAE  Carbon Assimilation Efficiency
     = assimilated C / uptaken C
     Fraction of uptaken C that is assimilated (not immediately respired for
     transport costs) (0–1).

CSV sources (all in each run directory)
----------------------------------------
CUE_system.csv       – 1 row × days; emergent system-level CUE
CAE_system.csv       – 1 row × days; emergent system-level CAE
CUE_TaxonSeries.csv  – n_taxa × days; per-taxon CUE (divided by taxon count)
CAE_TaxonSeries.csv  – n_taxa × days; per-taxon CAE (divided by taxon count)

Figures
-------
Fig 1 – System-level CUE and CAE: ctrl vs warm overlaid (1 × 2 panels)
Fig 2 – Community-mean taxon CUE and CAE: ctrl vs warm with ±1 SD (1 × 2 panels)
Fig 3 – Per-taxon mean CUE bar chart sorted by ctrl mean (1 × 2 panels ctrl/warm)
Fig 4 – Per-taxon CUE and CAE heatmaps: taxa × days (2 × 2 panels CUE/CAE × ctrl/warm)
Fig 5 – Warm − Control Δ CUE and Δ CAE per taxon (1 × 2 panels)
"""

import argparse
import sys
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from pathlib import Path

warnings.filterwarnings("ignore")

# Make scenario_utils importable (shared loaders live there)
_HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(_HERE))
try:
    import scenario_utils as su
except ImportError as _e:
    su = None
    _SU_ERR = _e

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
C_CTRL = "#1f77b4"   # blue — control
C_WARM = "#d62728"   # red  — warm
C_CUE  = "#2ca02c"   # green — CUE accent
C_CAE  = "#ff7f0e"   # orange — CAE accent

plt.rcParams.update({
    "font.size": 10, "axes.titlesize": 11, "axes.labelsize": 10,
    "legend.fontsize": 9, "figure.dpi": 150,
})

# ── Helpers ───────────────────────────────────────────────────────────────────
def load(run_dir: Path, name: str):
    f = run_dir / f"{name}.csv"
    if not f.exists():
        return None
    df = pd.read_csv(f, index_col=0)
    # Positional renaming: handles both correct (int) and broken ('C','C.1',...) column names
    df.columns = range(len(df.columns))        # 0 = init state, 1..n = days
    df = df.drop(columns=[0], errors="ignore") # drop init state column
    df.columns = range(1, len(df.columns) + 1) # renumber 1..n_days
    # CUE/CAE per-taxon series are divided by taxon_count in output.py.  Taxa
    # absent from the grid on a given day (count 0) yield +/-inf (nonzero/0) or
    # NaN (0/0).  These are absent-taxon artifacts, not real efficiencies, so we
    # convert inf to NaN; the skipna reductions then average only over present
    # taxa.  (Runs produced after the output.py inf-cleanup patch already store
    # NaN, so this is a no-op for them and a fix for older output.)
    df = df.replace([np.inf, -np.inf], np.nan)
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

def rep_mean(pool_list: list, key: str) -> pd.DataFrame | None:
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

def comm_ts(df: pd.DataFrame | None) -> tuple:
    """Community mean ± 1 SD across taxa (axis=0) at each day."""
    if df is None:
        return None, None
    mu  = df.mean(axis=0, skipna=True)
    sig = df.std(axis=0,  skipna=True)
    return mu, sig

def taxon_mean(df: pd.DataFrame | None) -> pd.Series | None:
    """Time-mean per taxon (axis=1)."""
    if df is None:
        return None
    return df.mean(axis=1, skipna=True)

def subsample(df: pd.DataFrame | None, n: int = 200) -> pd.DataFrame | None:
    if df is None or df.shape[1] <= n:
        return df
    idx = np.linspace(0, df.shape[1] - 1, n, dtype=int)
    return df.iloc[:, idx]

# ── Year / pulse helpers ───────────────────────────────────────────────────────
def read_runtime(run_dir) -> dict:
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

def add_year_labels_heatmap(ax, mat_ncols: int, col_start: int, col_end: int,
                             n_years: int, days_per_year: int) -> None:
    total = col_end - col_start + 1
    for yr in range(1, n_years + 1):
        boundary = yr * days_per_year
        if col_start < boundary < col_end:
            frac  = (boundary - col_start) / total
            x_pos = frac * (mat_ncols - 1)
            ax.axvline(x_pos, color="white", lw=0.8, ls=":", alpha=0.6, zorder=5)

# ── Per-run data loading ───────────────────────────────────────────────────────
def compute_run(run_dir: Path) -> dict | None:
    CUE_sys = load(run_dir, "CUE_system")
    CAE_sys = load(run_dir, "CAE_system")
    CUE_tax = load(run_dir, "CUE_TaxonSeries")
    CAE_tax = load(run_dir, "CAE_TaxonSeries")

    if all(x is None for x in [CUE_sys, CAE_sys, CUE_tax, CAE_tax]):
        return None

    # Clip efficiency values to [0, 1]; values outside this range reflect
    # transient numerical artefacts (e.g. taxa with zero count).
    def clip01(df):
        return df.clip(lower=0, upper=1) if df is not None else None

    return dict(
        CUE_sys=clip01(CUE_sys),
        CAE_sys=clip01(CAE_sys),
        CUE_tax=clip01(CUE_tax),
        CAE_tax=clip01(CAE_tax),
    )

# ── Parse CLI ──────────────────────────────────────────────────────────────────
_p = argparse.ArgumentParser(
    description="CUE / CAE figures for DEMENTpy",
    formatter_class=argparse.RawDescriptionHelpFormatter)
_p.add_argument("--base",        default=None)
_p.add_argument("--control-dir", default=None)
_p.add_argument("--warm-dir",    default=None)
_p.add_argument("--control",     nargs="*", default=None, metavar="FOLDER")
_p.add_argument("--warm",        nargs="*", default=None, metavar="FOLDER")
_p.add_argument("--fig-dir",     default=None)
_p.add_argument("--summary-dir", default=None, metavar="DIR",
                help="Read pre-computed summary CSVs (from export_summary_data.py) "
                     "instead of raw runs. --control-dir/--warm-dir then name the "
                     "scenario summary subdirs (by name under --summary-dir, or as "
                     "absolute paths).")
_p.add_argument("--days-per-year", type=int, default=365, metavar="N",
                help="Days per simulation year, used only in --summary-dir mode "
                     "for year-boundary labels (default: 365).")
_p.add_argument("--years",       nargs="+", type=int, default=None, metavar="YEAR",
                help="Year(s) to display. Single value = that year only. "
                     "Two values = inclusive range (e.g. --years 2 4). "
                     "Default: show all years.")
_args = _p.parse_args()

if _args.base:
    BASE = Path(_args.base)

FIG_DIR = Path(_args.fig_dir) if _args.fig_dir else BASE / "cue_figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

KEYS = ["CUE_sys", "CAE_sys", "CUE_tax", "CAE_tax"]

# Map this script's internal keys → export_summary_data.py stems
_SUMMARY_STEMS = {
    "CUE_sys": "system_CUE",
    "CAE_sys": "system_CAE",
    "CUE_tax": "cue_taxon",
    "CAE_tax": "cae_taxon",
}


def _summary_scenario_dir(spec, summary_root: Path, default_name: str) -> Path:
    """Resolve a ctrl/warm spec to a summary scenario directory."""
    if spec:
        p = Path(spec)
        if p.is_absolute():
            return p
        return summary_root / spec
    return summary_root / default_name


def _build_from_summary(sdir: Path) -> dict:
    """Build the ctrl/warm dict from one scenario's summary CSVs."""
    out = {}
    for key, stem in _SUMMARY_STEMS.items():
        mean, _sd = su.load_summary(sdir, stem)
        if mean is None:
            out[key] = None
        elif isinstance(mean, pd.Series):
            # 1-D system series → 1-row DataFrame to match raw-run shape (.iloc[0])
            out[key] = pd.DataFrame([mean.values], columns=list(mean.index))
        else:
            out[key] = mean
    # Match raw-run compute_run: clip efficiencies to [0, 1]
    for key in out:
        if out[key] is not None:
            out[key] = out[key].clip(lower=0, upper=1)
    return out


if _args.summary_dir:
    # ── Summary mode ───────────────────────────────────────────────────────────
    if su is None:
        sys.exit(f"ERROR: --summary-dir requires scenario_utils: {_SU_ERR}")
    summary_root = Path(_args.summary_dir)
    ctrl_sdir = _summary_scenario_dir(_args.control_dir, summary_root, "control")
    warm_sdir = _summary_scenario_dir(_args.warm_dir,    summary_root, "warm")
    print(f"Summary mode")
    print(f"  Control summary: {ctrl_sdir}")
    print(f"  Warm summary:    {warm_sdir}")

    ctrl = _build_from_summary(ctrl_sdir)
    warm = _build_from_summary(warm_sdir)

    # Determine day range directly from the loaded data (export already applied
    # any --years filtering at export time, so we display what is present).
    def _day_span(d):
        spans = []
        for v in d.values():
            if v is not None and v.shape[1] > 0:
                cols = [c for c in v.columns if isinstance(c, (int, np.integer))]
                if cols:
                    spans.append((min(cols), max(cols)))
        return spans

    spans = _day_span(ctrl) + _day_span(warm)
    if spans:
        COL_START = min(s for s, _ in spans)
        COL_END   = max(e for _, e in spans)
    else:
        COL_START, COL_END = 1, _args.days_per_year
    DAYS_PER_YEAR = int(_args.days_per_year)
    N_YEARS = max(1, -(-COL_END // DAYS_PER_YEAR))   # ceil division
    print(f"Simulation: inferred {N_YEARS} year(s), {DAYS_PER_YEAR} days/year  |  "
          f"Displaying: day {COL_START}–{COL_END}")

else:
    # ── Raw-run mode (original behaviour) ───────────────────────────────────────
    CTRL_BASE = Path(_args.control_dir) if _args.control_dir else BASE / "control"
    WARM_BASE = Path(_args.warm_dir)    if _args.warm_dir    else BASE / "warm"

    if _args.control is not None:
        ctrl_dirs = resolve_runs(_args.control, CTRL_BASE)
    elif CTRL_RUNS:
        ctrl_dirs = resolve_runs(CTRL_RUNS, CTRL_BASE)
    else:
        ctrl_dirs = pick_runs(CTRL_BASE)

    if _args.warm is not None:
        warm_dirs = resolve_runs(_args.warm, WARM_BASE)
    elif WARM_RUNS:
        warm_dirs = resolve_runs(WARM_RUNS, WARM_BASE)
    else:
        warm_dirs = pick_runs(WARM_BASE)

    print(f"Control: {[d.name for d in ctrl_dirs]}")
    print(f"Warm:    {[d.name for d in warm_dirs]}")

    ctrl_runs = [compute_run(d) for d in ctrl_dirs]
    warm_runs = [compute_run(d) for d in warm_dirs]

    ctrl = {k: rep_mean(ctrl_runs, k) for k in KEYS}
    warm = {k: rep_mean(warm_runs, k) for k in KEYS}

    # ── Runtime info and year filtering ───────────────────────────────────────
    _rt_dir = ctrl_dirs[0] if ctrl_dirs else (warm_dirs[0] if warm_dirs else None)
    _rt = read_runtime(_rt_dir)
    N_YEARS       = _rt["pulse"]
    DAYS_PER_YEAR = int(_rt["end_time"])
    COL_START, COL_END = year_range(_args.years, N_YEARS, DAYS_PER_YEAR)
    print(f"Simulation: {N_YEARS} year(s), {DAYS_PER_YEAR} days/year  |  "
          f"Displaying: day {COL_START}–{COL_END}")

    ctrl = {k: filter_cols(v, COL_START, COL_END) for k, v in ctrl.items()}
    warm = {k: filter_cols(v, COL_START, COL_END) for k, v in warm.items()}

# ────────────────────────────────────────────────────────────────────────────
# Fig 1: System-level CUE and CAE — ctrl vs warm overlaid
# ────────────────────────────────────────────────────────────────────────────
fig1, axes1 = plt.subplots(1, 2, figsize=(13, 5))
fig1.suptitle("System-Level CUE and CAE Over Time — Control vs Warm",
              fontsize=12, fontweight="bold")

for ax, key, ylabel, title in [
        (axes1[0], "CUE_sys", "Carbon Use Efficiency",        "CUE  (growth C / assimilated C)"),
        (axes1[1], "CAE_sys", "Carbon Assimilation Efficiency","CAE  (assimilated C / uptaken C)"),
]:
    ctrl_s = ctrl[key].iloc[0] if ctrl[key] is not None else None
    warm_s = warm[key].iloc[0] if warm[key] is not None else None

    if ctrl_s is None and warm_s is None:
        ax.set_title(f"{title} — no data"); continue

    if ctrl_s is not None:
        ax.plot(ctrl_s.index, ctrl_s.values, color=C_CTRL, lw=1.6,
                label="Control")
    if warm_s is not None:
        ax.plot(warm_s.index, warm_s.values, color=C_WARM, lw=1.6,
                label="Warm")

    ax.axhline(0.5, ls="--", lw=0.7, color="gray", alpha=0.5)
    ax.set_xlabel("Day")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.set_ylim(0, 1)
    days_ref = ctrl_s.index if ctrl_s is not None else warm_s.index
    ax.set_xlim(days_ref.min(), days_ref.max())
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.20),
              ncols=2, framealpha=0.85, fontsize=8)
    add_year_labels(ax, COL_START, COL_END, N_YEARS, DAYS_PER_YEAR)

fig1.tight_layout()
fig1.savefig(FIG_DIR / "fig1_system_cue_cae.png", dpi=150, bbox_inches="tight")
plt.close(fig1); print("Saved fig1")

# ────────────────────────────────────────────────────────────────────────────
# Fig 2: Community-mean taxon CUE and CAE — ctrl vs warm with ±1 SD
# ────────────────────────────────────────────────────────────────────────────
fig2, axes2 = plt.subplots(1, 2, figsize=(13, 5))
fig2.suptitle("Community-Mean Taxon CUE and CAE Over Time\n"
              "(mean ± 1 SD across alive taxa; ctrl vs warm overlaid)",
              fontsize=12, fontweight="bold")

for ax, key, ylabel, title in [
        (axes2[0], "CUE_tax", "CUE (per taxon)",  "Per-Taxon CUE"),
        (axes2[1], "CAE_tax", "CAE (per taxon)",  "Per-Taxon CAE"),
]:
    ctrl_mu, ctrl_sd = comm_ts(ctrl[key])
    warm_mu, warm_sd = comm_ts(warm[key])

    if ctrl_mu is None and warm_mu is None:
        ax.set_title(f"{title} — no data"); continue

    for mu, sd, color, label in [
            (ctrl_mu, ctrl_sd, C_CTRL, "Control"),
            (warm_mu, warm_sd, C_WARM, "Warm"),
    ]:
        if mu is None:
            continue
        days = mu.index
        ax.plot(days, mu.values, color=color, lw=1.6, label=label)
        ax.fill_between(days,
                        (mu - sd).clip(0).values,
                        (mu + sd).clip(upper=1).values,
                        color=color, alpha=0.18, linewidth=0)

    ax.axhline(0.5, ls="--", lw=0.7, color="gray", alpha=0.5)
    ax.set_xlabel("Day")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.set_ylim(0, 1)
    _ref = ctrl_mu if ctrl_mu is not None else warm_mu
    ax.set_xlim(_ref.index.min(), _ref.index.max())
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.20),
              ncols=2, framealpha=0.85, fontsize=8)
    add_year_labels(ax, COL_START, COL_END, N_YEARS, DAYS_PER_YEAR)

fig2.tight_layout()
fig2.savefig(FIG_DIR / "fig2_comm_mean_cue_cae.png", dpi=150, bbox_inches="tight")
plt.close(fig2); print("Saved fig2")

# ────────────────────────────────────────────────────────────────────────────
# Fig 3: Per-taxon mean CUE bar chart — sorted by ctrl mean
# ────────────────────────────────────────────────────────────────────────────
ctrl_cue_m = taxon_mean(ctrl["CUE_tax"])
warm_cue_m = taxon_mean(warm["CUE_tax"])
ctrl_cae_m = taxon_mean(ctrl["CAE_tax"])
warm_cae_m = taxon_mean(warm["CAE_tax"])

# Sort order by ctrl CUE mean; fall back to warm if ctrl unavailable
_ref_m = ctrl_cue_m if ctrl_cue_m is not None else warm_cue_m
order3 = _ref_m.sort_values(ascending=True).index if _ref_m is not None else None

if order3 is not None:
    fig3, axes3 = plt.subplots(1, 2, figsize=(14, 0.4 * len(order3) + 3),
                                sharey=True)
    fig3.suptitle("Per-Taxon Time-Mean CUE and CAE\n"
                  "(averaged over all days; taxa sorted by control CUE)",
                  fontsize=12, fontweight="bold")

    yticks3 = np.arange(len(order3))

    for ax, cue_m, cae_m, title in [
            (axes3[0], ctrl_cue_m, ctrl_cae_m, "Control"),
            (axes3[1], warm_cue_m, warm_cae_m, "Warm"),
    ]:
        if cue_m is None:
            ax.set_title(f"{title} — no data"); continue
        cue_vals = cue_m.reindex(order3).fillna(np.nan).values
        ax.barh(yticks3, cue_vals, color=C_CUE, alpha=0.80, label="CUE",
                edgecolor="white", linewidth=0.3)
        if cae_m is not None:
            cae_vals = cae_m.reindex(order3).fillna(np.nan).values
            ax.scatter(cae_vals, yticks3, color=C_CAE, s=18, zorder=4,
                       label="CAE", marker="D")
        ax.axvline(0.5, ls="--", lw=0.8, color="gray", alpha=0.5)
        ax.set_xlabel("Efficiency (0–1)")
        ax.set_xlim(0, 1)
        ax.set_title(title)
        ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.20),
                  ncols=2, framealpha=0.85, fontsize=8)
        ax.spines[["top", "right"]].set_visible(False)

    axes3[0].set_yticks(yticks3)
    axes3[0].set_yticklabels(order3, fontsize=7)

    fig3.tight_layout()
    fig3.savefig(FIG_DIR / "fig3_taxon_mean_cue_cae_bars.png", dpi=150, bbox_inches="tight")
    plt.close(fig3); print("Saved fig3")
else:
    print("Skipped fig3 — no taxon CUE data")

# ────────────────────────────────────────────────────────────────────────────
# Fig 4: Per-taxon CUE and CAE heatmaps (taxa × days)
# ────────────────────────────────────────────────────────────────────────────
# Sort rows by ctrl mean CUE descending (highest CUE at the top)
_sort_ref = taxon_mean(ctrl["CUE_tax"])
if _sort_ref is None:
    _sort_ref = taxon_mean(warm["CUE_tax"])
order4 = _sort_ref.sort_values(ascending=False).index if _sort_ref is not None else None

if order4 is not None:
    day_labels4 = np.linspace(COL_START, COL_END, 6, dtype=int)

    fig4, axes4 = plt.subplots(2, 2, figsize=(16, 0.35 * len(order4) + 4),
                                squeeze=False)
    fig4.suptitle(
        "Per-Taxon CUE and CAE Heatmaps  "
        "(rows = taxa sorted by ctrl CUE; cols = days)",
        fontsize=12, fontweight="bold",
    )

    panels = [
        (axes4[0, 0], ctrl["CUE_tax"], "Control — CUE", "YlGn",   "CUE"),
        (axes4[0, 1], warm["CUE_tax"], "Warm — CUE",    "YlGn",   "CUE"),
        (axes4[1, 0], ctrl["CAE_tax"], "Control — CAE", "YlOrBr", "CAE"),
        (axes4[1, 1], warm["CAE_tax"], "Warm — CAE",    "YlOrBr", "CAE"),
    ]

    for ax, df, title, cmap, cbar_label in panels:
        df_sub = subsample(df)
        if df_sub is None:
            ax.set_title(f"{title} — no data")
            ax.axis("off")
            continue
        mat = df_sub.reindex(order4).fillna(np.nan).values
        im = ax.imshow(mat, aspect="auto", cmap=cmap, vmin=0, vmax=1,
                       interpolation="nearest")
        ax.set_title(title, fontsize=10)
        ax.set_xlabel("Day (subsampled to 200)")
        ax.set_ylabel("Taxon")
        ax.set_yticks(np.arange(len(order4)))
        ax.set_yticklabels(order4, fontsize=5)
        xtick_pos = np.linspace(0, mat.shape[1] - 1, 6, dtype=int)
        ax.set_xticks(xtick_pos)
        ax.set_xticklabels(day_labels4)
        plt.colorbar(im, ax=ax, label=cbar_label, shrink=0.8)
        add_year_labels_heatmap(ax, mat.shape[1], COL_START, COL_END,
                                N_YEARS, DAYS_PER_YEAR)

    fig4.tight_layout()
    fig4.savefig(FIG_DIR / "fig4_taxon_cue_cae_heatmaps.png", dpi=150, bbox_inches="tight")
    plt.close(fig4); print("Saved fig4")
else:
    print("Skipped fig4 — no taxon efficiency data")

# ────────────────────────────────────────────────────────────────────────────
# Fig 5: Warm − Control Δ CUE and Δ CAE per taxon
# ────────────────────────────────────────────────────────────────────────────
delta_cue = None
delta_cae = None

if ctrl["CUE_tax"] is not None and warm["CUE_tax"] is not None:
    ctrl_m = taxon_mean(ctrl["CUE_tax"])
    warm_m = taxon_mean(warm["CUE_tax"])
    delta_cue = (warm_m - ctrl_m).dropna()

if ctrl["CAE_tax"] is not None and warm["CAE_tax"] is not None:
    ctrl_m = taxon_mean(ctrl["CAE_tax"])
    warm_m = taxon_mean(warm["CAE_tax"])
    delta_cae = (warm_m - ctrl_m).dropna()

if delta_cue is not None or delta_cae is not None:
    # Determine sort order from whichever delta is available
    _delta_ref = delta_cue if delta_cue is not None else delta_cae
    order5 = _delta_ref.sort_values(ascending=True).index
    yticks5 = np.arange(len(order5))

    n_panels = sum(x is not None for x in [delta_cue, delta_cae])
    fig5, axs5 = plt.subplots(1, n_panels,
                               figsize=(7 * n_panels, 0.4 * len(order5) + 3),
                               sharey=True)
    if n_panels == 1:
        axs5 = [axs5]
    fig5.suptitle("Warm − Control Difference in CUE and CAE Per Taxon\n"
                  "(positive = higher under warming; taxa sorted by Δ CUE)",
                  fontsize=12, fontweight="bold")

    panel_idx = 0
    for delta, xlabel, title in [
            (delta_cue, "Δ CUE  (Warm − Control)", "Δ Carbon Use Efficiency"),
            (delta_cae, "Δ CAE  (Warm − Control)", "Δ Carbon Assimilation Efficiency"),
    ]:
        if delta is None:
            continue
        ax = axs5[panel_idx]
        vals = delta.reindex(order5).fillna(0).values
        colors = [C_WARM if v > 0 else C_CTRL for v in vals]
        ax.barh(yticks5, vals, color=colors, alpha=0.85,
                edgecolor="white", linewidth=0.3)
        ax.axvline(0, color="black", lw=0.9)
        ax.set_xlabel(xlabel)
        ax.set_title(title)
        ax.spines[["top", "right"]].set_visible(False)
        panel_idx += 1

    axs5[0].set_yticks(yticks5)
    axs5[0].set_yticklabels(order5, fontsize=7)

    legend_handles = [Patch(facecolor=C_WARM, label="Higher under warming"),
                      Patch(facecolor=C_CTRL, label="Lower under warming")]
    axs5[-1].legend(handles=legend_handles, loc="upper center",
                    bbox_to_anchor=(0.5, -0.20), ncols=2, framealpha=0.85, fontsize=8)

    fig5.tight_layout()
    fig5.savefig(FIG_DIR / "fig5_warm_ctrl_delta_cue_cae.png", dpi=150, bbox_inches="tight")
    plt.close(fig5); print("Saved fig5")
else:
    print("Skipped fig5 — need both ctrl and warm taxon efficiency data")

print(f"\nAll figures saved to: {FIG_DIR}")
