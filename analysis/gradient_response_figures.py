#!/usr/bin/env python3
"""
gradient_response_figures.py

Temperature-gradient response curves for DeCA-DEMENT ensemble simulations.

Reads output from a gradient experiment structured as:
    <gradient-dir>/
        static_T{temp}_P{psi}/
            <community_run_1>/
            <community_run_2>/
            ...
        static_T{temp}_P{psi}/
            ...

For each temperature point the ensemble mean (± 1 SD across replicate community
runs) of each response variable is computed, then plotted as a function of
temperature.

Response variables
------------------
  C pools:
    - Total substrate C
    - Total monomer C
    - Total microbial biomass C
    - Cumulative respiration

  Carbon use efficiency:
    - System CUE
    - System CAE
    - Community-mean taxon CUE
    - Community-mean taxon CAE

  C pathways (community mean per alive microbe):
    - C uptake
    - Biosynthesis C
    - Enzyme C cost
    - HSP C cost
    - Osmolyte C cost

  Energy pathways (community mean per alive microbe):
    - E production
    - E gain (net)
    - E cost (total)
    - Transport E cost
    - Maintenance E cost
    - Enzyme E cost
    - HSP E cost
    - Osmolyte E cost

Figures produced
----------------
  Fig 1  – C pool response curves (substrate / monomer / microbial biomass)
  Fig 2  – Respiration response curve
  Fig 3  – CUE and CAE (system-level)
  Fig 4  – Community-mean taxon CUE and CAE
  Fig 5  – C pathway response curves
  Fig 6  – E pathway response curves
  Fig 7  – E cost breakdown response curves

Each figure shows ensemble mean ± 1 SD ribbons.

Usage
-----
  python analysis/gradient_response_figures.py \\
      --gradient-dir simulations/temp_gradient \\
      --fig-dir      figures/temp_gradient \\
      [--years 1 2]          # restrict to specific simulation year(s)
      [--summary last_year]  # how to collapse time into one point per temp:
                             #   mean (default), last_year, last_pulse
      [--min-cols 100]       # skip runs with fewer day-columns (default: 100)
"""

import argparse
import re
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm

warnings.filterwarnings("ignore")

# ── Style ──────────────────────────────────────────────────────────────────────
plt.rcParams.update({
    "font.size": 10, "axes.titlesize": 11, "axes.labelsize": 10,
    "legend.fontsize": 9, "figure.dpi": 150,
})

CMAP_GRADIENT = "plasma"

# Variable colours (match existing figure scripts for consistency)
C_SUB   = "#8c564b"
C_MON   = "#ff7f0e"
C_MIC   = "#1f77b4"
C_RESP  = "#d62728"
C_CUE   = "#2ca02c"
C_CAE   = "#ff7f0e"
C_UPTK  = "#2ca02c"
C_BIOSYN = "#1f77b4"
C_ENZ   = "#e377c2"
C_HSP   = "#8c564b"
C_OSMO  = "#bcbd22"
C_EPROD = "#2ca02c"
C_EGAIN = "#1f77b4"
C_ECOST = "#d62728"
C_TRANS = "#ff7f0e"
C_MAINT = "#9467bd"

# ── I/O helpers ───────────────────────────────────────────────────────────────

def _load_csv(run_dir: Path, name: str) -> pd.DataFrame | None:
    """Load a CSV output file from a run directory, normalising column indices."""
    f = run_dir / f"{name}.csv"
    if not f.exists():
        return None
    try:
        df = pd.read_csv(f, index_col=0)
    except Exception:
        return None
    # Normalise columns to integers 1..n (drop init-state column 0)
    df.columns = range(len(df.columns))
    df = df.drop(columns=[0], errors="ignore")
    df.columns = range(1, len(df.columns) + 1)
    # Per-taxon ratio outputs (CUE/CAE and pathway-per-count quantities) are
    # divided by taxon_count in output.py.  Taxa absent from the grid on a given
    # day have count 0, yielding +/-inf (nonzero numerator) or NaN (0/0).  These
    # are not real values — they are absent-taxon artifacts — so we convert inf
    # to NaN and let the nan-aware reductions average only over present taxa.
    df = df.replace([np.inf, -np.inf], np.nan)
    return df


def _pick_runs(scenario_dir: Path, min_cols: int) -> list[Path]:
    """
    Return one run directory per replicate in scenario_dir.
    Replicates are identified by the second underscore-delimited token of the
    folder name (e.g. community_0001_s42 → rep '0001').
    The most-recent / most-complete run per replicate is kept.
    """
    if not scenario_dir.exists():
        return []
    by_rep: dict[str, Path] = {}
    for d in sorted(scenario_dir.iterdir()):
        if not d.is_dir():
            continue
        # Probe a file that is always present for completed runs
        probe = d / "SubstratesSeries.csv"
        if not probe.exists():
            probe = d / "MicrobesSeries.csv"
        if not probe.exists():
            probe = d / "C_Uptake_Taxon.csv"
        if not probe.exists():
            continue
        try:
            nc = pd.read_csv(probe, index_col=0, nrows=1).shape[1]
        except Exception:
            continue
        if nc < min_cols:
            continue
        parts = d.name.split("_")
        rep = parts[1] if len(parts) > 1 else d.name
        if rep not in by_rep or d.name > by_rep[rep].name:
            by_rep[rep] = d
    return list(by_rep.values())


def _read_runtime(run_dir: Path) -> dict:
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


def _col_range(years_arg, n_years: int, days_per_year: int) -> tuple[int, int]:
    if years_arg is None:
        return 1, n_years * days_per_year
    if len(years_arg) == 1:
        yr = years_arg[0]
        return (yr - 1) * days_per_year + 1, yr * days_per_year
    return (years_arg[0] - 1) * days_per_year + 1, years_arg[-1] * days_per_year


def _filter_cols(df: pd.DataFrame | None, lo: int, hi: int) -> pd.DataFrame | None:
    if df is None:
        return None
    cols = [c for c in df.columns if lo <= c <= hi]
    return df[cols] if cols else df


def _last_year_cols(df: pd.DataFrame | None, n_years: int,
                    days_per_year: int) -> pd.DataFrame | None:
    lo = (n_years - 1) * days_per_year + 1
    hi = n_years * days_per_year
    return _filter_cols(df, lo, hi)


# ── Gradient directory parsing ─────────────────────────────────────────────────

# Matches both literal form  static_T3_P-0.5
# and encoded form           static_T3_Pn0p5  (n=minus, p=decimal point)
_STATIC_RE_LITERAL  = re.compile(
    r"static_T([+-]?\d+(?:\.\d+)?)_P([+-]?\d+(?:\.\d+)?)"
)
_STATIC_RE_ENCODED  = re.compile(
    r"static_T(\d+(?:p\d+)?)_P(n?)(\d+(?:p\d+)?)"
)


def _decode_encoded_value(sign: str, digits: str) -> float:
    """Convert n0p5 -> -0.5, 0p5 -> 0.5, 3 -> 3.0."""
    val = float(digits.replace("p", "."))
    return -val if sign == "n" else val


def _parse_treatment_dirs(gradient_dir: Path) -> list[tuple[float, Path]]:
    """
    Return list of (temperature, scenario_dir) tuples sorted by temperature,
    for all subdirectories matching the static_T{T}_P{P} naming convention.
    Handles both literal (static_T3_P-0.5) and encoded (static_T3_Pn0p5) forms.
    """
    results = []
    for d in gradient_dir.iterdir():
        if not d.is_dir():
            continue
        # Try literal form first
        m = _STATIC_RE_LITERAL.match(d.name)
        if m:
            temp = float(m.group(1))
            results.append((temp, d))
            continue
        # Try encoded form
        m = _STATIC_RE_ENCODED.match(d.name)
        if m:
            temp_str = m.group(1).replace("p", ".")
            temp = float(temp_str)
            results.append((temp, d))
    results.sort(key=lambda x: x[0])
    return results


# ── Per-run variable extraction ────────────────────────────────────────────────

def _extract_run(run_dir: Path, col_lo: int, col_hi: int) -> dict:
    """
    Load all tracked variables from a single run directory, restricted to
    [col_lo, col_hi], and return them as a dict of DataFrames/Series.
    """

    def lf(name):
        df = _load_csv(run_dir, name)
        return _filter_cols(df, col_lo, col_hi)

    return {
        # C pools (types × days)
        "SUB":        lf("SubstratesSeries"),
        "MON":        lf("MonomersSeries"),
        "MIC":        lf("MicrobesSeries"),
        # Respiration (1 × days)
        "RESP":       lf("RespSeries"),
        # CUE / CAE system (1 × days)
        "CUE_sys":    lf("CUE_system"),
        "CAE_sys":    lf("CAE_system"),
        # CUE / CAE taxon (taxa × days)
        "CUE_tax":    lf("CUE_TaxonSeries"),
        "CAE_tax":    lf("CAE_TaxonSeries"),
        # C pathway (taxa × days)
        "C_uptake":   lf("C_Uptake_Taxon"),
        "BC_uptake":  lf("BC_Uptake_Taxon"),
        "Enz_C":      lf("Enzyme_C_Taxon"),
        "HSP_C":      lf("HSP_C_Taxon"),
        "Osmo_C":     lf("Osmo_C_Taxon"),
        # E pathway (taxa × days)
        "E_prod":     lf("E_Prod_Taxon"),
        "E_gain":     lf("Microbe_E_Gain"),
        "E_cost":     lf("Microbe_E_Cost"),
        "Trans_E":    lf("Transport_E_Taxon"),
        "Maint_E":    lf("Maint_E_Taxon"),
        "Enz_E":      lf("Enzyme_E_Taxon"),
        "HSP_E":      lf("HSP_E_Taxon"),
        "Osmo_E":     lf("Osmo_E_Taxon"),
    }


# ── Summarisation ─────────────────────────────────────────────────────────────

def _community_total(df: pd.DataFrame | None) -> pd.Series | None:
    """Sum all rows → community total time series.

    min_count=1 keeps a column NaN if every taxon is missing there, rather than
    silently reporting 0; NaN cells from absent taxa contribute nothing.
    """
    if df is None or df.empty:
        return None
    return df.sum(axis=0, min_count=1)


def _community_mean_taxa(df: pd.DataFrame | None) -> pd.Series | None:
    """Mean across taxa rows → community-mean time series.

    pandas mean(axis=0) skips NaN, so the average is taken only over taxa that
    are actually present on each day (absent taxa were set to NaN at load time).
    """
    if df is None or df.empty:
        return None
    return df.mean(axis=0)


def _scalar(s: pd.Series | None, summary: str) -> float | None:
    """Collapse a time series to a single scalar according to summary mode."""
    if s is None or s.empty:
        return None
    # Drop NaN and any residual inf (e.g. from system-level series not routed
    # through _load_csv's cleaning step).
    s = s.replace([np.inf, -np.inf], np.nan)
    vals = s.dropna().values
    if len(vals) == 0:
        return None
    if summary == "last":
        return float(vals[-1])
    return float(np.mean(vals))      # default: mean over selected window


def _summarise_run(data: dict, summary: str) -> dict[str, float | None]:
    """Convert a run's dict of DataFrames into per-variable scalars."""
    out = {}
    # Community total pools
    out["SUB"]       = _scalar(_community_total(data["SUB"]),       summary)
    out["MON"]       = _scalar(_community_total(data["MON"]),       summary)
    out["MIC"]       = _scalar(_community_total(data["MIC"]),       summary)
    out["RESP"]      = _scalar(_community_total(data["RESP"]),      summary)
    # Mass-specific respiration: community respiration / total biomass (per day)
    _resp_ts = _community_total(data["RESP"])
    _bio_ts  = _community_total(data["MIC"])
    if _resp_ts is not None and _bio_ts is not None:
        cc = sorted(set(_resp_ts.index) & set(_bio_ts.index))
        _rpb = (_resp_ts[cc] / _bio_ts[cc].replace(0, np.nan)) if cc else None
        out["RESP_pb"] = _scalar(_rpb, summary) if _rpb is not None else None
    else:
        out["RESP_pb"] = None
    # System efficiency (already 1-row)
    out["CUE_sys"]   = _scalar(_community_total(data["CUE_sys"]),   summary)
    out["CAE_sys"]   = _scalar(_community_total(data["CAE_sys"]),   summary)
    # Community-mean taxon efficiency
    out["CUE_tax"]   = _scalar(_community_mean_taxa(data["CUE_tax"]), summary)
    out["CAE_tax"]   = _scalar(_community_mean_taxa(data["CAE_tax"]), summary)
    # C pathways — community mean across taxa
    for key in ("C_uptake", "BC_uptake", "Enz_C", "HSP_C", "Osmo_C"):
        out[key]     = _scalar(_community_mean_taxa(data[key]),     summary)
    # E pathways
    for key in ("E_prod", "E_gain", "E_cost",
                "Trans_E", "Maint_E", "Enz_E", "HSP_E", "Osmo_E"):
        out[key]     = _scalar(_community_mean_taxa(data[key]),     summary)
    return out


# ── Gradient aggregation ──────────────────────────────────────────────────────

def build_gradient_table(gradient_dir: Path, years_arg, summary: str,
                         min_cols: int) -> pd.DataFrame:
    """
    Build a DataFrame with one row per temperature point.
    Columns: temperature, then for each variable both mean and std across
    community replicates.

    Returns a DataFrame indexed by temperature with columns:
        <var>_mean, <var>_std
    """
    treatments = _parse_treatment_dirs(gradient_dir)
    if not treatments:
        raise ValueError(
            f"No static_T*_P* subdirectories found in {gradient_dir}.\n"
            "Run the temperature gradient simulations first."
        )

    rows = []
    for temp, scenario_dir in treatments:
        run_dirs = _pick_runs(scenario_dir, min_cols)
        if not run_dirs:
            print(f"  WARNING: no complete runs found in {scenario_dir.name} — skipping")
            continue

        # Determine column range from the first run's runtime
        rt = _read_runtime(run_dirs[0])
        n_years, days_per_year = rt["pulse"], rt["end_time"]

        if summary == "last_pulse":
            col_lo, col_hi = _last_year_cols(
                pd.DataFrame(columns=range(1, n_years * days_per_year + 1)),
                n_years, days_per_year
            )
            # last_year_cols returns a df; we just need the range
            col_lo = (n_years - 1) * days_per_year + 1
            col_hi = n_years * days_per_year
            _summary_mode = "mean"
        else:
            col_lo, col_hi = _col_range(years_arg, n_years, days_per_year)
            _summary_mode = "last" if summary == "last" else "mean"

        rep_scalars = []
        for run_dir in run_dirs:
            data = _extract_run(run_dir, col_lo, col_hi)
            s = _summarise_run(data, _summary_mode)
            rep_scalars.append(s)

        # Aggregate across replicates
        all_keys = rep_scalars[0].keys()
        row = {"temperature": temp}
        for k in all_keys:
            vals = [r[k] for r in rep_scalars if r[k] is not None]
            if vals:
                row[f"{k}_mean"] = float(np.mean(vals))
                row[f"{k}_std"]  = float(np.std(vals, ddof=1)) if len(vals) > 1 else np.nan
            else:
                row[f"{k}_mean"] = np.nan
                row[f"{k}_std"]  = np.nan

        rows.append(row)
        n_rep = len(run_dirs)
        print(f"  {scenario_dir.name}: {n_rep} replicate(s)")

    if not rows:
        raise RuntimeError("No data could be loaded from any temperature treatment.")

    df = pd.DataFrame(rows).set_index("temperature").sort_index()
    return df


# ── Summary-CSV mode ──────────────────────────────────────────────────────────
#
# When the gradient directory holds export_summary_data.py output (one
# static_T*_P* subdir per temperature, each with <stem>_mean.csv/<stem>_sd.csv)
# we read those pre-aggregated files instead of raw runs.  Each variable maps to
# an export stem; system pools come pre-summed as 1-D series, taxon quantities
# are per-taxon community means that we average across taxa.

# variable key → (export stem, kind)  where kind is "system" (1-D) or "taxon" (2-D)
_GRAD_SUMMARY_MAP = {
    "SUB":      ("system_substrate_C",        "system"),
    "MON":      ("system_monomer_C",          "system"),
    "MIC":      ("system_microbial_biomass_C","system"),
    "RESP":     ("system_resp",               "system"),
    "RESP_pb":  ("system_resp_per_biomass",   "system"),
    "CUE_sys":  ("system_CUE",                "system"),
    "CAE_sys":  ("system_CAE",                "system"),
    "CUE_tax":  ("cue_taxon",                 "taxon"),
    "CAE_tax":  ("cae_taxon",                 "taxon"),
    "C_uptake": ("c_uptake_taxon",            "taxon"),
    "BC_uptake":("c_biosynthesis_taxon",      "taxon"),
    "Enz_C":    ("c_enzyme_taxon",            "taxon"),
    "HSP_C":    ("c_hsp_taxon",               "taxon"),
    "Osmo_C":   ("c_osmolyte_taxon",          "taxon"),
    "E_prod":   ("e_production_taxon",        "taxon"),
    "E_gain":   ("e_gain_taxon",              "taxon"),
    "E_cost":   ("e_cost_taxon",              "taxon"),
    "Trans_E":  ("e_transport_taxon",         "taxon"),
    "Maint_E":  ("e_maintenance_taxon",       "taxon"),
    "Enz_E":    ("e_enzyme_taxon",            "taxon"),
    "HSP_E":    ("e_hsp_taxon",               "taxon"),
    "Osmo_E":   ("e_osmolyte_taxon",          "taxon"),
}


def _read_summary_pair(scenario_dir: Path, stem: str):
    """Read <stem>_mean.csv / <stem>_sd.csv → (mean, sd).

    Local reimplementation (no scenario_utils dependency) matching the
    export_summary_data.py layout: a 1-D series is written with a single column
    "mean"/"sd" and day index; a 2-D table has taxon index and day columns.
    """
    def _read(path, value_col):
        if not path.exists():
            return None
        df = pd.read_csv(path, index_col=0)
        if list(df.columns) == [value_col]:
            s = df[value_col]
            s.index = [int(i) if str(i).lstrip("-").isdigit() else i for i in s.index]
            return s
        df.columns = [int(c) if str(c).lstrip("-").isdigit() else c
                      for c in df.columns]
        return df

    mean = _read(scenario_dir / f"{stem}_mean.csv", "mean")
    sd   = _read(scenario_dir / f"{stem}_sd.csv",   "sd")
    return mean, sd


def _read_summary_meta_local(scenario_dir: Path) -> dict:
    """Read summary_meta.csv (gridsize/pulse/end_time/col range) → dict."""
    out = {"pulse": 1, "end_time": 365, "col_start": None, "col_end": None}
    f = scenario_dir / "summary_meta.csv"
    if not f.exists():
        return out
    try:
        s = pd.read_csv(f, index_col=0).iloc[:, 0]
        for k in out:
            if k in s.index and pd.notna(s.loc[k]):
                out[k] = int(float(s.loc[k]))
    except Exception:
        pass
    return out


def _slice_days(s, lo, hi):
    """Slice a day-indexed Series to [lo, hi]; pass through if no window set."""
    if s is None or lo is None:
        return s
    try:
        return s[(s.index >= lo) & (s.index <= hi)]
    except Exception:
        return s


def build_gradient_table_from_summary(summary_dir: Path,
                                      summary_mode: str,
                                      years_arg=None) -> pd.DataFrame:
    """Build the gradient response table from export_summary_data.py output.

    One row per temperature.  Each variable's time series (system series used
    directly; taxon tables averaged across taxa) is sliced to the chosen day
    window, then collapsed to a scalar:
      summary_mode='mean'        → mean over the window (default; full run unless
                                    --years is given)
      summary_mode='last'        → last day in the window
      summary_mode='last_pulse'  → mean over the final simulation year
      --years Y [Y2]             → restrict the window to year Y (or Y..Y2)
    Pulse/days-per-year come from each scenario's summary_meta.csv.
    The ribbon (<var>_std) is the same collapse of the *_sd.csv.
    """
    treatments = _parse_treatment_dirs(summary_dir)
    if not treatments:
        raise ValueError(
            f"No static_T*_P* summary subdirectories found in {summary_dir}.\n"
            "Run export_summary_data.py on the gradient simulations first."
        )

    rows = []
    for temp, scenario_dir in treatments:
        meta = _read_summary_meta_local(scenario_dir)
        n_years, dpy = meta["pulse"], meta["end_time"]

        # Determine day window and collapse mode
        if summary_mode == "last_pulse":
            lo, hi = (n_years - 1) * dpy + 1, n_years * dpy
            collapse = "mean"
        elif years_arg:
            if len(years_arg) == 1:
                yr = years_arg[0]
                lo, hi = (yr - 1) * dpy + 1, yr * dpy
            else:
                lo, hi = (years_arg[0] - 1) * dpy + 1, years_arg[-1] * dpy
            collapse = "last" if summary_mode == "last" else "mean"
        else:
            lo = hi = None
            collapse = "last" if summary_mode == "last" else "mean"

        row = {"temperature": temp}
        for key, (stem, kind) in _GRAD_SUMMARY_MAP.items():
            mean, sd = _read_summary_pair(scenario_dir, stem)
            if kind == "taxon":
                mean = _community_mean_taxa(mean) if mean is not None else None
                sd   = _community_mean_taxa(sd)   if sd   is not None else None
            if lo is not None:
                mean = _slice_days(mean, lo, hi)
                sd   = _slice_days(sd, lo, hi)
            row[f"{key}_mean"] = _scalar(mean, collapse) if mean is not None else np.nan
            sd_val = _scalar(sd, collapse) if sd is not None else None
            row[f"{key}_std"] = sd_val if sd_val is not None else np.nan
        rows.append(row)
        _win = f"days {lo}-{hi}" if lo is not None else "full run"
        print(f"  {scenario_dir.name}: summary loaded ({_win}, {collapse})")

    if not rows:
        raise RuntimeError("No summary data could be loaded from any temperature.")

    return pd.DataFrame(rows).set_index("temperature").sort_index()


# ── Plotting helpers ──────────────────────────────────────────────────────────

def _ribbon(ax, temps, mean_col, std_col, color, label, gt):
    """Plot a mean line with ±1 SD ribbon."""
    mu  = gt[mean_col]
    sd  = gt[std_col]
    ax.plot(temps, mu, color=color, lw=2, marker="o", ms=5, label=label)
    ax.fill_between(temps,
                    (mu - sd).values, (mu + sd).values,
                    color=color, alpha=0.20, linewidth=0)


def _finish(ax, xlabel="Temperature (°C)", ylabel="", title=""):
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(loc="best", framealpha=0.85, fontsize=8)


# ── Figure generators ─────────────────────────────────────────────────────────

def fig1_c_pools(gt: pd.DataFrame, temps, fig_dir: Path):
    """C pool response curves — substrate, monomer, microbial biomass."""
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    fig.suptitle("C Pool Response to Temperature\n(ensemble mean ± 1 SD)",
                 fontsize=12, fontweight="bold")

    specs = [
        (axes[0], "SUB", C_SUB,  "Substrate C",       "Total C (model units)"),
        (axes[1], "MON", C_MON,  "Monomer C",         "Total C (model units)"),
        (axes[2], "MIC", C_MIC,  "Microbial Biomass C","Total C (model units)"),
    ]
    for ax, key, color, label, ylabel in specs:
        if f"{key}_mean" not in gt or gt[f"{key}_mean"].isna().all():
            ax.text(0.5, 0.5, "No data", transform=ax.transAxes,
                    ha="center", va="center")
            ax.set_title(label); continue
        _ribbon(ax, temps, f"{key}_mean", f"{key}_std", color, label, gt)
        _finish(ax, ylabel=ylabel, title=label)

    fig.tight_layout()
    out = fig_dir / "fig1_c_pool_response.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {out.name}")


def fig2_respiration(gt: pd.DataFrame, temps, fig_dir: Path):
    """Respiration response curves: total and mass-specific (per unit biomass)."""
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    fig.suptitle("Respiration Response to Temperature\n(ensemble mean ± 1 SD)",
                 fontsize=12, fontweight="bold")

    # Left: total/community respiration
    ax = axes[0]
    if "RESP_mean" in gt and not gt["RESP_mean"].isna().all():
        _ribbon(ax, temps, "RESP_mean", "RESP_std", C_RESP, "Respiration", gt)
        _finish(ax, ylabel="Respiration (model units)", title="Community Respiration")
    else:
        ax.text(0.5, 0.5, "No data", transform=ax.transAxes, ha="center", va="center")
        ax.set_title("Community Respiration")

    # Right: mass-specific respiration (respiration per unit biomass C)
    ax = axes[1]
    if "RESP_pb_mean" in gt and not gt["RESP_pb_mean"].isna().all():
        _ribbon(ax, temps, "RESP_pb_mean", "RESP_pb_std", C_RESP,
                "Respiration / biomass", gt)
        _finish(ax, ylabel="Respiration per unit biomass C (d⁻¹)",
                title="Mass-Specific Respiration")
    else:
        ax.text(0.5, 0.5, "No data", transform=ax.transAxes, ha="center", va="center")
        ax.set_title("Mass-Specific Respiration")

    fig.tight_layout()
    out = fig_dir / "fig2_respiration_response.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {out.name}")


def fig3_system_efficiency(gt: pd.DataFrame, temps, fig_dir: Path):
    """System-level CUE and CAE response curves."""
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=True)
    fig.suptitle("System-Level Efficiency Response to Temperature\n(ensemble mean ± 1 SD)",
                 fontsize=12, fontweight="bold")

    specs = [
        (axes[0], "CUE_sys", C_CUE, "System CUE"),
        (axes[1], "CAE_sys", C_CAE, "System CAE"),
    ]
    for ax, key, color, label in specs:
        if f"{key}_mean" not in gt or gt[f"{key}_mean"].isna().all():
            ax.text(0.5, 0.5, "No data", transform=ax.transAxes,
                    ha="center", va="center")
            ax.set_title(label); continue
        _ribbon(ax, temps, f"{key}_mean", f"{key}_std", color, label, gt)
        ax.set_ylim(0, 1)
        ax.axhline(0.5, ls="--", lw=0.7, color="gray", alpha=0.5)
        _finish(ax, ylabel="Efficiency (0–1)", title=label)

    fig.tight_layout()
    out = fig_dir / "fig3_system_efficiency_response.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {out.name}")


def fig4_taxon_efficiency(gt: pd.DataFrame, temps, fig_dir: Path):
    """Community-mean taxon CUE and CAE response curves."""
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=True)
    fig.suptitle("Community-Mean Taxon Efficiency Response to Temperature\n(ensemble mean ± 1 SD)",
                 fontsize=12, fontweight="bold")

    specs = [
        (axes[0], "CUE_tax", C_CUE, "Taxon CUE (community mean)"),
        (axes[1], "CAE_tax", C_CAE, "Taxon CAE (community mean)"),
    ]
    for ax, key, color, label in specs:
        if f"{key}_mean" not in gt or gt[f"{key}_mean"].isna().all():
            ax.text(0.5, 0.5, "No data", transform=ax.transAxes,
                    ha="center", va="center")
            ax.set_title(label); continue
        _ribbon(ax, temps, f"{key}_mean", f"{key}_std", color, label, gt)
        ax.set_ylim(0, 1)
        ax.axhline(0.5, ls="--", lw=0.7, color="gray", alpha=0.5)
        _finish(ax, ylabel="Efficiency (0–1)", title=label)

    fig.tight_layout()
    out = fig_dir / "fig4_taxon_efficiency_response.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {out.name}")


def fig5_c_pathways(gt: pd.DataFrame, temps, fig_dir: Path):
    """C pathway response curves."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    fig.suptitle("C Pathway Response to Temperature\n"
                 "(community mean per taxon; ensemble mean ± 1 SD)",
                 fontsize=12, fontweight="bold")

    # Left: uptake / biosynthesis
    ax = axes[0]
    for key, color, label in [
        ("C_uptake",  C_UPTK,  "C Uptake"),
        ("BC_uptake", C_BIOSYN,"Biosynthesis C"),
    ]:
        if f"{key}_mean" in gt and not gt[f"{key}_mean"].isna().all():
            _ribbon(ax, temps, f"{key}_mean", f"{key}_std", color, label, gt)
    _finish(ax, ylabel="C (model units per taxon)", title="C Uptake and Biosynthesis")

    # Right: allocation costs
    ax = axes[1]
    for key, color, label in [
        ("Enz_C",  C_ENZ,  "Enzyme C"),
        ("HSP_C",  C_HSP,  "HSP C"),
        ("Osmo_C", C_OSMO, "Osmolyte C"),
    ]:
        if f"{key}_mean" in gt and not gt[f"{key}_mean"].isna().all():
            _ribbon(ax, temps, f"{key}_mean", f"{key}_std", color, label, gt)
    _finish(ax, ylabel="C (model units per taxon)", title="C Allocation Costs")

    fig.tight_layout()
    out = fig_dir / "fig5_c_pathway_response.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {out.name}")


def fig6_e_pathways(gt: pd.DataFrame, temps, fig_dir: Path):
    """Energy production and balance response curves."""
    fig, ax = plt.subplots(figsize=(7, 5))
    fig.suptitle("Energy Balance Response to Temperature\n"
                 "(community mean per taxon; ensemble mean ± 1 SD)",
                 fontsize=12, fontweight="bold")

    for key, color, label in [
        ("E_prod", C_EPROD, "E Production"),
        ("E_gain", C_EGAIN, "E Net Gain"),
        ("E_cost", C_ECOST, "E Total Cost"),
    ]:
        if f"{key}_mean" in gt and not gt[f"{key}_mean"].isna().all():
            _ribbon(ax, temps, f"{key}_mean", f"{key}_std", color, label, gt)
    ax.axhline(0, color="black", lw=0.8, ls="--", alpha=0.5)
    _finish(ax, ylabel="Energy (model units per taxon)", title="Energy Balance")

    fig.tight_layout()
    out = fig_dir / "fig6_e_balance_response.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {out.name}")


def fig7_e_cost_breakdown(gt: pd.DataFrame, temps, fig_dir: Path):
    """Energy cost pathway breakdown response curves."""
    fig, ax = plt.subplots(figsize=(8, 5))
    fig.suptitle("Energy Cost Breakdown Response to Temperature\n"
                 "(community mean per taxon; ensemble mean ± 1 SD)",
                 fontsize=12, fontweight="bold")

    for key, color, label in [
        ("Trans_E", C_TRANS, "Transport E"),
        ("Maint_E", C_MAINT, "Maintenance E"),
        ("Enz_E",   C_ENZ,   "Enzyme E"),
        ("HSP_E",   C_HSP,   "HSP E"),
        ("Osmo_E",  C_OSMO,  "Osmolyte E"),
    ]:
        if f"{key}_mean" in gt and not gt[f"{key}_mean"].isna().all():
            _ribbon(ax, temps, f"{key}_mean", f"{key}_std", color, label, gt)
    _finish(ax, ylabel="Energy (model units per taxon)", title="E Cost Breakdown")

    fig.tight_layout()
    out = fig_dir / "fig7_e_cost_breakdown_response.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {out.name}")


# ── CLI ───────────────────────────────────────────────────────────────────────

def _parse_args():
    p = argparse.ArgumentParser(
        description="Temperature-gradient response curve figures for DeCA-DEMENT.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    p.add_argument("--gradient-dir", default=None,
                   help="Directory containing static_T*_P* scenario subdirectories "
                        "of raw runs.")
    p.add_argument("--summary-dir", default=None, metavar="DIR",
                   help="Directory containing static_T*_P* subdirectories of "
                        "export_summary_data.py output (read pre-aggregated "
                        "summary CSVs instead of raw runs).")
    p.add_argument("--fig-dir", default=None,
                   help="Directory to save figures (default: <dir>/figures).")
    p.add_argument("--years", nargs="+", type=int, default=None, metavar="YEAR",
                   help="Restrict analysis to these simulation year(s). "
                        "E.g. --years 2 to use only year 2.")
    p.add_argument("--summary", default="mean",
                   choices=["mean", "last", "last_pulse"],
                   help="How to collapse the time window to a single scalar. "
                        "mean: average over the selected window (default). "
                        "last: last time point in the window. "
                        "last_pulse: average over the final simulation year.")
    p.add_argument("--min-cols", type=int, default=100,
                   help="Minimum number of day-columns for a run to be included "
                        "(default: 100).")
    return p.parse_args()


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    args = _parse_args()
    if not args.gradient_dir and not args.summary_dir:
        sys.exit("ERROR: supply --gradient-dir (raw runs) or --summary-dir (summary CSVs).")

    src_dir = Path(args.summary_dir if args.summary_dir else args.gradient_dir)
    fig_dir = Path(args.fig_dir) if args.fig_dir else src_dir / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    mode = "summary" if args.summary_dir else "raw runs"
    print(f"Source directory   : {src_dir}  ({mode})")
    print(f"Figure directory   : {fig_dir}")
    print(f"Summary collapse   : {args.summary}")
    if args.years and not args.summary_dir:
        print(f"Years              : {args.years}")
    print()

    print("Loading gradient data...")
    if args.summary_dir:
        gt = build_gradient_table_from_summary(src_dir, args.summary, args.years)
    else:
        gt = build_gradient_table(src_dir, args.years, args.summary, args.min_cols)
    temps = gt.index.values
    print(f"\nTemperature points : {list(temps)}")
    print(f"Replicates per point: see above\n")

    print("Generating figures...")
    fig1_c_pools(gt, temps, fig_dir)
    fig2_respiration(gt, temps, fig_dir)
    fig3_system_efficiency(gt, temps, fig_dir)
    fig4_taxon_efficiency(gt, temps, fig_dir)
    fig5_c_pathways(gt, temps, fig_dir)
    fig6_e_pathways(gt, temps, fig_dir)
    fig7_e_cost_breakdown(gt, temps, fig_dir)

    print(f"\nAll figures saved to: {fig_dir}")


if __name__ == "__main__":
    main()
