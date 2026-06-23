#!/usr/bin/env python3
"""
system_summary.py
=================
Compute community-averaged system-level outputs for each climate scenario
from DeCA-DEMENT simulations:

  • Respiration          (RespSeries)
  • Substrate degradation (SubstratesSeries — sum across substrate types)
  • Total microbial biomass C  (MicrobesSeries — sum across taxa)

For each variable, the script produces:
  1. A daily time-series CSV  (mean ± SD across communities, one column per scenario)
  2. A pulse-aggregate CSV    (total or end-of-pulse value per scenario)
  3. A plain-text summary printed to stdout

OUTPUT FILES
------------
  <out-dir>/daily_respiration.csv
  <out-dir>/daily_substrate_C.csv
  <out-dir>/daily_microbial_biomass_C.csv
  <out-dir>/pulse_summary.csv
  <out-dir>/pulse_summary_sd.csv

USAGE
-----
  # Auto-discover scenarios from a base output directory
  python system_summary.py --base /data/out/experiment1

  # Explicit scenario directories with ordered labels
  python system_summary.py \\
      --scenario-dirs /data/out/ambient /data/out/mod_warm /data/out/ext_warm \\
      --scenario-labels Ambient "Moderate warming" "Extreme warming"

  # Treatment-map from directory stems (most robust for machine-generated names)
  python system_summary.py \\
      --base /data/out/experiment1 \\
      --treatment-map static_T10_Pn0p5=Ambient \\
                      static_T15_Pn0p5="Moderate warming" \\
                      static_T20_Pn0p5="Extreme warming"

  # Restrict to years 2–4
  python system_summary.py --base /data/out/experiment1 --years 2 4

  # Write CSVs to a specific directory
  python system_summary.py --base /data/out/experiment1 --out-dir /data/summary
"""

from __future__ import annotations

import argparse
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

# ── Hardcoded fallback (edit if needed) ───────────────────────────────────────
BASE = Path("src/output")

# ── Import shared helpers ──────────────────────────────────────────────────────
_HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(_HERE))

try:
    from scenario_utils import (
        load_series,
        pick_runs,
        discover_scenarios,
        scenario_color_map,
        add_scenario_args,
        read_runtime,
        year_range,
        filter_cols,
    )
except ImportError as e:
    sys.exit(f"ERROR: could not import scenario_utils: {e}\n"
             f"Make sure scenario_utils.py is in the same directory as this script.")


# ── Series loaders ─────────────────────────────────────────────────────────────

def load_resp(run_dir: Path) -> pd.Series | None:
    """Load RespSeries → 1-D pandas Series (index = day)."""
    df = load_series(run_dir, "RespSeries")
    if df is None:
        return None
    # Could be 1-row DataFrame (taxon='RespSeries') or single-column DataFrame
    if isinstance(df, pd.DataFrame):
        if df.shape[0] == 1:
            s = df.iloc[0]
        elif df.shape[1] == 1:
            s = df.iloc[:, 0]
        else:
            # Multiple rows — try squeezing if all rows are identical or sum
            s = df.sum(axis=0)
    else:
        s = df
    s.index = [int(c) for c in s.index]
    return s.sort_index()


def load_substrate_total(run_dir: Path) -> pd.Series | None:
    """Load SubstratesSeries, sum across all substrate types → 1-D Series (day)."""
    df = load_series(run_dir, "SubstratesSeries")
    if df is None:
        return None
    # rows = substrate types, cols = days
    total = df.sum(axis=0)
    total.index = [int(c) for c in total.index]
    return total.sort_index()



def load_substrate_no_dead(run_dir: Path) -> pd.Series | None:
    """Load SubstratesSeries, exclude DeadMic and DeadEnz, sum -> 1-D Series (day)."""
    df = load_series(run_dir, "SubstratesSeries")
    if df is None:
        return None
    df = df.loc[~df.index.isin({"DeadMic", "DeadEnz"})]
    total = df.sum(axis=0)
    total.index = [int(c) for c in total.index]
    return total.sort_index()
def load_microbe_total(run_dir: Path) -> pd.Series | None:
    """Load MicrobesSeries, sum across all taxa → 1-D Series (day)."""
    df = load_series(run_dir, "MicrobesSeries")
    if df is None:
        return None
    # rows = taxa, cols = days
    total = df.sum(axis=0)
    total.index = [int(c) for c in total.index]
    return total.sort_index()


def load_scalar_series(run_dir: Path, name: str) -> pd.Series | None:
    """Load a 1-D scalar time series (CUE_system, CAE_system, etc.)."""
    df = load_series(run_dir, name)
    if df is None:
        return None
    if isinstance(df, pd.DataFrame):
        if df.shape[0] == 1:
            s = df.iloc[0]
        elif df.shape[1] == 1:
            s = df.iloc[:, 0]
        else:
            s = df.mean(axis=0)
    else:
        s = df
    s.index = [int(c) for c in s.index]
    return s.sort_index()


# ── Per-replicate loader ───────────────────────────────────────────────────────

def load_run(run_dir: Path) -> dict | None:
    """Load all series for one replicate run. Returns None on failure."""
    resp    = load_resp(run_dir)
    sub     = load_substrate_total(run_dir)
    sub_nd  = load_substrate_no_dead(run_dir)
    mic     = load_microbe_total(run_dir)
    cue     = load_scalar_series(run_dir, "CUE_system")
    cae     = load_scalar_series(run_dir, "CAE_system")
    if resp is None and sub is None and mic is None:
        print(f"    WARNING: no data found in {run_dir.name}")
        return None
    return {"resp": resp, "sub": sub, "sub_nd": sub_nd, "mic": mic, "cue": cue, "cae": cae}


# ── Aggregation across replicates ─────────────────────────────────────────────

def _stack_series(series_list: list[pd.Series | None]) -> np.ndarray | None:
    """Stack a list of 1-D Series into a (days, n_reps) array."""
    valid = [s for s in series_list if s is not None]
    if not valid:
        return None
    idx = valid[0].index
    arr = np.stack([s.reindex(idx).values for s in valid], axis=1)
    return arr


def aggregate_scenario(run_dicts: list[dict | None],
                       key: str) -> tuple[pd.Series, pd.Series]:
    """
    Returns (mean_series, std_series) across replicates for the given key.
    std is NaN when there is only one replicate.
    """
    series_list = [d[key] for d in run_dicts if d and d.get(key) is not None]
    if not series_list:
        empty = pd.Series(dtype=float)
        return empty, empty

    idx  = series_list[0].index
    arr  = np.stack([s.reindex(idx).fillna(np.nan).values for s in series_list], axis=1)
    mean = pd.Series(np.nanmean(arr, axis=1), index=idx)
    std  = pd.Series(np.nanstd(arr,  axis=1, ddof=1) if arr.shape[1] > 1
                     else np.full(len(idx), np.nan), index=idx)
    return mean, std


# ── Pulse-level aggregation ───────────────────────────────────────────────────

def pulse_aggregate(s: pd.Series, n_pulses: int, days_per_pulse: int,
                    mode: str = "sum") -> pd.Series:
    """
    Reduce a daily Series to one value per pulse.

    mode="sum"  → cumulative total respiration per pulse
    mode="end"  → value on the last day of each pulse (pool size)
    """
    values = []
    for p in range(1, n_pulses + 1):
        d_start = (p - 1) * days_per_pulse + 1
        d_end   = p * days_per_pulse
        window  = s.loc[d_start:d_end] if not s.empty else pd.Series(dtype=float)
        if window.empty:
            values.append(np.nan)
        elif mode == "sum":
            values.append(window.sum())
        else:
            values.append(window.iloc[-1] if len(window) else np.nan)
    return pd.Series(values, index=[f"Pulse {p}" for p in range(1, n_pulses + 1)])


# ── Daily CSV builder ─────────────────────────────────────────────────────────

def build_daily_df(scenario_means: dict[str, pd.Series],
                   scenario_sds:   dict[str, pd.Series],
                   labels: list[str]) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Assemble a wide DataFrame: columns = scenario labels, index = day.
    Column order is explicitly enforced from *labels*.
    Returns (mean_df, sd_df).
    """
    mean_df = pd.DataFrame(scenario_means).reindex(columns=labels)
    sd_df   = pd.DataFrame(scenario_sds).reindex(columns=labels)
    return mean_df, sd_df


# ── Scenario-loading utility ───────────────────────────────────────────────────

def _parse_treatment_map(entries: list[str] | None) -> dict[str, str]:
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


def resolve_scenario_list(args) -> list[tuple[str, list[Path]]]:
    """
    Build an ordered [(label, [run_dirs])] list respecting:
    1. --treatment-map (stem → label, order-independent)
    2. --scenario-labels (positional)
    3. Directory stem (fallback)

    Scenario order is:
    • When --scenario-dirs: the given order
    • Otherwise: alphabetical discovery order (then --scenario-labels re-labels)
    """
    tmap = _parse_treatment_map(getattr(args, "treatment_map", None))

    if args.scenario_dirs:
        paths = [Path(p) for p in args.scenario_dirs]
        raw_labels = (args.scenario_labels
                      if args.scenario_labels and len(args.scenario_labels) == len(paths)
                      else [p.stem for p in paths])
        result = []
        for raw_label, sdir in zip(raw_labels, paths):
            label = tmap.get(sdir.stem, raw_label)
            runs  = pick_runs(sdir, min_cols=100, sentinel="RespSeries")
            if not runs:
                print(f"  WARNING: no valid runs found in {sdir}")
            result.append((label, runs))
        return result

    base = Path(args.base) if args.base else BASE
    found = discover_scenarios(base)
    if not found:
        sys.exit(f"ERROR: no scenario subdirectories found under {base}")

    found_map = {stem: sdir for stem, sdir in found}

    # If treatment-map was given, emit scenarios in treatment-map key order
    if tmap:
        result = []
        for stem, label in tmap.items():
            if stem not in found_map:
                print(f"  WARNING: --treatment-map stem {stem!r} not found under {base}")
                continue
            runs = pick_runs(found_map[stem], min_cols=100, sentinel="RespSeries")
            result.append((label, runs))
        return result

    # Otherwise use alphabetical discovery order with positional labels
    pos_labels = args.scenario_labels or []
    result = []
    for i, (stem, sdir) in enumerate(found):
        label = pos_labels[i] if i < len(pos_labels) else stem
        runs = pick_runs(sdir, min_cols=100, sentinel="RespSeries")
        result.append((label, runs))
    return result


# ── Pretty printer ─────────────────────────────────────────────────────────────

def print_pulse_table(pulse_means: pd.DataFrame, pulse_sds: pd.DataFrame,
                      variable: str, unit: str) -> None:
    """Print a formatted pulse-level table to stdout."""
    print(f"\n{'─'*70}")
    print(f"  {variable}  [{unit}]")
    print(f"{'─'*70}")
    header = f"  {'Pulse':<12}" + "".join(f"{col:>22}" for col in pulse_means.columns)
    print(header)
    print(f"{'─'*70}")
    for idx in pulse_means.index:
        row = f"  {idx:<12}"
        for col in pulse_means.columns:
            m = pulse_means.loc[idx, col]
            s = pulse_sds.loc[idx, col] if col in pulse_sds.columns else np.nan
            if np.isnan(m):
                row += f"{'—':>22}"
            elif np.isnan(s):
                row += f"{m:>22.4g}"
            else:
                row += f"{f'{m:.4g} ± {s:.3g}':>22}"
        print(row)
    print(f"{'─'*70}")


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    add_scenario_args(parser)

    parser.add_argument(
        "--out-dir", default=None, metavar="DIR",
        help="Directory to write CSV summary files (default: <base>/system_summary).",
    )
    parser.add_argument(
        "--no-csv", action="store_true",
        help="Print to stdout only; do not write CSV files.",
    )

    args = parser.parse_args()

    # ── Resolve scenarios ──────────────────────────────────────────────────────
    scenario_list = resolve_scenario_list(args)
    if not scenario_list:
        sys.exit("ERROR: no scenarios resolved. Check --base / --scenario-dirs.")

    labels = [label for label, _ in scenario_list]
    print(f"\nScenarios ({len(labels)}): {labels}")

    # ── Output directory ───────────────────────────────────────────────────────
    if not args.no_csv:
        base = Path(args.base) if args.base else BASE
        out_dir = Path(args.out_dir) if args.out_dir else base / "system_summary"
        out_dir.mkdir(parents=True, exist_ok=True)
        print(f"Output directory: {out_dir}")

    # ── Detect runtime from first run dir found ────────────────────────────────
    first_run = None
    for _, runs in scenario_list:
        if runs:
            first_run = runs[0]
            break
    rt = read_runtime(first_run)
    n_pulses      = rt["pulse"]
    days_per_pulse = rt["end_time"]
    total_days    = n_pulses * days_per_pulse
    print(f"Runtime: {n_pulses} pulses × {days_per_pulse} days = {total_days} total days")

    # Year range filter
    col_start, col_end = year_range(args.years, n_pulses, days_per_pulse)
    if args.years:
        print(f"Filtering to days {col_start}–{col_end}")

    # ── Load data ──────────────────────────────────────────────────────────────
    print("\nLoading data...")

    resp_means,   resp_sds   = {}, {}
    sub_means,    sub_sds    = {}, {}
    sub_nd_means, sub_nd_sds = {}, {}
    mic_means,    mic_sds    = {}, {}
    cue_means,    cue_sds    = {}, {}
    cae_means,    cae_sds    = {}, {}

    for label, run_dirs in scenario_list:
        n = len(run_dirs)
        print(f"  {label}: {n} run(s)")

        run_data = []
        for rd in run_dirs:
            d = load_run(rd)
            run_data.append(d)

        loaded = sum(1 for d in run_data if d is not None)
        print(f"    loaded {loaded}/{n}")

        # Aggregate
        r_mean,  r_sd  = aggregate_scenario(run_data, "resp")
        s_mean,  s_sd  = aggregate_scenario(run_data, "sub")
        sn_mean, sn_sd = aggregate_scenario(run_data, "sub_nd")
        m_mean,  m_sd  = aggregate_scenario(run_data, "mic")
        c_mean,  c_sd  = aggregate_scenario(run_data, "cue")
        a_mean,  a_sd  = aggregate_scenario(run_data, "cae")

        def _clip(s):
            return s.loc[col_start:col_end] if (s is not None and not s.empty) else s

        resp_means[label]   = _clip(r_mean);  resp_sds[label]   = _clip(r_sd)
        sub_means[label]    = _clip(s_mean);  sub_sds[label]    = _clip(s_sd)
        sub_nd_means[label] = _clip(sn_mean); sub_nd_sds[label] = _clip(sn_sd)
        mic_means[label]    = _clip(m_mean);  mic_sds[label]    = _clip(m_sd)
        cue_means[label]    = _clip(c_mean);  cue_sds[label]    = _clip(c_sd)
        cae_means[label]    = _clip(a_mean);  cae_sds[label]    = _clip(a_sd)

    # ── Build daily DataFrames ─────────────────────────────────────────────────
    resp_mean_df,   resp_sd_df   = build_daily_df(resp_means,   resp_sds,   labels)
    sub_mean_df,    sub_sd_df    = build_daily_df(sub_means,    sub_sds,    labels)
    sub_nd_mean_df, sub_nd_sd_df = build_daily_df(sub_nd_means, sub_nd_sds, labels)
    mic_mean_df,    mic_sd_df    = build_daily_df(mic_means,    mic_sds,    labels)
    cue_mean_df,    cue_sd_df    = build_daily_df(cue_means,    cue_sds,    labels)
    cae_mean_df,    cae_sd_df    = build_daily_df(cae_means,    cae_sds,    labels)

    # ── Pulse-level aggregation ────────────────────────────────────────────────
    pulse_index = [f"Pulse {p}" for p in range(1, n_pulses + 1)]

    def pulse_agg_df(mean_dict: dict, sd_dict: dict, mode: str) -> tuple[pd.DataFrame, pd.DataFrame]:
        pm, ps = {}, {}
        for label in labels:
            s = mean_dict.get(label, pd.Series(dtype=float))
            sd_s = sd_dict.get(label, pd.Series(dtype=float))
            pm[label] = pulse_aggregate(s,    n_pulses, days_per_pulse, mode)
            ps[label] = pulse_aggregate(sd_s, n_pulses, days_per_pulse, mode)
        return pd.DataFrame(pm, index=pulse_index), pd.DataFrame(ps, index=pulse_index)

    resp_pulse_mean,   resp_pulse_sd   = pulse_agg_df(resp_means,   resp_sds,   mode="sum")
    sub_pulse_mean,    sub_pulse_sd    = pulse_agg_df(sub_means,    sub_sds,    mode="end")
    sub_nd_pulse_mean, sub_nd_pulse_sd = pulse_agg_df(sub_nd_means, sub_nd_sds, mode="end")
    mic_pulse_mean,    mic_pulse_sd    = pulse_agg_df(mic_means,    mic_sds,    mode="end")
    cue_pulse_mean,    cue_pulse_sd    = pulse_agg_df(cue_means,    cue_sds,    mode="end")
    cae_pulse_mean,    cae_pulse_sd    = pulse_agg_df(cae_means,    cae_sds,    mode="end")

    # ── Stdout tables ──────────────────────────────────────────────────────────
    print_pulse_table(resp_pulse_mean,   resp_pulse_sd,
                      "Cumulative Respiration per Pulse", "C units")
    print_pulse_table(sub_pulse_mean,    sub_pulse_sd,
                      "Total Substrate C — all (end of pulse)", "C units")
    print_pulse_table(sub_nd_pulse_mean, sub_nd_pulse_sd,
                      "Total Substrate C — excl. DeadMic & DeadEnz (end of pulse)", "C units")
    print_pulse_table(mic_pulse_mean,    mic_pulse_sd,
                      "Total Microbial Biomass C (end of pulse)", "C units")
    print_pulse_table(cue_pulse_mean,    cue_pulse_sd,
                      "System CUE (end of pulse)", "unitless")
    print_pulse_table(cae_pulse_mean,    cae_pulse_sd,
                      "System CAE (end of pulse)", "unitless")

    # ── Write CSVs ─────────────────────────────────────────────────────────────
    if not args.no_csv:
        # Daily mean time series
        resp_mean_df.index.name = "day"
        resp_mean_df.to_csv(out_dir / "daily_respiration_mean.csv")
        resp_sd_df.index.name = "day"
        resp_sd_df.to_csv(out_dir / "daily_respiration_sd.csv")

        sub_mean_df.index.name = "day"
        sub_mean_df.to_csv(out_dir / "daily_substrate_C_mean.csv")
        sub_sd_df.index.name = "day"
        sub_sd_df.to_csv(out_dir / "daily_substrate_C_sd.csv")

        sub_nd_mean_df.index.name = "day"
        sub_nd_mean_df.to_csv(out_dir / "daily_substrate_C_no_dead_mean.csv")
        sub_nd_sd_df.index.name = "day"
        sub_nd_sd_df.to_csv(out_dir / "daily_substrate_C_no_dead_sd.csv")

        mic_mean_df.index.name = "day"
        mic_mean_df.to_csv(out_dir / "daily_microbial_biomass_C_mean.csv")
        mic_sd_df.index.name = "day"
        mic_sd_df.to_csv(out_dir / "daily_microbial_biomass_C_sd.csv")

        cue_mean_df.index.name = "day"
        cue_mean_df.to_csv(out_dir / "daily_CUE_mean.csv")
        cue_sd_df.index.name = "day"
        cue_sd_df.to_csv(out_dir / "daily_CUE_sd.csv")

        cae_mean_df.index.name = "day"
        cae_mean_df.to_csv(out_dir / "daily_CAE_mean.csv")
        cae_sd_df.index.name = "day"
        cae_sd_df.to_csv(out_dir / "daily_CAE_sd.csv")

        # Pulse-level summary (mean)
        pulse_summary = pd.concat([
            resp_pulse_mean.add_prefix("Resp_"),
            sub_pulse_mean.add_prefix("SubC_"),
            sub_nd_pulse_mean.add_prefix("SubC_no_dead_"),
            mic_pulse_mean.add_prefix("MicC_"),
            cue_pulse_mean.add_prefix("CUE_"),
            cae_pulse_mean.add_prefix("CAE_"),
        ], axis=1)
        pulse_summary.index.name = "pulse"
        pulse_summary.to_csv(out_dir / "pulse_summary_mean.csv")

        # Pulse-level summary (SD)
        pulse_summary_sd = pd.concat([
            resp_pulse_sd.add_prefix("Resp_"),
            sub_pulse_sd.add_prefix("SubC_"),
            sub_nd_pulse_sd.add_prefix("SubC_no_dead_"),
            mic_pulse_sd.add_prefix("MicC_"),
            cue_pulse_sd.add_prefix("CUE_"),
            cae_pulse_sd.add_prefix("CAE_"),
        ], axis=1)
        pulse_summary_sd.index.name = "pulse"
        pulse_summary_sd.to_csv(out_dir / "pulse_summary_sd.csv")

        print(f"\nFiles written to {out_dir}:")
        for f in sorted(out_dir.glob("*.csv")):
            print(f"  {f.name}")

    print("\nDone.\n")


if __name__ == "__main__":
    main()
