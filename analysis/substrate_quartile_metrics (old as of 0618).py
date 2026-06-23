#!/usr/bin/env python3
"""
substrate_quartile_metrics.py

Extract per-pulse metrics at the time points where TOTAL SUBSTRATE crosses
75%, 50%, and 25% of its per-pulse maximum.

Two modes
---------
  --summary-dir DIR   Summary mode.  Reads export_summary_data.py output
                      (ensemble-mean trajectories) and produces a single point
                      estimate per scenario × pulse × quartile.  No
                      across-community error.

  --base DIR          Per-community mode (parallel).  Takes the same raw input
                      directory as export_summary_data.py (e.g. .../temp_gradient),
                      computes the crossing metrics for EACH community run, then
                      aggregates to mean / SD / n across communities.  Emits
                      '<metric>', '<metric>_sd', '<metric>_n' columns so
                      substrate_quartile_figures.py --error community can draw
                      ensemble error bands.

Reference variable
------------------
  Plant-litter substrate = total substrate EXCLUDING necromass (DeadMic, DeadEnz).
  Summary mode sums substrates_by_type over the litter rows; per-community mode
  sums SubstratesSeries over the litter rows.  The dead pools are excluded because
  they accumulate from microbial/enzyme death and would confound the plant-litter
  decomposition signal.  Per pulse:
      max       = litter-substrate max over the pulse window
      crossings = first day in the pulse where litter substrate <= q*max, q∈{.75,.5,.25}

Tracked variables
-----------------
  LEVEL / RATE variables (level + three rate definitions each):
    substrate monomer mbc enzymes        (level; substrate & monomer = plant-litter)
    respiration resp_per_biomass         (cumulative within pulse, reset per pulse)

    Per crossing day C:
      <var>_total          T(C)
      <var>_inst_rate      T(C) - T(C-1)
      <var>_cum_rate       (T(C) - T0) / (C - pulse_start)
                           T0 = 0 (cumulative) or pulse-start value (level)
      <var>_quartile_rate  (T(C) - T(prev crossing)) / (C - prev_day)
                           prev = pulse start for the 75% point

  INSTANT-only variables (value on the crossing day; taxon stems → community mean):
    cue cae
    c_uptake c_biosynthesis c_enzyme c_hsp c_osmolyte
    e_production e_gain e_cost
    e_transport e_maintenance e_hsp e_enzyme e_osmolyte

Output
------
  Tidy CSV, one row per (scenario, pulse, quartile).  Summary mode writes the
  metric columns directly; per-community mode writes mean + '_sd' + '_n'.

Usage
-----
  # ensemble error across communities (parallel), same input as the exporter:
  python analysis/substrate_quartile_metrics.py \\
      --base    /pub/jwalkup/.../temp_gradient \\
      --out-csv analysis_out/substrate_quartile_metrics.csv \\
      --jobs    0

  # quick point estimate from an existing summary directory:
  python analysis/substrate_quartile_metrics.py \\
      --summary-dir /pub/jwalkup/.../summary \\
      --out-csv     analysis_out/substrate_quartile_metrics_summary.csv

Notes
-----
- Per-community mode reuses compute_run() from export_summary_data.py so the
  per-community quantities match the summary definitions exactly; resp_per_biomass
  is derived here (community respiration / total biomass).
- Substrate and monomer are both reduced to their plant-litter components, using
  the same exclusions as export_summary_data.py: substrate drops the DeadMic/
  DeadEnz necromass rows; monomer drops the non-plant pools (NH4, PO4, FermProd,
  Mon_DeadMic, Mon_DeadEnz).  MBC, enzymes, and respiration are unchanged.
"""

import argparse
import os
import re
import sys
import warnings
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

# scenario_utils provides the summary loaders + metadata reader + run discovery
sys.path.insert(0, str(Path(__file__).resolve().parent))
try:
    from scenario_utils import (
        load_summary, discover_summary_scenarios, read_summary_meta,
        discover_scenarios, pick_runs, read_runtime,
    )
except ImportError as e:
    sys.exit(f"ERROR: could not import scenario_utils: {e}")

# compute_run (per-run derivation) is reused from the exporter so the
# per-community quantities match the summary definitions exactly.
try:
    from export_summary_data import compute_run as _compute_run
except ImportError:
    _compute_run = None


# ── Variable definitions ────────────────────────────────────────────────────────

# LEVEL / RATE variables: (output_name, stem, kind)  kind ∈ {"level", "cumulative"}
LEVEL_RATE_VARS = [
    ("substrate",        "system_substrate_C",         "level"),
    ("monomer",          "system_monomer_C",           "level"),
    ("mbc",              "system_microbial_biomass_C", "level"),
    ("enzymes",          "system_enzyme_C",            "level"),
    ("respiration",      "system_resp",                "cumulative"),
    ("resp_per_biomass", "system_resp_per_biomass",    "cumulative"),
]

# INSTANT-only variables: (output_name, stem, kind)  kind ∈ {"system", "taxon"}
INSTANT_VARS = [
    ("cue",            "system_CUE",          "system"),
    ("cae",            "system_CAE",          "system"),
    ("c_uptake",       "c_uptake_taxon",      "taxon"),
    ("c_biosynthesis", "c_biosynthesis_taxon","taxon"),
    ("c_enzyme",       "c_enzyme_taxon",      "taxon"),
    ("c_hsp",          "c_hsp_taxon",         "taxon"),
    ("c_osmolyte",     "c_osmolyte_taxon",    "taxon"),
    ("e_production",   "e_production_taxon",  "taxon"),
    ("e_gain",         "e_gain_taxon",        "taxon"),
    ("e_cost",         "e_cost_taxon",        "taxon"),
    ("e_transport",    "e_transport_taxon",   "taxon"),
    ("e_maintenance",  "e_maintenance_taxon", "taxon"),
    ("e_hsp",          "e_hsp_taxon",         "taxon"),
    ("e_enzyme",       "e_enzyme_taxon",      "taxon"),
    ("e_osmolyte",     "e_osmolyte_taxon",    "taxon"),
]

_TEMP_RE = re.compile(r"T([+-]?\d+(?:[._p]\d+)?)", re.IGNORECASE)


# ── Helpers ───────────────────────────────────────────────────────────────────

def _parse_temperature(name: str):
    """Pull a temperature out of a scenario name like 'static_T15_Pn0p5'."""
    m = _TEMP_RE.search(name)
    if not m:
        return np.nan
    try:
        return float(m.group(1).replace("p", ".").replace("_", "."))
    except ValueError:
        return np.nan


def _as_series(obj):
    """Coerce a loaded summary object to a day-indexed Series."""
    if obj is None:
        return None
    if isinstance(obj, pd.Series):
        return obj
    return obj.mean(axis=0)   # DataFrame (taxa × days) → community mean across taxa


# Necromass pools excluded from the "total substrate" reference: these accumulate
# from microbial/enzyme death and would confound plant-litter decomposition.
_DEAD_SUBSTRATES = ("DeadMic", "DeadEnz")

# Non-plant-derived monomers excluded from the plant-litter monomer total.
# Same definition as export_summary_data.py: nutrients (NH4, PO4), fermentation
# product, and necromass-derived monomers; candidate labels guard naming variants.
_NONPLANT_MONOMER_CANDS = {
    "nh4":      ("NH4",),
    "po4":      ("PO4",),
    "fermprod": ("FermProd",),
    "deadmic":  ("Mon_DeadMic", "DeadMic"),
    "deadenz":  ("Mon_DeadEnz", "DeadEnz"),
}


def _litter_total(df):
    """Sum a per-type substrate DataFrame across plant-litter types only.

    Drops the necromass rows (DeadMic, DeadEnz) before summing, so the result is
    the plant-litter substrate total used as the quartile reference.  Returns a
    day-indexed Series, or None.
    """
    if df is None:
        return None
    if isinstance(df, pd.Series):     # already a 1-D total; can't separate pools
        return df
    keep = [i for i in df.index if i not in _DEAD_SUBSTRATES]
    return df.loc[keep].sum(axis=0) if keep else None


def _litter_monomer_total(df):
    """Sum a per-type monomer DataFrame across plant-litter monomers only.

    Excludes the non-plant pools (NH4, PO4, FermProd, Mon_DeadMic, Mon_DeadEnz)
    by matching candidate labels, then sums the rest — identical in definition to
    export_summary_data.py's system_plant_litter_monomer_C.  Returns a day-indexed
    Series, or None.
    """
    if df is None:
        return None
    if isinstance(df, pd.Series):     # already a 1-D total; can't separate pools
        return df
    nonplant = set()
    for cands in _NONPLANT_MONOMER_CANDS.values():
        lbl = next((c for c in cands if c in df.index), None)
        if lbl is not None:
            nonplant.add(lbl)
    keep = [i for i in df.index if i not in nonplant]
    return df.loc[keep].sum(axis=0) if keep else None


def _pulse_windows(n_pulses: int, dpy: int):
    """Yield (pulse_index, start_global_day, end_global_day), 1-indexed days."""
    for p in range(1, n_pulses + 1):
        yield p, (p - 1) * dpy + 1, p * dpy


def _first_crossing(sub_window: pd.Series, threshold: float):
    """First day in the pulse window where substrate <= threshold (downward)."""
    below = sub_window[sub_window <= threshold]
    return int(below.index[0]) if len(below) else None


def _tracked_within_pulse(series: pd.Series, start: int, end: int, kind: str):
    """Build the tracked total T(d) over the pulse window.

    level      → the values themselves
    cumulative → running sum of the daily series within the pulse (resets to 0)
    """
    if series is None:
        return None
    win = series[(series.index >= start) & (series.index <= end)]
    if win.empty:
        return None
    return win.cumsum() if kind == "cumulative" else win


# ── Core extraction (shared by both modes) ──────────────────────────────────────

def _extract_from_series(substrate, lr_series, inst_series,
                         n_pulses, dpy, quartiles, base_row=None):
    """Per-pulse / per-quartile extraction from a set of day-indexed series.

    substrate   : Series (day-indexed) — the reference trajectory
    lr_series   : {name: (Series, kind)} for the level/rate variables
    inst_series : {name: Series} for the instant-only variables
    base_row    : optional dict of identifier fields copied into every row
    Returns a list of row dicts (one per pulse × quartile).
    """
    base_row = base_row or {}
    rows = []
    if substrate is None or substrate.empty:
        return rows

    for p, start, end in _pulse_windows(n_pulses, dpy):
        sub_win = substrate[(substrate.index >= start) & (substrate.index <= end)]
        if sub_win.empty:
            continue
        sub_max = float(sub_win.max())
        if not np.isfinite(sub_max) or sub_max <= 0:
            continue

        tracked = {name: _tracked_within_pulse(ser, start, end, kind)
                   for name, (ser, kind) in lr_series.items()}
        t0 = {}
        for name, (ser, kind) in lr_series.items():
            tv = tracked[name]
            if tv is None:
                t0[name] = np.nan
            elif kind == "cumulative":
                t0[name] = 0.0
            else:
                t0[name] = float(tv.iloc[0])

        prev_day = start
        prev_total = dict(t0)
        for q in quartiles:
            thr = q * sub_max
            C = _first_crossing(sub_win, thr)
            row = dict(base_row)
            row.update({
                "pulse": p,
                "quartile_pct": int(round(q * 100)),
                "threshold": thr,
                "substrate_max": sub_max,
                "crossing_day_global": C if C is not None else np.nan,
                "crossing_pulse_day": (C - start + 1) if C is not None else np.nan,
            })

            if C is None:
                for name, _stem, _k in LEVEL_RATE_VARS:
                    row[f"{name}_total"] = np.nan
                    row[f"{name}_inst_rate"] = np.nan
                    row[f"{name}_cum_rate"] = np.nan
                    row[f"{name}_quartile_rate"] = np.nan
                for name, _stem, _k in INSTANT_VARS:
                    row[f"{name}_value"] = np.nan
                rows.append(row)
                continue

            for name, _stem, kind in LEVEL_RATE_VARS:
                tv = tracked[name]
                if tv is None or C not in tv.index:
                    row[f"{name}_total"] = np.nan
                    row[f"{name}_inst_rate"] = np.nan
                    row[f"{name}_cum_rate"] = np.nan
                    row[f"{name}_quartile_rate"] = np.nan
                    continue

                t_now = float(tv.loc[C])
                row[f"{name}_total"] = t_now

                prev_idx = C - 1
                if prev_idx in tv.index:
                    row[f"{name}_inst_rate"] = t_now - float(tv.loc[prev_idx])
                elif kind == "cumulative":
                    row[f"{name}_inst_rate"] = t_now
                else:
                    row[f"{name}_inst_rate"] = np.nan

                denom_cum = C - start
                row[f"{name}_cum_rate"] = (
                    (t_now - t0[name]) / denom_cum if denom_cum > 0 else np.nan
                )
                denom_q = C - prev_day
                row[f"{name}_quartile_rate"] = (
                    (t_now - prev_total[name]) / denom_q if denom_q > 0 else np.nan
                )

            for name, _stem, _kind in INSTANT_VARS:
                ser = inst_series[name]
                row[f"{name}_value"] = (
                    float(ser.loc[C]) if (ser is not None and C in ser.index) else np.nan
                )

            rows.append(row)

            prev_day = C
            for name, _stem, _k in LEVEL_RATE_VARS:
                tv = tracked[name]
                prev_total[name] = (float(tv.loc[C])
                                    if (tv is not None and C in tv.index) else prev_total[name])

    return rows


# ── Summary mode ────────────────────────────────────────────────────────────────

def extract_scenario(scenario_dir: Path, label: str, quartiles) -> list[dict]:
    """Summary mode: crossings on the ensemble-mean trajectory (point estimate)."""
    meta = read_summary_meta(scenario_dir)
    n_pulses, dpy = meta["pulse"], int(meta["end_time"])

    # Reference = plant-litter substrate: sum substrates_by_type EXCLUDING the
    # necromass rows (DeadMic, DeadEnz).  Fall back to system_substrate_C only if
    # the per-type breakdown is unavailable (then dead pools can't be removed).
    sbt, _ = load_summary(scenario_dir, "substrates_by_type")
    substrate = _litter_total(sbt)
    if substrate is None or substrate.empty:
        s, _ = load_summary(scenario_dir, "system_substrate_C")
        substrate = _as_series(s)
        if substrate is not None and not substrate.empty:
            print(f"  WARNING: {label}: substrates_by_type missing; using "
                  f"system_substrate_C (includes DeadMic/DeadEnz).")
    if substrate is None or substrate.empty:
        print(f"  WARNING: no substrate series for {label}; skipping")
        return []

    # Plant-litter monomer: sum monomers_by_type EXCLUDING the non-plant pools,
    # handled exactly like substrate.  Fall back to the export's pre-computed
    # plant-litter stem, then the full monomer total, if the breakdown is absent.
    mbt, _ = load_summary(scenario_dir, "monomers_by_type")
    monomer = _litter_monomer_total(mbt)
    if monomer is None or (hasattr(monomer, "empty") and monomer.empty):
        s, _ = load_summary(scenario_dir, "system_plant_litter_monomer_C")
        monomer = _as_series(s)
    if monomer is None or (hasattr(monomer, "empty") and monomer.empty):
        s, _ = load_summary(scenario_dir, "system_monomer_C")
        monomer = _as_series(s)
        if monomer is not None and not monomer.empty:
            print(f"  WARNING: {label}: monomers_by_type and "
                  f"system_plant_litter_monomer_C missing; using full "
                  f"system_monomer_C (includes non-plant monomers).")

    lr_series = {}
    for out_name, stem, kind in LEVEL_RATE_VARS:
        if out_name == "substrate":
            lr_series["substrate"] = (substrate, kind)   # plant-litter substrate
            continue
        if out_name == "monomer":
            lr_series["monomer"] = (monomer, kind)        # plant-litter monomer
            continue
        s, _ = load_summary(scenario_dir, stem)
        lr_series[out_name] = (_as_series(s), kind)

    inst_series = {}
    for out_name, stem, _kind in INSTANT_VARS:
        s, _ = load_summary(scenario_dir, stem)
        inst_series[out_name] = _as_series(s)

    base = {"scenario": label, "temperature": _parse_temperature(label)}
    return _extract_from_series(substrate, lr_series, inst_series,
                                n_pulses, dpy, quartiles, base_row=base)


# ── Per-community mode (raw runs → mean/SD/n across communities) ────────────────

# compute_run dict key → instant-variable output name (taxon DataFrames)
_CR_INSTANT_KEYS = {
    "c_uptake": "CU", "c_biosynthesis": "BC", "c_enzyme": "EC",
    "c_hsp": "HC", "c_osmolyte": "OC",
    "e_production": "EP", "e_gain": "EGN", "e_cost": "ECT",
    "e_transport": "ET", "e_maintenance": "EM", "e_hsp": "EH",
    "e_enzyme": "EE", "e_osmolyte": "EO",
}


def _sum_rows(df):
    return df.sum(axis=0) if df is not None else None


def _mean_rows(df):
    return df.mean(axis=0) if df is not None else None


def _series_from_run(d: dict):
    """Build (substrate, lr_series, inst_series) from a compute_run() dict."""
    # Plant-litter substrate / monomer: exclude necromass (substrate) and the
    # non-plant pools (monomer), matching export_summary_data.py definitions.
    substrate = _litter_total(d.get("SUB"))
    monomer   = _litter_monomer_total(d.get("MON"))
    mbc       = _sum_rows(d.get("MIC"))
    enzymes   = _sum_rows(d.get("ENZ"))
    resp      = d.get("RESP")            # already a Series (squeezed)

    resp_pb = None
    if resp is not None and mbc is not None:
        cc = sorted(set(resp.index) & set(mbc.index))
        if cc:
            denom = mbc.loc[cc].replace(0, np.nan)
            resp_pb = resp.loc[cc] / denom

    lr_series = {
        "substrate":        (substrate, "level"),
        "monomer":          (monomer,   "level"),
        "mbc":              (mbc,        "level"),
        "enzymes":          (enzymes,    "level"),
        "respiration":      (resp,       "cumulative"),
        "resp_per_biomass": (resp_pb,    "cumulative"),
    }

    inst_series = {"cue": d.get("CUE_sys"), "cae": d.get("CAE_sys")}
    for out_name, key in _CR_INSTANT_KEYS.items():
        inst_series[out_name] = _mean_rows(d.get(key))

    return substrate, lr_series, inst_series


def _process_run(task):
    """Worker: compute one community's per-pulse/quartile rows."""
    scenario, temperature, run_dir_str, quartiles = task
    run_dir = Path(run_dir_str)
    try:
        d = _compute_run(run_dir)
        if d is None:
            return []
        substrate, lr_series, inst_series = _series_from_run(d)
        if substrate is None or substrate.empty:
            return []
        rt = read_runtime(run_dir)
        n_pulses, dpy = rt["pulse"], int(rt["end_time"])
        base = {"scenario": scenario, "temperature": temperature,
                "community": run_dir.name}
        return _extract_from_series(substrate, lr_series, inst_series,
                                    n_pulses, dpy, quartiles, base_row=base)
    except Exception as e:
        return [{"__error__": f"{run_dir.name}: {type(e).__name__}: {e}"}]


def _resolve_jobs(requested: int) -> int:
    """Resolve worker count: 0 → all available cores (SLURM-aware)."""
    cores = int(os.environ.get("SLURM_CPUS_PER_TASK", 0)) or (os.cpu_count() or 1)
    if requested and requested > 0:
        return max(1, min(requested, cores))
    return max(1, cores)


def _aggregate_across_communities(rows: list[dict]) -> pd.DataFrame:
    """Aggregate per-community rows → mean / SD / n across communities."""
    df = pd.DataFrame(rows)
    id_cols = ["scenario", "temperature", "pulse", "quartile_pct"]
    metric_cols = [c for c in df.columns if c not in id_cols + ["community"]]

    g = df.groupby(id_cols, dropna=False)
    mean = g[metric_cols].mean()
    sd   = g[metric_cols].std(ddof=1)
    n    = g[metric_cols].count()

    out = mean.copy()
    ordered = []
    for c in metric_cols:
        out[f"{c}_sd"] = sd[c]
        out[f"{c}_n"] = n[c].astype("Int64")
        ordered += [c, f"{c}_sd", f"{c}_n"]
    return out[ordered].reset_index()


def run_community_mode(base_dir: Path, quartiles, n_jobs: int) -> pd.DataFrame:
    """Compute crossings per community from raw runs, aggregate across communities."""
    if _compute_run is None:
        sys.exit("ERROR: could not import compute_run from export_summary_data.py "
                 "(needed for --base/per-community mode).")

    scenarios = discover_scenarios(base_dir)
    if not scenarios:
        sys.exit(f"No scenarios found under {base_dir}.")

    tasks = []
    per_scenario_counts = {}
    for label, sdir in scenarios:
        runs = pick_runs(sdir, min_cols=100, sentinel="SubstratesSeries")
        per_scenario_counts[label] = len(runs)
        temp = _parse_temperature(label)
        for rd in runs:
            tasks.append((label, temp, str(rd), quartiles))

    if not tasks:
        sys.exit("No community run directories found under any scenario.")

    workers = _resolve_jobs(n_jobs)
    print(f"Scenarios   : {len(scenarios)}")
    for label, c in per_scenario_counts.items():
        print(f"    {label}: {c} communities")
    print(f"Total runs  : {len(tasks)}")
    print(f"Workers     : {workers}\n")

    all_rows, errors = [], []
    done, total = 0, len(tasks)

    if workers == 1:
        for t in tasks:
            res = _process_run(t)
            done += 1
            for r in res:
                (errors if "__error__" in r else all_rows).append(r)
            if done % 25 == 0 or done == total:
                print(f"  [{done}/{total}] processed", flush=True)
    else:
        with ProcessPoolExecutor(max_workers=workers) as ex:
            futs = {ex.submit(_process_run, t): t for t in tasks}
            for fut in as_completed(futs):
                res = fut.result()
                done += 1
                for r in res:
                    (errors if "__error__" in r else all_rows).append(r)
                if done % 25 == 0 or done == total:
                    print(f"  [{done}/{total}] processed", flush=True)

    for e in errors[:20]:
        print(f"  WARNING: {e['__error__']}")
    if len(errors) > 20:
        print(f"  ... and {len(errors) - 20} more failures")

    if not all_rows:
        sys.exit("No crossing rows produced from any community.")

    return _aggregate_across_communities(all_rows)


# ── CLI / main ───────────────────────────────────────────────────────────────

def _parse_args():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    src = p.add_mutually_exclusive_group(required=True)
    src.add_argument("--base", default=None, metavar="DIR",
                     help="Raw gradient directory (same input as export_summary_data.py, "
                          "e.g. .../temp_gradient). Per-community mode with mean/SD/n.")
    src.add_argument("--summary-dir", default=None, metavar="DIR",
                     help="Root of export_summary_data.py output. Summary mode: single "
                          "point estimate per scenario (no across-community error).")
    p.add_argument("--out-csv", default=None, metavar="PATH",
                   help="Output CSV (default: <input-dir>/substrate_quartile_metrics.csv).")
    p.add_argument("--quartiles", nargs="+", type=float, default=[0.75, 0.50, 0.25],
                   metavar="Q",
                   help="Substrate thresholds as fractions of the per-pulse max "
                        "(default: 0.75 0.5 0.25).")
    p.add_argument("--jobs", "-j", type=int, default=0, metavar="N",
                   help="Worker processes for per-community mode "
                        "(0 = all cores; SLURM_CPUS_PER_TASK aware).")
    return p.parse_args()


def main():
    args = _parse_args()
    quartiles = sorted(args.quartiles, reverse=True)   # 0.75, 0.5, 0.25

    if args.base:
        base_dir = Path(args.base)
        out_csv = Path(args.out_csv) if args.out_csv else base_dir / "substrate_quartile_metrics.csv"
        print("Mode        : per-community (raw runs)")
        print(f"Base dir    : {base_dir}")
        print(f"Quartiles   : {quartiles}")
        df = run_community_mode(base_dir, quartiles, args.jobs)
    else:
        summary_dir = Path(args.summary_dir)
        out_csv = Path(args.out_csv) if args.out_csv else summary_dir / "substrate_quartile_metrics.csv"
        scenarios = discover_summary_scenarios(summary_dir)
        if not scenarios:
            sys.exit(f"No summary scenarios found under {summary_dir}.")
        print("Mode        : summary (ensemble-mean trajectory)")
        print(f"Summary dir : {summary_dir}")
        print(f"Scenarios   : {len(scenarios)}")
        print(f"Quartiles   : {quartiles}\n")
        all_rows = []
        for label, sdir in scenarios:
            print(f"  {label}")
            all_rows.extend(extract_scenario(sdir, label, quartiles))
        if not all_rows:
            sys.exit("No crossing rows produced — check the summary contents.")
        df = pd.DataFrame(all_rows)

    df = df.sort_values(["temperature", "scenario", "pulse", "quartile_pct"],
                        ascending=[True, True, True, False]).reset_index(drop=True)

    out_csv.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_csv, index=False)
    print(f"\nWrote {len(df)} rows × {df.shape[1]} columns to {out_csv}")


if __name__ == "__main__":
    main()
