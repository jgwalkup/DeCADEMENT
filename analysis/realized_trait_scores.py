#!/usr/bin/env python3
"""
realized_trait_scores.py

Extract per-taxon realized trait scores at substrate quartile crossing days
from DeCA-DEMENT simulation output. Produces a tidy CSV with one row per
taxon × community × pulse × quartile, covering all scenarios.

Trait axes
----------
  Acquisition:
    enzyme_C      – C invested in enzyme production at crossing day C
    transport_C   – C invested in transporter production at crossing day C

  Stress tolerance:
    hsp_C         – C invested in HSP production at crossing day C

Each trait is recorded in three normalizations:
    _raw          – raw C cost (model units)
    _per_biomass  – C cost / biomass C at day C
    _per_uptake   – C cost / C uptake at day C

Dead taxa (biomass = 0 at day C) are recorded with NaN trait values.

Taxon ID correction
-------------------
Simulation outputs produced before the output.py fix carry positional taxon
labels (Tax1…TaxN). This script always remaps those labels to global IDs by
reading Initialization/Microbes_pp.csv from each run directory, which carries
the correct global IDs in community order regardless of the output.py version.

Output columns
--------------
  scenario         – scenario label (e.g. static_T15_Pn0p5)
  temperature      – numeric temperature parsed from scenario name
  community        – run directory name
  taxon            – global taxon ID (e.g. Tax312)
  pulse            – pulse number (1-indexed)
  quartile_pct     – substrate threshold as integer percent (75, 50, 25)
  crossing_day     – global simulation day of the crossing
  biomass_C        – taxon biomass C at crossing day (NaN if dead)
  uptake_C         – taxon C uptake at crossing day (NaN if dead)
  alive            – True if taxon was alive (biomass > 0) at crossing day
  enzyme_C_raw          – enzyme C investment, raw
  enzyme_C_per_biomass  – enzyme C / biomass C
  enzyme_C_per_uptake   – enzyme C / uptake C
  transport_C_raw
  transport_C_per_biomass
  transport_C_per_uptake
  hsp_C_raw
  hsp_C_per_biomass
  hsp_C_per_uptake

Usage
-----
    # Single scenario directory
    python analysis/realized_trait_scores.py \\
        --base    /path/to/temp_gradient \\
        --out-csv realized_traits.csv \\
        --jobs    0

    # Specific scenarios only
    python analysis/realized_trait_scores.py \\
        --base    /path/to/temp_gradient \\
        --scenarios static_T15_Pn0p5 static_T27_Pn0p5 \\
        --out-csv realized_traits.csv

Notes
-----
- Crossing days are computed using the same plant-litter substrate definition
  as substrate_quartile_metrics.py (DeadMic and DeadEnz excluded).
- Taxon-level matrices are (taxa × days); values at the crossing day are
  extracted by column label after global ID remapping.
- All divisions use safe division (NaN where denominator is 0 or NaN).
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

_HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(_HERE))
try:
    from scenario_utils import discover_scenarios, pick_runs, read_runtime
except ImportError as e:
    sys.exit(f"ERROR: could not import scenario_utils: {e}")

try:
    from export_summary_data import compute_run as _compute_run
except ImportError:
    _compute_run = None

# Necromass pools excluded from the plant-litter substrate reference
_DEAD_SUBSTRATES = {"DeadMic", "DeadEnz"}

# Temperature parser (matches static_T15_Pn0p5 etc.)
_TEMP_RE = re.compile(r"T([+-]?\d+(?:[._p]\d+)?)", re.IGNORECASE)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _parse_temperature(name: str):
    m = _TEMP_RE.search(name)
    if not m:
        return np.nan
    try:
        return float(m.group(1).replace("p", ".").replace("_", "."))
    except ValueError:
        return np.nan


def _safe_div(num: pd.Series, denom: pd.Series) -> pd.Series:
    """Element-wise division; NaN where denominator is 0 or NaN."""
    d = denom.replace(0, np.nan)
    return num / d


def _load_taxon_csv(run_dir: Path, stem: str) -> pd.DataFrame | None:
    """Load a taxa × days CSV from the run directory (CSV or NetCDF)."""
    csv_path = run_dir / f"{stem}.csv"
    nc_path  = run_dir / f"{stem}.nc"
    if csv_path.exists():
        df = pd.read_csv(csv_path, index_col=0)
        # Columns are day numbers — coerce to int
        df.columns = [int(c) for c in df.columns]
        return df
    if nc_path.exists():
        try:
            import xarray as xr
            ds = xr.open_dataset(nc_path)
            df = ds.to_dataframe().unstack()
            ds.close()
            return df
        except Exception:
            return None
    return None


def _global_taxon_ids(run_dir: Path) -> list[str] | None:
    """
    Read the correct global taxon IDs from Initialization/Microbes_pp.csv.
    Returns an ordered list matching the positional Tax1…TaxN labels used in
    output CSVs, or None if the file cannot be read.
    """
    init_path = run_dir / "Initialization" / "Microbes_pp.csv"
    if not init_path.exists():
        return None
    try:
        df = pd.read_csv(init_path, index_col=0)
        return list(df.index)
    except Exception:
        return None


def _remap_index(df: pd.DataFrame, global_ids: list[str]) -> pd.DataFrame:
    """
    Replace positional Tax1…TaxN row labels with global IDs.
    If the number of rows doesn't match, returns df unchanged with a warning.
    """
    if df is None:
        return None
    if len(df) != len(global_ids):
        warnings.warn(
            f"Row count mismatch: DataFrame has {len(df)} rows but "
            f"{len(global_ids)} global IDs. Index not remapped."
        )
        return df
    df = df.copy()
    df.index = global_ids
    return df


def _litter_substrate(sub_df: pd.DataFrame | pd.Series | None) -> pd.Series | None:
    """Sum substrate DataFrame across plant-litter types (exclude necromass)."""
    if sub_df is None:
        return None
    if isinstance(sub_df, pd.Series):
        return sub_df
    keep = [i for i in sub_df.index if i not in _DEAD_SUBSTRATES]
    return sub_df.loc[keep].sum(axis=0) if keep else None


def _pulse_windows(n_pulses: int, dpy: int):
    for p in range(1, n_pulses + 1):
        yield p, (p - 1) * dpy + 1, p * dpy


def _first_crossing(sub_window: pd.Series, threshold: float) -> int | None:
    below = sub_window[sub_window <= threshold]
    return int(below.index[0]) if len(below) else None


# ---------------------------------------------------------------------------
# Per-run processing
# ---------------------------------------------------------------------------

def _extract_run(run_dir: Path, scenario: str, temperature: float,
                 quartiles: list[float]) -> list[dict]:
    """
    Extract per-taxon trait scores at each substrate crossing day for one run.
    Returns a list of row dicts (one per taxon × pulse × quartile).
    """
    # ── Taxon ID remapping ────────────────────────────────────────────────────
    global_ids = _global_taxon_ids(run_dir)
    if global_ids is None:
        warnings.warn(f"{run_dir.name}: Initialization/Microbes_pp.csv not found; "
                      f"positional taxon IDs will be used.")

    # ── Runtime ───────────────────────────────────────────────────────────────
    rt = read_runtime(run_dir)
    n_pulses = int(rt["pulse"])
    dpy      = int(rt["end_time"])

    # ── Load taxon-level matrices ─────────────────────────────────────────────
    def load_remap(stem):
        df = _load_taxon_csv(run_dir, stem)
        if df is not None and global_ids is not None:
            df = _remap_index(df, global_ids)
        return df

    MIC  = load_remap("MicrobesSeries")     # biomass C
    ENZ  = load_remap("Enzyme_C_Taxon")     # enzyme C investment
    TRN  = load_remap("Transport_C_Taxon")  # transporter C investment
    HSP  = load_remap("HSP_C_Taxon")        # HSP C investment
    UPT  = load_remap("C_Uptake_Taxon")     # C uptake

    if MIC is None:
        warnings.warn(f"{run_dir.name}: MicrobesSeries not found; skipping.")
        return []

    # ── Load substrate for crossing detection ─────────────────────────────────
    # Use compute_run() SUB matrix if available; fall back to SubstratesSeries
    d = None
    if _compute_run is not None:
        try:
            d = _compute_run(run_dir)
        except Exception:
            pass

    if d is not None and d.get("SUB") is not None:
        substrate = _litter_substrate(d["SUB"])
    else:
        sub_df = _load_taxon_csv(run_dir, "SubstratesSeries")
        substrate = _litter_substrate(sub_df)

    if substrate is None or substrate.empty:
        warnings.warn(f"{run_dir.name}: no substrate series; skipping.")
        return []

    # ── Per-pulse / per-quartile extraction ───────────────────────────────────
    rows = []
    taxon_ids = list(MIC.index)

    for p, start, end in _pulse_windows(n_pulses, dpy):
        sub_win = substrate[(substrate.index >= start) & (substrate.index <= end)]
        if sub_win.empty:
            continue
        sub_max = float(sub_win.max())
        if not np.isfinite(sub_max) or sub_max <= 0:
            continue

        for q in quartiles:
            thr = q * sub_max
            C   = _first_crossing(sub_win, thr)

            for taxon in taxon_ids:
                row = {
                    "scenario":    scenario,
                    "temperature": temperature,
                    "community":   run_dir.name,
                    "taxon":       taxon,
                    "pulse":       p,
                    "quartile_pct": int(round(q * 100)),
                    "crossing_day": C if C is not None else np.nan,
                }

                if C is None:
                    # No crossing — record NaN for all trait values
                    row.update({
                        "biomass_C": np.nan, "uptake_C": np.nan,
                        "alive": False,
                        "enzyme_C_raw": np.nan,
                        "enzyme_C_per_biomass": np.nan,
                        "enzyme_C_per_uptake": np.nan,
                        "transport_C_raw": np.nan,
                        "transport_C_per_biomass": np.nan,
                        "transport_C_per_uptake": np.nan,
                        "hsp_C_raw": np.nan,
                        "hsp_C_per_biomass": np.nan,
                        "hsp_C_per_uptake": np.nan,
                    })
                    rows.append(row)
                    continue

                # Extract scalar values at crossing day C
                def _val(mat, tid=taxon, day=C):
                    if mat is None or tid not in mat.index or day not in mat.columns:
                        return np.nan
                    return float(mat.loc[tid, day])

                biomass = _val(MIC)
                uptake  = _val(UPT)
                alive   = (biomass > 0) if np.isfinite(biomass) else False

                row["biomass_C"] = biomass
                row["uptake_C"]  = uptake
                row["alive"]     = alive

                if not alive:
                    # Dead taxon — NaN for all trait columns
                    row.update({
                        "enzyme_C_raw": np.nan,
                        "enzyme_C_per_biomass": np.nan,
                        "enzyme_C_per_uptake": np.nan,
                        "transport_C_raw": np.nan,
                        "transport_C_per_biomass": np.nan,
                        "transport_C_per_uptake": np.nan,
                        "hsp_C_raw": np.nan,
                        "hsp_C_per_biomass": np.nan,
                        "hsp_C_per_uptake": np.nan,
                    })
                    rows.append(row)
                    continue

                # Alive taxon — compute all three normalizations
                for prefix, mat in [("enzyme_C",    ENZ),
                                     ("transport_C", TRN),
                                     ("hsp_C",       HSP)]:
                    raw = _val(mat)
                    row[f"{prefix}_raw"] = raw

                    if np.isfinite(raw) and np.isfinite(biomass) and biomass > 0:
                        row[f"{prefix}_per_biomass"] = raw / biomass
                    else:
                        row[f"{prefix}_per_biomass"] = np.nan

                    if np.isfinite(raw) and np.isfinite(uptake) and uptake > 0:
                        row[f"{prefix}_per_uptake"] = raw / uptake
                    else:
                        row[f"{prefix}_per_uptake"] = np.nan

                rows.append(row)

    return rows


def _process_run(task: tuple) -> list[dict]:
    """Worker function — unpacks task tuple and calls _extract_run."""
    scenario, temperature, run_dir_str, quartiles = task
    try:
        return _extract_run(Path(run_dir_str), scenario, temperature, quartiles)
    except Exception as e:
        return [{"__error__": f"{Path(run_dir_str).name}: {type(e).__name__}: {e}"}]


# ---------------------------------------------------------------------------
# Job resolution
# ---------------------------------------------------------------------------

def _resolve_jobs(requested: int) -> int:
    cores = int(os.environ.get("SLURM_CPUS_PER_TASK", 0)) or (os.cpu_count() or 1)
    if requested and requested > 0:
        return max(1, min(requested, cores))
    return max(1, cores)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--base", required=True, metavar="DIR",
                   help="Root directory of raw simulation output "
                        "(same input layout as export_summary_data.py)")
    p.add_argument("--out-csv", default=None, metavar="PATH",
                   help="Output CSV path "
                        "(default: <base>/realized_trait_scores.csv)")
    p.add_argument("--quartiles", nargs="+", type=float,
                   default=[0.75, 0.50, 0.25], metavar="Q",
                   help="Substrate thresholds as fractions of per-pulse max "
                        "(default: 0.75 0.5 0.25)")
    p.add_argument("--scenarios", nargs="+", default=None, metavar="NAME",
                   help="Subset of scenario names to process "
                        "(default: all discovered scenarios)")
    p.add_argument("--jobs", "-j", type=int, default=0, metavar="N",
                   help="Parallel worker processes (0 = all cores; "
                        "SLURM_CPUS_PER_TASK aware)")
    args = p.parse_args()

    base_dir = Path(args.base)
    out_csv  = (Path(args.out_csv) if args.out_csv
                else base_dir / "realized_trait_scores.csv")
    quartiles = sorted(args.quartiles, reverse=True)
    workers   = _resolve_jobs(args.jobs)

    # ── Discover scenarios ────────────────────────────────────────────────────
    scenarios = discover_scenarios(base_dir)
    if not scenarios:
        sys.exit(f"ERROR: no scenarios found under {base_dir}")

    if args.scenarios:
        requested = set(args.scenarios)
        scenarios = [(label, sdir) for label, sdir in scenarios
                     if label in requested]
        if not scenarios:
            sys.exit(f"ERROR: none of the requested scenarios found. "
                     f"Available: {[l for l, _ in discover_scenarios(base_dir)]}")

    # ── Build task list ───────────────────────────────────────────────────────
    tasks = []
    per_scenario_counts = {}
    for label, sdir in scenarios:
        runs = pick_runs(sdir, min_cols=100, sentinel="MicrobesSeries")
        per_scenario_counts[label] = len(runs)
        temp = _parse_temperature(label)
        for rd in runs:
            tasks.append((label, temp, str(rd), quartiles))

    if not tasks:
        sys.exit("ERROR: no community run directories found.")

    print(f"Base dir    : {base_dir}")
    print(f"Scenarios   : {len(scenarios)}")
    for label, c in per_scenario_counts.items():
        print(f"    {label}: {c} communities")
    print(f"Total runs  : {len(tasks)}")
    print(f"Quartiles   : {quartiles}")
    print(f"Workers     : {workers}\n")

    # ── Process runs ──────────────────────────────────────────────────────────
    all_rows, errors = [], []
    done, total = 0, len(tasks)

    if workers == 1:
        for t in tasks:
            res = _process_run(t)
            done += 1
            for r in res:
                (errors if "__error__" in r else all_rows).append(r)
            if done % 10 == 0 or done == total:
                print(f"  [{done}/{total}] processed", flush=True)
    else:
        with ProcessPoolExecutor(max_workers=workers) as ex:
            futs = {ex.submit(_process_run, t): t for t in tasks}
            for fut in as_completed(futs):
                res = fut.result()
                done += 1
                for r in res:
                    (errors if "__error__" in r else all_rows).append(r)
                if done % 10 == 0 or done == total:
                    print(f"  [{done}/{total}] processed", flush=True)

    # ── Report errors ─────────────────────────────────────────────────────────
    if errors:
        print(f"\n  {len(errors)} run(s) failed:")
        for e in errors[:20]:
            print(f"    {e['__error__']}")
        if len(errors) > 20:
            print(f"    ... and {len(errors) - 20} more")

    if not all_rows:
        sys.exit("ERROR: no trait rows produced from any community.")

    # ── Assemble and sort output ──────────────────────────────────────────────
    df = pd.DataFrame(all_rows)
    df = df.sort_values(
        ["temperature", "scenario", "community", "taxon", "pulse", "quartile_pct"],
        ascending=[True, True, True, True, True, False],
    ).reset_index(drop=True)

    out_csv.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_csv, index=False)
    print(f"\nWrote {len(df)} rows × {df.shape[1]} columns → {out_csv}")

    # ── Summary stats ─────────────────────────────────────────────────────────
    n_alive = df["alive"].sum()
    n_dead  = (~df["alive"]).sum()
    print(f"Taxon-crossing records: {len(df)} "
          f"({n_alive} alive, {n_dead} dead/no-crossing)")
    print(f"Unique taxa observed  : {df['taxon'].nunique()}")
    print(f"Communities           : {df['community'].nunique()}")


if __name__ == "__main__":
    main()
