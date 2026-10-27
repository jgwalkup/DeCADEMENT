"""
mass_balance.py  –  Carbon mass balance checker for DEMENTpy simulation outputs
================================================================================

USAGE
-----
    # Check individual replicate directories (or parent dirs):
    python mass_balance.py <run_dir> [<run_dir2> ...]

    # Multi-scenario convenience: walk <output-dir>/<scenario>/<rep>/ layout
    python mass_balance.py --scenario-dir <output-dir>

    Each <run_dir> should be the output directory of a single DEMENTpy replicate,
    i.e. the folder that contains SubstratesSeries.csv/.nc, MonomersSeries.csv/.nc,
    etc.  Both CSV and NetCDF output formats are supported and auto-detected.

    You can also pass a parent directory that contains multiple replicate sub-dirs:
        python mass_balance.py output/control/

WHAT IT CHECKS  (Carbon only)
------------------------------
For each pulse p independently:

    Initial C(p) =  [Sub + Mon + Mic + Enz]  at t=0
    Final C(p)   =  [Sub + Mon + Mic + Enz]  at the last day of pulse p
                 +  CumRespByPulse[p]   (cumulative respiration within pulse p)

    Error (%)    =  100 × (Final − Initial) / Initial

The initial state is always the t=0 snapshot for every pulse.  This is correct
because repopulation() resets Substrates, Monomers, and Enzymes back to their
original initial values (Substrates_init etc.) at the end of each pulse, before
the next pulse begins.  The t=0 snapshot captures exactly those reset values.

FORMAT DETECTION
----------------
CSV runs    (output_format=csv):    SubstratesSeries.csv etc.
NetCDF runs (output_format=netcdf): SubstratesSeries.nc etc. via xarray.
Zarr runs   (output_format=zarr):   <outname>.zarr store with groups
                                    per variable, read via xarray+zarr.
Auto-detection checks for each format in that order.

REQUIREMENTS
------------
Requires output.py to contain the CumRespByPulse series (added in the
accompanying output.py patch).  NetCDF runs additionally require xarray:
    pip install xarray

NOTES ON N AND P
----------------
N and P are not fully checkable from standard outputs because only the C
dimension of each pool is saved as a time series.  To add N/P balance tracking
you would need to enable the commented-out NH4Series/PO4Series outputs and add
a cumulative leaching tracker in output.py.

OUTPUT FILES
------------
  mass_balance_summary.csv  –  one row per (replicate, pulse) with all pool
                                values and error %
  mass_balance_report.txt   –  human-readable summary; also printed to stdout
"""

import sys
import textwrap
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from branch_output import read_table  # noqa: E402


# ── helpers ───────────────────────────────────────────────────────────────────

def _is_csv_run(run_dir: Path) -> bool:
    return (run_dir / "SubstratesSeries.csv").exists()


def _is_nc_run(run_dir: Path) -> bool:
    return (run_dir / "SubstratesSeries.nc").exists()

def _is_zarr_run(run_dir: Path) -> bool:
    # Case 1: run_dir itself is the .zarr store
    if run_dir.suffix == ".zarr" and run_dir.is_dir():
        return True
    # Case 2: run_dir contains a .zarr store (non-init)
    return any(p for p in run_dir.glob("*.zarr")
               if not p.name.endswith(".init.zarr"))

# ── CSV loaders ───────────────────────────────────────────────────────────────

def _load_series_csv(run_dir: Path, name: str) -> pd.DataFrame | None:
    """Load a pool-series DataFrame from CSV.  Columns normalised to int.

    branch_output.read_table joins split branched output (spin-up days +
    branch days) into the full series; unsplit output is read as is."""
    p = run_dir / f"{name}.csv"
    if not p.exists():
        return None
    df = read_table(run_dir, name)
    df.columns = [int(c) if str(c).lstrip("-").isdigit() else c for c in df.columns]
    return df


# ── NetCDF loaders ────────────────────────────────────────────────────────────

def _load_series_nc(run_dir: Path, name: str) -> pd.DataFrame | None:
    """Load a named series from NetCDF, returning a DataFrame.

    DataFrames are saved by export_to_netcdf() as xr.Datasets (pool names as
    the index dimension, day numbers as string-named variables).  pd.Series are
    saved as individual xr.DataArrays (dims=["index"]).  This function handles
    both by trying open_dataset first and falling back to open_dataarray, so
    that load_series and load_cum_resp_by_pulse can share the same code path.
    """
    try:
        import xarray as xr
    except ImportError:
        raise ImportError("xarray is required to read NetCDF output. "
                          "Install it with: pip install xarray")
    p = run_dir / f"{name}.nc"
    if not p.exists():
        return None
    try:
        ds = xr.open_dataset(p)
        df = ds.to_dataframe()
    except Exception:
        # File was written as a DataArray (e.g. a pd.Series export)
        da = xr.open_dataarray(p)
        df = da.to_series().to_frame()
    # Normalise column labels to integers where possible
    df.columns = [int(c) if str(c).lstrip("-").isdigit() else c for c in df.columns]
    return df

# ── Zarr loaders ──────────────────────────────────────────────────────────────

def _get_zarr_store(run_dir: Path) -> str | None:
    """Return the path to the main Zarr store for a run directory."""
    # Case 1: run_dir itself is the store
    if run_dir.suffix == ".zarr" and run_dir.is_dir():
        return str(run_dir)
    # Case 2: store is inside run_dir
    stores = [p for p in run_dir.glob("*.zarr")
              if not p.name.endswith(".init.zarr")]
    return str(stores[0]) if stores else None


def _load_series_zarr(run_dir: Path, name: str) -> pd.DataFrame | None:
    try:
        import xarray as xr
    except ImportError:
        raise ImportError("xarray and zarr are required to read Zarr output. "
                          "Install with: pip install xarray zarr")

    store = _get_zarr_store(run_dir)
    if store is None:
        return None

    try:
        ds = xr.open_zarr(store, group=name)
        df = ds.to_dataframe()
    except Exception:
        return None

    df.columns = [int(c) if str(c).lstrip("-").isdigit() else c for c in df.columns]
    return df

# ── format-agnostic loaders ───────────────────────────────────────────────────

def load_series(run_dir: Path, name: str) -> pd.DataFrame | None:
    """Load a pool-series DataFrame from whichever format is present (CSV, NC, or Zarr)."""
    if _is_csv_run(run_dir):
        return _load_series_csv(run_dir, name)
    if _is_zarr_run(run_dir):
        return _load_series_zarr(run_dir, name)
    return _load_series_nc(run_dir, name)


def _to_pulse_series(df: pd.DataFrame | None) -> pd.Series | None:
    """Extract a pulse-indexed Series from a single-column DataFrame.

    Pulse-indexed outputs (CumRespByPulse, MicC_PulseInit) are loaded by
    load_series like any other file.  This converts the resulting single-column
    DataFrame to a Series with an integer pulse-number index.
    """
    if df is None:
        return None
    s = df.iloc[:, 0]
    s.index = s.index.astype(int)
    return s


# ── remaining helpers ─────────────────────────────────────────────────────────

def infer_pulse_structure(sub_df: pd.DataFrame,
                          cum_resp: pd.Series) -> tuple[int, int]:
    """Infer n_pulses and days_per_pulse from available data.

    n_pulses is taken from the length of CumRespByPulse.
    days_per_pulse = (max day column) / n_pulses.
    """
    n_pulses = len(cum_resp)
    int_cols = sorted(c for c in sub_df.columns if isinstance(c, int) and c > 0)
    if not int_cols:
        raise ValueError("SubstratesSeries has no integer day columns > 0")
    max_day = max(int_cols)
    cycle = max_day // n_pulses
    return n_pulses, cycle


def pool_total_at(series_df: pd.DataFrame, day: int) -> float:
    """Sum all pool rows at a given day column.  Returns nan if the column is
    absent (never substitutes a neighbouring day: a missing day means the run
    does not hold that state, and a silent substitute would hide it)."""
    if series_df is None or day not in series_df.columns:
        return np.nan
    return float(series_df[day].sum())

def find_run_dirs(paths: list[str]) -> list[Path]:
    candidates = []
    for raw in paths:
        p = Path(raw)
        if not p.exists():
            print(f"  [WARN] Path not found: {p}")
            continue
        if _is_csv_run(p) or _is_nc_run(p) or _is_zarr_run(p):
            candidates.append(p)
        else:
            subdirs = [d for d in sorted(p.iterdir())
                       if d.is_dir() and (_is_csv_run(d) or _is_nc_run(d) or _is_zarr_run(d))]
            if subdirs:
                candidates.extend(subdirs)
            else:
                print(f"  [WARN] No SubstratesSeries.csv/.nc or *.zarr found in {p} "
                      f"or its sub-dirs")
    return candidates


def find_run_dirs_multi_scenario(scenario_dir: Path) -> list[tuple[str, Path]]:
    """Walk <scenario_dir>/<scenario_stem>/<rep_dir>/ and return (scenario, run_dir) pairs.

    Falls back to treating scenario_dir itself as a single scenario when no
    qualifying scenario subdirectories are found (flat layout).
    """
    if not scenario_dir.exists():
        print(f"  [WARN] --scenario-dir not found: {scenario_dir}")
        return []

    results: list[tuple[str, Path]] = []

    for scen_d in sorted(scenario_dir.iterdir()):
        if not scen_d.is_dir() or scen_d.name == "spinup":   # spin-up states of branched runs are not a scenario
            continue
        rep_dirs = [d for d in sorted(scen_d.iterdir())
                    if d.is_dir() and (_is_csv_run(d) or _is_nc_run(d) or _is_zarr_run(d))]
        if rep_dirs:
            for rd in rep_dirs:
                results.append((scen_d.name, rd))

    if not results:
        # Flat layout: scenario_dir itself contains replicate dirs
        flat = [d for d in sorted(scenario_dir.iterdir())
                if d.is_dir() and (_is_csv_run(d) or _is_nc_run(d) or _is_zarr_run(d))]
        label = scenario_dir.name
        results = [(label, rd) for rd in flat]

    if not results:
        print(f"  [WARN] No run directories found under {scenario_dir}")
    return results


# ── per-run mass balance ───────────────────────────────────────────────────────

def check_run(run_dir: Path, threshold_pct: float = 0.1,
              scenario: str | None = None) -> list[dict]:
    """Compute per-pulse carbon mass balance for one replicate directory.

    Returns a list of dicts, one per pulse.
    """
    base = {"scenario": scenario or "", "run": str(run_dir)}

    sub        = load_series(run_dir, "SubstratesSeries")
    mon        = load_series(run_dir, "MonomersSeries")
    mic        = load_series(run_dir, "MicrobesSeries")
    enz        = load_series(run_dir, "EnzymesSeries")
    cum_resp   = _to_pulse_series(load_series(run_dir, "CumRespByPulse"))
    mic_c_init = _to_pulse_series(load_series(run_dir, "MicC_PulseInit"))

    # Check for missing files
    missing = [name for name, obj in [
        ("SubstratesSeries", sub), ("MonomersSeries", mon),
        ("MicrobesSeries",   mic), ("EnzymesSeries",  enz),
    ] if obj is None]
    if missing:
        return [{**base, "error": f"Missing pool files: {', '.join(missing)}"}]

    if cum_resp is None:
        return [{**base, "error": (
            "CumRespByPulse not found. "
            "Re-run the simulation with the updated output.py to generate this file."
        )}]

    if mic_c_init is None:
        return [{**base, "error": (
            "MicC_PulseInit not found. "
            "Re-run the simulation with the updated output.py to generate this file."
        )}]

    # Infer pulse structure
    try:
        n_pulses, cycle = infer_pulse_structure(sub, cum_resp)
    except ValueError as e:
        return [{**base, "error": str(e)}]

    results = []
    for p in range(1, n_pulses + 1):
        # repopulation() resets Substrates, Monomers, and Enzymes back to their
        # initial values (Substrates_init etc.) at the end of every pulse.  The
        # t=0 snapshot captures exactly those reset values, so it is the correct
        # initial state for Sub/Mon/Enz in every pulse.
        # Microbes are different: repopulation() resets to Microbes_pp (pre-
        # placement) and re-places stochastically, so the actual post-repopulation
        # microbial C is recorded in MicC_PulseInit and used here instead.
        day_start = 0             # post-reset initial state for Sub/Mon/Enz
        day_end   = p * cycle     # last recorded day of this pulse

        sub_0 = pool_total_at(sub, day_start)
        mon_0 = pool_total_at(mon, day_start)
        mic_0 = float(mic_c_init.get(p, np.nan))   # actual post-repopulation microbial C
        enz_0 = pool_total_at(enz, day_start)
        total_0 = sub_0 + mon_0 + mic_0 + enz_0

        sub_T = pool_total_at(sub, day_end)
        mon_T = pool_total_at(mon, day_end)
        mic_T = pool_total_at(mic, day_end)
        enz_T = pool_total_at(enz, day_end)

        resp_p = float(cum_resp.get(p, np.nan))
        total_T = sub_T + mon_T + mic_T + enz_T + resp_p

        error_abs = total_T - total_0
        error_pct = 100.0 * error_abs / total_0 if total_0 != 0 else np.nan
        flag = "OK" if abs(error_pct) <= threshold_pct else "FAIL"

        results.append({
            **base,
            "pulse":     p,
            "n_pulses":  n_pulses,
            "cycle":     cycle,
            "day_start": day_start,
            "day_end":   day_end,
            # Initial pools
            "Sub_C_0":   sub_0,
            "Mon_C_0":   mon_0,
            "Mic_C_0":   mic_0,
            "Enz_C_0":   enz_0,
            "Total_C_0": total_0,
            # Final pools
            "Sub_C_T":   sub_T,
            "Mon_C_T":   mon_T,
            "Mic_C_T":   mic_T,
            "Enz_C_T":   enz_T,
            "Cum_Resp":  resp_p,
            "Total_C_T": total_T,
            # Balance
            "Error_abs": error_abs,
            "Error_pct": error_pct,
            "Flag":      flag,
        })

    return results


# ── report formatting ──────────────────────────────────────────────────────────

HEADER = textwrap.dedent("""
    ╔══════════════════════════════════════════════════════════╗
    ║         DEMENTpy Carbon Mass Balance Checker             ║
    ╚══════════════════════════════════════════════════════════╝
""").strip()


def format_run_results(run_results: list[dict]) -> str:
    """Format results for one run (potentially multiple pulses) as a string."""
    if not run_results:
        return ""

    lines = [f"\n  Run : {run_results[0]['run']}"]

    if "error" in run_results[0]:
        lines.append(f"  !! {run_results[0]['error']}")
        return "\n".join(lines)

    n_pulses = run_results[0]["n_pulses"]
    cycle    = run_results[0]["cycle"]
    lines.append(f"  Pulses: {n_pulses}   Days per pulse: {cycle}")

    for r in run_results:
        p = r["pulse"]
        pulse_day_start = (p - 1) * cycle + 1
        lines.append(f"\n  ── Pulse {p}/{n_pulses}  "
                     f"(days {pulse_day_start}–{r['day_end']}) {'─'*30}")
        lines.append(f"     Initial C (t=0 reset, Mic=post-repop):  "
                     f"Sub={r['Sub_C_0']:>12.4f}  "
                     f"Mon={r['Mon_C_0']:>12.4f}  "
                     f"Mic={r['Mic_C_0']:>12.4f}  "
                     f"Enz={r['Enz_C_0']:>12.4f}")
        lines.append(f"       TOTAL = {r['Total_C_0']:>14.4f}")
        lines.append(f"     Final C   (t={r['day_end']:>4d}):        Sub={r['Sub_C_T']:>12.4f}  "
                     f"Mon={r['Mon_C_T']:>12.4f}  "
                     f"Mic={r['Mic_C_T']:>12.4f}  "
                     f"Enz={r['Enz_C_T']:>12.4f}")
        lines.append(f"     Cum. Resp (pulse {p:>2d}):          {r['Cum_Resp']:>14.4f}")
        lines.append(f"       TOTAL = {r['Total_C_T']:>14.4f}")
        lines.append(f"     Balance:  abs={r['Error_abs']:>+14.6f}   "
                     f"rel={r['Error_pct']:>+10.6f} %   [{r['Flag']}]")

    return "\n".join(lines)


# ── main ──────────────────────────────────────────────────────────────────────

def main():
    import argparse

    parser = argparse.ArgumentParser(
        description="DEMENTpy carbon mass balance checker.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "run_dirs", nargs="*", metavar="run_dir",
        help="One or more replicate output directories (or parent directories).",
    )
    parser.add_argument(
        "--scenario-dir", metavar="DIR", default=None,
        help="Root output directory produced by run_multi_climate_simulations.py. "
             "Each subdirectory is treated as a named scenario; its replicate "
             "subdirectories are checked individually.  Adds a 'scenario' column "
             "to the CSV output.  Mutually exclusive with positional run_dirs.",
    )
    parser.add_argument(
        "--treatment-map", nargs="+", default=None, metavar="DIR=LABEL",
        help="Map directory stems to display labels in the scenario column. "
             "Each entry is DIR_STEM=LABEL, e.g.: "
             "--treatment-map static_T10_Pn0p5=Control static_T15_Pn0p5='Warm +2C'. "
             "Unmatched stems are used as-is.",
    )
    parser.add_argument(
        "--report-dir", metavar="DIR", default=None,
        help="Directory in which to save report files (default: current directory).",
    )
    parser.add_argument(
        "--report-id", metavar="ID", default=None,
        help="Identifier appended to output filenames, e.g. 'ctrl' produces "
             "mass_balance_summary_ctrl.csv and mass_balance_report_ctrl.txt.",
    )
    args = parser.parse_args()

    if args.scenario_dir and args.run_dirs:
        print("ERROR: --scenario-dir and positional run_dirs are mutually exclusive.")
        sys.exit(1)
    if not args.scenario_dir and not args.run_dirs:
        parser.print_help()
        sys.exit(1)

    threshold = 0.1   # % – tweak as needed

    # Resolve output directory
    out_dir = Path(args.report_dir) if args.report_dir else Path(".")
    out_dir.mkdir(parents=True, exist_ok=True)

    # Build file stems
    suffix   = f"_{args.report_id}" if args.report_id else ""
    csv_path = out_dir / f"mass_balance_summary{suffix}.csv"
    txt_path = out_dir / f"mass_balance_report{suffix}.txt"

    # Build treatment label map
    tmap = {}
    if args.treatment_map:
        for entry in args.treatment_map:
            if "=" not in entry:
                print(f"  WARNING: --treatment-map entry has no '=': {entry!r}  (ignored)")
                continue
            stem, label = entry.split("=", 1)
            tmap[stem.strip()] = label.strip()

    # Build (scenario, run_dir) pairs
    if args.scenario_dir:
        scenario_run_pairs = find_run_dirs_multi_scenario(Path(args.scenario_dir))
        if not scenario_run_pairs:
            print("No valid run directories found under --scenario-dir. Exiting.")
            sys.exit(1)
        # Apply treatment-map labels
        scenario_run_pairs = [(tmap.get(s, s), rd) for s, rd in scenario_run_pairs]
        run_dirs    = [rd   for _, rd   in scenario_run_pairs]
        scenarios_m = [scen for scen, _ in scenario_run_pairs]
    else:
        run_dirs   = find_run_dirs(args.run_dirs)
        scenarios_m = [None] * len(run_dirs)
        if not run_dirs:
            print("No valid run directories found. Exiting.")
            sys.exit(1)

    print(HEADER)
    print(f"\n  Threshold for FAIL: ±{threshold} %")
    n_scenarios = len(set(s for s in scenarios_m if s))
    if n_scenarios > 1:
        print(f"  Scenarios detected: {sorted(set(s for s in scenarios_m if s))}")
    print(f"  Checking {len(run_dirs)} run(s)...\n")

    all_results: list[dict] = []
    report_lines = [HEADER, f"\nThreshold: ±{threshold} %\n"]

    for rd, scen in zip(run_dirs, scenarios_m):
        run_results = check_run(rd, threshold_pct=threshold, scenario=scen)
        all_results.extend(run_results)
        block = format_run_results(run_results)
        print(block)
        report_lines.append(block)

    # ── summary table ──────────────────────────────────────────────────────────
    divider = "─" * 70
    print(f"\n{divider}")
    print("  SUMMARY")
    print(divider)
    report_lines += [f"\n{divider}", "SUMMARY", divider]

    ok_count   = sum(1 for r in all_results if r.get("Flag") == "OK")
    fail_count = sum(1 for r in all_results if r.get("Flag") == "FAIL")
    err_count  = sum(1 for r in all_results if "error" in r)
    total_pulse_checks = ok_count + fail_count

    summary = (f"  Pulse-level checks: {total_pulse_checks}   "
               f"OK: {ok_count}   FAIL: {fail_count}   "
               f"Errors (missing files): {err_count}")
    print(summary)
    report_lines.append(summary)

    # Per-scenario breakdown when multi-scenario
    if n_scenarios > 1:
        scen_names = sorted(set(r.get("scenario", "") for r in all_results
                                if "error" not in r))
        for scen in scen_names:
            s_ok   = sum(1 for r in all_results
                         if r.get("scenario") == scen and r.get("Flag") == "OK")
            s_fail = sum(1 for r in all_results
                         if r.get("scenario") == scen and r.get("Flag") == "FAIL")
            line   = f"    {scen:30s}  OK: {s_ok}   FAIL: {s_fail}"
            print(line);  report_lines.append(line)

    if fail_count:
        print("\n  Failing pulses:")
        report_lines.append("\n  Failing pulses:")
        for r in all_results:
            if r.get("Flag") == "FAIL":
                scen_tag = f"[{r['scenario']}] " if r.get("scenario") else ""
                line = (f"    {scen_tag}{r['run']}  pulse {r['pulse']}  "
                        f"({r['Error_pct']:+.4f} %)")
                print(line)
                report_lines.append(line)

    # ── write CSV ──────────────────────────────────────────────────────────────
    good = [r for r in all_results if "error" not in r]
    if good:
        df_out = pd.DataFrame(good)
        # Put scenario first if present
        cols = df_out.columns.tolist()
        if "scenario" in cols and df_out["scenario"].str.len().max() > 0:
            cols = ["scenario"] + [c for c in cols if c != "scenario"]
            df_out = df_out[cols]
        df_out.to_csv(csv_path, index=False)
        print(f"\n  CSV saved  → {csv_path.resolve()}")
        report_lines.append(f"\nCSV saved → {csv_path.resolve()}")

    # ── write report ───────────────────────────────────────────────────────────
    txt_path.write_text("\n".join(report_lines))
    print(f"  Report saved → {txt_path.resolve()}\n")


if __name__ == "__main__":
    main()
