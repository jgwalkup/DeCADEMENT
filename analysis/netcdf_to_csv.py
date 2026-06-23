#!/usr/bin/env python3
"""
netcdf_to_csv.py  —  Convert DEMENTpy NetCDF output files to CSV
=================================================================

Converts all .nc files in a DEMENTpy replicate output directory to .csv,
matching the format produced by output_format=csv in runtime.txt.

USAGE
-----
    # Convert a single replicate
    python netcdf_to_csv.py path/to/control/control_00001/

    # Convert all replicates in a scenario directory
    python netcdf_to_csv.py path/to/control/ --recursive

    # Convert and write CSVs to a different location
    python netcdf_to_csv.py path/to/control/control_00001/ --out-dir path/to/csvs/

    # Convert only specific variables
    python netcdf_to_csv.py path/to/control/control_00001/ --vars SubstratesSeries MicrobesSeries

    # Preview what would be converted without writing anything
    python netcdf_to_csv.py path/to/control/ --recursive --dry-run

OPTIONS
-------
    run_dir             Path to a single replicate folder, or a parent folder
                        containing multiple replicate folders (with --recursive)
    --recursive, -r     Recurse into subdirectories and convert all replicates found
    --out-dir DIR       Write CSVs to DIR instead of alongside the .nc files.
                        In recursive mode, mirrors the replicate subfolder structure.
    --vars VAR [VAR]    Only convert the named variables (e.g. SubstratesSeries RespSeries).
                        By default all .nc files are converted.
    --skip-existing     Skip conversion if a .csv already exists for a given .nc file
    --dry-run           Print what would be converted without writing any files
    --quiet, -q         Suppress per-file progress messages; only print summary
"""

import argparse
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr


# ── Conversion helpers ────────────────────────────────────────────────────────

def nc_to_dataframe(nc_path: Path) -> pd.DataFrame | pd.Series | None:
    """
    Load a single .nc file and return a DataFrame or Series.

    DEMENTpy NetCDF files fall into three shapes:
      - DataFrame variables (e.g. SubstratesSeries, MicrobesSeries):
            dimensions are (index, columns) or similar 2-D layout
      - Series variables (e.g. RespSeries, CUE_system):
            single dimension 'index'
      - Scalars (scalars.nc):
            xarray Dataset with scalar data variables, no dimensions
    """
    try:
        ds = xr.open_dataset(nc_path, engine="netcdf4")
    except Exception as e:
        warnings.warn(f"Could not open {nc_path.name}: {e}")
        return None

    # scalars.nc: Dataset with scalar (0-d) variables
    if all(ds[v].ndim == 0 for v in ds.data_vars):
        scalars = {v: float(ds[v].values) for v in ds.data_vars}
        ds.close()
        return pd.Series(scalars)

    # Single-variable DataArray (most Series and DataFrame outputs)
    if len(ds.data_vars) == 1:
        var_name = list(ds.data_vars)[0]
        da = ds[var_name]
        ds.close()

        if da.ndim == 1:
            # 1-D → pandas Series
            coord_name = da.dims[0]
            return pd.Series(da.values, index=da.coords[coord_name].values,
                             name=var_name)

        elif da.ndim == 2:
            # 2-D → pandas DataFrame
            dim0, dim1 = da.dims
            return pd.DataFrame(da.values,
                                index=da.coords[dim0].values,
                                columns=da.coords[dim1].values)

        elif da.ndim > 2:
            # Higher-dimensional (e.g. MultiIndex Series unstack to 3-D)
            # Flatten to 2-D by stacking all but the last dimension
            df = da.to_series().unstack()
            ds.close()
            return df

    # Multi-variable Dataset → return as DataFrame (variables become columns)
    try:
        df = ds.to_dataframe()
        ds.close()
        return df
    except Exception as e:
        warnings.warn(f"Could not convert {nc_path.name} to DataFrame: {e}")
        ds.close()
        return None


def write_csv(data: pd.DataFrame | pd.Series, csv_path: Path) -> None:
    """Write a DataFrame or Series to CSV, matching DEMENTpy's CSV output style."""
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    if isinstance(data, pd.Series):
        data.to_csv(csv_path, header=[data.name or "value"])
    else:
        data.to_csv(csv_path)


def convert_file(nc_path: Path,
                 out_dir: Path | None,
                 skip_existing: bool,
                 dry_run: bool,
                 quiet: bool) -> bool:
    """
    Convert a single .nc file to .csv.
    Returns True if conversion was performed (or would be in dry-run), False if skipped.
    """
    # Determine output path
    if out_dir is not None:
        csv_path = out_dir / nc_path.with_suffix(".csv").name
    else:
        csv_path = nc_path.with_suffix(".csv")

    if skip_existing and csv_path.exists():
        if not quiet:
            print(f"  SKIP  {nc_path.name}  (CSV already exists)")
        return False

    if dry_run:
        print(f"  DRY   {nc_path}  →  {csv_path}")
        return True

    data = nc_to_dataframe(nc_path)
    if data is None:
        print(f"  ERROR {nc_path.name}  (could not load)")
        return False

    write_csv(data, csv_path)
    if not quiet:
        shape = data.shape if hasattr(data, "shape") else f"({len(data)},)"
        print(f"  OK    {nc_path.name}  →  {csv_path.name}  {shape}")
    return True


def find_nc_files(run_dir: Path, var_names: list[str] | None) -> list[Path]:
    """Return all .nc files in run_dir, optionally filtered by variable name."""
    nc_files = sorted(run_dir.glob("*.nc"))
    if var_names:
        nc_files = [f for f in nc_files if f.stem in var_names]
    return nc_files


def is_replicate_dir(path: Path) -> bool:
    """Heuristic: a replicate directory contains at least one .nc or .csv file."""
    return path.is_dir() and (
        any(path.glob("*.nc")) or any(path.glob("*.csv"))
    )


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Convert DEMENTpy NetCDF outputs to CSV.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("run_dir", type=Path,
                        help="Replicate folder, or parent folder (with --recursive)")
    parser.add_argument("--recursive", "-r", action="store_true",
                        help="Recurse into subdirectories and convert all replicates found")
    parser.add_argument("--out-dir", type=Path, default=None,
                        help="Write CSVs here instead of alongside the .nc files")
    parser.add_argument("--vars", nargs="+", default=None, metavar="VAR",
                        help="Only convert these variable names (e.g. SubstratesSeries)")
    parser.add_argument("--skip-existing", action="store_true",
                        help="Skip if a .csv already exists for a given .nc file")
    parser.add_argument("--dry-run", action="store_true",
                        help="Print what would be converted without writing files")
    parser.add_argument("--quiet", "-q", action="store_true",
                        help="Suppress per-file messages; only print summary")
    args = parser.parse_args()

    run_dir: Path = args.run_dir.resolve()
    if not run_dir.exists():
        print(f"ERROR: path does not exist: {run_dir}", file=sys.stderr)
        sys.exit(1)

    # Build list of (replicate_dir, out_dir) pairs to process
    tasks: list[tuple[Path, Path | None]] = []

    if args.recursive:
        # Find all subdirectories that look like replicate folders
        subdirs = [d for d in sorted(run_dir.iterdir()) if is_replicate_dir(d)]
        if not subdirs:
            # run_dir itself might be the replicate
            subdirs = [run_dir]
        for sub in subdirs:
            if args.out_dir:
                out = args.out_dir / sub.name
            else:
                out = None
            tasks.append((sub, out))
    else:
        tasks.append((run_dir, args.out_dir))

    # Process
    total_converted = 0
    total_skipped = 0
    total_errors = 0

    for rep_dir, out_dir in tasks:
        nc_files = find_nc_files(rep_dir, args.vars)
        if not nc_files:
            if not args.quiet:
                print(f"  (no .nc files found in {rep_dir.name})")
            continue

        if not args.quiet or args.recursive:
            print(f"\n{'[DRY RUN] ' if args.dry_run else ''}Converting: {rep_dir.name}")

        for nc_path in nc_files:
            result = convert_file(nc_path, out_dir, args.skip_existing,
                                  args.dry_run, args.quiet)
            if result is True:
                total_converted += 1
            elif result is False and args.skip_existing:
                total_skipped += 1
            elif result is False:
                total_errors += 1

    # Summary
    action = "Would convert" if args.dry_run else "Converted"
    print(f"\n{action} {total_converted} file(s)", end="")
    if total_skipped:
        print(f", skipped {total_skipped}", end="")
    if total_errors:
        print(f", {total_errors} error(s)", end="")
    print()


if __name__ == "__main__":
    main()
