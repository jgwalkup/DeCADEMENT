"""
run_community_simulations.py
============================
Run DeCA-DEMENT simulations for a set of pre-assembled communities produced
by assemble_communities.py. Each community directory is used in place of the
random microbial initialization that dementpy.py normally performs.

This script is a drop-in replacement for calling dementpy.py directly when
running community-based experiments. It handles:

  * Iterating over community directories in a specified folder
  * Passing each community directory to initialize_data() via the
    community_dir argument added to initialization.py
  * Naming output directories to reflect the community used:
      <output-base>/<community_name>_<seed>/
  * Optionally running control and warming climate scenarios for each
    community (paired design)

Expected project layout
-----------------------
    Decoupled-Energy-DEMENTpy/   <- run from here
    ├── src/                     <- --src-dir
    ├── input/                   <- --input-dir  (runtime.txt, CSVs)
    ├── community_pool_03/
    │   ├── communities/         <- --community-dir
    │   │   ├── community_0001/
    │   │   ├── community_0002/
    │   │   └── ...
    │   └── ReqEnz.csv
    └── analysis/
        └── run_community_simulations.py

Usage
-----
    # Single climate scenario
    python analysis/run_community_simulations.py \\
        --src-dir       src \\
        --input-dir     input \\
        --community-dir community_pool_03/communities \\
        --output-dir    simulations/control \\
        --climate       climate.csv \\
        --seed          42

    # Paired control + warming
    python analysis/run_community_simulations.py \\
        --src-dir           src \\
        --input-dir         input \\
        --community-dir     community_pool_03/communities \\
        --output-dir        simulations \\
        --climate           climate.csv \\
        --climate-warm      climate_warm.csv \\
        --seed              42

Options
-------
  --src-dir         Directory containing DEMENTpy .py source files (default: src)
  --input-dir       Directory containing runtime.txt and input CSVs (default: input)
  --community-dir   Directory of community subdirectories to simulate
  --output-dir      Base directory for simulation outputs
  --climate         Climate CSV filename inside --input-dir (default: climate.csv)
  --climate-warm    Optional warming climate CSV for paired design. When provided,
                    each community is run twice: once with --climate (written to
                    <output-dir>/control/) and once with --climate-warm (written
                    to <output-dir>/warm/)
  --seed            Base random seed (default: 0). Each community gets a unique
                    seed derived from this value + community index, ensuring
                    reproducibility while keeping runs independent
  --communities     Optional comma-separated list of community folder names to
                    run (e.g. community_0001,community_0005). If omitted, all
                    subdirectories in --community-dir are run in sorted order
  --start-from      Optional community folder name to resume from (inclusive).
                    Skips all communities that sort before this name
  --quiet           Suppress per-day progress output from the model

Output directory structure
--------------------------
  Single climate:
    <output-dir>/<community_name>_s<seed>/

  Paired (control + warm):
    <output-dir>/control/<community_name>_s<seed>/
    <output-dir>/warm/<community_name>_s<seed>/
"""

import argparse
import os
import sys
import warnings
import traceback
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore", category=UserWarning, module="zarr")


# ---------------------------------------------------------------------------
# Runtime loading  (mirrors dementpy.py exactly)
# ---------------------------------------------------------------------------

def load_runtime(input_dir, runtime_file='runtime.txt'):
    path = os.path.join(input_dir, runtime_file)
    runtime = pd.read_csv(path, header=None, index_col=0, sep='\t')
    output_format = str(runtime.iloc[runtime.index.get_loc('output_format'), 0]).strip().lower()
    runtime = runtime.drop('output_format')
    runtime = runtime.astype(float)
    # Rename column to integer 1 explicitly — this is what all DEMENTpy source
    # files expect (substrate.py, microbe.py, grid.py all use runtime.loc[key, 1])
    runtime.columns = [1]
    return runtime, output_format


# ---------------------------------------------------------------------------
# Single simulation run
# ---------------------------------------------------------------------------

def run_simulation(runtime, output_format, community_dir, climate_path,
                   output_path, outname, src_dir, quiet=False):
    """
    Run one complete DEMENTpy simulation for a single community.

    Parameters
    ----------
    runtime        : DataFrame — numeric runtime parameters
    output_format  : str — 'csv', 'netcdf', or 'zarr'
    community_dir  : Path — absolute path to the community directory
    climate_path   : str — absolute path to the climate CSV
    output_path    : Path — absolute path to the output directory
    outname        : str — output file label (used as filename prefix)
    src_dir        : str — absolute path to DEMENTpy source directory
    quiet          : bool — suppress per-day print statements
    """
    # Imports must be done after sys.path is set
    sys.path.insert(0, src_dir)
    from initialization import initialize_data
    from grid import Grid
    from output import Output

    # Use the actual column label rather than hardcoding 1 — it may be
    # integer or string depending on pandas version / load path.
    col = runtime.columns[0]

    pulse    = int(runtime.loc['pulse',    col])
    cycle    = int(runtime.loc['end_time', col])
    interval = int(runtime.loc['interval', col])
    mic_reinit = True

    # Capture project root once here so both chdir calls restore correctly.
    project_root = os.getcwd()
    input_dir = os.path.dirname(climate_path)

    os.chdir(input_dir)
    try:
        data = initialize_data(runtime, climate_file=climate_path,
                               community_dir=str(community_dir))
    finally:
        os.chdir(project_root)

    Output_init = Output(runtime, data)
    Ecosystem   = Grid(runtime, data)

    for p in range(pulse):
        for i in range(p * cycle, (p + 1) * cycle):
            if not quiet:
                print(f"Day {i} starting", flush=True)

            Ecosystem.degradation(p, i)
            Ecosystem.uptake(p, i)

            try:
                Ecosystem.metabolism(i)
            except Exception as e:
                print(f"CRASH in metabolism() on day {i}:", flush=True)
                traceback.print_exc()
                raise

            Ecosystem.mortality(i)
            Ecosystem.reproduction(i)

            if i == 0:
                Output_init.output(Ecosystem, i)
            elif i % interval == interval - 1:
                Output_init.output(Ecosystem, i)

            Output_init.microbes_abundance(Ecosystem, i)

            if i == (p + 1) * cycle - 1:
                Ecosystem.repopulation(Output_init, i, mic_reinit)

    # Save output
    output_path.mkdir(parents=True, exist_ok=True)
    os.chdir(output_path)
    try:
        if output_format == 'csv':
            Output_init.export_to_csv(outname)
        elif output_format == 'netcdf':
            Output_init.export_to_netcdf(outname)
        elif output_format == 'zarr':
            Output_init.export_to_zarr(outname)
        else:
            warnings.warn(f"Unknown output_format '{output_format}', defaulting to NetCDF.")
            Output_init.export_to_netcdf(outname)
    finally:
        os.chdir(project_root)


# ---------------------------------------------------------------------------
# Community directory detection  (same logic as community_summary.py)
# ---------------------------------------------------------------------------

def get_community_dirs(community_dir, selection=None, start_from=None):
    """
    Return sorted list of community subdirectories (those containing
    EnzGenes.csv). Optionally filter to a specific selection or resume
    from a given community name.
    """
    all_dirs = sorted(
        [d for d in Path(community_dir).iterdir()
         if d.is_dir() and (d / 'EnzGenes.csv').exists()],
        key=lambda p: p.name
    )

    if selection:
        names = set(selection)
        all_dirs = [d for d in all_dirs if d.name in names]

    if start_from:
        all_dirs = [d for d in all_dirs if d.name >= start_from]

    return all_dirs


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Run DeCA-DEMENT community simulations.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument('--src-dir',        type=str, default='src',
                        help='Directory containing DEMENTpy .py source files')
    parser.add_argument('--input-dir',      type=str, default='input',
                        help='Directory containing runtime.txt and input CSVs')
    parser.add_argument('--community-dir',  type=str, required=True,
                        help='Directory of community subdirectories to simulate')
    parser.add_argument('--output-dir',     type=str, required=True,
                        help='Base directory for simulation outputs')
    parser.add_argument('--climate',        type=str, default='climate.csv',
                        help='Climate CSV filename inside --input-dir')
    parser.add_argument('--climate-warm',   type=str, default=None,
                        help='Warming climate CSV filename inside --input-dir '
                             '(enables paired control/warm design)')
    parser.add_argument('--seed',           type=int, default=0,
                        help='Base random seed; each community gets seed + index')
    parser.add_argument('--communities',    type=str, default=None,
                        help='Comma-separated list of community folder names to '
                             'run (default: all)')
    parser.add_argument('--range',          type=int, nargs=2,
                        metavar=('START', 'END'),
                        help='1-based inclusive range of communities to run by '
                             'position in the sorted list, e.g. --range 1 50 '
                             'runs the first 50 communities. Cannot be combined '
                             'with --communities or --start-from.')
    parser.add_argument('--start-from',     type=str, default=None,
                        help='Resume from this community name (inclusive). '
                             'Cannot be combined with --range.')
    parser.add_argument('--quiet',          action='store_true',
                        help='Suppress per-day progress output')
    args = parser.parse_args()

    # ── Resolve paths ────────────────────────────────────────────────────────
    src_dir       = os.path.abspath(args.src_dir)
    input_dir     = os.path.abspath(args.input_dir)
    community_dir = os.path.abspath(args.community_dir)
    output_dir    = Path(os.path.abspath(args.output_dir))
    paired        = args.climate_warm is not None

    for label, path in [('--src-dir', src_dir),
                         ('--input-dir', input_dir),
                         ('--community-dir', community_dir)]:
        if not os.path.isdir(path):
            print(f"ERROR: {label} '{path}' does not exist.")
            sys.exit(1)

    climate_path      = os.path.join(input_dir, args.climate)
    climate_warm_path = (os.path.join(input_dir, args.climate_warm)
                         if paired else None)

    for label, path in [('--climate', climate_path)] + (
            [('--climate-warm', climate_warm_path)] if paired else []):
        if not os.path.isfile(path):
            print(f"ERROR: {label} file '{path}' not found.")
            sys.exit(1)

    # ── Load runtime once — shared across all communities ────────────────────
    # Note: dementpy.py normally os.chdir into input_dir before reading
    # runtime.txt. We load it here with an explicit path instead, so we can
    # stay in the project root throughout.
    runtime, output_format = load_runtime(input_dir)

    # ── Validate mutually exclusive range options ────────────────────────────
    if args.range and (args.communities or args.start_from):
        print("ERROR: --range cannot be combined with --communities or --start-from.")
        sys.exit(1)

    if args.range:
        start, end = args.range
        if start < 1:
            print("ERROR: --range START must be >= 1.")
            sys.exit(1)
        if end < start:
            print("ERROR: --range END must be >= START.")
            sys.exit(1)

    # ── Collect community directories ─────────────────────────────────────────
    selection = ([s.strip() for s in args.communities.split(',')]
                 if args.communities else None)
    comm_dirs = get_community_dirs(community_dir,
                                   selection=selection,
                                   start_from=args.start_from)

    # Apply --range as a slice over the sorted list (1-based, inclusive)
    if args.range:
        start, end = args.range
        comm_dirs = comm_dirs[start - 1 : end]
        if not comm_dirs:
            print(f"ERROR: --range {start} {end} produced no communities "
                  f"(pool has {len(get_community_dirs(community_dir))} communities).")
            sys.exit(1)

    if not comm_dirs:
        print("ERROR: no community directories found.")
        sys.exit(1)

    n = len(comm_dirs)
    scenario_label = "paired control+warm" if paired else "single climate"
    print(f"Communities : {n}")
    print(f"Scenario    : {scenario_label}")
    print(f"Base seed   : {args.seed}")
    print(f"Output dir  : {output_dir}")
    print()

    # ── Run simulations ───────────────────────────────────────────────────────
    failed = []

    for idx, cdir in enumerate(comm_dirs):
        comm_name   = cdir.name                        # e.g. community_0042
        comm_seed   = (args.seed + idx) % (2**32)     # unique seed per community
        outname     = f"{comm_name}_s{comm_seed}"      # e.g. community_0042_s42
        run_label   = f"[{idx+1}/{n}] {comm_name}  (seed={comm_seed})"

        # Set NumPy seed for this community
        np.random.seed(comm_seed)

        if paired:
            scenarios = [
                ('control', climate_path,     output_dir / 'control'),
                ('warm',    climate_warm_path, output_dir / 'warm'),
            ]
        else:
            scenarios = [
                ('run', climate_path, output_dir),
            ]

        for scenario_name, clim_path, out_path in scenarios:
            label = f"{run_label}  [{scenario_name}]" if paired else run_label
            print(f"{label}", flush=True)

            try:
                run_simulation(
                    runtime       = runtime.copy(),  # copy to prevent mutation across calls
                    output_format = output_format,
                    community_dir = cdir,
                    climate_path  = clim_path,
                    output_path   = out_path,
                    outname       = outname,
                    src_dir       = src_dir,
                    quiet         = args.quiet,
                )
                print(f"  -> saved to {out_path / outname}", flush=True)

            except Exception as e:
                import traceback as _tb
                print(f"  !! FAILED: {e}", flush=True)
                _tb.print_exc()
                failed.append((comm_name, scenario_name, str(e)))
                continue   # move to next scenario / community

    # ── Summary ───────────────────────────────────────────────────────────────
    print()
    print(f"Done. {n - len(failed)}/{n} communities completed successfully.")
    if failed:
        print(f"\nFailed runs ({len(failed)}):")
        for comm, scenario, err in failed:
            print(f"  {comm} [{scenario}]: {err}")


if __name__ == '__main__':
    main()
