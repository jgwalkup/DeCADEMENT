"""
run_multi_climate_simulations.py
=================================
Run DeCA-DEMENT simulations for a set of pre-assembled communities across one
or more climate scenarios. This script supersedes run_community_simulations.py
and can reproduce all of its behaviour while also supporting an arbitrary number
of climate treatments.

Climate source (--climate)
--------------------------
  Omitted          → climate.csv from --input-dir is used (warning printed)
  A .csv file      → single scenario; output written flat into --output-dir
  A directory      → every .csv inside is a scenario; each is validated before
                     any simulations run; output written to
                     <output-dir>/<stem>/<community_outname>/

Static overrides (--static-temp / --static-psi)
------------------------------------------------
  Mutually exclusive with a climate directory.
  Accept one or more values (must be equal in count):
    --static-temp 10 15 20  --static-psi -0.5 -0.5 -0.3
  Each value-pair generates a synthetic climate CSV by replacing Temp and Psi
  in the FM-column template (climate.csv from --input-dir) with the constant
  values supplied.  Output subdirectory stem: static_T{temp}_P{psi}.
  When used with a single --climate CSV the FM template is that CSV instead of
  the default.  When --climate is omitted or is a directory, the FM template is
  always climate.csv from --input-dir.

Output directory structure
--------------------------
  Single CSV, no overrides:
    <output-dir>/<community_name>_s<seed>/

  Multiple scenarios (directory or static):
    <output-dir>/<scenario_stem>/<community_name>_s<seed>/

Expected project layout
-----------------------
    Decoupled-Energy-DEMENTpy/   ← run from here
    ├── src/                     ← --src-dir
    ├── input/                   ← --input-dir  (runtime.txt, CSVs)
    ├── community_pool_03/
    │   └── communities/         ← --community-dir
    │       ├── community_0001/
    │       └── ...
    └── analysis/
        └── run_multi_climate_simulations.py

Usage examples
--------------
  # Single climate (same as run_community_simulations.py, single-scenario mode)
  python analysis/run_multi_climate_simulations.py \\
      --src-dir src --input-dir input \\
      --community-dir community_pool_03/communities \\
      --output-dir simulations \\
      --climate climate.csv --seed 42

  # Climate directory (arbitrary number of treatments)
  python analysis/run_multi_climate_simulations.py \\
      --src-dir src --input-dir input \\
      --community-dir community_pool_03/communities \\
      --output-dir simulations \\
      --climate climate_scenarios/ --seed 42

  # Static override — three temperature treatments, shared psi
  python analysis/run_multi_climate_simulations.py \\
      --src-dir src --input-dir input \\
      --community-dir community_pool_03/communities \\
      --output-dir simulations \\
      --static-temp 10 15 20 --static-psi -0.5 -0.5 -0.5 \\
      --seed 42

Options
-------
  --src-dir         Directory containing DEMENTpy .py source files (default: src)
  --input-dir       Directory containing runtime.txt and input CSVs (default: input)
  --community-dir   Directory of community subdirectories to simulate (required)
  --output-dir      Base directory for simulation outputs (required)
  --climate         Climate CSV filename or directory of climate CSVs inside
                    --input-dir.  Defaults to climate.csv from --input-dir
                    (a warning is printed when omitted).
  --static-temp     One or more constant temperature values (°C) for synthetic
                    climate generation.  Must be paired with --static-psi.
  --static-psi      One or more constant water-potential values to pair with
                    --static-temp.  Count must match --static-temp.
  --seed            Base random seed (default: 0)
  --communities     Comma-separated community folder names to run (default: all)
  --start-from      Resume from this community name (inclusive, sorted order)
  --range           1-based inclusive range: --range START END
  --quiet           Suppress per-day progress output
"""

import argparse
import os
import sys
import tempfile
import time
import warnings
import traceback
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore", category=UserWarning, module="zarr")


# ---------------------------------------------------------------------------
# Constants — expected climate CSV schema
# ---------------------------------------------------------------------------

CLIMATE_REQUIRED_COLS = {'Temp', 'FM', 'Psi', 'Date'}


# ---------------------------------------------------------------------------
# Runtime loading  (mirrors dementpy.py / run_community_simulations.py)
# ---------------------------------------------------------------------------

def load_runtime(input_dir, runtime_file='runtime.txt'):
    path = os.path.join(input_dir, runtime_file)
    runtime = pd.read_csv(path, header=None, index_col=0, sep='\t')
    output_format = str(runtime.iloc[runtime.index.get_loc('output_format'), 0]).strip().lower()
    runtime = runtime.drop('output_format')
    runtime = runtime.astype(float)
    runtime.columns = [1]
    return runtime, output_format


# ---------------------------------------------------------------------------
# Climate file validation
# ---------------------------------------------------------------------------

def validate_climate_file(path, reference_nrows, label=None):
    """
    Validate a single climate CSV against the expected schema.

    Checks:
      - Required columns present: Temp, FM, Psi, Date
      - Same number of data rows as the reference climate file
      - Temp, FM, Psi are numeric
      - Date column is non-empty strings

    Parameters
    ----------
    path          : str or Path — file to validate
    reference_nrows : int — expected number of data rows
    label         : str — display name for error messages (defaults to path)

    Raises
    ------
    ValueError with a descriptive message if any check fails.
    """
    label = label or str(path)
    try:
        df = pd.read_csv(path, header=0, index_col=0)
    except Exception as e:
        raise ValueError(f"{label}: could not be read as CSV — {e}")

    missing = CLIMATE_REQUIRED_COLS - set(df.columns)
    if missing:
        raise ValueError(
            f"{label}: missing required column(s): {sorted(missing)}. "
            f"Expected columns: {sorted(CLIMATE_REQUIRED_COLS)}"
        )

    if len(df) != reference_nrows:
        raise ValueError(
            f"{label}: has {len(df)} data rows but reference climate has "
            f"{reference_nrows}. All climate files must have the same length."
        )

    for col in ('Temp', 'FM', 'Psi'):
        if not pd.api.types.is_numeric_dtype(df[col]):
            raise ValueError(
                f"{label}: column '{col}' is not numeric "
                f"(dtype={df[col].dtype})."
            )

    if df['Date'].isnull().any():
        raise ValueError(f"{label}: column 'Date' contains null values.")

    return df   # return parsed df so callers can reuse it if needed


def load_reference_climate(input_dir):
    """Load climate.csv from input_dir as the reference / FM template."""
    ref_path = os.path.join(input_dir, 'climate.csv')
    if not os.path.isfile(ref_path):
        raise FileNotFoundError(
            f"Reference climate file not found: {ref_path}\n"
            "  climate.csv must be present in --input-dir."
        )
    df = pd.read_csv(ref_path, header=0, index_col=0)
    return ref_path, df


# ---------------------------------------------------------------------------
# Scenario building
# ---------------------------------------------------------------------------

def build_scenarios_from_file(climate_path):
    """Single CSV → one scenario, output written flat (no subdirectory)."""
    return [
        {
            'name':   None,         # None = flat output, no subdirectory
            'stem':   None,
            'path':   climate_path,
            'tmpfile': None,
        }
    ]


def build_scenarios_from_directory(climate_dir, input_dir):
    """
    Scan climate_dir for .csv files, validate each against the reference
    climate, and return a list of scenario dicts sorted alphabetically by stem.

    Validation is performed in full before any dicts are returned — a single
    bad file aborts the whole run with a clear message.
    """
    climate_dir = Path(climate_dir)
    csv_files = sorted(climate_dir.glob('*.csv'), key=lambda p: p.stem)

    if not csv_files:
        raise ValueError(
            f"No .csv files found in climate directory: {climate_dir}"
        )

    # Load reference for row-count check
    ref_path, ref_df = load_reference_climate(input_dir)
    reference_nrows = len(ref_df)

    print(f"Validating {len(csv_files)} climate file(s) in {climate_dir} …")
    errors = []
    for f in csv_files:
        try:
            validate_climate_file(f, reference_nrows, label=f.name)
            print(f"  ✓  {f.name}")
        except ValueError as e:
            errors.append(str(e))
            print(f"  ✗  {e}")

    if errors:
        print(f"\n{len(errors)} climate file(s) failed validation. Aborting.")
        sys.exit(1)

    print()

    return [
        {
            'name':    f.stem,
            'stem':    f.stem,
            'path':    str(f),
            'tmpfile': None,
        }
        for f in csv_files
    ]


def build_scenarios_from_static(static_temps, static_psis, input_dir,
                                 fm_template_path=None):
    """
    Generate one synthetic climate CSV per (temp, psi) pair.

    The FM column is taken from fm_template_path if provided, otherwise from
    climate.csv in input_dir.  Temp and Psi columns are replaced with the
    supplied constant values.  Synthetic files are written to a temporary
    directory that the caller is responsible for cleaning up.

    Returns
    -------
    scenarios : list of dicts with keys name/stem/path/tmpfile
    tmpdir    : str — path to the temporary directory (caller must delete)
    """
    if fm_template_path:
        template_df = pd.read_csv(fm_template_path, header=0, index_col=0)
    else:
        _, template_df = load_reference_climate(input_dir)

    tmpdir = tempfile.mkdtemp(prefix='dement_static_climate_')
    scenarios = []

    for temp, psi in zip(static_temps, static_psis):
        # Format stem: replace '.' with 'p' and '-' with 'n' for filesystem safety
        def _fmt(v):
            s = f"{v:g}"
            s = s.replace('-', 'n').replace('.', 'p')
            return s

        stem = f"static_T{_fmt(temp)}_P{_fmt(psi)}"
        out_path = os.path.join(tmpdir, f"{stem}.csv")

        synthetic = template_df.copy()
        synthetic['Temp'] = float(temp)
        synthetic['Psi']  = float(psi)
        synthetic.to_csv(out_path)

        scenarios.append({
            'name':    stem,
            'stem':    stem,
            'path':    out_path,
            'tmpfile': out_path,
        })

    return scenarios, tmpdir


# ---------------------------------------------------------------------------
# Single simulation run  (identical to run_community_simulations.py)
# ---------------------------------------------------------------------------

def run_simulation(runtime, output_format, community_dir, climate_path,
                   output_path, outname, src_dir, input_dir, quiet=False):
    """
    Run one complete DEMENTpy simulation for a single community.

    Parameters
    ----------
    runtime        : DataFrame — numeric runtime parameters
    output_format  : str — 'csv', 'netcdf', or 'zarr'
    community_dir  : Path — absolute path to the community directory
    climate_path   : str — absolute path to the climate CSV
    output_path    : Path — absolute path to the output directory
    outname        : str — output file label (used as filename/directory prefix)
    src_dir        : str — absolute path to DEMENTpy source directory
    quiet          : bool — suppress per-day print statements
    """
    sys.path.insert(0, src_dir)
    from initialization import initialize_data
    from grid import Grid
    from output import Output

    col = runtime.columns[0]

    pulse      = int(runtime.loc['pulse',    col])
    cycle      = int(runtime.loc['end_time', col])
    interval   = int(runtime.loc['interval', col])
    mic_reinit = True

    project_root = os.getcwd()
    # chdir to the model input directory (where parameters.csv lives),
    # NOT os.path.dirname(climate_path) — which breaks when climate files
    # are in a subdirectory of input_dir (e.g. input/climate_data/).
    os.chdir(input_dir)
    try:
        data = initialize_data(runtime, climate_file=climate_path,
                               community_dir=str(community_dir))
    finally:
        os.chdir(project_root)

    Output_init = Output(runtime, data)
    Ecosystem   = Grid(runtime, data)

    # ── Timing accumulators ───────────────────────────────────────────────────
    PHASES = ["degradation", "uptake", "metabolism", "mortality",
              "reproduction", "output", "microbes_abundance", "repopulation"]
    t_total   = {ph: 0.0 for ph in PHASES}
    run_start = time.perf_counter()

    def _fmt_time(secs):
        m, s = divmod(secs, 60)
        return f"{int(m):3d}m {s:05.2f}s"

    def _print_timing(label, timings, n_days):
        phase_total = sum(timings.values())
        print(f"\n  ── Timing: {label} ({n_days} days) ─────────────────────────────")
        for ph in PHASES:
            t = timings[ph]
            if t < 0.001:
                continue
            pct = 100.0 * t / phase_total if phase_total > 0 else 0.0
            print(f"    {ph:22s}  {_fmt_time(t)}  ({pct:5.1f}%)")
        print(f"    {'TOTAL':22s}  {_fmt_time(phase_total)}")
        print(f"    {'Per-day average':22s}  {phase_total / n_days * 1000:.1f} ms/day",
              flush=True)

    for p in range(pulse):
        t_pulse = {ph: 0.0 for ph in PHASES}

        for i in range(p * cycle, (p + 1) * cycle):
            if not quiet:
                print(f"Day {i} starting", flush=True)

            t0 = time.perf_counter()
            Ecosystem.degradation(p, i)
            t_pulse["degradation"] += time.perf_counter() - t0

            t0 = time.perf_counter()
            Ecosystem.uptake(p, i)
            t_pulse["uptake"] += time.perf_counter() - t0

            t0 = time.perf_counter()
            try:
                Ecosystem.metabolism(i)
            except Exception as e:
                print(f"CRASH in metabolism() on day {i}:", flush=True)
                traceback.print_exc()
                raise
            t_pulse["metabolism"] += time.perf_counter() - t0

            t0 = time.perf_counter()
            Ecosystem.mortality(i)
            t_pulse["mortality"] += time.perf_counter() - t0

            t0 = time.perf_counter()
            Ecosystem.reproduction(i)
            t_pulse["reproduction"] += time.perf_counter() - t0

            t0 = time.perf_counter()
            if i == 0:
                Output_init.output(Ecosystem, i)
            elif i % interval == interval - 1:
                Output_init.output(Ecosystem, i)
            t_pulse["output"] += time.perf_counter() - t0

            t0 = time.perf_counter()
            Output_init.microbes_abundance(Ecosystem, i)
            t_pulse["microbes_abundance"] += time.perf_counter() - t0

            t0 = time.perf_counter()
            if i == (p + 1) * cycle - 1:
                # Flush Taxon_count_repop buffer before repopulation() reads it
                Output_init.flush_pulse_buffers()
                Ecosystem.repopulation(Output_init, i, mic_reinit)
            t_pulse["repopulation"] += time.perf_counter() - t0

        for ph in PHASES:
            t_total[ph] += t_pulse[ph]
        _print_timing(f"Pulse {p + 1}/{pulse}", t_pulse, cycle)

    # ── Full-run timing summary ───────────────────────────────────────────────
    _print_timing(f"Full run ({pulse} pulse(s))", t_total, pulse * cycle)
    print(f"    {'Wall time (incl. init)':22s}  {_fmt_time(time.perf_counter() - run_start)}",
          flush=True)

    # Save output
    t0 = time.perf_counter()
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
            warnings.warn(
                f"Unknown output_format '{output_format}', defaulting to NetCDF."
            )
            Output_init.export_to_netcdf(outname)
    finally:
        os.chdir(project_root)
    print(f"    {'Export':22s}  {_fmt_time(time.perf_counter() - t0)}", flush=True)


# ---------------------------------------------------------------------------
# Community directory detection  (same as run_community_simulations.py)
# ---------------------------------------------------------------------------

def get_community_dirs(community_dir, selection=None, start_from=None):
    """
    Return sorted list of community subdirectories (those containing
    EnzGenes.csv).  Optionally filter to a specific selection or resume
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
        description=(
            "Run DeCA-DEMENT community simulations across one or more climate "
            "scenarios."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # ── Paths ─────────────────────────────────────────────────────────────────
    parser.add_argument('--src-dir', type=str, default='src',
                        help='Directory containing DEMENTpy .py source files')
    parser.add_argument('--input-dir', type=str, default='input',
                        help='Directory containing runtime.txt and input CSVs')
    parser.add_argument('--community-dir', type=str, required=True,
                        help='Directory of community subdirectories to simulate')
    parser.add_argument('--output-dir', type=str, required=True,
                        help='Base directory for simulation outputs')

    # ── Climate source ────────────────────────────────────────────────────────
    parser.add_argument('--climate', type=str, default=None,
                        help=(
                            'Climate CSV filename or directory of climate CSVs '
                            'inside --input-dir.  If omitted, climate.csv from '
                            '--input-dir is used and a warning is printed.'
                        ))

    # ── Static overrides ──────────────────────────────────────────────────────
    parser.add_argument('--static-temp', type=float, nargs='+', default=None,
                        metavar='TEMP',
                        help=(
                            'One or more constant temperature values (°C) for '
                            'synthetic climate generation.  Must be paired with '
                            '--static-psi and the same number of values.'
                        ))
    parser.add_argument('--static-psi', type=float, nargs='+', default=None,
                        metavar='PSI',
                        help=(
                            'One or more constant water-potential values to pair '
                            'with --static-temp.  Count must equal --static-temp.'
                        ))

    # ── Community selection ───────────────────────────────────────────────────
    parser.add_argument('--seed', type=int, default=0,
                        help='Base random seed; each community gets seed + index')
    parser.add_argument('--communities', type=str, default=None,
                        help='Comma-separated community folder names to run '
                             '(default: all)')
    parser.add_argument('--taxon-list', type=str, default=None,
                        metavar='CSV',
                        help='Path to taxon_summary.csv produced by '
                             'extract_community_taxa.py.  Selects the '
                             'single-taxon community directory for every '
                             'taxon_id in the file.  --community-dir should '
                             'point to the single_taxon/ directory.  '
                             'Cannot be combined with --communities or --range.')
    parser.add_argument('--start-from', type=str, default=None,
                        help='Resume from this community name (inclusive, sorted '
                             'order).  Cannot be combined with --range.')
    parser.add_argument('--range', type=int, nargs=2, metavar=('START', 'END'),
                        help='1-based inclusive range of communities by sorted '
                             'position, e.g. --range 1 50.  Cannot be combined '
                             'with --communities or --start-from.')
    parser.add_argument('--quiet', action='store_true',
                        help='Suppress per-day progress output')

    args = parser.parse_args()

    # ── Resolve absolute paths ────────────────────────────────────────────────
    src_dir       = os.path.abspath(args.src_dir)
    input_dir     = os.path.abspath(args.input_dir)
    community_dir = os.path.abspath(args.community_dir)
    output_dir    = Path(os.path.abspath(args.output_dir))

    for label, path in [('--src-dir', src_dir),
                         ('--input-dir', input_dir),
                         ('--community-dir', community_dir)]:
        if not os.path.isdir(path):
            print(f"ERROR: {label} '{path}' does not exist.")
            sys.exit(1)

    # ── Validate static override arguments ───────────────────────────────────
    using_static = args.static_temp is not None or args.static_psi is not None

    if using_static:
        if args.static_temp is None or args.static_psi is None:
            print("ERROR: --static-temp and --static-psi must be supplied together.")
            sys.exit(1)
        if len(args.static_temp) != len(args.static_psi):
            print(
                f"ERROR: --static-temp has {len(args.static_temp)} value(s) but "
                f"--static-psi has {len(args.static_psi)}.  Counts must match."
            )
            sys.exit(1)

    # ── Resolve climate source ────────────────────────────────────────────────
    # Determine climate_source_path (absolute) and whether it is a directory.
    if args.climate is None:
        default_climate = os.path.join(input_dir, 'climate.csv')
        if not os.path.isfile(default_climate):
            print(
                f"ERROR: --climate not specified and default climate.csv not "
                f"found in --input-dir ({input_dir})."
            )
            sys.exit(1)
        print(
            f"WARNING: --climate not specified. Defaulting to climate.csv "
            f"from --input-dir ({input_dir})."
        )
        climate_source_path = default_climate
        climate_is_dir = False
    else:
        candidate = os.path.join(input_dir, args.climate)
        if os.path.isdir(candidate):
            climate_source_path = candidate
            climate_is_dir = True
        elif os.path.isfile(candidate):
            climate_source_path = candidate
            climate_is_dir = False
        else:
            # Also try as a bare absolute/relative path outside input_dir
            candidate2 = os.path.abspath(args.climate)
            if os.path.isdir(candidate2):
                climate_source_path = candidate2
                climate_is_dir = True
            elif os.path.isfile(candidate2):
                climate_source_path = candidate2
                climate_is_dir = False
            else:
                print(
                    f"ERROR: --climate '{args.climate}' not found as a file or "
                    f"directory (looked in --input-dir and as an absolute/relative path)."
                )
                sys.exit(1)

    # Static overrides and a climate directory are mutually exclusive.
    if using_static and climate_is_dir:
        print(
            "ERROR: --static-temp/--static-psi cannot be combined with a climate "
            "directory.  Use a single climate CSV as the FM template instead, or "
            "omit --climate to use climate.csv from --input-dir."
        )
        sys.exit(1)

    # ── Build scenario list ───────────────────────────────────────────────────
    tmpdir = None   # set if static scenarios allocate a temp directory

    if using_static:
        fm_template = (climate_source_path
                       if not climate_is_dir else None)  # None → uses default
        scenarios, tmpdir = build_scenarios_from_static(
            args.static_temp, args.static_psi, input_dir,
            fm_template_path=fm_template,
        )
    elif climate_is_dir:
        scenarios = build_scenarios_from_directory(climate_source_path, input_dir)
    else:
        scenarios = build_scenarios_from_file(climate_source_path)

    multi_scenario = len(scenarios) > 1 or scenarios[0]['name'] is not None

    # ── Validate community selection args ─────────────────────────────────────
    if args.taxon_list and (args.range or args.communities or args.start_from):
        print("ERROR: --taxon-list cannot be combined with --range, "
              "--communities, or --start-from.")
        sys.exit(1)
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
    if args.taxon_list:
        # Read taxon_summary.csv and convert taxon IDs → folder names
        taxon_csv = Path(os.path.abspath(args.taxon_list))
        if not taxon_csv.exists():
            print(f"ERROR: --taxon-list file not found: {taxon_csv}")
            sys.exit(1)
        try:
            taxon_df  = pd.read_csv(taxon_csv)
            taxon_ids = taxon_df['taxon_id'].dropna().unique().tolist()
        except Exception as e:
            print(f"ERROR: could not read --taxon-list '{taxon_csv}': {e}")
            sys.exit(1)
        # Map Tax47 → taxon_Tax47 (folder name convention from assemble_communities.py)
        selection = [f"taxon_{tid}" for tid in taxon_ids]
        print(f"Taxon list    : {len(selection)} taxa from {taxon_csv.name}")
    else:
        selection = ([s.strip() for s in args.communities.split(',')]
                     if args.communities else None)

    comm_dirs = get_community_dirs(community_dir,
                                   selection=selection,
                                   start_from=args.start_from)

    if args.range:
        start, end = args.range
        comm_dirs = comm_dirs[start - 1 : end]
        if not comm_dirs:
            print(
                f"ERROR: --range {start} {end} produced no communities "
                f"(pool has {len(get_community_dirs(community_dir))} communities)."
            )
            sys.exit(1)

    if not comm_dirs:
        print("ERROR: no community directories found.")
        sys.exit(1)

    # ── Load runtime — shared across all communities and scenarios ────────────
    runtime, output_format = load_runtime(input_dir)

    # ── Startup summary ───────────────────────────────────────────────────────
    n = len(comm_dirs)
    scenario_names = [s['stem'] or '(single)' for s in scenarios]
    print(f"Communities : {n}")
    print(f"Scenarios   : {len(scenarios)}  {scenario_names}")
    print(f"Base seed   : {args.seed}")
    print(f"Output dir  : {output_dir}")
    print()

    # ── Main loop ─────────────────────────────────────────────────────────────
    failed = []

    try:
        for idx, cdir in enumerate(comm_dirs):
            comm_name  = cdir.name
            comm_seed  = (args.seed + idx) % (2**32)
            outname    = f"{comm_name}_s{comm_seed}"
            run_label  = f"[{idx+1}/{n}] {comm_name}  (seed={comm_seed})"

            for scenario in scenarios:
                scen_stem  = scenario['stem']
                clim_path  = scenario['path']

                # Reseed for every scenario, not once per community: all
                # randomness (placement, mortality, dispersal, ReqEnz) uses
                # numpy's global generator, so each scenario of a community
                # then starts from the same random state and the runs stay
                # identical until their climate inputs differ (paired design).
                np.random.seed(comm_seed)

                # Output path: flat when single unnamed scenario, else subdir.
                if scen_stem is None:
                    out_path = output_dir
                else:
                    out_path = output_dir / scen_stem

                label = (f"{run_label}  [{scen_stem}]"
                         if multi_scenario else run_label)
                print(f"{label}", flush=True)

                try:
                    run_simulation(
                        runtime       = runtime.copy(),
                        output_format = output_format,
                        community_dir = cdir,
                        climate_path  = clim_path,
                        output_path   = out_path,
                        outname       = outname,
                        src_dir       = src_dir,
                        input_dir     = input_dir,
                        quiet         = args.quiet,
                    )
                    print(f"  -> saved to {out_path / outname}", flush=True)

                except Exception as e:
                    print(f"  !! FAILED: {e}", flush=True)
                    traceback.print_exc()
                    failed.append((comm_name, scen_stem or 'run', str(e)))
                    continue

    finally:
        # Always clean up synthetic climate temp files
        if tmpdir is not None:
            import shutil
            shutil.rmtree(tmpdir, ignore_errors=True)

    # ── Summary ───────────────────────────────────────────────────────────────
    total_runs = n * len(scenarios)
    n_failed   = len(failed)
    print()
    print(
        f"Done. {total_runs - n_failed}/{total_runs} run(s) completed successfully "
        f"({n} communities × {len(scenarios)} scenario(s))."
    )
    if failed:
        print(f"\nFailed runs ({n_failed}):")
        for comm, scenario, err in failed:
            print(f"  {comm} [{scenario}]: {err}")


if __name__ == '__main__':
    main()
