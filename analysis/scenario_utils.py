"""
scenario_utils.py
=================
Shared helpers used by all DEMENTpy analysis scripts for multi-scenario output
produced by run_multi_climate_simulations.py.

Output layout from run_multi_climate_simulations.py (multi-scenario):
    <output-dir>/<scenario_stem>/<community_name>_s<seed>/

Single-scenario (flat) layout:
    <output-dir>/<community_name>_s<seed>/

Public API
----------
load_series(run_dir, name)      – format-agnostic series loader (CSV/NC/Zarr)
pick_runs(scenario_dir, ...)    – select one run-dir per replicate from a scenario dir
resolve_runs(names, base_dir)   – resolve explicit folder names to Paths
discover_scenarios(base)        – auto-detect scenario subdirs under base
scenario_color_map(labels)      – assign a colour to each scenario label
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from pathlib import Path


# ── Format detection ──────────────────────────────────────────────────────────

def _sentinel_csv(d: Path) -> Path | None:
    for name in ("SubstratesSeries.csv", "MicrobesSeries.csv",
                 "C_Uptake_Taxon.csv", "CUE_system.csv"):
        f = d / name
        if f.exists():
            return f
    return None


def _sentinel_nc(d: Path) -> Path | None:
    for name in ("SubstratesSeries.nc", "MicrobesSeries.nc",
                 "C_Uptake_Taxon.nc", "CUE_system.nc"):
        f = d / name
        if f.exists():
            return f
    return None


def _is_csv_run(d: Path) -> bool:
    return _sentinel_csv(d) is not None


def _is_nc_run(d: Path) -> bool:
    return _sentinel_nc(d) is not None


def _is_zarr_run(d: Path) -> bool:
    if d.suffix == ".zarr" and d.is_dir():
        return True
    return any(p for p in d.glob("*.zarr") if not p.name.endswith(".init.zarr"))


def _get_zarr_store(d: Path) -> str | None:
    if d.suffix == ".zarr" and d.is_dir():
        return str(d)
    stores = [p for p in d.glob("*.zarr") if not p.name.endswith(".init.zarr")]
    return str(stores[0]) if stores else None


def is_run_dir(d: Path) -> bool:
    return _is_csv_run(d) or _is_nc_run(d) or _is_zarr_run(d)


# ── Per-format loaders ────────────────────────────────────────────────────────

def _normalise_cols(df: pd.DataFrame) -> pd.DataFrame:
    """Normalise day columns to int where possible."""
    df.columns = [int(c) if str(c).lstrip("-").isdigit() else c for c in df.columns]
    return df


def _load_csv(run_dir: Path, name: str) -> pd.DataFrame | None:
    p = run_dir / f"{name}.csv"
    if not p.exists():
        return None
    df = pd.read_csv(p, index_col=0)
    # Positional renaming: handles both correct (int) and broken ('C','C.1',...) column names
    df.columns = range(len(df.columns))        # 0 = init state, 1..n = days
    df = df.drop(columns=[0], errors="ignore") # drop init state column
    df.columns = range(1, len(df.columns) + 1) # renumber 1..n_days
    return df


def _load_nc(run_dir: Path, name: str) -> pd.DataFrame | None:
    try:
        import xarray as xr
    except ImportError:
        raise ImportError("xarray is required for NetCDF output. pip install xarray")
    p = run_dir / f"{name}.nc"
    if not p.exists():
        return None
    try:
        ds = xr.open_dataset(p)
        df = ds.to_dataframe()
    except Exception:
        da = xr.open_dataarray(p)
        df = da.to_series().to_frame()
    return _normalise_cols(df)


def _load_zarr(run_dir: Path, name: str) -> pd.DataFrame | None:
    try:
        import xarray as xr
    except ImportError:
        raise ImportError("xarray/zarr are required for Zarr output. pip install xarray zarr")
    store = _get_zarr_store(run_dir)
    if store is None:
        return None
    try:
        ds = xr.open_zarr(store, group=name)
        df = ds.to_dataframe()
    except Exception:
        return None
    return _normalise_cols(df)


def load_series(run_dir: Path, name: str) -> pd.DataFrame | None:
    """Load a named series from whichever output format is present (CSV > NC > Zarr).

    Columns are normalised: day columns become ints, the init-state column (0 or
    'C') is dropped, and columns are renumbered 1..N for CSV runs.  For NC/Zarr
    runs the raw integer column labels from the file are preserved.
    """
    if _is_csv_run(run_dir):
        return _load_csv(run_dir, name)
    if _is_nc_run(run_dir):
        return _load_nc(run_dir, name)
    return _load_zarr(run_dir, name)


# ── Run-directory selection ───────────────────────────────────────────────────

def pick_runs(scenario_dir: Path, min_cols: int = 100,
              sentinel: str | None = None) -> list[Path]:
    """Select the best run directory per replicate from a scenario directory.

    A run directory is valid if:
      - it is a directory
      - it contains any recognised output sentinel file
      - the sentinel has at least *min_cols* columns

    When multiple directories share the same replicate id (the second
    underscore-delimited token of the folder name), the lexicographically
    latest one wins.

    Parameters
    ----------
    scenario_dir : Path
        Parent directory containing <community_name>_s<seed>/ subdirs.
    min_cols : int
        Minimum number of day-columns a sentinel file must have.
    sentinel : str or None
        Preferred sentinel filename stem (no extension) to check for
        column count.  Falls back to any recognised format sentinel.
    """
    if not scenario_dir.exists():
        return []

    by_rep: dict[str, Path] = {}

    for d in sorted(scenario_dir.iterdir()):
        if not d.is_dir():
            continue
        if not is_run_dir(d):
            continue

        # Column-count filter ─────────────────────────────────────────────────
        nc = _count_cols(d, sentinel)
        if nc is not None and nc < min_cols:
            continue

        parts = d.name.split("_")
        rep = parts[1] if len(parts) > 1 else "000"
        if rep not in by_rep or d.name > by_rep[rep].name:
            by_rep[rep] = d

    return list(by_rep.values())


def _count_cols(run_dir: Path, preferred: str | None) -> int | None:
    """Return column count for the preferred sentinel (or any available one)."""
    candidates: list[str] = []
    if preferred:
        candidates.append(preferred)
    candidates += ["SubstratesSeries", "MicrobesSeries",
                   "C_Uptake_Taxon", "CUE_system"]

    for stem in candidates:
        # CSV
        f = run_dir / f"{stem}.csv"
        if f.exists():
            try:
                return pd.read_csv(f, index_col=0, nrows=1).shape[1]
            except Exception:
                pass
        # NC
        f = run_dir / f"{stem}.nc"
        if f.exists():
            try:
                import xarray as xr
                ds = xr.open_dataset(f)
                return len(ds.data_vars)
            except Exception:
                pass
    return None


def resolve_runs(names: list[str], base_dir: Path) -> list[Path]:
    """Resolve a list of folder names to Paths, relative to base_dir or absolute."""
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


# ── Scenario discovery ────────────────────────────────────────────────────────

def discover_scenarios(base: Path) -> list[tuple[str, Path]]:
    """Auto-detect scenario subdirectories under *base*.

    A subdirectory is counted as a scenario if it contains at least one
    recognised run directory (i.e. a dir that has CSV/NC/Zarr output).
    Returns a list of (label, path) tuples sorted alphabetically by label.

    Falls back to treating *base* itself as a single unnamed scenario when
    no qualifying subdirectories are found (flat single-scenario layout).
    """
    if not base.exists():
        return []

    scenarios: list[tuple[str, Path]] = []
    for d in sorted(base.iterdir()):
        if not d.is_dir():
            continue
        # Does this subdir contain any run dirs?
        has_runs = any(is_run_dir(sub) for sub in d.iterdir() if sub.is_dir())
        if has_runs:
            scenarios.append((d.name, d))

    if not scenarios:
        # Flat layout: base itself contains run dirs
        if any(is_run_dir(d) for d in base.iterdir() if d.is_dir()):
            scenarios = [(base.name, base)]

    return scenarios


# ── Summary-CSV loading (output of export_summary_data.py) ────────────────────
#
# export_summary_data.py writes, per scenario, a directory of paired files:
#     <stem>_mean.csv   and   <stem>_sd.csv
# A 1-D system series is written with a single column "mean"/"sd" and day index.
# A 2-D quantity (taxa×days or type×days) has the taxon/type as index and day
# numbers as columns.  These loaders return the same shapes the figure scripts
# get from their own per-run aggregation, so a script can swap raw-run loading
# for summary loading without changing plotting code.

def _read_summary_file(path: Path, value_col: str):
    """Read one <stem>_{mean,sd}.csv → Series (1-D) or DataFrame (2-D)."""
    if not path.exists():
        return None
    df = pd.read_csv(path, index_col=0)
    if list(df.columns) == [value_col]:
        s = df[value_col]
        s.index = [int(i) if str(i).lstrip("-").isdigit() else i for i in s.index]
        return s
    return _normalise_cols(df)


def load_summary(scenario_dir: Path, stem: str):
    """Load a (mean, sd) pair for one variable from a summary scenario dir.

    Returns (mean, sd); each is a Series (1-D) or DataFrame (taxa/type × days),
    or None if the file is absent.
    """
    scenario_dir = Path(scenario_dir)
    mean = _read_summary_file(scenario_dir / f"{stem}_mean.csv", "mean")
    sd   = _read_summary_file(scenario_dir / f"{stem}_sd.csv",   "sd")
    return mean, sd


def has_summary(scenario_dir: Path) -> bool:
    """True if *scenario_dir* looks like an export_summary_data.py scenario dir."""
    scenario_dir = Path(scenario_dir)
    return scenario_dir.is_dir() and any(scenario_dir.glob("*_mean.csv"))


def _natural_key(name: str):
    """Sort key that orders embedded numbers numerically.

    'static_T3_Pn0p5' < 'static_T15_Pn0p5' < 'static_T27_Pn0p5' (not alphabetical),
    which is what a temperature gradient expects.  Named scenarios with no digits
    fall back to plain alphabetical order.
    """
    import re as _re
    return [int(t) if t.isdigit() else t.lower()
            for t in _re.split(r"(\d+)", str(name))]


def discover_summary_scenarios(root: Path) -> list[tuple[str, Path]]:
    """Auto-detect summary scenario subdirectories under *root*.

    A subdirectory qualifies if it contains at least one <stem>_mean.csv.
    Returns [(label, path)] in natural-sorted order (numeric-aware, so a
    temperature gradient comes out T3, T9, T15, ...).  Falls back to treating
    *root* itself as a single scenario when it directly contains summary files.
    """
    root = Path(root)
    if not root.exists():
        return []
    found = [(d.name, d) for d in root.iterdir() if d.is_dir() and has_summary(d)]
    found.sort(key=lambda lp: _natural_key(lp[0]))
    if not found and has_summary(root):
        found = [(root.name, root)]
    return found


def resolve_summary_scenarios(summary_dir,
                              treatment_map: dict | None = None,
                              label_order: list[str] | None = None
                              ) -> list[tuple[str, Path]]:
    """Return [(label, scenario_summary_dir)] for a summary root.

    treatment_map : optional {dir_name: display_label} to rename discovered dirs.
    label_order   : optional desired left-to-right ordering of display labels;
                    scenarios not listed are appended in discovered order.
    """
    found = discover_summary_scenarios(Path(summary_dir))
    tmap = treatment_map or {}
    relabelled = [(tmap.get(name, name), path) for name, path in found]
    if label_order:
        rank = {lab: i for i, lab in enumerate(label_order)}
        relabelled.sort(key=lambda lp: rank.get(lp[0], len(rank)))
    return relabelled


def add_summary_arg(parser) -> None:
    """Add the --summary-dir flag that switches a figure script to summary mode."""
    parser.add_argument(
        "--summary-dir", default=None, metavar="DIR",
        help="Read pre-computed summary CSVs (from export_summary_data.py) under "
             "this root instead of raw runs. Scenario subdirectories are "
             "auto-discovered; all are plotted (full gradient or N-scenario set).",
    )


def read_summary_meta(scenario_dir: Path) -> dict:
    """Read summary_meta.csv from a summary scenario dir → metadata dict.

    Returns {'gridsize', 'pulse', 'end_time', 'col_start', 'col_end'} with sane
    defaults (gridsize=1, pulse=1, end_time=365) when the file is absent — so
    callers can rely on the keys existing.
    """
    defaults = {"gridsize": 1, "pulse": 1, "end_time": 365,
                "col_start": None, "col_end": None}
    f = Path(scenario_dir) / "summary_meta.csv"
    if not f.exists():
        return defaults
    try:
        s = pd.read_csv(f, index_col=0).iloc[:, 0]
        out = dict(defaults)
        for k in defaults:
            if k in s.index and pd.notna(s.loc[k]):
                out[k] = int(float(s.loc[k]))
        return out
    except Exception:
        return defaults


# ── Colour palette ────────────────────────────────────────────────────────────

# Fallback discrete palette for non-temperature contexts (2–8 scenarios).
_PALETTE = [
    "#1f77b4",  # blue
    "#d62728",  # red
    "#2ca02c",  # green
    "#ff7f0e",  # orange
    "#9467bd",  # purple
    "#8c564b",  # brown
    "#e377c2",  # pink
    "#17becf",  # cyan
]

# Temperature gradient: blue (cool) → yellow (mid) → red (hot).
# Uses matplotlib's 'coolwarm' colormap which is perceptually uniform,
# diverging, and works well for both screen and print.
# Custom blue → orange → red gradient for temperature scenarios.
# Blue = cool/ambient, orange = moderate warming, red = extreme warming.
# No white midpoint; all colours are vivid on white backgrounds.
# Registered as 'temperature' so --palette temperature uses it automatically.
from matplotlib.colors import LinearSegmentedColormap as _LSC
_TEMP_CMAP_OBJ = _LSC.from_list(
    "blue_orange_red",
    ["#2166ac", "#f4a836", "#d6172a"],  # blue, orange, red
    N=256,
)
# Register with matplotlib (API differs by version)
try:
    import matplotlib
    matplotlib.colormaps.register(_TEMP_CMAP_OBJ, name="blue_orange_red", force=True)
except AttributeError:
    import matplotlib.cm as _cm
    _cm.register_cmap(name="blue_orange_red", cmap=_TEMP_CMAP_OBJ)
_TEMP_CMAP = "blue_orange_red"


def scenario_color_map(labels: list[str],
                       palette: str = "temperature") -> dict[str, str]:
    """Return a mapping from scenario label to hex colour.

    Parameters
    ----------
    labels : list of str
        Scenario labels in display order (first = coolest, last = warmest).
    palette : str
        'temperature' (default) — evenly-spaced colours from a blue→yellow→red
        gradient, scaled to the number of scenarios.  Any number of treatments
        is supported; the colours are always drawn from the full gradient range.
        'discrete' — the original fixed discrete palette (up to 8 scenarios).
        Any valid matplotlib colormap name — used directly as a gradient.

    Returns
    -------
    dict mapping each label to a hex colour string.
    """
    import matplotlib.pyplot as plt
    import matplotlib.colors as mcolors
    n = len(labels)
    if n == 0:
        return {}

    if palette == "discrete":
        return {lab: _PALETTE[i % len(_PALETTE)] for i, lab in enumerate(labels)}

    # Resolve colormap name
    cmap_name = _TEMP_CMAP if palette == "temperature" else palette
    try:
        cmap = plt.get_cmap(cmap_name)
    except ValueError:
        # Unknown colormap — fall back to discrete
        return {lab: _PALETTE[i % len(_PALETTE)] for i, lab in enumerate(labels)}

    if n == 1:
        # Single scenario — use the mid-point of the gradient
        positions = [0.5]
    else:
        # Space evenly across [0.05, 0.95] to avoid the extreme white ends
        # of symmetric colormaps like coolwarm
        positions = [0.05 + 0.90 * i / (n - 1) for i in range(n)]

    result = {}
    for lab, pos in zip(labels, positions):
        rgba = cmap(pos)
        result[lab] = mcolors.to_hex(rgba)
    return result


# ── Aggregation helpers ───────────────────────────────────────────────────────

def rep_mean(pool_list: list[dict | None], key: str) -> pd.DataFrame | None:
    """Average a (types × days) DataFrame or Series across replicate run dicts."""
    frames = [a[key] for a in pool_list if a and a.get(key) is not None]
    if not frames:
        return None
    if isinstance(frames[0], pd.Series):
        try:
            arr = np.stack([f.reindex(frames[0].index).values for f in frames], axis=1)
            return pd.Series(np.nanmean(arr, axis=1), index=frames[0].index)
        except Exception:
            return frames[0]
    idx, cols = frames[0].index, frames[0].columns
    try:
        arr = np.stack([f.reindex(index=idx, columns=cols).values
                        for f in frames], axis=2)
        return pd.DataFrame(np.nanmean(arr, axis=2), index=idx, columns=cols)
    except Exception:
        return frames[0]


def rep_std(pool_list: list[dict | None], key: str) -> pd.DataFrame | None:
    """Std-dev of a (types × days) DataFrame or Series across replicate run dicts."""
    frames = [a[key] for a in pool_list if a and a.get(key) is not None]
    if len(frames) < 2:
        return None
    if isinstance(frames[0], pd.Series):
        try:
            arr = np.stack([f.reindex(frames[0].index).values for f in frames], axis=1)
            return pd.Series(np.nanstd(arr, axis=1, ddof=1), index=frames[0].index)
        except Exception:
            return None
    idx, cols = frames[0].index, frames[0].columns
    try:
        arr = np.stack([f.reindex(index=idx, columns=cols).values
                        for f in frames], axis=2)
        return pd.DataFrame(np.nanstd(arr, axis=2, ddof=1), index=idx, columns=cols)
    except Exception:
        return None


# ── Runtime / year helpers ────────────────────────────────────────────────────

def read_runtime(run_dir: Path | None) -> dict:
    """Read Runtime.csv → {'pulse': n_years, 'end_time': days_per_year}."""
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
    """Convert --years argument to (col_start, col_end) inclusive day range."""
    if years_arg is None:
        return 1, n_years * days_per_year
    if len(years_arg) == 1:
        yr = years_arg[0]
        return (yr - 1) * days_per_year + 1, yr * days_per_year
    yr_start, yr_end = years_arg[0], years_arg[-1]
    return (yr_start - 1) * days_per_year + 1, yr_end * days_per_year


def filter_cols(df, col_start: int, col_end: int):
    """Slice a DataFrame (rows × days) or Series (days,) to [col_start, col_end]."""
    if df is None:
        return None
    if isinstance(df, pd.Series):
        return df.loc[col_start:col_end]
    cols = [c for c in df.columns if col_start <= c <= col_end]
    return df[cols] if cols else df


def add_year_labels(ax, col_start: int, col_end: int,
                    n_years: int, days_per_year: int) -> None:
    """Draw dashed year-boundary lines and 'Year N' labels on a time-series axes."""
    trans = ax.get_xaxis_transform()
    for yr in range(1, n_years + 1):
        boundary   = yr * days_per_year
        yr_start_d = (yr - 1) * days_per_year + 1
        yr_end_d   = yr * days_per_year
        if col_start < boundary < col_end:
            ax.axvline(boundary, color="gray", lw=0.6, ls=":", alpha=0.45, zorder=0)
        label_day = (max(yr_start_d, col_start) + min(yr_end_d, col_end)) / 2
        if col_start <= label_day <= col_end:
            ax.text(label_day, 0.98, f"Year {yr}",
                    transform=trans, ha="center", va="top",
                    fontsize=6.5, color="dimgray", alpha=0.75)


def add_year_labels_heatmap(ax, mat_ncols: int, col_start: int, col_end: int,
                             n_years: int, days_per_year: int) -> None:
    """Draw year-boundary markers on a subsampled imshow axes."""
    total = col_end - col_start + 1
    for yr in range(1, n_years + 1):
        boundary = yr * days_per_year
        if col_start < boundary < col_end:
            frac  = (boundary - col_start) / total
            x_pos = frac * (mat_ncols - 1)
            ax.axvline(x_pos, color="white", lw=0.8, ls=":", alpha=0.6, zorder=5)


# ── CLI helpers ───────────────────────────────────────────────────────────────

def add_scenario_args(parser) -> None:
    """Add the standard multi-scenario CLI flags to an argparse parser."""
    parser.add_argument(
        "--base", default=None,
        help="Root output directory.  Scenario subdirectories are auto-discovered "
             "from here unless --scenario-dirs is also given.",
    )
    parser.add_argument(
        "--scenario-dirs", nargs="+", default=None, metavar="PATH",
        help="Explicit list of scenario directories (each contains replicate "
             "sub-dirs).  Labels default to directory stems.  Overrides "
             "auto-discovery from --base.",
    )
    parser.add_argument(
        "--scenario-labels", nargs="+", default=None, metavar="LABEL",
        help="Human-readable display labels for scenarios, in the same order as "
             "the scenarios being loaded.  Works with both --scenario-dirs "
             "(positional match) and --base auto-discovery (positional match "
             "against alphabetically sorted subdirectories).  Use quotes for "
             "labels that contain spaces: --scenario-labels Control 'Warm +2C'.",
    )
    parser.add_argument(
        "--treatment-map", nargs="+", default=None, metavar="DIR=LABEL",
        help="Map directory stems to display labels without re-listing paths. "
             "Each entry is DIR_STEM=LABEL, e.g.: "
             "--treatment-map static_T10_Pn0p5=Control static_T15_Pn0p5='Warm +2C'. "
             "Applies during both --base auto-discovery and --scenario-dirs. "
             "Takes precedence over --scenario-labels when both are given for "
             "the same directory.",
    )
    parser.add_argument(
        "--palette", default="temperature", metavar="PALETTE",
        help="Colour palette for scenarios. 'temperature' (default): blue→orange→red "
             "gradient scaled to treatment order. 'discrete': fixed distinguishable "
             "colours. Any matplotlib colormap name is also accepted (e.g. 'viridis', "
             "'plasma'). Treatment order is set by --treatment-map or --scenario-labels.",
    )
    parser.add_argument(
        "--fig-dir", default=None,
        help="Directory to write figures (default: <base>/<script>_figures).",
    )
    parser.add_argument(
        "--years", nargs="+", type=int, default=None, metavar="YEAR",
        help="Year(s) to display.  Single value = that year only. "
             "Two values = inclusive range (e.g. --years 2 4). "
             "Default: all years.",
    )


def _parse_treatment_map(entries: list[str] | None) -> dict[str, str]:
    """Parse --treatment-map DIR=LABEL entries into a {stem: label} dict.

    Entries are split on the *first* '=' only, so labels may contain '='.
    Entries that contain no '=' are ignored with a warning.
    """
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


def resolve_scenario_list(
    args,
    base_fallback: Path,
    min_cols: int = 100,
    sentinel: str | None = None,
) -> list[tuple[str, list[Path]]]:
    """Return [(label, [run_dirs])] from parsed CLI args.

    Label resolution priority (highest wins):
    1. --treatment-map  DIR_STEM=LABEL  (overrides everything for matched stems)
    2. --scenario-labels  (positional, matched in discovery/dir order)
    3. Directory stem (fallback)

    Path resolution priority:
    1. --scenario-dirs  (explicit list, with optional --scenario-labels or --treatment-map)
    2. Auto-discovery from --base  (subdirs of base that contain run dirs)
    3. base_fallback  (hardcoded BASE path in the calling script)
    """
    base = Path(args.base) if args.base else base_fallback
    tmap = _parse_treatment_map(getattr(args, "treatment_map", None))

    def _apply_labels(pairs: list[tuple[str, Path]],
                      positional_labels: list[str] | None) -> list[tuple[str, Path]]:
        """Apply label overrides to a [(stem, path)] list."""
        result = []
        for i, (stem, path) in enumerate(pairs):
            # treatment-map wins; then positional label; then stem
            if stem in tmap:
                label = tmap[stem]
            elif positional_labels and i < len(positional_labels):
                label = positional_labels[i]
            else:
                label = stem
            result.append((label, path))
        return result

    if args.scenario_dirs:
        paths = [Path(p) for p in args.scenario_dirs]
        pairs = [(p.stem, p) for p in paths]
        labeled = _apply_labels(pairs, getattr(args, "scenario_labels", None))
        result = []
        for label, sdir in labeled:
            runs = pick_runs(sdir, min_cols=min_cols, sentinel=sentinel)
            if not runs:
                print(f"  WARNING: no valid runs found in {sdir}")
            result.append((label, runs))
        return result

    # Auto-discover from base
    found = discover_scenarios(base)   # [(stem, path), ...]
    if not found:
        print(f"  WARNING: no scenario subdirs found under {base}")
        return []

    labeled = _apply_labels(found, getattr(args, "scenario_labels", None))

    # Re-order to match --treatment-map entry order when provided.
    # tmap is an ordered dict (Python 3.7+) so tmap.values() preserves
    # the order the user typed the entries on the command line.
    if tmap:
        tmap_labels  = list(tmap.values())   # desired display-label order
        label_to_pair = {lbl: (lbl, sdir) for lbl, sdir in labeled}
        ordered  = [label_to_pair[lbl] for lbl in tmap_labels
                    if lbl in label_to_pair]
        # Append any scenarios not covered by treatment-map at the end
        covered  = set(tmap_labels)
        ordered += [(lbl, sdir) for lbl, sdir in labeled if lbl not in covered]
        labeled  = ordered

    result = []
    for label, sdir in labeled:
        runs = pick_runs(sdir, min_cols=min_cols, sentinel=sentinel)
        result.append((label, runs))
    return result


# ── Gridsize helpers ──────────────────────────────────────────────────────────

def read_gridsize(run_dir: Path | None) -> int:
    """Read gridsize from Runtime.csv in a run directory.

    Falls back to reading x*y from runtime parameters if gridsize row is absent.
    Returns 1 if the file cannot be read (so division is a no-op).
    """
    if run_dir is None:
        return 1
    f = run_dir / "Runtime.csv"
    if not f.exists():
        return 1
    try:
        rt  = pd.read_csv(f, index_col=0, header=0)
        col = rt.columns[0]
        if "gridsize" in rt.index:
            return int(float(rt.loc["gridsize", col]))
        # Fall back to x * y
        if "x" in rt.index and "y" in rt.index:
            return int(float(rt.loc["x", col])) * int(float(rt.loc["y", col]))
    except Exception:
        pass
    return 1


def normalise_by_gridsize(df, gridsize: int):
    """Divide a DataFrame or Series by gridsize (in-place copy).

    Returns None unchanged.  gridsize=1 is a no-op (safe fallback).
    """
    if df is None or gridsize <= 1:
        return df
    return df / gridsize


def get_scenario_gridsize(scenario_list: list, run_index: int = 0) -> int:
    """Get gridsize from the first run of the first scenario in scenario_list.

    Parameters
    ----------
    scenario_list : list of (label, [run_dirs])
        As returned by resolve_scenario_list.
    run_index : int
        Which run within the first scenario to read from (default 0).
    """
    if not scenario_list or not scenario_list[0][1]:
        return 1
    run_dirs = scenario_list[0][1]
    if run_index >= len(run_dirs):
        run_index = 0
    return read_gridsize(run_dirs[run_index])

def get_palette_from_args(args) -> str:
    """Extract --palette value from parsed args, defaulting to 'temperature'."""
    return getattr(args, "palette", "temperature")
