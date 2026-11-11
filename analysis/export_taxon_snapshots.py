#!/usr/bin/env python3
"""
export_taxon_snapshots.py
=========================
Per-taxon snapshots on chosen days, for every community and climate scenario,
kept per community (no averaging across communities), so scenarios can be
compared taxon by taxon within a community (paired design).

One row per community x scenario x taxon x day.  Taxa that have died out stay
in the table with cells = 0.

COLUMNS
-------
  community, scenario, taxon, day      day k = state after day k-1 of the
                                       climate file (output column label k)

  State at the end of the day (grid totals)
    cells           living cells
    biomass_C       biomass C (HSP and osmolyte C included)
    energy          energy pool
    hsp_pool_C      HSP pool C (summed over HSP types)
    osmo_pool_C     osmolyte pool C (summed over osmolytes)
    drought_tol     drought tolerance, mean over living cells (NaN when cells = 0)

  Fluxes during the day (grid totals = per-living-microbe value in the model
  output x living cells; 0 when cells = 0, see NOTES)
    uptake_C, biosynthesis_C          C taken up; the part sent to biosynthesis
    enzyme_C, hsp_C, osmolyte_C, transport_C    C to each production pathway
    respiration_C   all CO2: uptake metabolism + endogenous + rebalancing
    endo_resp_C     biomass C respired to pay an energy deficit (R_endo)
    rebalance_C     C respired in stoichiometric rebalancing (R_stoich E / yield)
    retained_C      biosynthesis_C - enzyme_C - endo_resp_C - rebalance_C
    production_E    energy from uptake metabolism (endogenous respiration and
                    rebalancing not included)
    stoich_E        energy credited from rebalanced C
    transport_E, maintenance_E, enzyme_E, hsp_E, osmolyte_E    energy costs
    unpaid_E        energy demand left unpaid after endogenous respiration

  Sums over the WINDOW days ending on the day (default 7)
    uptake_C_7d, retained_C_7d

Derived at analysis time, e.g.:
    daily CUE   = retained_C / uptake_C        7-day CUE = retained_C_7d / uptake_C_7d
    energy fraction of uptake = (uptake_C - biosynthesis_C) / uptake_C
    per biomass = x / biomass_C                per cell = x / cells
    energy from endogenous respiration = endo_resp_C x E_yield (in the meta file)

OUTPUT  (<out-dir>/)
------
  parts/<community>.parquet        all scenarios of one community (work files;
  parts/<community>.traits.csv     a rerun with the same settings skips
                                   communities that have them)
  taxon_snapshots_<scenario>.parquet   one table per scenario
  taxon_traits.csv                 one row per taxon: trait values from the
                                   community inputs
  taxon_snapshots_meta.csv         days, window, E_yield, version, git hash

  --format csv writes .csv.gz instead of .parquet (no pyarrow needed).

NOTES
-----
* Per-taxon labels are used as the model wrote them, so the runs must come
  from an output.py that writes the pool's taxon IDs (all runs since Sep 2026).
* The model divides each flux by the living cells at the end of the day, so
  flux x cells is the grid total exactly while a taxon has living cells.  On
  the day its last cell dies the flux of that day is not recoverable (written
  as 0 here).
* Only the scenario folders are read when they hold every day needed (split
  branched output and days after the branch); otherwise the full series is
  read through branch_output.read_table().

USAGE
-----
  python analysis/export_taxon_snapshots.py --base /path/to/sim_output \\
      --out-dir /path/to/taxon_snapshots -j 16

  # other days / window, a subset of scenarios or communities
  python analysis/export_taxon_snapshots.py --base SIM --out-dir OUT \\
      --days 2013 2020 2027 --window 7 --scenarios climate_ambient_temp climate_June_heatwave
"""

from __future__ import annotations

import argparse
import gzip
import json
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd

_HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(_HERE))

from branch_output import initialization_dir, read_table          # noqa: E402
from export_summary_data import fs_retry, _git_short_hash, _resolve_jobs  # noqa: E402
from scenario_utils import discover_scenarios, pick_runs          # noqa: E402

SNAPSHOT_VERSION = "2026-11-10"

# Summer heatwave-timing experiment (year 6): every 4 weeks from day 1845, the
# day before / first / last day of each 7-day heatwave (June 2021-2027, July
# 2049-2055, August 2077-2083), weekly between and for 3 weeks after them, then
# every 3 weeks to the end of the year.
DEFAULT_DAYS = [1845, 1873, 1901, 1929, 1957, 1985,
                2013, 2020, 2021, 2027, 2034, 2041,
                2048, 2049, 2055, 2062, 2069,
                2076, 2077, 2083, 2090, 2097, 2104,
                2118, 2139, 2160, 2181, 2190]

# column -> model output table, per living microbe (x cells = grid total)
FLUXES = {
    "uptake_C":       "C_Uptake_Taxon",
    "biosynthesis_C": "BC_Uptake_Taxon",
    "enzyme_C":       "Enzyme_C_Taxon",
    "hsp_C":          "HSP_C_Taxon",
    "osmolyte_C":     "Osmo_C_Taxon",
    "transport_C":    "Transport_C_Taxon",
    "respiration_C":  "Respiration_Taxon",
    "endo_resp_C":    "R_endo_C_Taxon",
    "production_E":   "E_Prod_Taxon",
    "stoich_E":       "R_stoich_E_Taxon",
    "transport_E":    "Transport_E_Taxon",
    "maintenance_E":  "Maint_E_Taxon",
    "enzyme_E":       "Enzyme_E_Taxon",
    "hsp_E":          "HSP_E_Taxon",
    "osmolyte_E":     "Osmo_E_Taxon",
    "unpaid_E":       "E_unpaid_Taxon",
}
# column -> model output table, grid totals per taxon
POOLS = {"biomass_C": "MicrobesSeries", "energy": "EnergySeries"}
# column -> model output table with a (taxon, type) index, summed over types
TYPED_POOLS = {"hsp_pool_C": "Taxon_HSP_Series", "osmo_pool_C": "Taxon_Osmo_Series"}
WINDOW_SUMS = ("uptake_C", "retained_C")

COLUMNS = (["community", "scenario", "taxon", "day", "cells", "biomass_C", "energy",
            "hsp_pool_C", "osmo_pool_C", "drought_tol",
            "uptake_C", "biosynthesis_C", "enzyme_C", "hsp_C", "osmolyte_C", "transport_C",
            "respiration_C", "endo_resp_C", "rebalance_C", "retained_C",
            "production_E", "stoich_E",
            "transport_E", "maintenance_E", "enzyme_E", "hsp_E", "osmolyte_E", "unpaid_E"]
           + [f"{c}_{{w}}d" for c in WINDOW_SUMS])

# taxon traits: column -> (Initialization table, how to reduce it to one value per taxon)
TRAITS = {
    "enz_genes":        ("EnzGenes", "sum"),
    "enz_prod_consti":  ("EnzProdConsti_trait", "value"),
    "enz_prod_induci":  ("EnzProdInduci_trait", "value"),
    "uptake_genes":     ("UptakeGenes", "sum"),
    "uptake_trait":     ("UptakeGenes_trait", "value"),
    "hsp_genes":        ("HSPGenes", "sum"),
    "hsp_prod_consti":  ("HSPProdConsti_trait", "value"),
    "hsp_prod_induci":  ("HSPProdInduci_trait", "value"),
    "osmo_genes":       ("OsmoGenes", "sum"),
    "osmo_prod_consti": ("OsmoProdConsti_trait", "value"),
    "osmo_prod_induci": ("OsmoProdInduci_trait", "value"),
}


def _columns(window: int) -> list[str]:
    return [c.format(w=window) for c in COLUMNS]


# ── Reading one run ───────────────────────────────────────────────────────────

def _table(run_dir: Path, name: str, need: list[int], index_cols: int = 1) -> pd.DataFrame:
    """Output table `name` restricted to the day columns in `need` (int labels)."""
    def _read():
        return pd.read_csv(run_dir / f"{name}.csv", index_col=list(range(index_cols)))
    df = fs_retry(_read, f"reading {run_dir.name}/{name}")
    df.columns = df.columns.astype(int)
    if not set(need) <= set(df.columns):               # days before the branch: full series
        df = fs_retry(lambda: read_table(run_dir, name), f"reading {run_dir.name}/{name} (full)")
        df.columns = df.columns.astype(int)
    missing = sorted(set(need) - set(df.columns))
    if missing:
        raise ValueError(f"{run_dir.name}/{name}: day(s) {missing[:5]} not in the output")
    df.index = df.index.map(str) if index_cols == 1 else df.index
    return df[need]


def _e_yield(run_dir: Path) -> float:
    """Energy per C respired endogenously or in rebalancing (E_resp of Mon_DeadMic)."""
    p = initialization_dir(run_dir) / "MonEnergy_Attrib.csv"
    return float(fs_retry(lambda: pd.read_csv(p, index_col=0), f"reading {p.name}")
                 .loc["Mon_DeadMic", "E_resp"])


def snapshot_run(run_dir: Path, days: list[int], window: int, e_yield: float) -> pd.DataFrame:
    """Snapshot table of one run: one row per taxon x day, columns as in COLUMNS
    (without community and scenario)."""
    need = sorted({d - k for d in days for k in range(window)})
    cells = _table(run_dir, "Taxon_count", need)
    taxa = cells.index
    if taxa.has_duplicates:
        raise ValueError(f"{run_dir.name}: duplicate taxon labels in the output")
    alive = cells.to_numpy() > 0
    n = cells.to_numpy(np.float64)

    def total(name):                                   # per living microbe -> grid total
        x = _table(run_dir, name, need).reindex(taxa).to_numpy(np.float64)
        return np.where(alive, x, 0.0) * n

    flux = {col: total(name) for col, name in FLUXES.items()}
    flux["rebalance_C"] = flux["stoich_E"] / e_yield
    flux["retained_C"] = (flux["biosynthesis_C"] - flux["enzyme_C"]
                          - flux["endo_resp_C"] - flux["rebalance_C"])

    pos = {d: i for i, d in enumerate(need)}
    at = [pos[d] for d in days]                        # columns of the snapshot days
    out = {"cells": cells.to_numpy()[:, at]}
    for col, name in POOLS.items():
        out[col] = _table(run_dir, name, days).reindex(taxa).to_numpy()
    for col, name in TYPED_POOLS.items():
        out[col] = (_table(run_dir, name, days, index_cols=2)
                    .groupby(level=0, sort=False).sum().rename(index=str).reindex(taxa).to_numpy())
    tol = _table(run_dir, "Drought_tol_Taxon", days).reindex(taxa).to_numpy(np.float64)
    out["drought_tol"] = np.where(alive[:, at], tol, np.nan)
    for col, x in flux.items():
        out[col] = x[:, at]
    for col in WINDOW_SUMS:                            # trailing sums over `window` days
        out[f"{col}_{window}d"] = np.stack(
            [flux[col][:, pos[d] - window + 1: pos[d] + 1].sum(axis=1) for d in days], axis=1)

    df = pd.DataFrame({"taxon": np.repeat(taxa.to_numpy(), len(days)),
                       "day": np.tile(np.asarray(days, dtype=np.int16), len(taxa))})
    for col, x in out.items():
        df[col] = np.asarray(x).reshape(-1).astype(np.uint32 if col == "cells" else np.float32)
    return df


def taxon_traits(run_dir: Path) -> pd.DataFrame:
    """Trait values of a run's taxa, from the community inputs in Initialization."""
    init = initialization_dir(run_dir)
    cols = {}
    for col, (name, how) in TRAITS.items():
        p = init / f"{name}.csv"
        if not p.exists():
            continue
        t = fs_retry(lambda: pd.read_csv(p, index_col=0), f"reading {p.name}")
        cols[col] = t.sum(axis=1) if how == "sum" else t.iloc[:, 0]
    df = pd.DataFrame(cols)
    df.index = df.index.map(str)
    return df.rename_axis("taxon")


# ── One community (worker) ────────────────────────────────────────────────────

def _write_table(df: pd.DataFrame, path: Path, fmt: str) -> None:
    tmp = path.with_name(path.name + ".tmp")

    def _write():
        if fmt == "parquet":
            df.to_parquet(tmp, index=False, compression="zstd")
        else:
            df.to_csv(tmp, index=False, compression="gzip")
        os.replace(tmp, path)
    fs_retry(_write, f"writing {path.name}")


def _ext(fmt: str) -> str:
    return "parquet" if fmt == "parquet" else "csv.gz"


def export_community(name: str, runs: list[tuple[str, str]], parts: str,
                     days: list[int], window: int, fmt: str) -> tuple[str, int]:
    """Write parts/<name>.<ext> (all scenarios) and parts/<name>.traits.csv."""
    first = Path(runs[0][1])
    e_yield = _e_yield(first)
    frames = []
    for scenario, run_dir in runs:
        df = snapshot_run(Path(run_dir), days, window, e_yield)
        df.insert(0, "scenario", scenario)
        df.insert(0, "community", name)
        frames.append(df)
    table = pd.concat(frames, ignore_index=True)[_columns(window)]
    traits = taxon_traits(first)
    fs_retry(lambda: traits.to_csv(Path(parts) / f"{name}.traits.csv"), f"writing {name} traits")
    _write_table(table, Path(parts) / f"{name}.{_ext(fmt)}", fmt)   # last: marks the community done
    return name, len(table)


# ── Combining the parts ───────────────────────────────────────────────────────

def _read_part(path: Path, fmt: str) -> pd.DataFrame:
    read = pd.read_parquet if fmt == "parquet" else pd.read_csv
    return fs_retry(lambda: read(path), f"reading {path.name}")


def combine(parts: Path, out_dir: Path, names: list[str], scenarios: list[str],
            fmt: str, batch: int = 50) -> None:
    """Write taxon_snapshots_<scenario>.<ext> from the community parts (streamed,
    `batch` communities at a time) and taxon_traits.csv."""
    ext = _ext(fmt)
    finals = {s: out_dir / f"taxon_snapshots_{s}.{ext}" for s in scenarios}
    tmps = {s: p.with_name(p.name + ".tmp") for s, p in finals.items()}
    writers = {}
    try:
        for i in range(0, len(names), batch):
            chunk = pd.concat([_read_part(parts / f"{n}.{ext}", fmt) for n in names[i:i + batch]],
                              ignore_index=True)
            for s, sub in chunk.groupby("scenario", sort=False):
                if fmt == "parquet":
                    import pyarrow as pa
                    import pyarrow.parquet as pq
                    tbl = pa.Table.from_pandas(sub, preserve_index=False)
                    if s not in writers:
                        writers[s] = pq.ParquetWriter(tmps[s], tbl.schema, compression="zstd")
                    writers[s].write_table(tbl)
                else:
                    if s not in writers:
                        writers[s] = gzip.open(tmps[s], "wt", newline="")
                    sub.to_csv(writers[s], index=False, header=(i == 0))
            print(f"  combined {min(i + batch, len(names))}/{len(names)}", flush=True)
    finally:
        for w in writers.values():
            w.close()
    for s in scenarios:
        os.replace(tmps[s], finals[s])

    traits = pd.concat([fs_retry(lambda n=n: pd.read_csv(parts / f"{n}.traits.csv", index_col=0),
                                 f"reading {n} traits") for n in names])
    traits = traits[~traits.index.duplicated()].sort_index(
        key=lambda ix: ix.str.extract(r"(\d+)", expand=False).astype(float))
    traits.to_csv(out_dir / "taxon_traits.csv")


# ── Main ──────────────────────────────────────────────────────────────────────

def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("COLUMNS")[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--base", required=True, help="simulation output directory (scenario folders inside)")
    ap.add_argument("--out-dir", required=True, help="where the tables are written")
    ap.add_argument("--days", type=int, nargs="+", default=DEFAULT_DAYS,
                    help="snapshot days (output column labels); default: the summer heatwave-timing set")
    ap.add_argument("--window", type=int, default=7, help="days in the trailing sums (default 7)")
    ap.add_argument("--scenarios", nargs="+", default=None, help="scenario folder names to include (default all)")
    ap.add_argument("--communities", nargs="+", default=None, help="run folder names to include (default all)")
    ap.add_argument("-j", "--jobs", type=int, default=0, help="workers (0 = $SLURM_CPUS_PER_TASK or all cores)")
    ap.add_argument("--format", choices=("parquet", "csv"), default="parquet")
    ap.add_argument("--overwrite", action="store_true", help="redo communities that already have a part file")
    args = ap.parse_args()

    if args.format == "parquet":
        try:
            import pyarrow  # noqa: F401
        except ImportError:
            sys.exit("ERROR: pyarrow is not installed (conda install -c conda-forge pyarrow), "
                     "or use --format csv")
    days = sorted(set(args.days))
    base, out_dir = Path(args.base), Path(args.out_dir)
    parts = out_dir / "parts"
    parts.mkdir(parents=True, exist_ok=True)

    found = [(s, d) for s, d in discover_scenarios(base) if not args.scenarios or s in args.scenarios]
    if not found or (args.scenarios and len(found) != len(set(args.scenarios))):
        sys.exit(f"ERROR: scenarios not found under {base}: "
                 f"{sorted(set(args.scenarios or []) - {s for s, _ in found}) or 'none at all'}")
    runs = {s: {r.name: r for r in pick_runs(d, min_cols=100)} for s, d in found}
    scenarios = [s for s, _ in found]
    names = sorted(set.intersection(*(set(r) for r in runs.values())))
    if args.communities:
        names = [n for n in names if n in set(args.communities)]
    if not names:
        sys.exit("ERROR: no community is complete in every scenario")
    print(f"{len(scenarios)} scenarios x {len(names)} communities x {len(days)} days "
          f"(window {args.window}) -> {out_dir}", flush=True)

    # Part files are only reused when they were made with the same settings.
    settings = {"snapshot_version": SNAPSHOT_VERSION, "days": days, "window": args.window,
                "scenarios": scenarios, "format": args.format}
    stamp = parts / "settings.json"
    if stamp.exists() and not args.overwrite and json.loads(stamp.read_text()) != settings:
        sys.exit(f"ERROR: {parts} holds part files made with other settings (see {stamp.name}); "
                 f"use a new --out-dir, or --overwrite to redo them")
    stamp.write_text(json.dumps(settings, indent=1))

    ext = _ext(args.format)
    todo = [n for n in names if args.overwrite or not (parts / f"{n}.{ext}").exists()]
    if len(todo) < len(names):
        print(f"  {len(names) - len(todo)} communities already done, skipped", flush=True)
    failed, t0 = [], time.perf_counter()
    with ProcessPoolExecutor(max_workers=_resolve_jobs(args.jobs)) as ex:
        futures = {ex.submit(export_community, n, [(s, str(runs[s][n])) for s in scenarios],
                             str(parts), days, args.window, args.format): n for n in todo}
        for k, fut in enumerate(as_completed(futures), start=1):
            n = futures[fut]
            try:
                fut.result()
                print(f"  [{k}/{len(todo)}] {n}", flush=True)
            except Exception as e:
                failed.append(n)
                print(f"  [{k}/{len(todo)}] {n} FAILED: {e}", flush=True)
    print(f"  extraction: {time.perf_counter() - t0:.0f} s", flush=True)
    if failed:
        sys.exit(f"ERROR: {len(failed)} communities failed ({', '.join(failed[:10])}"
                 f"{' ...' if len(failed) > 10 else ''}); rerun the same command to retry them")

    combine(parts, out_dir, names, scenarios, args.format)
    pd.Series({
        "snapshot_version": SNAPSHOT_VERSION,
        "export_git": _git_short_hash(),
        "export_time": time.strftime("%Y-%m-%d %H:%M:%S"),
        "base": str(base),
        "n_communities": len(names),
        "scenarios": " ".join(scenarios),
        "days": " ".join(map(str, days)),
        "window": args.window,
        "E_yield": _e_yield(runs[scenarios[0]][names[0]]),
        "format": args.format,
    }, name="value").to_csv(out_dir / "taxon_snapshots_meta.csv")
    print(f"Done.  All output written to {out_dir}", flush=True)


if __name__ == "__main__":
    main()
