#!/usr/bin/env python3
"""
diagnose_nan.py
===============
Find where a NaN first enters the model state of one community.

Loads the community's saved spin-up state (<output-dir>/spinup/<community>_s<seed>.pkl,
written by run_multi_climate_simulations.py), continues it day by day under one
climate file exactly as a branch run does, and after every phase of every day
(degradation, uptake, metabolism, mortality, reproduction) checks the state
tables for NaN or inf. At the first one it reports the day, the phase, the
tables and rows affected (grid cell and taxon / substrate / monomer / enzyme),
and the values of that grid cell before and after the phase, then stops.

Run it from the project root, with the same --src-dir / --input-dir /
--community-dir / --output-dir as the simulation (a compute node: it needs the
memory of one simulation task):

  python analysis/diagnose_nan.py --community community_0308 \\
      --src-dir src --input-dir input_summer_warming_20261113 \\
      --community-dir community_pool_20261113/communities \\
      --output-dir /pub/.../sim_summer_warm_20261113 \\
      --climate climate_data/climate_ambient_temp.csv [--max-days 120] [--report FILE]
"""

from __future__ import annotations

import argparse
import os
import pickle
import sys
from pathlib import Path

import numpy as np
import pandas as pd

_HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(_HERE))
import run_multi_climate_simulations as runner      # noqa: E402

STATE = ["Microbes", "Substrates", "Monomers", "Enzymes", "_hsp_pool", "_osmo_pool"]
PHASE_NAMES = ["degradation", "uptake", "metabolism", "mortality", "reproduction"]
_out = []


def say(*a):
    line = " ".join(str(x) for x in a)
    print(line, flush=True)
    _out.append(line)


def values(x):
    return x.to_numpy() if hasattr(x, "to_numpy") else np.asarray(x)


def bad_mask(x):
    v = values(x)
    return ~np.isfinite(v.astype(np.float64)) if v.dtype.kind in "fc" else np.zeros(v.shape, bool)


def snapshot(eco):
    return {n: (getattr(eco, n).copy() if hasattr(eco, n) else None) for n in STATE}


def n_bad(eco):
    return {n: int(bad_mask(getattr(eco, n)).sum()) for n in STATE if hasattr(eco, n)}


def restore(state_file, cdir, spinup_climate, seed, src_dir, input_dir):
    """(Ecosystem, Output_init, numpy random state) from a spin-up file of either format."""
    with open(state_file, "rb") as fh:
        meta = pickle.load(fh)
        payload = fh.read()
    say(f"Spin-up file: {state_file}  ({len(payload) / 1e6:.0f} MB, format {meta.get('format')})")
    runner._import_model(src_dir)
    if hasattr(runner, "SpinupState") and meta.get("format", 1) >= 2:
        st = runner.SpinupState(payload, cdir, spinup_climate, seed, src_dir, input_dir)
        return runner.restore_spinup(st)
    return pickle.loads(payload)                     # format 1: the three objects, pickled


def rows_of(obj, gridsize):
    """Rows per grid cell of a state table, and its row labels (or None)."""
    n = values(obj).shape[0]
    per = n // gridsize
    labels = list(obj.index[:per]) if hasattr(obj, "index") else None
    return per, labels


def describe(eco, name, before, gridsize, limit=15):
    """List the non-finite entries of one state table; return the grid cells affected."""
    obj = getattr(eco, name)
    mask = bad_mask(obj)
    if mask.ndim == 1:
        mask = mask[:, None]
    per, labels = rows_of(obj, gridsize)
    cols = list(obj.columns) if hasattr(obj, "columns") else list(range(mask.shape[1]))
    rc = np.argwhere(mask)
    cells = sorted({int(r) // per for r, _ in rc})
    say(f"  {name}: {len(rc)} non-finite value(s) in {len(cells)} grid cell(s); "
        f"{len({int(r) % per for r, _ in rc})} of {per} rows per cell affected")
    v_new = values(obj).reshape(mask.shape)
    v_old = values(before[name]).reshape(mask.shape) if before.get(name) is not None else None
    for r, c in rc[:limit]:
        lab = labels[int(r) % per] if labels else int(r) % per
        old = "" if v_old is None else f"   before the phase: {v_old[r, c]!r}"
        say(f"    cell {int(r) // per:5d}  {str(lab):16s} {str(cols[c]):8s} = {v_new[r, c]!r}{old}")
    if len(rc) > limit:
        say(f"    ... {len(rc) - limit} more")
    return cells


def dump_cell(eco, before, cell, gridsize):
    """State of one grid cell before and after the phase (rows that are non-zero)."""
    for name in STATE:
        if not hasattr(eco, name):
            continue
        obj = getattr(eco, name)
        per, labels = rows_of(obj, gridsize)
        sl = slice(cell * per, (cell + 1) * per)
        new = values(obj)[sl]
        old = values(before[name])[sl]
        new2, old2 = new.reshape(per, -1), old.reshape(per, -1)
        keep = np.flatnonzero((old2 != 0).any(axis=1) | (new2 != 0).any(axis=1) | ~np.isfinite(new2).all(axis=1))
        cols = list(obj.columns) if hasattr(obj, "columns") else [f"c{i}" for i in range(new2.shape[1])]
        say(f"\n  {name}, grid cell {cell}: {len(keep)} of {per} rows non-zero   (columns: {cols})")
        for i in keep[:40]:
            lab = labels[i] if labels else i
            say(f"    {str(lab):16s} before {np.array2string(old2[i], precision=6, max_line_width=200)}"
                f"   after {np.array2string(new2[i], precision=6, max_line_width=200)}")
        if len(keep) > 40:
            say(f"    ... {len(keep) - 40} more rows")


def other_bad_attributes(eco, known):
    """Names of every other Grid attribute that holds a non-finite number."""
    found = []
    for k, v in vars(eco).items():
        if k in known:
            continue
        try:
            if isinstance(v, (pd.DataFrame, pd.Series, np.ndarray)) and bad_mask(v).any():
                found.append(f"{k} ({int(bad_mask(v).sum())})")
            elif isinstance(v, (float, np.floating)) and not np.isfinite(v):
                found.append(f"{k} (scalar)")
        except (TypeError, ValueError):
            pass
    return found


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--community", required=True, help="community folder name, e.g. community_0308")
    ap.add_argument("--src-dir", default="src")
    ap.add_argument("--input-dir", required=True)
    ap.add_argument("--community-dir", required=True, help="the pool's communities folder")
    ap.add_argument("--output-dir", required=True, help="simulation output folder (holds spinup/)")
    ap.add_argument("--climate", required=True, help="climate file to continue with (inside --input-dir or a path)")
    ap.add_argument("--seed", type=int, default=0, help="base seed of the simulation (default 0)")
    ap.add_argument("--max-days", type=int, default=0, help="stop after this many days (default: to the end)")
    ap.add_argument("--report", default=None, help="also write the report to this file "
                                                   "(default diagnose_nan_<community>.txt in the working directory)")
    a = ap.parse_args()

    src_dir, input_dir = os.path.abspath(a.src_dir), os.path.abspath(a.input_dir)
    pool = Path(os.path.abspath(a.community_dir))
    cdir = pool / a.community
    if not cdir.is_dir():
        sys.exit(f"ERROR: {cdir} not found")
    clim = next((p for p in (os.path.join(input_dir, a.climate), os.path.abspath(a.climate))
                 if os.path.isfile(p)), None)
    if clim is None:
        sys.exit(f"ERROR: climate file {a.climate!r} not found")
    position = {d.name: k for k, d in enumerate(runner.get_community_dirs(pool))}
    seed = (a.seed + position[a.community]) % (2 ** 32)
    state_file = Path(a.output_dir) / "spinup" / f"{a.community}_s{seed}.pkl"
    if not state_file.exists():
        sys.exit(f"ERROR: {state_file} not found")

    runtime, _ = runner.load_runtime(input_dir)
    col = runtime.columns[0]
    cycle, pulse = int(runtime.loc["end_time", col]), int(runtime.loc["pulse", col])
    branch_pulse = int(runtime.loc["branch_pulse", col])
    interval = int(runtime.loc["interval", col])
    gridsize = int(runtime.loc["gridsize", col])

    say(f"Community {a.community} (seed {seed}); climate {Path(clim).name}; "
        f"days {branch_pulse * cycle}-{pulse * cycle - 1}")
    eco, out, rng = restore(state_file, cdir, clim, seed, src_dir, input_dir)
    temp, psi = runner.read_climate_series(clim)
    eco.temp, eco.psi = temp, psi
    np.random.set_state(rng)

    start = n_bad(eco)
    if any(start.values()):
        say(f"Non-finite values are already in the spin-up state: {start}")
        for name, n in start.items():
            if n:
                describe(eco, name, {}, gridsize)
        _finish(a)
        return

    steps = {"degradation": lambda p, i: eco.degradation(p, i), "uptake": lambda p, i: eco.uptake(p, i),
             "metabolism": lambda p, i: eco.metabolism(i), "mortality": lambda p, i: eco.mortality(i),
             "reproduction": lambda p, i: eco.reproduction(i)}
    n_run = 0
    for p in range(branch_pulse, pulse):
        for i in range(p * cycle, (p + 1) * cycle):
            for phase in PHASE_NAMES:
                before = snapshot(eco)
                with np.errstate(all="ignore"):
                    steps[phase](p, i)
                counts = n_bad(eco)
                if any(counts.values()):
                    say(f"\nFIRST NON-FINITE VALUE: day {i} (output column {i + 1}), after {phase}()")
                    say(f"  Temp {float(np.asarray(temp)[i]):.3f}   Psi {float(np.asarray(psi)[i]):.3f}")
                    cells = []
                    for name, n in counts.items():
                        if n:
                            cells += describe(eco, name, before, gridsize)
                    others = other_bad_attributes(eco, set(STATE))
                    if others:
                        say("  Other Grid attributes holding non-finite values: " + ", ".join(others[:40]))
                    for cell in sorted(set(cells))[:2]:
                        dump_cell(eco, before, cell, gridsize)
                    _finish(a)
                    return
            if i == 0 or i % interval == interval - 1:
                out.output(eco, i)
            out.microbes_abundance(eco, i)
            if i == (p + 1) * cycle - 1:
                out.flush_pulse_buffers()
                eco.repopulation(out, i, True)
            n_run += 1
            if n_run % 10 == 0:
                print(f"  day {i} done, state finite", flush=True)
            if a.max_days and n_run >= a.max_days:
                say(f"\nNo non-finite value in the state tables in {n_run} days (through day {i}).")
                _finish(a)
                return
    say(f"\nNo non-finite value in the state tables through day {pulse * cycle - 1}.")
    _finish(a)


def _finish(a):
    path = a.report or f"diagnose_nan_{a.community}.txt"
    with open(path, "w") as fh:
        fh.write("\n".join(_out) + "\n")
    print(f"\nReport written to {path}", flush=True)


if __name__ == "__main__":
    main()
