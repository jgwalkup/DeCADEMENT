"""Check that the climate scenarios of each community are paired with the baseline run.

    python analysis/check_pairing.py OUTPUT_DIR [--baseline climate_ambient_temp] [--jobs N]

OUTPUT_DIR is the --output-dir of run_multi_climate_simulations.py (one subfolder per scenario, each holding
<community>_s<seed>/ run folders). For every community and every non-baseline scenario, the script finds:

  - the first day on which the scenario's temperature or water potential differs from the baseline's
    (Initialization/Temp.csv and Psi.csv, which each run writes);
  - the first day on which the model output differs (SubstratesSeries, MicrobesSeries, RespSeries).

Series columns are the initial state followed by one column per day, so column k is the state after day k - 1.
For split branched output (branch_output.json in the run folders) the scenario folders hold only the days from the
branch on; the spin-up days are shared by construction, so only the branch days are compared.
A paired run's output matches the baseline exactly until its climate differs. Output that differs before the
climate does means the runs did not start from the same random state, and paired comparisons between them are
not valid.
"""
import argparse
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

SERIES = ["SubstratesSeries", "MicrobesSeries", "RespSeries"]


def _load_climate(run):
    out = {}
    for f in ("Temp", "Psi"):
        pa = run / "Initialization" / f"{f}.csv"
        if pa.exists():
            out[f] = pd.read_csv(pa, index_col=0).iloc[:, 0].to_numpy()
    return out


def _load_series(run):
    out = {}
    for f in SERIES:
        pa = run / f"{f}.csv"
        if pa.exists():
            x = pd.read_csv(pa, index_col=0).sort_index()
            out[f] = x[sorted(x.columns, key=int)]
    return out


def first_climate_diff(a_clim, b_clim):
    """First day on which Temp or Psi differs, or None if identical."""
    days = []
    for f in ("Temp", "Psi"):
        if f not in a_clim or f not in b_clim:
            continue
        a, b = a_clim[f], b_clim[f]
        n = min(len(a), len(b))
        d = np.flatnonzero(a[:n] != b[:n])
        if d.size:
            days.append(int(d[0]))
    return min(days) if days else None


def first_output_diff(a_ser, b_ser):
    """First day whose end state differs in any series (-1 = initial state), or None."""
    first = None
    for f in SERIES:
        if f not in a_ser or f not in b_ser:
            continue
        a, b = a_ser[f].align(b_ser[f], fill_value=0.0)
        a = a[sorted(a.columns, key=int)]
        b = b[a.columns]
        d = np.flatnonzero((a.to_numpy() != b.to_numpy()).any(axis=0))
        if d.size:
            # column labels are absolute output days (label k = after day k-1),
            # also in split branched output, which starts at the branch day
            day = int(a.columns[d[0]]) - 1
            first = day if first is None else min(first, day)
    return first


def check_community(args):
    """All scenarios of one community against its baseline run (one worker task).

    The baseline run is read once; in split branched output the Initialization
    folder (Temp/Psi) lives in the shared spin-up folder, so the climate files
    are looked up there when the run folder has none.
    """
    run, scenarios = args
    lines, n_ok, n_bad = [], 0, 0
    base_clim = _load_climate(run)
    base_ser = _load_series(run)
    for sc in scenarios:
        other = sc / run.name
        if not other.is_dir():
            lines.append(f"{run.name:24s} {sc.name:24s} missing")
            continue
        c = first_climate_diff(base_clim, _load_climate(other))
        o = first_output_diff(base_ser, _load_series(other))
        if o is None:
            ok, msg = True, "output identical throughout"
        elif c is None:
            ok, msg = False, f"output differs from day {o} but the climate is identical"
        elif o < c:
            ok, msg = False, f"output differs from day {o}, before the climate does (day {c})"
        else:
            ok, msg = True, f"identical until the climate diverges on day {c}; output from day {o}"
        n_ok += ok
        n_bad += not ok
        lines.append(f"{run.name:24s} {sc.name:24s} {'PAIRED  ' if ok else 'UNPAIRED'} {msg}")
    return lines, n_ok, n_bad


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("output_dir")
    ap.add_argument("--baseline", default="climate_ambient_temp")
    ap.add_argument("--jobs", type=int, default=0,
                    help="worker processes (default: $SLURM_CPUS_PER_TASK, else all cores)")
    a = ap.parse_args()
    root = Path(a.output_dir)
    base = root / a.baseline
    if not base.is_dir():
        sys.exit(f"baseline folder {base} not found")
    # scenario folders; skip the baseline and the spin-up states of branched runs
    scenarios = sorted(p for p in root.iterdir()
                       if p.is_dir() and p.name not in (a.baseline, "spinup"))
    runs = sorted(p for p in base.iterdir() if p.is_dir())
    jobs = a.jobs or int(os.environ.get("SLURM_CPUS_PER_TASK", 0)) or os.cpu_count() or 1
    tasks = [(run, scenarios) for run in runs]
    n_bad = n_ok = 0
    if jobs > 1 and len(tasks) > 1:
        from multiprocessing import Pool
        it = Pool(jobs).imap(check_community, tasks, chunksize=1)
    else:
        it = map(check_community, tasks)
    for lines, ok, bad in it:
        n_ok += ok
        n_bad += bad
        print("\n".join(lines), flush=True)
    print(f"\n{n_ok} paired, {n_bad} not paired")
    sys.exit(1 if n_bad else 0)


if __name__ == "__main__":
    main()
