"""Check that the climate scenarios of each community are paired with the baseline run.

    python analysis/check_pairing.py OUTPUT_DIR [--baseline climate_ambient_temp]

OUTPUT_DIR is the --output-dir of run_multi_climate_simulations.py (one subfolder per scenario, each holding
<community>_s<seed>/ run folders). For every community and every non-baseline scenario, the script finds:

  - the first day on which the scenario's temperature or water potential differs from the baseline's
    (Initialization/Temp.csv and Psi.csv, which each run writes);
  - the first day on which the model output differs (SubstratesSeries, MicrobesSeries, RespSeries).

Series columns are the initial state followed by one column per day, so column k is the state after day k - 1.
A paired run's output matches the baseline exactly until its climate differs. Output that differs before the
climate does means the runs did not start from the same random state, and paired comparisons between them are
not valid.
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

SERIES = ["SubstratesSeries", "MicrobesSeries", "RespSeries"]


def first_climate_diff(a_run, b_run):
    days = []
    for f in ("Temp", "Psi"):
        pa, pb = a_run / "Initialization" / f"{f}.csv", b_run / "Initialization" / f"{f}.csv"
        if not (pa.exists() and pb.exists()):
            continue
        a = pd.read_csv(pa, index_col=0).iloc[:, 0].to_numpy()
        b = pd.read_csv(pb, index_col=0).iloc[:, 0].to_numpy()
        n = min(len(a), len(b))
        d = np.flatnonzero(a[:n] != b[:n])
        if d.size:
            days.append(int(d[0]))
    return min(days) if days else None


def first_output_diff(a_run, b_run):
    """First day whose end state differs in any series (-1 = initial state), or None."""
    first = None
    for f in SERIES:
        pa, pb = a_run / f"{f}.csv", b_run / f"{f}.csv"
        if not (pa.exists() and pb.exists()):
            continue
        a = pd.read_csv(pa, index_col=0).sort_index()
        b = pd.read_csv(pb, index_col=0).sort_index()
        a, b = a.align(b, fill_value=0.0)
        d = np.flatnonzero((a.to_numpy() != b.to_numpy()).any(axis=0))
        if d.size:
            day = int(d[0]) - 1
            first = day if first is None else min(first, day)
    return first


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("output_dir")
    ap.add_argument("--baseline", default="climate_ambient_temp")
    a = ap.parse_args()
    root = Path(a.output_dir)
    base = root / a.baseline
    if not base.is_dir():
        sys.exit(f"baseline folder {base} not found")
    scenarios = sorted(p for p in root.iterdir() if p.is_dir() and p.name != a.baseline)
    n_bad = n_ok = 0
    for run in sorted(p for p in base.iterdir() if p.is_dir()):
        for sc in scenarios:
            other = sc / run.name
            if not other.is_dir():
                print(f"{run.name:24s} {sc.name:24s} missing")
                continue
            c, o = first_climate_diff(run, other), first_output_diff(run, other)
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
            print(f"{run.name:24s} {sc.name:24s} {'PAIRED  ' if ok else 'UNPAIRED'} {msg}")
    print(f"\n{n_ok} paired, {n_bad} not paired")
    sys.exit(1 if n_bad else 0)


if __name__ == "__main__":
    main()
