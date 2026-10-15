"""
Net effect of HSP production vs temperature across a parameter grid.

For every (temperature, grid point) this runs the minimal scenario twice -- once with HSP
production on, once with HSP production switched off (control: HSP_Consti_Prod = HSP_Induci_Prod = 0)
-- and tabulates the difference. The control only depends on the grid parameters that matter without
HSP (e.g. Maint_E_Rate_Const, Maint_E_HSP_Ea_Intercept), so it is cached per such combination.

Usage:
    python sweep_hsp_effect.py ../input_minimal --temps 4 8 12 16 20 24 28 32 \\
        --grid HSP_prod_E_cost 200 600 1200 --grid Maint_E_HSP_Ea_Slope -30000 -60000 -90000 \\
        --set HSP_Induci_Prod=0.3 --set HSP_Consti_Prod=0.0003 --csv out.csv

Output columns: Temp, grid params, <metric> (HSP run), ref_<metric> (no-HSP control),
d_<metric> (HSP - control), rel_biomC (% change in biomass C vs control).
Runs use a process pool (run_one chdirs into a temp copy, so processes, not threads).
"""
import argparse
import itertools
import os
import sys
from multiprocessing import Pool

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sweep_minimal import run_one  # noqa: E402

# Parameters that have no effect when HSP production is off (so the control need not be re-run for them).
HSP_ONLY = {"HSP_prod_E_cost", "Maint_E_HSP_Ea_Slope", "Maint_E_HSP_Ea_n",
            "HSP_Induci_Prod", "HSP_Consti_Prod", "HSP_Consti_Ea", "HSP_Induci_Ea"}
OFF = {"HSP_Induci_Prod": 0.0, "HSP_Consti_Prod": 0.0}
METRICS = ["deg_pct", "biomC", "EC_end", "net_E", "CUE", "HSP_prop", "maint_Ea", "alive", "enz_pool"]


def _job(args):
    src, changes, side, days, seed, temp = args
    summary, _ = run_one(src, changes, side, days, seed, temp=temp)
    return summary


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("input_folder")
    ap.add_argument("--temps", type=float, nargs="+", required=True)
    ap.add_argument("--grid", nargs="+", action="append", default=[], metavar="NAME V1 V2 ...",
                    help="parameter followed by its values (repeatable; full product is run)")
    ap.add_argument("--set", action="append", default=[], metavar="NAME=VALUE")
    ap.add_argument("--side", type=int, default=20)
    ap.add_argument("--days", type=int, default=365)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--procs", type=int, default=os.cpu_count())
    ap.add_argument("--csv", required=True)
    a = ap.parse_args()

    src = os.path.abspath(a.input_folder)
    fixed = {k: float(v) for k, v in (s.split("=") for s in a.set)}
    names = [g[0] for g in a.grid]
    values = [[float(v) for v in g[1:]] for g in a.grid]
    combos = list(itertools.product(*values)) if names else [()]
    ref_names = [n for n in names if n not in HSP_ONLY]

    jobs, keys = [], []
    for t in a.temps:
        for c in combos:
            ch = {**fixed, **dict(zip(names, c))}
            jobs.append((src, ch, a.side, a.days, a.seed, t))
            keys.append(("hsp", t, c))
    ref_keys = {}
    for t in a.temps:
        for c in combos:
            rc = tuple(v for n, v in zip(names, c) if n in ref_names)
            k = (t, rc)
            if k not in ref_keys:
                ch = {**fixed, **{n: v for n, v in zip(names, c) if n in ref_names}, **OFF}
                ref_keys[k] = len(jobs)
                jobs.append((src, ch, a.side, a.days, a.seed, t))
                keys.append(("ref", t, rc))
    print(f"# {len(jobs)} runs ({len(ref_keys)} controls) on {a.procs} processes", file=sys.stderr, flush=True)

    with Pool(a.procs) as pool:
        results = pool.map(_job, jobs, chunksize=1)

    rows = []
    for i, (kind, t, c) in enumerate(keys):
        if kind != "hsp":
            continue
        rc = tuple(v for n, v in zip(names, c) if n in ref_names)
        h, r = results[i], results[ref_keys[(t, rc)]]
        row = {"Temp": t, **dict(zip(names, c))}
        for m in METRICS:
            row[m] = h[m]
            row["ref_" + m] = r[m]
            row["d_" + m] = h[m] - r[m]
        row["rel_biomC"] = 100 * (h["biomC"] - r["biomC"]) / r["biomC"]
        rows.append(row)
    pd.DataFrame(rows).to_csv(a.csv, index=False)
    print(f"wrote {a.csv} ({len(rows)} rows)", file=sys.stderr)


if __name__ == "__main__":
    main()
