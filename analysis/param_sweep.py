"""One-at-a-time parameter sweep of the minimal model under the scenario climates.

    python analysis/param_sweep.py OUT_DIR --vary NAME=v1,v2,... [NAME=...] [--input DIR] [--climates CSV ...]
        [--seeds 0 1] [--days 365] [--side 20] [--set NAME=VALUE ...] [--runtime NAME=VALUE ...] [--procs 2]

Each parameter is varied alone, all others at the input values (NAME may join several parameters with "+", e.g.
Enz_Prod_min+Enz_Prod_max, all set to the same value); the baseline (input values) is run once and stored
under every varied parameter at its input value. R_endo is on unless --set R_endo=0. Writes OUT_DIR/summary.csv and
OUT_DIR/timeseries.csv.gz (columns param, value, climate, seed); restartable.
"""
import argparse
import itertools
import os
import sys
from multiprocessing import Pool

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from sweep_minimal import run_one  # noqa: E402

CLIM_DIR = os.path.join(HERE, "..", "input", "climate_data")
DEFAULT_CLIMATES = [os.path.join(CLIM_DIR, f) for f in
                    ("climate_ambient_temp.csv", "climate_moderate_warm.csv", "climate_extreme_warm.csv")]


def _name(path):
    return os.path.splitext(os.path.basename(path))[0].replace("climate_", "")


def _job(args):
    src, clim, param, value, seed, side, days, fixed, rt = args
    changes = dict(fixed) if param == "baseline" else {**fixed, **{q: value for q in param.split("+")}}
    summary, series = run_one(src, changes, side, days, seed, climate=clim, runtime=rt)
    key = dict(param=param, value=value, climate=_name(clim), seed=seed)
    return {**key, **summary}, pd.DataFrame(series).assign(**key)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("out_dir")
    ap.add_argument("--vary", nargs="+", required=True, metavar="NAME=v1,v2,...")
    ap.add_argument("--input", default=os.path.join(HERE, "..", "input_minimal_e1500_crich_km25_in"))
    ap.add_argument("--climates", nargs="+", default=DEFAULT_CLIMATES)
    ap.add_argument("--seeds", type=int, nargs="+", default=[0, 1])
    ap.add_argument("--days", type=int, default=365)
    ap.add_argument("--side", type=int, default=20)
    ap.add_argument("--set", nargs="*", default=["R_endo=1"], metavar="NAME=VALUE")
    ap.add_argument("--runtime", nargs="*", default=[], metavar="NAME=VALUE",
                    help="runtime.txt rows to set for every run (appended if absent), e.g. Decay_MM_Form=1")
    ap.add_argument("--procs", type=int, default=2)
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    src = os.path.abspath(a.input)
    params = pd.read_csv(os.path.join(src, "parameters.csv"), header=None, index_col=0)[1]
    fixed = {k: float(v) for k, v in (x.split("=") for x in a.set)}
    vary = {k: [float(x) for x in v.split(",")] for k, v in (x.split("=") for x in a.vary)}
    base = {k: float(fixed.get(k.split("+")[0], params[k.split("+")[0]])) for k in vary}
    pv = [("baseline", 0.0)] + [(k, v) for k, vals in vary.items() for v in vals if v != base[k]]
    sp, tp = os.path.join(a.out_dir, "summary.csv"), os.path.join(a.out_dir, "timeseries.csv.gz")
    done = set()
    if os.path.exists(sp):
        d = pd.read_csv(sp)
        done = set(zip(d["param"], d["value"].astype(float), d["climate"], d["seed"]))
    rt = dict(x.split("=") for x in a.runtime)
    jobs = [(src, os.path.abspath(c), p, float(v), s, a.side, a.days, fixed, rt)
            for (p, v), c, s in itertools.product(pv, a.climates, a.seeds) if (p, float(v), _name(c), s) not in done]
    pd.Series(base).to_csv(os.path.join(a.out_dir, "baseline_values.csv"), header=False)
    print(f"{len(jobs)} runs ({len(done)} already done) -> {a.out_dir}; baseline {base}; fixed {fixed}; runtime {rt}",
          flush=True)
    with Pool(a.procs) as pool:
        for k, (row, ts) in enumerate(pool.imap_unordered(_job, jobs), 1):
            pd.DataFrame([row]).to_csv(sp, mode="a", header=not os.path.exists(sp), index=False)
            ts.to_csv(tp, mode="a", header=not os.path.exists(tp), index=False, compression="gzip")
            print(f"  {k}/{len(jobs)} {row['param']}={row['value']} {row['climate']} seed={row['seed']}", flush=True)


if __name__ == "__main__":
    main()
