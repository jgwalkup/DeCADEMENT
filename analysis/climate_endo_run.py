"""
Minimal model under daily climate forcing, endogenous respiration (R_endo, the E_deficit_erosion switch) off and on.

Each run is one pulse driven by the first --days days of a climate file (Temp and Psi), with the input folder otherwise
unchanged. Default: the three scenario files in ../input/climate_data (ambient, moderate warming +5 C, extreme warming
+10 C), bacterial minimal folder, 20x20, 365 days, seeds 0-2.

Usage:
    python climate_endo_run.py OUT_DIR [--input ../input_minimal_e1500_crich_km25_in]
        [--climates ../input/climate_data/climate_ambient_temp.csv ...] [--seeds 0 1 2] [--days 365] [--side 20]
        [--endo 0 1] [--set NAME=VALUE ...] [--procs 2]
Writes OUT_DIR/summary.csv and OUT_DIR/timeseries.csv.gz (columns climate, endo, seed + sweep_minimal.run_one output);
restartable.
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
    src, clim, endo, seed, side, days, fixed = args
    summary, series = run_one(src, {**fixed, "E_deficit_erosion": float(endo)}, side, days, seed, climate=clim)
    key = dict(climate=_name(clim), endo=endo, seed=seed)
    return {**key, **summary}, pd.DataFrame(series).assign(**key)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("out_dir")
    ap.add_argument("--input", default=os.path.join(HERE, "..", "input_minimal_e1500_crich_km25_in"))
    ap.add_argument("--climates", nargs="+", default=DEFAULT_CLIMATES)
    ap.add_argument("--seeds", type=int, nargs="+", default=[0, 1, 2])
    ap.add_argument("--endo", type=int, nargs="+", default=[0, 1])
    ap.add_argument("--days", type=int, default=365)
    ap.add_argument("--side", type=int, default=20)
    ap.add_argument("--set", nargs="*", default=[], metavar="NAME=VALUE")
    ap.add_argument("--procs", type=int, default=2)
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    sp, tp = os.path.join(a.out_dir, "summary.csv"), os.path.join(a.out_dir, "timeseries.csv.gz")
    done = set()
    if os.path.exists(sp):
        d = pd.read_csv(sp)
        done = set(zip(d["climate"], d["endo"], d["seed"]))
    fixed = {k: float(v) for k, v in (x.split("=") for x in a.set)}
    src = os.path.abspath(a.input)
    jobs = [(src, os.path.abspath(c), e, s, a.side, a.days, fixed) for c, e, s in itertools.product(a.climates, a.endo, a.seeds)
            if (_name(c), e, s) not in done]
    print(f"{len(jobs)} runs ({len(done)} already done) -> {a.out_dir}; fixed {fixed}", flush=True)
    with Pool(a.procs) as pool:
        for k, (row, ts) in enumerate(pool.imap_unordered(_job, jobs), 1):
            pd.DataFrame([row]).to_csv(sp, mode="a", header=not os.path.exists(sp), index=False)
            ts.to_csv(tp, mode="a", header=not os.path.exists(tp), index=False, compression="gzip")
            print(f"  {k}/{len(jobs)} {row['climate']} endo={row['endo']} seed={row['seed']}", flush=True)


if __name__ == "__main__":
    main()
