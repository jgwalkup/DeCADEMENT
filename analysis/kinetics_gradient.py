"""MMRT vs Arrhenius kinetics across a constant-temperature gradient (minimal model, flat climate file).

    python analysis/kinetics_gradient.py OUT_DIR [--input DIR] [--temps 0 5 ... 40] [--kinetics 1 0] [--seeds 0 1]
        [--days 365] [--side 20] [--set NAME=VALUE ...] [--procs 2]

The temperature of the input's flat climate.csv is replaced by each value (Psi as in the file); Kinetics_Model is set
in runtime.txt (1 = MMRT, 0 = Arrhenius). R_endo on unless --set R_endo=0. Writes OUT_DIR/summary.csv and
OUT_DIR/timeseries.csv.gz (columns kinetics, Temp_set, seed); restartable.
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

KIN = {1: "MMRT", 0: "Arrhenius"}


def _job(args):
    src, T, k, seed, side, days, fixed = args
    summary, series = run_one(src, fixed, side, days, seed, temp=T, runtime={"Kinetics_Model": k})
    key = dict(kinetics=KIN[k], Temp_set=T, seed=seed)
    return {**key, **summary}, pd.DataFrame(series).assign(**key)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("out_dir")
    ap.add_argument("--input", default=os.path.join(HERE, "..", "input_minimal_e1500_crich_km25_in"))
    ap.add_argument("--temps", type=float, nargs="+", default=[0, 5, 10, 15, 20, 25, 30, 35, 40])
    ap.add_argument("--kinetics", type=int, nargs="+", default=[1, 0])
    ap.add_argument("--seeds", type=int, nargs="+", default=[0, 1])
    ap.add_argument("--days", type=int, default=365)
    ap.add_argument("--side", type=int, default=20)
    ap.add_argument("--set", nargs="*", default=["R_endo=1"], metavar="NAME=VALUE")
    ap.add_argument("--procs", type=int, default=2)
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    sp, tp = os.path.join(a.out_dir, "summary.csv"), os.path.join(a.out_dir, "timeseries.csv.gz")
    done = set()
    if os.path.exists(sp):
        d = pd.read_csv(sp)
        done = set(zip(d["kinetics"], d["Temp_set"].astype(float), d["seed"]))
    fixed = {k: float(v) for k, v in (x.split("=") for x in a.set)}
    src = os.path.abspath(a.input)
    jobs = [(src, float(T), k, s, a.side, a.days, fixed) for k, T, s in itertools.product(a.kinetics, a.temps, a.seeds)
            if (KIN[k], float(T), s) not in done]
    print(f"{len(jobs)} runs ({len(done)} already done) -> {a.out_dir}; fixed {fixed}", flush=True)
    with Pool(a.procs) as pool:
        for n, (row, ts) in enumerate(pool.imap_unordered(_job, jobs), 1):
            pd.DataFrame([row]).to_csv(sp, mode="a", header=not os.path.exists(sp), index=False)
            ts.to_csv(tp, mode="a", header=not os.path.exists(tp), index=False, compression="gzip")
            print(f"  {n}/{len(jobs)} {row['kinetics']} {row['Temp_set']:g} C seed={row['seed']}", flush=True)


if __name__ == "__main__":
    main()
