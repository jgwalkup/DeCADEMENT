"""
Maintenance-energy sensitivity with the energy-deficit erosion switch off and on (minimal scenario).

Sweeps Maint_E_Rate_Const as multiples of its baseline at flat temperatures, for E_deficit_erosion 0 and 1, several
seeds; same design as the stage-A re-run (analysis_output/stageA_rerun: 20x20, 365 days, 4/15/28/36 C, 15 multiples,
3 seeds) but HSP production as in the input folder (no no-HSP control arm).

Usage:
    python maint_erosion_sweep.py OUT_DIR [--input ../input_minimal_e1500_crich_km25_in] [--param Maint_E_Rate_Const]
        [--base 2.5] [--mult 0.2 0.4 ...] [--temps 4 15 28 36] [--seeds 0 1 2] [--side 20] [--days 365] [--procs 2]
Writes OUT_DIR/summary.csv (one row per run: end-of-run summary from sweep_minimal.run_one) and
OUT_DIR/timeseries.csv.gz (daily series of every run). Runs already in summary.csv are skipped (restartable).
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

MULT = [0.2, 0.4, 0.6, 0.8, 1, 1.25, 1.5, 1.75, 2, 2.5, 3, 3.5, 4, 4.5, 5]


def _job(args):
    src, param, base, mult, temp, erosion, seed, side, days = args
    summary, series = run_one(src, {param: base * mult, "E_deficit_erosion": float(erosion)}, side, days, seed, temp=temp)
    key = dict(param=param, mult=mult, value=base * mult, Temp=temp, erosion=erosion, seed=seed)
    return {**key, **summary}, pd.DataFrame(series).assign(**key)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("out_dir")
    ap.add_argument("--input", default=os.path.join(HERE, "..", "input_minimal_e1500_crich_km25_in"))
    ap.add_argument("--param", default="Maint_E_Rate_Const")
    ap.add_argument("--base", type=float, default=2.5)
    ap.add_argument("--mult", type=float, nargs="+", default=MULT)
    ap.add_argument("--temps", type=float, nargs="+", default=[4, 15, 28, 36])
    ap.add_argument("--seeds", type=int, nargs="+", default=[0, 1, 2])
    ap.add_argument("--erosion", type=int, nargs="+", default=[0, 1])
    ap.add_argument("--side", type=int, default=20)
    ap.add_argument("--days", type=int, default=365)
    ap.add_argument("--procs", type=int, default=2)
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    sp, tp = os.path.join(a.out_dir, "summary.csv"), os.path.join(a.out_dir, "timeseries.csv.gz")
    done = set()
    if os.path.exists(sp):
        d = pd.read_csv(sp)
        done = set(zip(d["mult"].round(6), d["Temp"], d["erosion"], d["seed"]))
    src = os.path.abspath(a.input)
    jobs = [(src, a.param, a.base, m, t, e, s, a.side, a.days)
            for t, e, s, m in itertools.product(a.temps, a.erosion, a.seeds, a.mult) if (round(m, 6), t, e, s) not in done]
    print(f"{len(jobs)} runs ({len(done)} already done) -> {a.out_dir}", flush=True)
    with Pool(a.procs) as pool:
        for k, (row, ts) in enumerate(pool.imap_unordered(_job, jobs), 1):
            pd.DataFrame([row]).to_csv(sp, mode="a", header=not os.path.exists(sp), index=False)
            ts.to_csv(tp, mode="a", header=not os.path.exists(tp), index=False, compression="gzip")
            print(f"  {k}/{len(jobs)} T={row['Temp']:g} erosion={row['erosion']} mult={row['mult']:g} seed={row['seed']}", flush=True)


if __name__ == "__main__":
    main()
