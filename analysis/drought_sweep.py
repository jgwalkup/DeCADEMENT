"""Drought-mortality sensitivity of the minimal model under the scenario climates.

One-at-a-time changes of the drought-mortality parameters (beta_bac/beta_fun: slope of mortality on water potential
below wp_fc; Osmo_tol_K: osmolyte share of biomass giving drought tolerance 0.5), first year, R_endo on by default.

    python analysis/drought_sweep.py OUT_DIR [--input DIR] [--climates CSV ...] [--seeds 0 1] [--procs 2]
        [--beta 10 5 2.5 1 0] [--tolk 0.05 0.02 0.01 0.005] [--set NAME=VALUE ...]

Writes OUT_DIR/summary.csv and OUT_DIR/timeseries.csv.gz (columns param, value, climate, seed); restartable.
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
DEFAULT_CLIMATES = [os.path.join(CLIM_DIR, f) for f in ("climate_ambient_temp.csv", "climate_extreme_warm.csv")]
BASE = dict(beta=10.0, Osmo_tol_K=0.05)


def _name(path):
    return os.path.splitext(os.path.basename(path))[0].replace("climate_", "")


def _changes(param, value):
    if param == "beta":
        return dict(beta_bac=value, beta_fun=value)
    return {param: value}


def _job(args):
    src, clim, param, value, seed, side, days, fixed = args
    summary, series = run_one(src, {**fixed, **_changes(param, value)}, side, days, seed, climate=clim)
    key = dict(param=param, value=value, climate=_name(clim), seed=seed)
    return {**key, **summary}, pd.DataFrame(series).assign(**key)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("out_dir")
    ap.add_argument("--input", default=os.path.join(HERE, "..", "input_minimal_e1500_crich_km25_in"))
    ap.add_argument("--climates", nargs="+", default=DEFAULT_CLIMATES)
    ap.add_argument("--seeds", type=int, nargs="+", default=[0, 1])
    ap.add_argument("--beta", type=float, nargs="+", default=[10, 5, 2.5, 1, 0])
    ap.add_argument("--tolk", type=float, nargs="+", default=[0.05, 0.02, 0.01, 0.005])
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
        done = set(zip(d["param"], d["value"].astype(float), d["climate"], d["seed"]))
    fixed = {k: float(v) for k, v in (x.split("=") for x in a.set)}
    src = os.path.abspath(a.input)
    # baseline once (as param "beta" at its baseline value), then the one-at-a-time changes
    pv = [("beta", v) for v in a.beta] + [("Osmo_tol_K", v) for v in a.tolk if v != BASE["Osmo_tol_K"]]
    jobs = [(src, os.path.abspath(c), p, float(v), s, a.side, a.days, fixed)
            for (p, v), c, s in itertools.product(pv, a.climates, a.seeds) if (p, float(v), _name(c), s) not in done]
    print(f"{len(jobs)} runs ({len(done)} already done) -> {a.out_dir}; fixed {fixed}", flush=True)
    with Pool(a.procs) as pool:
        for k, (row, ts) in enumerate(pool.imap_unordered(_job, jobs), 1):
            pd.DataFrame([row]).to_csv(sp, mode="a", header=not os.path.exists(sp), index=False)
            ts.to_csv(tp, mode="a", header=not os.path.exists(tp), index=False, compression="gzip")
            print(f"  {k}/{len(jobs)} {row['param']}={row['value']} {row['climate']} seed={row['seed']}", flush=True)


if __name__ == "__main__":
    main()
