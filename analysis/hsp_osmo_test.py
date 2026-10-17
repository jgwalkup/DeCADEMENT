"""HSP and osmolyte effects and production on/off (minimal model under the scenario climates).

    python analysis/hsp_osmo_test.py OUT_DIR [--input DIR] [--climates CSV ...] [--seeds 0 1] [--days 365] [--procs 2]

Variants: baseline; HSP effect off (Maint_E_HSP_Ea_Slope 0: HSP no longer lowers the maintenance Ea; production and
costs kept); osmolyte effect off (Osmo_tol_K 1e9: drought tolerance about 0; production and costs kept); both effects
off; HSP production off (HSP_Consti_Prod = HSP_Induci_Prod = 0); osmolyte production off (Osmo_Consti_Prod =
Osmo_Induci_Prod = 0); both production off. R_endo on. Writes OUT_DIR/summary.csv and timeseries.csv.gz (columns
variant, climate, seed); restartable.
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
VARIANTS = {
    "baseline": {},
    "HSP effect off": {"Maint_E_HSP_Ea_Slope": 0.0},
    "osmolyte effect off": {"Osmo_tol_K": 1e9},
    "both effects off": {"Maint_E_HSP_Ea_Slope": 0.0, "Osmo_tol_K": 1e9},
    "HSP production off": {"HSP_Consti_Prod": 0.0, "HSP_Induci_Prod": 0.0},
    "osmolyte production off": {"Osmo_Consti_Prod": 0.0, "Osmo_Induci_Prod": 0.0},
    "both production off": {"HSP_Consti_Prod": 0.0, "HSP_Induci_Prod": 0.0,
                            "Osmo_Consti_Prod": 0.0, "Osmo_Induci_Prod": 0.0},
}


def _name(path):
    return os.path.splitext(os.path.basename(path))[0].replace("climate_", "")


def _job(args):
    src, clim, var, seed, days = args
    summary, series = run_one(src, {"R_endo": 1.0, **VARIANTS[var]}, 20, days, seed, climate=clim)
    key = dict(variant=var, climate=_name(clim), seed=seed)
    return {**key, **summary}, pd.DataFrame(series).assign(**key)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("out_dir")
    ap.add_argument("--input", default=os.path.join(HERE, "..", "input_minimal_e1500_crich_km25_in"))
    ap.add_argument("--climates", nargs="+", default=DEFAULT_CLIMATES)
    ap.add_argument("--seeds", type=int, nargs="+", default=[0, 1])
    ap.add_argument("--days", type=int, default=365)
    ap.add_argument("--procs", type=int, default=2)
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    sp, tp = os.path.join(a.out_dir, "summary.csv"), os.path.join(a.out_dir, "timeseries.csv.gz")
    done = set()
    if os.path.exists(sp):
        d = pd.read_csv(sp)
        done = set(zip(d["variant"], d["climate"], d["seed"]))
    src = os.path.abspath(a.input)
    jobs = [(src, os.path.abspath(c), v, s, a.days) for v, c, s in itertools.product(VARIANTS, a.climates, a.seeds)
            if (v, _name(c), s) not in done]
    print(f"{len(jobs)} runs ({len(done)} already done) -> {a.out_dir}", flush=True)
    with Pool(a.procs) as pool:
        for k, (row, ts) in enumerate(pool.imap_unordered(_job, jobs), 1):
            pd.DataFrame([row]).to_csv(sp, mode="a", header=not os.path.exists(sp), index=False)
            ts.to_csv(tp, mode="a", header=not os.path.exists(tp), index=False, compression="gzip")
            print(f"  {k}/{len(jobs)} {row['variant']} {row['climate']} seed={row['seed']}", flush=True)


if __name__ == "__main__":
    main()
