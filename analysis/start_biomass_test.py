"""
Single-pulse tests of what drives the end-of-pulse energy:C: starting biomass and its composition (real climate, 365 days).

Conditions (each applied to the same communities, 20x20 by default):
  base       as initialised (pulse-1 start: about 4 bacterial rows per cell, fungi 17 % of the C)
  x3, x10    every microbe's C, N, P and E multiplied by 3 or 10 (same composition, more biomass)
  bac1       one bacterial taxon per cell (frequency-weighted uniformly over bacteria), initial C = 1   (like a reseeded cell)
  bac25      same but every cell scaled to C = 25 (as large as a fungal cell)
  fun25      one fungal taxon per cell, initial C = 25 (as fungal cells start in a reseeded pulse)

Usage: python start_biomass_test.py OUT.csv ../input/climate.csv [--input ../input_e1500_enz2_km25_osmorange2] [--seeds 0 1 2 3 4] [--side 20] [--procs 6]
"""
import argparse
import os
import shutil
import sys
import tempfile
from multiprocessing import Pool

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "src"))
sys.path.insert(0, HERE)
from grid import Grid  # noqa: E402
from initialization import initialize_data  # noqa: E402
from sweep_minimal import set_runtime  # noqa: E402

CONDS = ["base", "x3", "x10", "bac1", "bac25", "fun25"]
DAYS = 365


def run(job):
    src, climate, cond, seed, side = job
    tmp = tempfile.mkdtemp()
    w = os.path.join(tmp, "in")
    shutil.copytree(src, w)
    shutil.copy(climate, os.path.join(w, "climate.csv"))
    set_runtime(os.path.join(w, "runtime.txt"), dict(x=side, y=side, gridsize=side * side, pulse=1, end_time=DAYS))
    cwd = os.getcwd()
    os.chdir(w)
    rows = []
    try:
        np.random.seed(seed)
        rt = pd.read_csv("runtime.txt", header=None, index_col=0, sep="\t").drop("output_format").astype(float)
        eco = Grid(rt, initialize_data(rt))
        G, T = side * side, int(rt.loc["n_taxa", 1])
        fb = np.asarray(eco.fb[:T])
        M = eco.Microbes_init.to_numpy(dtype=np.float32).reshape(G, T, -1)
        cur = eco.Microbes.to_numpy(dtype=np.float32).reshape(G, T, -1).copy()
        rng = np.random.default_rng(1000 + seed)
        if cond == "x3":
            cur *= 3
        elif cond == "x10":
            cur *= 10
        elif cond in ("bac1", "bac25", "fun25"):
            pool = np.where(fb == (1 if cond == "fun25" else 0))[0]
            cur = np.zeros_like(cur)
            pick = rng.choice(pool, G)
            for g in range(G):
                cur[g, pick[g]] = M[g, pick[g]]
            if cond == "bac25":
                cur *= 25
        eco.Microbes[:] = cur.reshape(G * T, -1)
        subs = [s for s in eco.Substrates.index.get_level_values(0).unique() if s not in ("DeadMic", "DeadEnz")]
        tot = lambda: float(sum(np.asarray(eco.Substrates.loc[s, "C"]).sum() for s in subs))
        sub0 = tot()
        for i in range(DAYS):
            eco.degradation(0, i)
            eco.uptake(0, i)
            eco.metabolism(i)
            eco.mortality(i)
            eco.reproduction(i)
            C = np.asarray(eco.Microbes["C"])
            E = np.asarray(eco.Microbes["E"])
            alive = C > 0
            rows.append(dict(cond=cond, seed=seed, day=i + 1, biomass_C=C.sum() / G, alive_rows=int(alive.sum()),
                             energy_to_C=E[alive].sum() / C[alive].sum() if alive.any() else np.nan,
                             remaining_pct=100 * tot() / sub0, CUE=float(eco.CUE_system),
                             fungal_C_share=float(C.reshape(G, T)[:, fb == 1].sum() / max(C.sum(), 1e-9))))
    finally:
        os.chdir(cwd)
        shutil.rmtree(tmp, ignore_errors=True)
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("out_csv")
    ap.add_argument("climate")
    ap.add_argument("--input", default=os.path.join(HERE, "..", "input_e1500_enz2_km25_osmorange2"))
    ap.add_argument("--seeds", type=int, nargs="+", default=[0, 1, 2, 3, 4])
    ap.add_argument("--side", type=int, default=20)
    ap.add_argument("--conds", nargs="+", default=CONDS)
    ap.add_argument("--procs", type=int, default=6)
    a = ap.parse_args()
    jobs = [(os.path.abspath(a.input), os.path.abspath(a.climate), c, s, a.side) for c in a.conds for s in a.seeds]
    with Pool(min(len(jobs), a.procs)) as p:
        parts = p.map(run, jobs, chunksize=1)
    pd.concat(parts, ignore_index=True).to_csv(a.out_csv, index=False)
    print("wrote", a.out_csv)


if __name__ == "__main__":
    main()
