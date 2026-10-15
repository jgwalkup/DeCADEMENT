"""
Decay of each substrate pool, including the necromass pools DeadMic and DeadEnz, and the matching monomer pools, day by day.

Usage: python necromass_check.py OUT.csv ../analysis_output/climate_12_years.csv [--input ../input_e1500_enz2_km25_osmorange2] [--pulses 2]
    [--side 20] [--seeds 0 1 2] [--sub-n-mult 1] [--procs 3]
Per seed, pulse, day and substrate: mass (C+N+P per grid cell), decay (mass degraded per cell per day, from Grid.DecayRates), the pool of its
monomer (Mon_<substrate>), and Psi. The first-order rate constant is decay / mass.
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
from full_climate_run import CYCLE, _StubOutput  # noqa: E402
from grid import Grid  # noqa: E402
from initialization import initialize_data  # noqa: E402
from sweep_minimal import set_params, set_runtime  # noqa: E402


def run(job):
    src, climate, pulses, seed, side, nmult = job
    tmp = tempfile.mkdtemp()
    w = os.path.join(tmp, "in")
    shutil.copytree(src, w)
    shutil.copy(climate, os.path.join(w, "climate.csv"))
    set_runtime(os.path.join(w, "runtime.txt"), dict(x=side, y=side, gridsize=side * side, pulse=pulses, end_time=CYCLE))
    set_params(os.path.join(w, "parameters.csv"), {"Osmo_Induci_Prod_min": 0.04, "Osmo_Induci_Prod_max": 0.08,
                                                   "Osmo_Consti_Prod_min": 0.002, "Osmo_Consti_Prod_max": 0.005})
    if nmult != 1.0:
        sp = os.path.join(w, "initial_substrates.csv")
        sdf = pd.read_csv(sp, index_col=0)
        sdf["N"] = sdf["N"] * nmult
        sdf.to_csv(sp)
    cwd = os.getcwd()
    os.chdir(w)
    rows = []
    try:
        np.random.seed(seed)
        rt = pd.read_csv("runtime.txt", header=None, index_col=0, sep="\t").drop("output_format").astype(float)
        eco = Grid(rt, initialize_data(rt))
        G, T = side * side, int(rt.loc["n_taxa", 1])
        stub = _StubOutput(T, pulses * CYCLE, pd.unique(eco.Microbes.index))
        subs = list(pd.unique(eco.Substrates.index))
        sub_idx = np.asarray(eco.Substrates.index)
        mon_idx = np.asarray(eco.Monomers.index)
        for pl in range(pulses):
            for i in range(pl * CYCLE, (pl + 1) * CYCLE):
                eco.degradation(pl, i)
                dec = np.asarray(eco.DecayRates, dtype=np.float64)
                S = eco.Substrates.to_numpy(dtype=np.float64).sum(axis=1)
                eco.uptake(pl, i)
                M = eco.Monomers.to_numpy(dtype=np.float64).sum(axis=1)
                eco.metabolism(i)
                eco.mortality(i)
                eco.reproduction(i)
                for s in subs:
                    ms = sub_idx == s
                    mm = mon_idx == "Mon_" + s
                    rows.append(dict(seed=seed, pulse=pl + 1, doy=i - pl * CYCLE + 1, substrate=s, mass=S[ms].sum() / G,
                                     decay=dec[ms].sum() / G, monomer_after_uptake=M[mm].sum() / G if mm.any() else np.nan,
                                     psi=float(eco.psi[i])))
                stub.Taxon_count_repop.iloc[:, i + 1] = (eco.Microbes["C"].to_numpy().reshape(G, T) > 0).sum(axis=0)
                if i == (pl + 1) * CYCLE - 1 and pl < pulses - 1:
                    eco.repopulation(stub, i, True)
    finally:
        os.chdir(cwd)
        shutil.rmtree(tmp, ignore_errors=True)
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("out_csv")
    ap.add_argument("climate")
    ap.add_argument("--input", default=os.path.join(HERE, "..", "input_e1500_enz2_km25_osmorange2"))
    ap.add_argument("--pulses", type=int, default=2)
    ap.add_argument("--side", type=int, default=20)
    ap.add_argument("--seeds", type=int, nargs="+", default=[0, 1, 2])
    ap.add_argument("--sub-n-mult", type=float, default=1.0)
    ap.add_argument("--procs", type=int, default=3)
    a = ap.parse_args()
    jobs = [(os.path.abspath(a.input), os.path.abspath(a.climate), a.pulses, s, a.side, a.sub_n_mult) for s in a.seeds]
    with Pool(min(len(jobs), a.procs)) as p:
        parts = p.map(run, jobs, chunksize=1)
    pd.concat(parts, ignore_index=True).to_csv(a.out_csv, index=False)
    print("wrote", a.out_csv)


if __name__ == "__main__":
    main()
