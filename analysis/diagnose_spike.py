"""
What drives the burst of energy production at the rewetting around day 295-300 of the real-climate minimal-model runs?
Runs the minimal model (default: bacterium, C 100 substrate) for 365 days under the real climate and under counterfactual climates, and records per day:
  Psi, T, the moisture factors on enzyme Vmax (exp(0.25 (Psi - wp_fc))) and on uptake (exp(0.5 (Psi - wp_fc))), the C released by decay, the monomer C pool that
  uptake sees (left after uptake + taken), the C and energy taken up, the biomass C and the energy:C of the cells.
Climates: real | psi_flat (Psi -0.5, real T) | t_flat (T 15, real Psi) | no_rewet (Psi held at the day-285 value from day 285, real T) | t_flat_late (real Psi, T 15 from day 285)
Usage: python diagnose_spike.py OUT.csv [--input ../input_minimal_e1500_km25] [--seeds 0 1 2]
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

CLIMATES = ["real", "psi_flat", "t_flat", "no_rewet", "t_flat_late"]


def make_climate(kind):
    c = pd.read_csv(os.path.join(HERE, "..", "input", "climate.csv"), index_col=0).iloc[:365].reset_index(drop=True)
    if kind == "psi_flat":
        c["Psi"] = -0.5
    elif kind == "t_flat":
        c["Temp"] = 15.0
    elif kind == "no_rewet":
        c.loc[284:, "Psi"] = float(c.loc[284, "Psi"])
    elif kind == "t_flat_late":
        c.loc[284:, "Temp"] = 15.0
    return c


def run(job):
    src, kind, seed = job
    tmp = tempfile.mkdtemp()
    w = os.path.join(tmp, "in")
    shutil.copytree(src, w)
    clim = make_climate(kind)
    pd.concat([clim, clim.iloc[:5]], ignore_index=True).to_csv(os.path.join(w, "climate.csv"))
    set_runtime(os.path.join(w, "runtime.txt"), dict(x=20, y=20, gridsize=400, pulse=1, end_time=365))
    cwd = os.getcwd()
    os.chdir(w)
    rows = []
    try:
        np.random.seed(seed)
        rt = pd.read_csv("runtime.txt", header=None, index_col=0, sep="\t").drop("output_format").astype(float)
        eco = Grid(rt, initialize_data(rt))
        G, T = 400, 1
        mon_idx = np.asarray(eco.Monomers.index)
        is_sub = mon_idx == "Mon_IdealSub"
        for i in range(365):
            eco.degradation(0, i)
            decay_C = float(np.asarray(eco.DecayRates, dtype=np.float64).sum() / G)
            eco.uptake(0, i)
            Mc = eco.Monomers.to_numpy(dtype=np.float64)[:, 0]
            upC = float(np.asarray(eco.Taxon_Uptake_C, dtype=np.float64).sum() / G)
            upE = float(np.asarray(eco.Taxon_Uptake_E, dtype=np.float64).sum() / G)
            mon_left = float(Mc[is_sub].sum() / G)
            C0 = float(eco.Microbes["C"].sum() / G)
            eco.metabolism(i)
            eco.mortality(i)
            eco.reproduction(i)
            C = eco.Microbes["C"].to_numpy(dtype=np.float64)
            E = eco.Microbes["E"].to_numpy(dtype=np.float64)
            psi, temp = float(eco.psi[i]), float(eco.temp[i])
            fe = np.exp(0.25 * (psi - eco.wp_fc)) if psi < eco.wp_fc else 1.0
            fu = np.exp(0.5 * (psi - eco.wp_fc)) if psi < eco.wp_fc else 1.0
            rows.append(dict(kind=kind, seed=seed, doy=i + 1, Psi=psi, Temp=temp, f_enzyme=fe, f_uptake=fu, decay_C=decay_C, mon_pool_before=mon_left + upC, mon_left=mon_left,
                             uptake_C=upC, uptake_E=upE, biomass_C=C0, E_per_biomass=upE / C0 if C0 > 0 else np.nan, uptake_C_per_biomass=upC / C0 if C0 > 0 else np.nan,
                             frac_pool_taken=upC / (mon_left + upC) if (mon_left + upC) > 0 else np.nan, EC=float(E.sum() / C.sum()) if C.sum() > 0 else np.nan))
    finally:
        os.chdir(cwd)
        shutil.rmtree(tmp, ignore_errors=True)
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("out_csv")
    ap.add_argument("--input", default=os.path.join(HERE, "..", "input_minimal_e1500_km25"))
    ap.add_argument("--seeds", type=int, nargs="+", default=[0, 1, 2])
    a = ap.parse_args()
    jobs = [(os.path.abspath(a.input), k, s) for k in CLIMATES for s in a.seeds]
    with Pool(min(len(jobs), 10)) as p:
        parts = p.map(run, jobs, chunksize=1)
    pd.concat(parts, ignore_index=True).to_csv(a.out_csv, index=False)
    print("wrote", a.out_csv)


if __name__ == "__main__":
    main()
