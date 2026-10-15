"""
Short test runs of the full multi-taxon, multi-substrate model on a small grid under a flat climate.

Usage:
    python full_model_test.py OUT.csv --inputs ../input ../input_e1500 --temps 15 28 [--side 20] [--days 365] [--seed 0]
Each (input folder, temperature) is run in its own process; the daily summary of every run is written to OUT.csv.
Flat climate: Temp as given, Psi -0.5, FM 50 (taken from the minimal scenario's climate.csv structure).
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
from sweep_minimal import set_params, set_runtime  # noqa: E402


def run(job):
    src, temp, side, days, seed, mult, smult = job
    tmp = tempfile.mkdtemp()
    w = os.path.join(tmp, "in")
    shutil.copytree(src, w)
    set_runtime(os.path.join(w, "runtime.txt"), dict(x=side, y=side, gridsize=side * side, pulse=1, end_time=days))
    if smult != 1.0:  # scale the initial amount of every substrate (C, N, P)
        sp = os.path.join(w, "initial_substrates.csv")
        sdf = pd.read_csv(sp, index_col=0)
        (sdf * smult).to_csv(sp)
    if mult != 1.0:  # scale inducible and constitutive enzyme production rate ranges together
        pr = pd.read_csv(os.path.join(w, "parameters.csv"), header=None, index_col=0)[1].astype(float)
        names = ["Enz_Prod_min", "Enz_Prod_max", "Constit_Prod_min", "Constit_Prod_max"]
        set_params(os.path.join(w, "parameters.csv"), {n: float(pr[n]) * mult for n in names})
    clim = pd.read_csv(os.path.join(w, "climate.csv"), index_col=0)
    n = max(len(clim), days + 5)
    flat = pd.DataFrame({"Temp": float(temp), "FM": 50.0, "Psi": -0.5, "Date": "1/1/00"}, index=range(n))
    flat.to_csv(os.path.join(w, "climate.csv"))
    cwd = os.getcwd()
    os.chdir(w)
    rows = []
    try:
        np.random.seed(seed)
        rt = pd.read_csv("runtime.txt", header=None, index_col=0, sep="\t").drop("output_format").astype(float)
        data = initialize_data(rt)
        eco = Grid(rt, data)
        G = side * side
        subs = [s for s in data["Substrate_index"] if s not in ("DeadMic", "DeadEnz")]

        def sub_total():
            S = eco.Substrates
            return float(sum(np.asarray(S.loc[s, "C"]).sum() for s in subs))

        sub0 = sub_total()
        sub_c = lambda s: float(np.asarray(eco.Substrates.loc[s, "C"]).sum())
        sub0_each = {s: sub_c(s) for s in subs}
        ps = lambda x: float(np.asarray(x).sum())
        for i in range(days):
            eco.degradation(0, i)
            eco.uptake(0, i)
            eco.metabolism(i)
            eco.mortality(i)
            eco.reproduction(i)
            C = np.asarray(eco.Microbes["C"])
            E = np.asarray(eco.Microbes["E"])
            alive = C > 0
            n_alive = int(alive.sum())
            tax = np.asarray(eco.Microbes.index)[alive]
            rows.append(dict(
                input=os.path.basename(src), Temp=temp, seed=seed, enz_mult=mult, sub_mult=smult, day=i + 1,
                biomass_C=C.sum() / G, alive_rows=n_alive, taxa_alive=len(set(tax)),
                energy_to_C=E[alive].sum() / C[alive].sum() if n_alive else np.nan,
                CUE=float(eco.CUE_system), CAE=float(eco.CAE_system),
                substrate_remaining_pct=100 * sub_total() / sub0,
                HSP_prop=float(np.asarray(eco.HSP_prop)[alive].mean()) if n_alive else np.nan,
                maint_Ea=float(np.asarray(eco.maint_E_cost_Ea)[alive].mean()) if n_alive else np.nan,
                net_energy=ps(eco.Microbe_E_Gain) / max(n_alive, 1),
                E_cost_total=ps(eco.Microbe_E_Cost) / max(n_alive, 1),
                E_cost_maint=ps(eco.Tax_cell_maint_E_cost) / max(n_alive, 1),
                **{f"sub_{s}": 100 * sub_c(s) / sub0_each[s] if sub0_each[s] > 0 else np.nan for s in subs},
            ))
    finally:
        os.chdir(cwd)
        shutil.rmtree(tmp, ignore_errors=True)
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("out_csv")
    ap.add_argument("--inputs", nargs="+", required=True)
    ap.add_argument("--temps", type=float, nargs="+", default=[15.0])
    ap.add_argument("--side", type=int, default=20)
    ap.add_argument("--days", type=int, default=365)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--seeds", type=int, nargs="+", default=None, help="several seeds (overrides --seed)")
    ap.add_argument("--sub-mult", type=float, nargs="+", default=[1.0], help="multipliers for the initial substrate amounts")
    ap.add_argument("--enz-mult", type=float, nargs="+", default=[1.0], help="multipliers for Enz_Prod and Constit_Prod ranges")
    a = ap.parse_args()
    seeds = a.seeds if a.seeds else [a.seed]
    jobs = [(os.path.abspath(i), t, a.side, a.days, s, m, sm) for i in a.inputs for t in a.temps for m in a.enz_mult for sm in a.sub_mult for s in seeds]
    with Pool(min(len(jobs), os.cpu_count())) as p:
        parts = p.map(run, jobs, chunksize=1)
    pd.concat(parts, ignore_index=True).to_csv(a.out_csv, index=False)
    print("wrote", a.out_csv)


if __name__ == "__main__":
    main()
