"""
Full model under the real climate with osmolyte production varied across taxa: original (relative) drought tolerance
vs absolute tolerance (Osmo_tol_ref), on identical communities.

Usage:
    python full_climate_osmo_range.py OUT_DIR ../input/climate.csv [--inputs ../input_e1500_enz2] [--pulses 6]
        [--seeds 0 1 2] [--side 20] [--induci MIN MAX] [--consti MIN MAX] [--conds 0:-1 0:0 0.055:0 0.095:0]
A condition is REF:LOWER = Osmo_tol_ref : Osmo_tol_lower. REF > 0 is absolute tolerance min(1, allocation / REF); REF 0 is
relative tolerance with lower bound LOWER (0 = default, negative = community minimum, the original method).
Each (condition, seed) is one run; all runs of a seed share the same random stream, so the communities and their
osmolyte traits are identical and only the tolerance mapping differs.
Writes OUT_DIR/timeseries.csv (daily community summary) and OUT_DIR/taxa.csv (per taxon, per pulse: allocation, tolerance,
cell count on day 1 and day 365).
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
from full_climate_run import CYCLE, _StubOutput  # noqa: E402


def run(job):
    src, climate, pulses, seed, side, induci, consti, (ref, lower) = job
    tmp = tempfile.mkdtemp()
    w = os.path.join(tmp, "in")
    shutil.copytree(src, w)
    shutil.copy(climate, os.path.join(w, "climate.csv"))
    set_runtime(os.path.join(w, "runtime.txt"), dict(x=side, y=side, gridsize=side * side, pulse=pulses, end_time=CYCLE))
    set_params(os.path.join(w, "parameters.csv"), {
        "Osmo_Induci_Prod_min": induci[0], "Osmo_Induci_Prod_max": induci[1],
        "Osmo_Consti_Prod_min": consti[0], "Osmo_Consti_Prod_max": consti[1]})
    pp = os.path.join(w, "parameters.csv")  # overwrite the tolerance rows if the folder has them, else append them
    txt = open(pp).read().rstrip("\n")
    have = {l.split(",")[0] for l in txt.split("\n") if "," in l}
    if {"Osmo_tol_ref", "Osmo_tol_lower"} <= have:
        set_params(pp, {"Osmo_tol_ref": ref, "Osmo_tol_lower": lower})
    else:
        assert not ({"Osmo_tol_ref", "Osmo_tol_lower"} & have), "only one of the tolerance rows is present"
        open(pp, "w").write(txt + f"\nOsmo_tol_ref,{ref}\nOsmo_tol_lower,{lower}\n")
    cond = f"abs{ref:g}" if ref > 0 else ("rel_min" if lower < 0 else f"rel_lower{lower:g}")
    clim = pd.read_csv(os.path.join(w, "climate.csv"), index_col=0)
    cwd = os.getcwd()
    os.chdir(w)
    rows, trows = [], []
    try:
        np.random.seed(seed)
        rt = pd.read_csv("runtime.txt", header=None, index_col=0, sep="\t").drop("output_format").astype(float)
        data = initialize_data(rt)
        eco = Grid(rt, data)
        G = side * side
        T = int(rt.loc["n_taxa", 1])
        taxa = pd.unique(eco.Microbes.index)
        stub = _StubOutput(T, pulses * CYCLE, taxa)
        # per-taxon allocation (constitutive + inducible osmolyte C fraction) and tolerance actually used by the model
        alloc = (np.asarray(data["OsmoProdConsti"].sum(axis=1)) + np.asarray(data["OsmoProdInduci"].sum(axis=1))).reshape(G, T)[0]
        tol = np.asarray(eco.tolerance, dtype=float).reshape(G, T)[0]
        subs = [s for s in data["Substrate_index"] if s not in ("DeadMic", "DeadEnz")]
        sub_total = lambda: float(sum(np.asarray(eco.Substrates.loc[s, "C"]).sum() for s in subs))
        sub0 = sub_total()
        ps = lambda x: float(np.asarray(x).sum())
        count0 = None
        for pl in range(pulses):
            for i in range(pl * CYCLE, (pl + 1) * CYCLE):
                eco.degradation(pl, i)
                eco.uptake(pl, i)
                eco.metabolism(i)
                eco.mortality(i)
                eco.reproduction(i)
                C = np.asarray(eco.Microbes["C"])
                E = np.asarray(eco.Microbes["E"])
                alive = C > 0
                n = int(alive.sum())
                cnt = alive.reshape(G, T).sum(axis=0)
                stub.Taxon_count_repop.iloc[:, i + 1] = cnt
                if i == pl * CYCLE:
                    count0 = cnt.copy()
                if i == (pl + 1) * CYCLE - 1:
                    for t in range(T):
                        trows.append(dict(cond=cond, ref=ref, lower=lower, seed=seed, pulse=pl + 1, taxon=taxa[t], alloc=alloc[t], tol=tol[t],
                                          cells_day1=int(count0[t]), cells_day365=int(cnt[t])))
                nn = max(n, 1)
                rows.append(dict(
                    cond=cond, ref=ref, lower=lower, seed=seed, pulse=pl + 1, day=i + 1, doy=i - pl * CYCLE + 1,
                    Temp=clim["Temp"].iloc[i], Psi=clim["Psi"].iloc[i],
                    biomass_C=C.sum() / G, alive_rows=n, taxa_alive=int((cnt > 0).sum()),
                    energy_to_C=E[alive].sum() / C[alive].sum() if n else np.nan, CUE=float(eco.CUE_system),
                    substrate_remaining_pct=100 * sub_total() / sub0,
                    net_energy=ps(eco.Microbe_E_Gain) / nn, E_cost_total=ps(eco.Microbe_E_Cost) / nn,
                    osmo_E_cost=ps(eco.Taxon_Osmo_E_Cost_total) / nn,
                ))
                if i == (pl + 1) * CYCLE - 1 and pl < pulses - 1:
                    eco.repopulation(stub, i, True)
    finally:
        os.chdir(cwd)
        shutil.rmtree(tmp, ignore_errors=True)
    return pd.DataFrame(rows), pd.DataFrame(trows)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("out_dir")
    ap.add_argument("climate")
    ap.add_argument("--inputs", default=os.path.join(HERE, "..", "input_e1500_enz2"))
    ap.add_argument("--pulses", type=int, default=6)
    ap.add_argument("--seeds", type=int, nargs="+", default=[0, 1, 2])
    ap.add_argument("--side", type=int, default=20)
    ap.add_argument("--induci", type=float, nargs=2, default=[0.01, 0.09])
    ap.add_argument("--consti", type=float, nargs=2, default=[0.005, 0.005])
    ap.add_argument("--conds", nargs="+", default=["0:-1", "0:0", "0.055:0", "0.095:0"], help="REF:LOWER pairs")
    ap.add_argument("--procs", type=int, default=None, help="maximum number of simultaneous runs (default: all CPUs)")
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    src, clim = os.path.abspath(a.inputs), os.path.abspath(a.climate)
    conds = [tuple(float(x) for x in c.split(":")) for c in a.conds]
    jobs = [(src, clim, a.pulses, s, a.side, a.induci, a.consti, c) for c in conds for s in a.seeds]
    with Pool(min(len(jobs), a.procs or os.cpu_count())) as p:
        parts = p.map(run, jobs, chunksize=1)
    pd.concat([p[0] for p in parts], ignore_index=True).to_csv(os.path.join(a.out_dir, "timeseries.csv"), index=False)
    pd.concat([p[1] for p in parts], ignore_index=True).to_csv(os.path.join(a.out_dir, "taxa.csv"), index=False)
    print("wrote", a.out_dir)


if __name__ == "__main__":
    main()
