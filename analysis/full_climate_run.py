"""
Multi-pulse runs of the full multi-taxon, multi-substrate model under a climate file (real Psi unless the file says otherwise).

Mirrors src/dementpy.py's pulse loop: at the end of each pulse Grid.repopulation() resets substrate, monomers and enzymes and
re-seeds the community from the taxa's cumulative abundance (a stub provides the per-taxon counts that Output normally tracks).

Usage:
    python full_climate_run.py OUT.csv ../input/climate.csv --inputs ../input_e1500_enz2 ../input_e1500_enz2_osmo \\
        [--pulses 6] [--seeds 0 1 2] [--side 20]
Writes one row per (input, seed, day).
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

CYCLE = 365


class _StubOutput:
    def __init__(self, n_taxa, total_days, index):
        self.Taxon_count_repop = pd.DataFrame(np.zeros((n_taxa, total_days + 2)), index=index)

    def record_pulse_init(self, grid, next_pulse):
        pass


def run(job):
    src, climate, pulses, seed, side = job
    tmp = tempfile.mkdtemp()
    w = os.path.join(tmp, "in")
    shutil.copytree(src, w)
    shutil.copy(climate, os.path.join(w, "climate.csv"))
    set_runtime(os.path.join(w, "runtime.txt"), dict(x=side, y=side, gridsize=side * side, pulse=pulses, end_time=CYCLE))
    clim = pd.read_csv(os.path.join(w, "climate.csv"), index_col=0)
    cwd = os.getcwd()
    os.chdir(w)
    rows = []
    try:
        np.random.seed(seed)
        rt = pd.read_csv("runtime.txt", header=None, index_col=0, sep="\t").drop("output_format").astype(float)
        data = initialize_data(rt)
        eco = Grid(rt, data)
        G = side * side
        T = int(rt.loc["n_taxa", 1])
        taxa = pd.unique(eco.Microbes.index)
        stub = _StubOutput(T, pulses * CYCLE, taxa)
        subs = [s for s in data["Substrate_index"] if s not in ("DeadMic", "DeadEnz")]
        sub_total = lambda: float(sum(np.asarray(eco.Substrates.loc[s, "C"]).sum() for s in subs))
        sub0 = sub_total()
        ps = lambda x: float(np.asarray(x).sum())
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
                stub.Taxon_count_repop.iloc[:, i + 1] = alive.reshape(G, T).sum(axis=0)
                nn = max(n, 1)
                rows.append(dict(
                    input=os.path.basename(src), seed=seed, pulse=pl + 1, day=i + 1, doy=i - pl * CYCLE + 1,
                    Temp=clim["Temp"].iloc[i], Psi=clim["Psi"].iloc[i],
                    biomass_C=C.sum() / G, alive_rows=n, taxa_alive=len(set(np.asarray(eco.Microbes.index)[alive])),
                    energy_to_C=E[alive].sum() / C[alive].sum() if n else np.nan, CUE=float(eco.CUE_system),
                    substrate_remaining_pct=100 * sub_total() / sub0,
                    HSP_prop=float(np.asarray(eco.HSP_prop)[alive].mean()) if n else np.nan,
                    net_energy=ps(eco.Microbe_E_Gain) / nn, E_cost_total=ps(eco.Microbe_E_Cost) / nn,
                    osmo_E_cost=ps(eco.Taxon_Osmo_E_Cost_total) / nn,
                ))
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
    ap.add_argument("--inputs", nargs="+", required=True)
    ap.add_argument("--pulses", type=int, default=6)
    ap.add_argument("--seeds", type=int, nargs="+", default=[0, 1, 2])
    ap.add_argument("--side", type=int, default=20)
    a = ap.parse_args()
    clim = os.path.abspath(a.climate)
    jobs = [(os.path.abspath(i), clim, a.pulses, s, a.side) for i in a.inputs for s in a.seeds]
    with Pool(min(len(jobs), os.cpu_count())) as p:
        parts = p.map(run, jobs, chunksize=1)
    pd.concat(parts, ignore_index=True).to_csv(a.out_csv, index=False)
    print("wrote", a.out_csv)


if __name__ == "__main__":
    main()
