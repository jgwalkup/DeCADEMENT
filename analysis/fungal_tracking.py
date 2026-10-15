"""
Track bacterial vs fungal abundance, biomass, births, deaths and stoichiometry day by day across pulses (real climate).

Same runs as full_climate_osmo_range.py (identical random stream; only reads the model state), with extra records:
  groups.csv      per seed, pulse, day, group (bacteria / fungi): cell-taxa alive, biomass C, mean C per cell, N:C and P:C of the
                  group's biomass, cells that start the day below 2 x C_min, new cells (colonisation), starvation deaths,
                  stochastic (drought / basal) deaths
  taxa_track.csv  per seed, pulse, taxon: fungal flag, cells and biomass C on day 1 and day 365, cumulative daily cell counts over the
                  pulse and the reseeding frequency derived from it (the probability with which the taxon is placed in a cell next pulse)
Usage:
    python fungal_tracking.py OUT_DIR ../analysis_output/climate_12_years.csv [--input ../input_e1500_enz2_km25_osmorange2]
        [--induci 0.04 0.08] [--consti 0.002 0.005] [--pulses 12] [--side 20] [--seeds 0 1 2 3 4] [--procs 5]
        [--sub-n-mult K] [--set NAME=VALUE ...]   override parameters.csv rows (e.g. Cfrac_f=0.825 Nfrac_f=0.16 Pfrac_f=0.015 max_size_f=2)
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
    src, climate, pulses, seed, side, induci, consti, extra = job
    tmp = tempfile.mkdtemp()
    w = os.path.join(tmp, "in")
    shutil.copytree(src, w)
    shutil.copy(climate, os.path.join(w, "climate.csv"))
    set_runtime(os.path.join(w, "runtime.txt"), dict(x=side, y=side, gridsize=side * side, pulse=pulses, end_time=CYCLE))
    set_params(os.path.join(w, "parameters.csv"), {"Osmo_Induci_Prod_min": induci[0], "Osmo_Induci_Prod_max": induci[1],
                                                   "Osmo_Consti_Prod_min": consti[0], "Osmo_Consti_Prod_max": consti[1]})
    if extra.get("params"):
        set_params(os.path.join(w, "parameters.csv"), extra["params"])
    if extra.get("runtime"):
        set_runtime(os.path.join(w, "runtime.txt"), extra["runtime"])
    if extra.get("sub_n_mult", 1.0) != 1.0:  # scale the N of every initial substrate (raises the N supply, C unchanged)
        sp = os.path.join(w, "initial_substrates.csv")
        sdf = pd.read_csv(sp, index_col=0)
        sdf["N"] = sdf["N"] * extra["sub_n_mult"]
        sdf.to_csv(sp)
    cwd = os.getcwd()
    os.chdir(w)
    grows, trows, drows = [], [], []
    try:
        np.random.seed(seed)
        rt = pd.read_csv("runtime.txt", header=None, index_col=0, sep="\t").drop("output_format").astype(float)
        eco = Grid(rt, initialize_data(rt))
        G, T = side * side, int(rt.loc["n_taxa", 1])
        fb = np.asarray(eco.fb[:T]).astype(int)
        names = list(pd.unique(eco.Microbes.index))
        stub = _StubOutput(T, pulses * CYCLE, pd.unique(eco.Microbes.index))
        subs = [x for x in initialize_data.__globals__.get("_unused", [])] or [x for x in eco.Substrates.index.get_level_values(0).unique()
                                                                              if x not in ("DeadMic", "DeadEnz")]
        sub_c = lambda: float(sum(np.asarray(eco.Substrates.loc[x, "C"]).sum() for x in subs))
        sub_n = lambda: float(sum(np.asarray(eco.Substrates.loc[x, "N"]).sum() for x in subs))
        sub0 = sub_c()
        Cmin, Nmin, Pmin = float(eco.C_min), float(eco.N_min), float(eco.P_min)
        groups = {"bacteria": fb == 0, "fungi": fb == 1}
        for pl in range(pulses):
            cum = np.zeros(T)
            for i in range(pl * CYCLE, (pl + 1) * CYCLE):
                doy = i - pl * CYCLE + 1
                eco.degradation(pl, i)
                eco.uptake(pl, i)
                eco.metabolism(i)
                M0 = eco.Microbes.to_numpy(dtype=np.float32).reshape(G, T, 4)           # before mortality (start of mortality step)
                eco.mortality(i)
                M1 = eco.Microbes.to_numpy(dtype=np.float32).reshape(G, T, 4)           # after mortality
                eco.reproduction(i)
                M2 = eco.Microbes.to_numpy(dtype=np.float32).reshape(G, T, 4)           # after reproduction (end of day)
                a0, a1, a2 = M0[..., 0] > 0, M1[..., 0] > 0, M2[..., 0] > 0
                starve = a0 & ((M0[..., 0] < Cmin) | (M0[..., 1] < Nmin) | (M0[..., 2] < Pmin))
                stoch = a0 & ~starve & ~a1
                born = a2 & ~a1
                cnt = a2.sum(axis=0)
                E2 = M2[..., 3]
                drows.append(dict(seed=seed, pulse=pl + 1, doy=doy, remaining_pct=100 * sub_c() / sub0, substrate_N=sub_n() / G,
                                  energy_to_C=float(E2[a2].sum() / M2[..., 0][a2].sum()) if a2.any() else np.nan,
                                  NH4=float(np.asarray(eco.Monomers.loc["NH4", "N"]).sum()) / G, biomass_C=float(M2[..., 0].sum()) / G))
                stub.Taxon_count_repop.iloc[:, i + 1] = cnt
                cum += cnt
                for gname, mask in groups.items():
                    sel = np.zeros(T, dtype=bool) if not mask.any() else mask
                    C, N, P = (M2[:, sel, k] for k in range(3))
                    nC = int(a2[:, sel].sum())
                    sC = float(C.sum())
                    grows.append(dict(seed=seed, pulse=pl + 1, doy=doy, group=gname, cells=nC, biomass_C=sC / G,
                                      mean_C_per_cell=sC / nC if nC else np.nan, N_to_C=float(N.sum()) / sC if sC else np.nan,
                                      P_to_C=float(P.sum()) / sC if sC else np.nan,
                                      near_starve=int((a2[:, sel] & (C < 2 * Cmin)).sum()),
                                      new_cells=int(born[:, sel].sum()), starve_deaths=int(starve[:, sel].sum()),
                                      stochastic_deaths=int(stoch[:, sel].sum())))
                if doy == 1:
                    cells1, C1 = a2.sum(axis=0).copy(), M2[..., 0].sum(axis=0).copy()
                if doy == CYCLE:
                    cells2, C2 = a2.sum(axis=0).copy(), M2[..., 0].sum(axis=0).copy()
                if i == (pl + 1) * CYCLE - 1:
                    freq = cum / cum.sum() if cum.sum() else np.zeros(T)
                    for t in range(T):
                        trows.append(dict(seed=seed, pulse=pl + 1, taxon=names[t], fungal=int(fb[t]), cells_day1=int(cells1[t]),
                                          C_day1=float(C1[t]) / G, cells_day365=int(cells2[t]), C_day365=float(C2[t]) / G,
                                          cum_cells=float(cum[t]), reseed_freq=float(freq[t])))
                    if pl < pulses - 1:
                        eco.repopulation(stub, i, True)
    finally:
        os.chdir(cwd)
        shutil.rmtree(tmp, ignore_errors=True)
    return pd.DataFrame(grows), pd.DataFrame(trows), pd.DataFrame(drows)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("out_dir")
    ap.add_argument("climate")
    ap.add_argument("--input", default=os.path.join(HERE, "..", "input_e1500_enz2_km25_osmorange2"))
    ap.add_argument("--induci", type=float, nargs=2, default=[0.04, 0.08])
    ap.add_argument("--consti", type=float, nargs=2, default=[0.002, 0.005])
    ap.add_argument("--pulses", type=int, default=12)
    ap.add_argument("--side", type=int, default=20)
    ap.add_argument("--seeds", type=int, nargs="+", default=[0, 1, 2, 3, 4])
    ap.add_argument("--procs", type=int, default=5)
    ap.add_argument("--set", nargs="*", default=[], metavar="NAME=VALUE")
    ap.add_argument("--sub-n-mult", type=float, default=1.0, help="multiply the N of all initial substrates (N supply)")
    a = ap.parse_args()
    extra = {"params": {k: float(v) for k, v in (x.split("=") for x in a.set)}, "sub_n_mult": a.sub_n_mult}
    os.makedirs(a.out_dir, exist_ok=True)
    jobs = [(os.path.abspath(a.input), os.path.abspath(a.climate), a.pulses, s, a.side, a.induci, a.consti, extra) for s in a.seeds]
    with Pool(min(len(jobs), a.procs)) as p:
        parts = p.map(run, jobs, chunksize=1)
    pd.concat([x[0] for x in parts], ignore_index=True).to_csv(os.path.join(a.out_dir, "groups.csv"), index=False)
    pd.concat([x[1] for x in parts], ignore_index=True).to_csv(os.path.join(a.out_dir, "taxa_track.csv"), index=False)
    pd.concat([x[2] for x in parts], ignore_index=True).to_csv(os.path.join(a.out_dir, "daily.csv"), index=False)
    print("wrote", a.out_dir)


if __name__ == "__main__":
    main()
