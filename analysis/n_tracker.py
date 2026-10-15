"""
N pools, N fluxes and a daily N budget check for the full model (real climate, several pulses). Reads the model state around each
process step; the model itself is not changed and the random stream is identical to full_climate_osmo_range.py / fungal_tracking.py.

Pools (per grid cell, mean over cells, N units):
  sub_N_litter    substrates other than DeadMic / DeadEnz            sub_N_dead   DeadMic + DeadEnz substrates
  mon_N_org       N in organic (and dead) monomers                   mon_N_nh4    NH4
  mic_N_bac/fun   microbial N by group                                enz_N        N in the enzyme pools (implicit: enzyme C x N_cost)
Fluxes (per grid cell per day):
  decay_N          substrates -> monomers (degradation)
  uptake_N_bac/fun monomers -> microbes (Taxon_Uptake_N)
  enz_prod_N       microbes -> enzymes (inducible + constitutive enzyme N cost), = uptake_N - change in microbial N in metabolism
  deadenz_N        enzymes -> DeadEnz substrate (enzyme turnover)
  death_N_bac/fun  killed / starved cells -> DeadMic substrate
  release_N_bac/fun  N released from living cells to NH4 when their N:C is above the minimum (MicLoss N)
  leach_N          NH4 leaching (the only loss of N from the system)
  overflow_C_bac/fun  C released from living cells to respiration in the same step (MicLoss C: C in excess of the limiting element);
  uptake_C_bac/fun, death_C_bac/fun for context
  resid_N          change in total N over the day + leach_N (should be 0 up to float32 rounding)
HSP and osmolyte N costs are not deducted from any pool in the model (they only enter the N limitation of inducible production).

Usage:
    python n_tracker.py OUT_DIR ../analysis_output/climate_12_years.csv [--input ../input_e1500_enz2_km25_osmorange2]
        [--induci 0.04 0.08] [--consti 0.002 0.005] [--pulses 6] [--side 20] [--seeds 0 1 2 3 4] [--procs 5] [--sub-n-mult K] [--keep-osmo]
Writes OUT_DIR/n_daily.csv (one row per seed, pulse, day).
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
    src, climate, pulses, seed, side, induci, consti, sub_n_mult, keep_osmo = job[:9]
    params = job[9] if len(job) > 9 else {}  # extra parameters.csv overrides, e.g. Input_NH4
    tmp = tempfile.mkdtemp()
    w = os.path.join(tmp, "in")
    shutil.copytree(src, w)
    shutil.copy(climate, os.path.join(w, "climate.csv"))
    set_runtime(os.path.join(w, "runtime.txt"), dict(x=side, y=side, gridsize=side * side, pulse=pulses, end_time=CYCLE))
    if not keep_osmo:
        set_params(os.path.join(w, "parameters.csv"), {"Osmo_Induci_Prod_min": induci[0], "Osmo_Induci_Prod_max": induci[1],
                                                       "Osmo_Consti_Prod_min": consti[0], "Osmo_Consti_Prod_max": consti[1]})
    if params:
        set_params(os.path.join(w, "parameters.csv"), params)
    if sub_n_mult != 1.0:
        sp = os.path.join(w, "initial_substrates.csv")
        sdf = pd.read_csv(sp, index_col=0)
        sdf["N"] = sdf["N"] * sub_n_mult
        sdf.to_csv(sp)
    cwd = os.getcwd()
    os.chdir(w)
    rows = []
    try:
        np.random.seed(seed)
        rt = pd.read_csv("runtime.txt", header=None, index_col=0, sep="\t").drop("output_format").astype(float)
        eco = Grid(rt, initialize_data(rt))
        G, T = side * side, int(rt.loc["n_taxa", 1])
        NE = int(rt.loc["n_enzymes", 1])
        fb = np.asarray(eco.fb[:T]).astype(int)
        grp = {"bac": fb == 0, "fun": fb == 1}
        stub = _StubOutput(T, pulses * CYCLE, pd.unique(eco.Microbes.index))
        sub_names = np.asarray(eco.Substrates.index)
        is_dead = np.isin(sub_names, ["DeadMic", "DeadEnz"])
        mon_names = np.asarray(eco.Monomers.index)
        is_nh4 = mon_names == "NH4"
        n_cost = eco.Enz_Attrib["N_cost"].reindex(eco._enz_index).to_numpy(dtype=np.float64)
        p_cost = eco.Enz_Attrib["P_cost"].reindex(eco._enz_index).to_numpy(dtype=np.float64)
        is_po4 = mon_names == "PO4"
        mr = eco.MinRatios.to_numpy(dtype=np.float64).reshape(G, T, 4)[0]          # (T,4) minimum mass fractions C, N, P, E
        minC, minN, minP = mr[:, 0], mr[:, 1], mr[:, 2]
        in_N = float(getattr(eco, "_nh4_input", 0.0))
        in_P = float(getattr(eco, "_po4_input", 0.0))

        def pools():
            S = eco.Substrates.to_numpy(dtype=np.float64)[:, 1]
            Mo = eco.Monomers.to_numpy(dtype=np.float64)[:, 1]
            Mi = eco.Microbes.to_numpy(dtype=np.float64).reshape(G, T, 4)
            Ez = eco.Enzymes.to_numpy(dtype=np.float64).reshape(G, NE)
            Sc = eco.Substrates.to_numpy(dtype=np.float64)[:, 0]
            Sp = eco.Substrates.to_numpy(dtype=np.float64)[:, 2]
            Mp = eco.Monomers.to_numpy(dtype=np.float64)[:, 2]
            return dict(P_sub=Sp.sum() / G, P_mon_org=Mp[~is_po4].sum() / G, P_mon_po4=Mp[is_po4].sum() / G, P_mic_bac=Mi[:, grp["bac"], 2].sum() / G,
                        P_mic_fun=Mi[:, grp["fun"], 2].sum() / G, P_enz=(Ez * p_cost[None, :]).sum() / G,
                        sub_C_litter=Sc[~is_dead].sum() / G, mic_C_bac=Mi[:, grp["bac"], 0].sum() / G, mic_C_fun=Mi[:, grp["fun"], 0].sum() / G,
                        sub_litter=S[~is_dead].sum() / G, sub_dead=S[is_dead].sum() / G, mon_org=Mo[~is_nh4].sum() / G,
                        mon_nh4=Mo[is_nh4].sum() / G, mic_bac=Mi[:, grp["bac"], 1].sum() / G, mic_fun=Mi[:, grp["fun"], 1].sum() / G,
                        enz=(Ez * n_cost[None, :]).sum() / G), Mi

        for pl in range(pulses):
            for i in range(pl * CYCLE, (pl + 1) * CYCLE):
                doy = i - pl * CYCLE + 1
                p0, _ = pools()
                eco.degradation(pl, i)
                p1, _ = pools()
                eco.uptake(pl, i)
                p2, _ = pools()
                upN = np.asarray(eco.Taxon_Uptake_N, dtype=np.float64).reshape(G, T)
                upC = np.asarray(eco.Taxon_Uptake_C, dtype=np.float64).reshape(G, T)
                eco.metabolism(i)
                p3, M0 = pools()
                eco.mortality(i)
                p4, M1 = pools()
                eco.reproduction(i)
                p5, _ = pools()
                a0, a1 = M0[..., 0] > 0, M1[..., 0] > 0
                dead = a0 & ~a1
                row = dict(seed=seed, pulse=pl + 1, doy=doy)
                for k in ("sub_litter", "sub_dead", "mon_org", "mon_nh4", "mic_bac", "mic_fun", "enz"):
                    row[{"sub_litter": "sub_N_litter", "sub_dead": "sub_N_dead", "mon_org": "mon_N_org", "mon_nh4": "mon_N_nh4",
                         "mic_bac": "mic_N_bac", "mic_fun": "mic_N_fun", "enz": "enz_N"}[k]] = p5[k]
                tot = lambda p: sum(v for k, v in p.items() if not k.startswith(("sub_C", "mic_C", "P_")))
                totP = lambda p: sum(v for k, v in p.items() if k.startswith("P_"))
                row["total_N"] = tot(p5)
                row["input_N"], row["input_P"] = in_N, in_P
                row["total_P"] = totP(p5)
                row["mon_P_po4"], row["sub_P"], row["mon_P_org"] = p5["P_mon_po4"], p5["P_sub"], p5["P_mon_org"]
                row["mic_P_bac"], row["mic_P_fun"] = p5["P_mic_bac"], p5["P_mic_fun"]
                alive_end = eco.Microbes["C"].to_numpy().reshape(G, T) > 0
                row["cells_bac"], row["cells_fun"] = int(alive_end[:, grp["bac"]].sum()), int(alive_end[:, grp["fun"]].sum())
                row["mic_C_bac"], row["mic_C_fun"], row["sub_C_litter"] = p5["mic_C_bac"], p5["mic_C_fun"], p5["sub_C_litter"]
                Ee = M1[..., 3]
                row["energy_to_C"] = float(Ee[a1].sum() / M1[..., 0][a1].sum()) if a1.any() else np.nan
                row["CUE"] = float(eco.CUE_system)
                row["decay_N"] = (p0["sub_litter"] + p0["sub_dead"]) - (p1["sub_litter"] + p1["sub_dead"])
                row["deadenz_N"] = (p3["sub_litter"] + p3["sub_dead"]) - (p2["sub_litter"] + p2["sub_dead"])
                tot_up = 0.0
                upP = np.asarray(eco.Taxon_Uptake_P, dtype=np.float64).reshape(G, T)
                with np.errstate(invalid="ignore", divide="ignore"):
                    S3 = M0[..., 0] + M0[..., 1] + M0[..., 2]
                    fC, fN, fP = M0[..., 0] / S3, M0[..., 1] / S3, M0[..., 2] / S3
                below = a1 & ((fC < minC[None, :]) | (fN < minN[None, :]) | (fP < minP[None, :]))
                with np.errstate(invalid="ignore", divide="ignore"):
                    amounts = np.stack([M0[..., 0] / minC[None, :], M0[..., 1] / minN[None, :], M0[..., 2] / minP[None, :]], axis=-1)
                limiter = np.argmin(np.where(np.isnan(amounts), np.inf, amounts), axis=-1)       # 0 = C, 1 = N, 2 = P
                loss_C = (M0[..., 0] - M1[..., 0]) * a1
                for g, m in grp.items():
                    row[f"uptake_N_{g}"] = upN[:, m].sum() / G
                    row[f"uptake_C_{g}"] = upC[:, m].sum() / G
                    row[f"death_N_{g}"] = (M0[:, m, 1] * dead[:, m]).sum() / G
                    row[f"death_C_{g}"] = (M0[:, m, 0] * dead[:, m]).sum() / G
                    live = a1[:, m]
                    row[f"release_N_{g}"] = ((M0[:, m, 1] - M1[:, m, 1]) * live).sum() / G
                    row[f"overflow_C_{g}"] = ((M0[:, m, 0] - M1[:, m, 0]) * live).sum() / G
                    row[f"uptake_P_{g}"] = upP[:, m].sum() / G
                    for k, nm in ((0, "C"), (1, "N"), (2, "P")):
                        sel = below[:, m] & (limiter[:, m] == k)
                        row[f"overflow_C_lim{nm}_{g}"] = (loss_C[:, m] * sel).sum() / G
                        row[f"cells_lim{nm}_{g}"] = int(sel.sum())
                    row[f"cells_belowN_{g}"] = int((a1[:, m] & (fN[:, m] < minN[m][None, :])).sum())
                    row[f"cells_belowP_{g}"] = int((a1[:, m] & (fP[:, m] < minP[m][None, :])).sum())
                    tot_up += row[f"uptake_N_{g}"]
                row["enz_prod_N"] = tot_up - (p3["mic_bac"] + p3["mic_fun"] - p2["mic_bac"] - p2["mic_fun"])
                r = 0.1 * np.exp(0.5 * (eco.psi[i] - eco.wp_fc))
                row["leach_N"] = p4["mon_nh4"] * r / (1 - r)
                row["leach_P"] = p4["P_mon_po4"] * r / (1 - r)
                row["resid_N"] = (tot(p5) - tot(p0)) + row["leach_N"] - in_N
                row["resid_P"] = (totP(p5) - totP(p0)) + row["leach_P"] - in_P
                row["mon_N_check_uptake"] = (p2["mon_org"] + p2["mon_nh4"]) - (p1["mon_org"] + p1["mon_nh4"]) + tot_up - row["decay_N"] - in_N
                rows.append(row)
                cnt = (M1[..., 0] > 0).sum(axis=0)
                stub.Taxon_count_repop.iloc[:, i + 1] = (eco.Microbes["C"].to_numpy().reshape(G, T) > 0).sum(axis=0)
                if i == (pl + 1) * CYCLE - 1 and pl < pulses - 1:
                    eco.repopulation(stub, i, True)
    finally:
        os.chdir(cwd)
        shutil.rmtree(tmp, ignore_errors=True)
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("out_dir")
    ap.add_argument("climate")
    ap.add_argument("--input", default=os.path.join(HERE, "..", "input_e1500_enz2_km25_osmorange2"))
    ap.add_argument("--induci", type=float, nargs=2, default=[0.04, 0.08])
    ap.add_argument("--consti", type=float, nargs=2, default=[0.002, 0.005])
    ap.add_argument("--pulses", type=int, default=6)
    ap.add_argument("--side", type=int, default=20)
    ap.add_argument("--seeds", type=int, nargs="+", default=[0, 1, 2, 3, 4])
    ap.add_argument("--procs", type=int, default=5)
    ap.add_argument("--sub-n-mult", type=float, default=1.0)
    ap.add_argument("--keep-osmo", action="store_true", help="do not overwrite the osmolyte production parameters of the input folder")
    ap.add_argument("--set", nargs="*", default=[], metavar="NAME=VALUE", help="override parameters.csv rows, e.g. Input_NH4=0.002 Input_PO4=0.00018")
    a = ap.parse_args()
    params = {k: float(v) for k, v in (x.split("=") for x in a.set)}
    os.makedirs(a.out_dir, exist_ok=True)
    jobs = [(os.path.abspath(a.input), os.path.abspath(a.climate), a.pulses, s, a.side, a.induci, a.consti, a.sub_n_mult, a.keep_osmo, params) for s in a.seeds]
    with Pool(min(len(jobs), a.procs)) as p:
        parts = p.map(run, jobs, chunksize=1)
    pd.concat(parts, ignore_index=True).to_csv(os.path.join(a.out_dir, "n_daily.csv"), index=False)
    print("wrote", a.out_dir)


if __name__ == "__main__":
    main()
