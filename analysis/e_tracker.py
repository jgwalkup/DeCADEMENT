"""
Energy (E) and C-allocation tracker for the full / minimal model (real or flat climate, several pulses). Reads the model state around each
process step; the model is not changed and the random stream is identical to n_tracker.py / full_climate_osmo_range.py.

Per grid cell and day, per group g (bac = bacteria, fun = fungi):
  E_uptake_g      energy produced from the C routed to energy (Taxon_Uptake_E)
  E_maint_g, E_transport_g (active transport + transporter production), E_enzyme_g, E_hsp_g, E_osmo_g, E_cost_g (= Microbe_E_Cost)
  E_pool_g        energy in the cells at the end of the day;   biomass_C_g;   EC_g = E_pool / biomass C
  E_demand_g      mean fraction of uptake C that the cells route to energy (Tax_E_demand, from E:C and E_max), computed on the state used by uptake
  uptake_C_g, toE_C_g (C routed to energy = uptake_C - uptake_BC), uptake_BC_g (C routed to biomass), enzyme_C_g (C of enzyme production, taken from
                  the biomass C), overflow_C_g (C released by cells short of N or P in mortality), death_C_g
  E_clip_g        energy shortfall (costs above the energy available) that brings a negative E pool back to zero;  cells_Eclip_g = number of cells with a shortfall
  E_unpaid_g      part of the shortfall not paid for (all of it unless R_endo = 1, which pays it by respiring biomass C)
  R_endo_C_g     biomass C respired endogenously to pay shortfalls (R_endo = 1)
  E_cap_g         energy lost by capping E:C at E_max;   E_dead_g = energy in cells killed that day
  resid_E         (E at end - E at start of the day) - (uptake - costs + clip - cap - dead), summed over groups (should be ~0)
  CUE (model), cells_g (occupied cells)

Usage:
    python e_tracker.py OUT_DIR CLIMATE.csv [--input ../input_minimal_e1500_km25_fb1] [--pulses 1] [--side 20] [--seeds 0 1 2] [--procs 3]
        [--sub-n-mult K] [--osmo-override]
Writes OUT_DIR/e_daily.csv
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
    src, climate, pulses, seed, side, nmult, osmo = job[:7]
    params = job[7] if len(job) > 7 else {}
    tmp = tempfile.mkdtemp()
    w = os.path.join(tmp, "in")
    shutil.copytree(src, w)
    shutil.copy(climate, os.path.join(w, "climate.csv"))
    set_runtime(os.path.join(w, "runtime.txt"), dict(x=side, y=side, gridsize=side * side, pulse=pulses, end_time=CYCLE))
    if osmo:
        set_params(os.path.join(w, "parameters.csv"), {"Osmo_Induci_Prod_min": 0.04, "Osmo_Induci_Prod_max": 0.08,
                                                       "Osmo_Consti_Prod_min": 0.002, "Osmo_Consti_Prod_max": 0.005})
    if params:  # parameters.csv overrides; a row that is not in the file yet is appended (e.g. Uptake_psi_slope)
        pp = os.path.join(w, "parameters.csv")
        have = {l.split(",")[0] for l in open(pp).read().split("\n") if "," in l}
        old = {k: v for k, v in params.items() if k in have}
        new = {k: v for k, v in params.items() if k not in have}
        if old:
            set_params(pp, old)
        if new:
            txt = open(pp).read().rstrip("\n")
            open(pp, "w").write(txt + "".join(f"\n{k},{v}" for k, v in new.items()) + "\n")
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
        fb = np.asarray(eco.fb[:T]).astype(int)
        grp = {"bac": fb == 0, "fun": fb == 1}
        stub = _StubOutput(T, pulses * CYCLE, pd.unique(eco.Microbes.index))
        f32 = lambda x: np.asarray(x, dtype=np.float64).reshape(G, T)
        sub_names = np.asarray(eco.Substrates.index)
        sub_litter = ~np.isin(sub_names, ["DeadMic", "DeadEnz"])
        for pl in range(pulses):
            for i in range(pl * CYCLE, (pl + 1) * CYCLE):
                M_start = eco.Microbes.to_numpy(dtype=np.float64).reshape(G, T, 4)
                eco.degradation(pl, i)
                eco.uptake(pl, i)
                Mu = eco.Microbes.to_numpy(dtype=np.float64).reshape(G, T, 4)          # state used by uptake and metabolism
                upC, upBC, upE = f32(eco.Taxon_Uptake_C), f32(eco.Taxon_Uptake_BC), f32(eco.Taxon_Uptake_E)
                demand = np.asarray(eco.Tax_E_demand, dtype=np.float64).reshape(G, T)  # published by Grid.uptake()
                eco.metabolism(i)
                Mm = eco.Microbes.to_numpy(dtype=np.float64).reshape(G, T, 4)
                cost = {"maint": f32(eco.Taxon_Maint_E_Cost), "transport": f32(eco.Taxon_Transport_E_Cost),
                        "enzyme": f32(eco.Taxon_Enzyme_E_Cost), "hsp": f32(eco.Taxon_HSP_E_Cost), "osmo": f32(eco.Taxon_Osmo_E_Cost_total)}
                cost_tot = f32(eco.Microbe_E_Cost)
                enzC = f32(eco.Taxon_Enzyme_Cost_C)
                pre = Mu[..., 3] + f32(eco.Microbe_E_Gain)                              # E pool before the floor at zero
                clip = f32(eco.Taxon_E_deficit)                                         # shortfall (= -pre where pre < 0)
                unpaid, endo_C = f32(eco.Taxon_E_unpaid), f32(eco.Taxon_R_endo_C)
                eco.mortality(i)
                M1 = eco.Microbes.to_numpy(dtype=np.float64).reshape(G, T, 4)
                eco.reproduction(i)
                M2 = eco.Microbes.to_numpy(dtype=np.float64).reshape(G, T, 4)
                a0, a1, a2 = Mm[..., 0] > 0, M1[..., 0] > 0, M2[..., 0] > 0
                dead = a0 & ~a1
                row = dict(seed=seed, pulse=pl + 1, doy=i - pl * CYCLE + 1, CUE=float(eco.CUE_system),
                           sub_C_litter=float(eco.Substrates.to_numpy(dtype=np.float64)[:, 0][sub_litter].sum() / G))
                tot_up = tot_cost = tot_clip = tot_cap = tot_dead = 0.0
                for g, m in grp.items():
                    live = a1[:, m]
                    row[f"E_uptake_{g}"] = upE[:, m].sum() / G
                    for k, v in cost.items():
                        row[f"E_{k}_{g}"] = v[:, m].sum() / G
                    row[f"E_cost_{g}"] = cost_tot[:, m].sum() / G
                    row[f"E_pool_{g}"] = (M2[:, m, 3] * a2[:, m]).sum() / G
                    row[f"biomass_C_{g}"] = (M2[:, m, 0] * a2[:, m]).sum() / G
                    with np.errstate(invalid="ignore", divide="ignore"):
                        row[f"EC_{g}"] = row[f"E_pool_{g}"] / row[f"biomass_C_{g}"] if row[f"biomass_C_{g}"] > 0 else np.nan
                    w_ = Mu[:, m, 0] * (Mu[:, m, 0] > 0)
                    row[f"E_demand_{g}"] = float(np.nansum(demand[:, m] * w_) / w_.sum()) if w_.sum() > 0 else np.nan
                    row[f"uptake_C_{g}"] = upC[:, m].sum() / G
                    row[f"toE_C_{g}"] = (upC[:, m] - upBC[:, m]).sum() / G
                    row[f"uptake_BC_{g}"] = upBC[:, m].sum() / G
                    row[f"enzyme_C_{g}"] = enzC[:, m].sum() / G
                    row[f"overflow_C_{g}"] = ((Mm[:, m, 0] - M1[:, m, 0]) * live).sum() / G
                    row[f"death_C_{g}"] = (Mm[:, m, 0] * dead[:, m]).sum() / G
                    row[f"E_clip_{g}"] = clip[:, m].sum() / G
                    row[f"E_unpaid_{g}"] = unpaid[:, m].sum() / G
                    row[f"R_endo_C_{g}"] = endo_C[:, m].sum() / G
                    row[f"cells_Eclip_{g}"] = int(((pre[:, m] < 0) & (Mu[:, m, 0] > 0)).sum())
                    row[f"E_cap_{g}"] = ((Mm[:, m, 3] - M1[:, m, 3]) * live).sum() / G
                    row[f"E_dead_{g}"] = (Mm[:, m, 3] * dead[:, m]).sum() / G
                    row[f"cells_{g}"] = int(a2[:, m].sum())
                    tot_up += row[f"E_uptake_{g}"]
                    tot_cost += row[f"E_cost_{g}"]
                    tot_clip += row[f"E_clip_{g}"]
                    tot_cap += row[f"E_cap_{g}"]
                    tot_dead += row[f"E_dead_{g}"]
                e0 = (M_start[..., 3] * (M_start[..., 0] > 0)).sum() / G
                e1 = (M2[..., 3] * a2).sum() / G
                row["resid_E"] = (e1 - e0) - (tot_up - tot_cost + tot_clip - tot_cap - tot_dead)
                rows.append(row)
                stub.Taxon_count_repop.iloc[:, i + 1] = a2.sum(axis=0)
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
    ap.add_argument("--input", default=os.path.join(HERE, "..", "input_minimal_e1500_km25_fb1"))
    ap.add_argument("--pulses", type=int, default=1)
    ap.add_argument("--side", type=int, default=20)
    ap.add_argument("--seeds", type=int, nargs="+", default=[0, 1, 2])
    ap.add_argument("--procs", type=int, default=3)
    ap.add_argument("--sub-n-mult", type=float, default=1.0)
    ap.add_argument("--osmo-override", action="store_true", help="set the osmorange2 osmolyte production ranges (multi-taxon inputs)")
    ap.add_argument("--set", nargs="*", default=[], metavar="NAME=VALUE", help="parameters.csv overrides, e.g. Uptake_psi_slope=0.25")
    a = ap.parse_args()
    params = {k: float(v) for k, v in (x.split("=") for x in a.set)}
    os.makedirs(a.out_dir, exist_ok=True)
    jobs = [(os.path.abspath(a.input), os.path.abspath(a.climate), a.pulses, s, a.side, a.sub_n_mult, a.osmo_override, params) for s in a.seeds]
    with Pool(min(len(jobs), a.procs)) as p:
        parts = p.map(run, jobs, chunksize=1)
    pd.concat(parts, ignore_index=True).to_csv(os.path.join(a.out_dir, "e_daily.csv"), index=False)
    print("wrote", a.out_dir)


if __name__ == "__main__":
    main()
