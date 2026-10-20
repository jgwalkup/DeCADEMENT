"""Full-model runs (12 substrates) with per-substrate and per-monomer bookkeeping, for the single-taxon port check
(SENSITIVITY_PLAN 2.2) and small community runs.

    python analysis/full_single_taxon.py OUT_DIR --label NAME [--input input_full_port_bac] [--seeds 0 1 2]
        [--climates CSV ...] [--kinetics 1] [--side 20] [--days 365] [--procs 2]
        [--set NAME=VALUE ...] [--runtime NAME=VALUE ...] [--litter-C TOTAL]

Defaults give the single-taxon check: runtime n_taxa 1, taxa_per_box 1, fb 0, one pulse; parameters
Enz_per_taxon_min = max = 40 (a generalist with 40 of the 50 enzyme genes), R_endo 1. --set and --runtime override
any of these (NAME may be a NAME_min / NAME_max stem, or monomer_energy:<row>:<column>). --litter-C scales every litter
substrate's C, N and P in initial_substrates.csv (and sub_mon_inputs.csv) so the litter C totals TOTAL (composition and
stoichiometry kept). Appends to OUT_DIR/timeseries.csv.gz and OUT_DIR/summary.csv (rows labelled NAME; existing rows with
the same label, climate, kinetics and seed are skipped).

Daily series: litter remaining (% of initial litter C, all substrates except DeadMic and DeadEnz), each substrate's C
(per grid cell), biomass, CUE (model), C uptake, energy production, endogenous respiration C, respiration, and the C
taken up per monomer (from the monomer pools before and after uptake; FermProd excluded). Summary: % litter degraded,
day 50 % / 95 %, uptake-weighted CUE, endogenous-respiration share, and the uptake-weighted maximum energy per unit C taken
up (sum_m U_m f_m (0.9 E_resp + 0.1 E_ferm) / sum_m U_m, f_m the monomer's energy routing share; the minimal IdealSub
gives 1373).
"""
import argparse
import itertools
import os
import shutil
import sys
import tempfile
from multiprocessing import Pool

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "src"))
from sweep_minimal import set_params, set_runtime, _RepopStub  # noqa: E402

CLIM_DIR = os.path.join(HERE, "..", "input", "climate_data")
NECRO = ("DeadMic", "DeadEnz")
DEFAULT_SET = {"Enz_per_taxon": 40, "R_endo": 1}
DEFAULT_RUNTIME = {"n_taxa": 1, "taxa_per_box": 1, "fb": 0}


def _parse(items):
    out = {}
    for it in items or []:
        k, v = it.split("=", 1)
        out[k] = v
    return out


def _prepare(src, work, changes, runtime, side, days, climate, litter_C):
    shutil.copytree(src, work)
    me = {k: v for k, v in changes.items() if k.startswith("monomer_energy:")}
    set_params(os.path.join(work, "parameters.csv"), {k: float(v) for k, v in changes.items() if k not in me})
    if me:
        mp = os.path.join(work, "monomer_energy.csv")
        mdf = pd.read_csv(mp, index_col=0).astype(float)
        for k, v in me.items():
            _, row, col = k.split(":")
            mdf.loc[row, col] = float(v)
        mdf.to_csv(mp)
    if litter_C:
        sp = os.path.join(work, "initial_substrates.csv")
        s = pd.read_csv(sp, index_col=0)
        lit = [i for i in s.index if i not in NECRO]
        f = float(litter_C) / s.loc[lit, "C"].sum()
        s.loc[lit, ["C", "N", "P"]] *= f
        s.to_csv(sp)
        ip = os.path.join(work, "sub_mon_inputs.csv")
        si = pd.read_csv(ip, index_col=0)
        si.loc[[i for i in si.index if i not in NECRO]] *= f
        si.to_csv(ip)
    set_runtime(os.path.join(work, "runtime.txt"),
                dict(x=side, y=side, gridsize=side * side, pulse=1, end_time=days, **runtime))
    if climate:
        shutil.copy(climate, os.path.join(work, "climate.csv"))


def _job(args):
    label, src, changes, runtime, side, days, climate, kin, seed, litter_C = args
    from initialization import initialize_data
    from grid import Grid
    tmp = tempfile.mkdtemp()
    work = os.path.join(tmp, "in")
    _prepare(src, work, changes, dict(runtime, Kinetics_Model=kin), side, days, climate, litter_C)
    cwd = os.getcwd()
    os.chdir(work)
    try:
        np.random.seed(seed)
        rt = pd.read_csv("runtime.txt", header=None, index_col=0, sep="\t").drop("output_format").astype(float)
        data = initialize_data(rt)
        eco = Grid(rt, data)
        G = side * side
        subs = list(data["Substrate_index"])
        lit = [s for s in subs if s not in NECRO]
        mons = list(data["Monomer_index"])
        org = [m for m in mons if m not in ("NH4", "PO4", "FermProd")]
        sub_idx = np.asarray(eco.Substrates.index)
        mon_idx = np.asarray(eco.Monomers.index)

        def by(index, df, names):
            c = np.asarray(df["C"], dtype=float)
            return {n: float(c[index == n].sum()) / G for n in names}

        lit0 = sum(by(sub_idx, eco.Substrates, lit).values())
        me = data["MonEnergy_Attrib"].astype(float)
        fc = float(eco.ferm_cond)
        route = me.E_value ** ((2 / np.log10(me.E_ferm)) ** fc)
        e_max = (route * ((1 - fc) * me.E_resp + fc * me.E_ferm)).to_dict()
        T = int(rt.loc["n_taxa", 1])
        stub = _RepopStub(pd.unique(eco.Microbes.index), days + 2)
        rows = []
        reb = {"C": 0.0}
        _orig_update_cue = eco._update_cue

        def _update_cue(rebalance_C_loss=None):            # record the C respired by stoichiometric rebalancing
            reb["C"] = 0.0 if rebalance_C_loss is None else float(np.sum(rebalance_C_loss))
            return _orig_update_cue(rebalance_C_loss)
        eco._update_cue = _update_cue
        is_nh4, is_po4 = mon_idx == "NH4", mon_idx == "PO4"
        for i in range(days):
            reb["C"] = 0.0
            eco.degradation(0, i)
            m_before = by(mon_idx, eco.Monomers, org)
            eco.uptake(0, i)
            m_after = by(mon_idx, eco.Monomers, org)
            eco.metabolism(i)
            eco.mortality(i)
            eco.reproduction(i)
            C = np.asarray(eco.Microbes["C"])
            alive = C > 0
            sc = by(sub_idx, eco.Substrates, subs)
            # uptake() first adds the decay products of this step to the monomer pools, then removes the uptake
            mr = eco.Monomers.index
            is_org_m = (mr != "NH4") & (mr != "PO4") & (mr != "FermProd")
            dec_c = np.asarray(eco.Monomer_ratios["C"], dtype=float)[is_org_m] * np.asarray(eco.DecayRates, dtype=float)
            dec = {m: float(dec_c[np.asarray(mr)[is_org_m] == m].sum()) / G for m in org}
            upt = {m: max(m_before[m] + dec[m] - m_after[m], 0.0) for m in org}
            row = dict(doy=i + 1, litter_remaining_pct=100 * sum(sc[s] for s in lit) / lit0,
                       biomass_C=float(C.sum()) / G, alive_cells=int(alive.sum()), CUE=float(eco.CUE_system),
                       C_uptake=float(np.asarray(eco.Taxon_Uptake_C).sum()) / G,
                       E_production=float(np.asarray(eco.Taxon_Uptake_E).sum()) / G,
                       R_endo_C=float(np.sum(getattr(eco, "Taxon_R_endo_C", 0.0))) / G,
                       respiration=float(np.asarray(eco.Taxon_Resp).sum()) / G,
                       enzyme_pool=float(np.asarray(eco.Enzymes).sum()) / G,
                       killed=int(eco.Kill), Temp=float(eco.temp[i]), Psi=float(eco.psi[i]),
                       rebalance_C=reb["C"] / G,
                       R_stoich_E=float(np.sum(getattr(eco, "Taxon_R_stoich_E", 0.0))) / G,
                       N_uptake=float(np.asarray(eco.Taxon_Uptake_N).sum()) / G,
                       P_uptake=float(np.asarray(eco.Taxon_Uptake_P).sum()) / G,
                       NH4_N=float(np.asarray(eco.Monomers["N"])[is_nh4].sum()) / G,
                       PO4_P=float(np.asarray(eco.Monomers["P"])[is_po4].sum()) / G,
                       leached_N=float(np.sum(eco.Leached_N)) / G,
                       biomass_NC=float(np.asarray(eco.Microbes["N"])[alive].sum() / C[alive].sum()) if alive.any() else np.nan,
                       biomass_PC=float(np.asarray(eco.Microbes["P"])[alive].sum() / C[alive].sum()) if alive.any() else np.nan,
                       induci_limitation=float(np.nanmean(eco.Induci_limitation)) if np.size(eco.Induci_limitation) else np.nan)
            row.update({f"S_{s}": v for s, v in sc.items()})
            row.update({f"U_{m}": v for m, v in upt.items()})
            rows.append(row)
            stub.Taxon_count_repop.iloc[:, i + 1] = alive.reshape(G, T).sum(axis=0)
    finally:
        os.chdir(cwd)
        shutil.rmtree(tmp, ignore_errors=True)
    s = pd.DataFrame(rows)
    deg = 100 - s.litter_remaining_pct.to_numpy()
    w = s.C_uptake.clip(lower=0).to_numpy()
    U = s[[f"U_{m}" for m in org]].sum()
    eyield = sum(U[f"U_{m}"] * e_max[m] for m in org) / U.sum() if U.sum() > 0 else np.nan
    summ = dict(label=label, climate=os.path.basename(climate or "climate.csv").replace(".csv", "").replace("climate_", ""),
                kinetics="MMRT" if kin == 1 else "Arrhenius", seed=seed,
                litter_degraded=float(deg[-1]),
                day_50pct=float(np.argmax(deg >= 50) + 1) if (deg >= 50).any() else 366.0,
                day_95pct=float(np.argmax(deg >= 95) + 1) if (deg >= 95).any() else 366.0,
                CUE_annual=float(np.average(s.CUE.fillna(0), weights=w)) if w.sum() > 0 else 0.0,
                endo_share_pct=100 * s.R_endo_C.sum() / s.respiration.sum() if s.respiration.sum() > 0 else 0.0,
                rebalance_share_of_uptake_pct=100 * s.rebalance_C.sum() / s.C_uptake.sum() if s.C_uptake.sum() > 0 else 0.0,
                E_max_per_uptake_C=float(eyield), biomass_max=float(s.biomass_C.max()),
                **{f"uptake_share_{m}": float(U[f"U_{m}"] / U.sum()) for m in org},
                **{f"left_{x}": float(100 * s[f"S_{x}"].iloc[-1] / max(s[f"S_{x}"].iloc[0], 1e-12)) for x in lit})
    for k, v in (("label", label), ("climate", summ["climate"]), ("kinetics", summ["kinetics"]), ("seed", seed)):
        s.insert(0, k, v)
    return summ, s


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("out_dir")
    ap.add_argument("--label", required=True)
    ap.add_argument("--input", default=os.path.join(HERE, "..", "input_full_port_bac"))
    ap.add_argument("--seeds", type=int, nargs="+", default=[0, 1, 2])
    ap.add_argument("--climates", nargs="+", default=[os.path.join(CLIM_DIR, "climate_ambient_temp.csv")])
    ap.add_argument("--kinetics", type=int, nargs="+", default=[1])
    ap.add_argument("--side", type=int, default=20)
    ap.add_argument("--days", type=int, default=365)
    ap.add_argument("--procs", type=int, default=2)
    ap.add_argument("--set", nargs="*", default=[])
    ap.add_argument("--runtime", nargs="*", default=[])
    ap.add_argument("--litter-C", type=float, default=None)
    a = ap.parse_args()
    changes = dict(DEFAULT_SET, **_parse(a.set))
    runtime = dict(DEFAULT_RUNTIME, **_parse(a.runtime))
    os.makedirs(a.out_dir, exist_ok=True)
    sp, tp = os.path.join(a.out_dir, "summary.csv"), os.path.join(a.out_dir, "timeseries.csv.gz")
    done = set()
    if os.path.exists(sp):
        prev = pd.read_csv(sp)
        done = set(zip(prev.label, prev.climate, prev.kinetics, prev.seed))
    jobs = []
    for c, k, sd in itertools.product(a.climates, a.kinetics, a.seeds):
        cn = os.path.basename(c).replace(".csv", "").replace("climate_", "")
        if (a.label, cn, "MMRT" if k == 1 else "Arrhenius", sd) not in done:
            jobs.append((a.label, os.path.abspath(a.input), changes, runtime, a.side, a.days, os.path.abspath(c), k, sd,
                         a.litter_C))
    print(f"{len(jobs)} runs -> {a.out_dir} [{a.label}]", flush=True)
    with Pool(a.procs) as pool:
        for n, (summ, s) in enumerate(pool.imap_unordered(_job, jobs), 1):
            prev = pd.read_csv(sp) if os.path.exists(sp) else None   # rewritten whole: columns differ between inputs
            pd.concat([prev, pd.DataFrame([summ])] if prev is not None else [pd.DataFrame([summ])]).to_csv(sp, index=False)
            old = pd.read_csv(tp) if os.path.exists(tp) else None
            pd.concat([old, s] if old is not None else [s]).to_csv(tp, index=False)
            print(f"  {n}/{len(jobs)} {summ['climate']} {summ['kinetics']} seed={summ['seed']}: "
                  f"{summ['litter_degraded']:.1f} % degraded, CUE {summ['CUE_annual']:.3f}, "
                  f"E/uptake C {summ['E_max_per_uptake_C']:.0f}", flush=True)


if __name__ == "__main__":
    main()
