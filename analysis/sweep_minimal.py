"""
Sweep one parameter of the minimal scenario and tabulate outcomes (in-process, temp copies only).

Usage (from anywhere):
    python sweep_minimal.py ../input_minimal Maint_E_Rate_Const 0 1 2 4 8 16
    python sweep_minimal.py ../input_minimal Enz_Prod 0.001 0.002 --set Uptake_C_cost=0.0375 --csv out.csv

PARAM may be a parameters.csv row name, or the stem of a NAME_min / NAME_max pair (both are set).
Every run: single pulse, side x side grid, fixed seed, input folder copied to a temp dir.
"""
import argparse
import os
import shutil
import sys
import tempfile

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "src"))
from initialization import initialize_data  # noqa: E402
from grid import Grid  # noqa: E402


def set_params(path, changes):
    lines = open(path).read().split("\n")
    names = {l.split(",")[0] for l in lines if "," in l}
    expanded = {}
    for k, v in changes.items():
        if k in names:
            expanded[k] = v
        elif k + "_min" in names and k + "_max" in names:
            expanded[k + "_min"] = expanded[k + "_max"] = v
        else:
            raise KeyError(f"{k} not found in parameters.csv (nor {k}_min/{k}_max)")
    out = []
    for l in lines:
        n = l.split(",")[0]
        out.append(f"{n},{float(expanded[n])!r}" if "," in l and n in expanded else l)
    open(path, "w").write("\n".join(out))


def set_runtime(path, changes):
    lines = open(path).read().split("\n")
    out = []
    for l in lines:
        n = l.split("\t")[0]
        out.append(f"{n}\t{changes[n]}" if "\t" in l and n in changes else l)
    open(path, "w").write("\n".join(out))


def set_climate_temp(path, temp):
    """Overwrite the Temp column of a flat climate.csv with a new constant value,
    leaving Psi/FM/Date (row count too) untouched -- mirrors run_multi_climate_simulations.py's
    --static-temp, but in place on an already-flat template rather than building a new file."""
    df = pd.read_csv(path, index_col=0)
    df["Temp"] = float(temp)
    df.to_csv(path)


def run_one(src_input, changes, side, days, seed, temp=None):
    tmp = tempfile.mkdtemp()
    work = os.path.join(tmp, "in")
    shutil.copytree(src_input, work)
    set_params(os.path.join(work, "parameters.csv"), changes)
    set_runtime(os.path.join(work, "runtime.txt"),
                dict(x=side, y=side, gridsize=side * side, pulse=1))
    if temp is not None:
        set_climate_temp(os.path.join(work, "climate.csv"), temp)
    cwd = os.getcwd()
    os.chdir(work)
    try:
        np.random.seed(seed)
        rt = pd.read_csv("runtime.txt", header=None, index_col=0, sep="\t").drop("output_format").astype(float)
        data = initialize_data(rt)
        eco = Grid(rt, data)
        G = side * side
        sub = [s for s in data["Substrate_index"] if s not in ("DeadMic", "DeadEnz")][0]
        sub0 = float(np.asarray(eco.Substrates.loc[sub, "C"]).sum())
        resp_cum, ec_min, zero_prod = 0.0, np.inf, 0
        tail = {k: [] for k in ("E_prod", "E_cost", "maint", "CUE", "HSP_prop", "maint_Ea", "maint_rate_percap")}
        series = []
        mon = "Mon_" + sub
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
            n = int(alive.sum())
            resp_cum += float(np.asarray(eco.Taxon_Resp).sum()) / G
            zero_prod += float(np.asarray(eco.Enzyme_Production).sum()) <= 0
            nn = max(n, 1)
            series.append(dict(
                day=i + 1,
                substrate_remaining_pct=100 * ps(eco.Substrates.loc[sub, "C"]) / sub0,
                decay_rate=ps(eco.DecayRates) / G,
                biomass_C=C.sum() / G,
                alive_cells=n,
                biomass_NC=np.asarray(eco.Microbes["N"])[alive].sum() / C[alive].sum() if n else np.nan,
                biomass_PC=np.asarray(eco.Microbes["P"])[alive].sum() / C[alive].sum() if n else np.nan,
                energy_to_C=E[alive].sum() / C[alive].sum() if n else np.nan,
                enzyme_pool=ps(eco.Enzymes) / G,
                enzyme_production=ps(eco.Enzyme_Production) / G,
                monomer_pool=ps(eco.Monomers.loc[mon, "C"]) / G,
                C_uptake=ps(eco.Taxon_Uptake_C) / nn,
                C_to_biomass=ps(eco.Taxon_Uptake_BC) / nn,
                E_production=ps(eco.Taxon_Uptake_E) / nn,
                E_cost_total=ps(eco.Microbe_E_Cost) / nn,
                E_cost_maintenance=ps(eco.Tax_cell_maint_E_cost) / nn,
                net_energy=ps(eco.Microbe_E_Gain) / nn,
                CUE=float(eco.CUE_system),
                CAE=float(eco.CAE_system),
                respiration=float(np.asarray(eco.Taxon_Resp).sum()) / G,
                HSP_prop=float(np.asarray(eco.HSP_prop)[alive].mean()) if n else np.nan,
                maint_Ea_effective=float(np.asarray(eco.maint_E_cost_Ea)[alive].mean()) if n else np.nan,
                maint_rate_percap=float(np.asarray(eco.tax_maint_E_cost)[alive].mean()) if n else np.nan,
            ))
            if n:
                ec_min = min(ec_min, E[alive].sum() / C[alive].sum())
                if i >= days - 30:
                    tail["E_prod"].append(float(np.asarray(eco.Taxon_Uptake_E).sum()) / n)
                    tail["E_cost"].append(float(np.asarray(eco.Microbe_E_Cost).sum()) / n)
                    tail["maint"].append(float(np.asarray(eco.Tax_cell_maint_E_cost).sum()) / n)
                    tail["CUE"].append(float(eco.CUE_system))
                    tail["HSP_prop"].append(float(np.asarray(eco.HSP_prop)[alive].mean()))
                    tail["maint_Ea"].append(float(np.asarray(eco.maint_E_cost_Ea)[alive].mean()))
                    tail["maint_rate_percap"].append(float(np.asarray(eco.tax_maint_E_cost)[alive].mean()))
        C = np.asarray(eco.Microbes["C"])
        E = np.asarray(eco.Microbes["E"])
        alive = C > 0
        sub_end = float(np.asarray(eco.Substrates.loc[sub, "C"]).sum())
        m = {k: (np.nanmean(v) if v else np.nan) for k, v in tail.items()}
        summary = dict(
            deg_pct=100 * (1 - sub_end / sub0),
            alive=int(alive.sum()),
            biomC=float(C.sum()) / G,
            EC_end=float(E[alive].sum() / C[alive].sum()) if alive.any() else np.nan,
            EC_min=ec_min,
            enz_pool=float(np.asarray(eco.Enzymes).sum()) / G,
            E_prod=m["E_prod"],
            E_cost=m["E_cost"],
            maint=m["maint"],
            net_E=m["E_prod"] - m["E_cost"],
            CUE=m["CUE"],
            resp_cum=resp_cum,
            zero_enz_prod_days=zero_prod,
            HSP_prop=m["HSP_prop"],
            maint_Ea=m["maint_Ea"],
            maint_rate_percap=m["maint_rate_percap"],
        )
        return summary, series
    finally:
        os.chdir(cwd)
        shutil.rmtree(tmp, ignore_errors=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("input_folder")
    ap.add_argument("param")
    ap.add_argument("values", type=float, nargs="+")
    ap.add_argument("--set", action="append", default=[], metavar="NAME=VALUE",
                    help="fixed override applied to every run (repeatable)")
    ap.add_argument("--side", type=int, default=20)
    ap.add_argument("--days", type=int, default=365)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--csv", default=None, help="save the end-of-run summary table")
    ap.add_argument("--timeseries", default=None,
                    help="save daily time series of every run (long format: one row per value x day)")
    a = ap.parse_args()

    src = os.path.abspath(a.input_folder)
    fixed = {k: float(v) for k, v in (s.split("=") for s in a.set)}
    print(f"# {a.param} sweep | input={src} | grid={a.side}x{a.side} | days={a.days} | seed={a.seed} | fixed={fixed}")
    print("# E_prod/E_cost/maint/net_E: per living cell per day, mean of last 30 days; CUE: mean of last 30 days")
    rows, ts = [], []
    for v in a.values:
        summary, series = run_one(src, {**fixed, a.param: v}, a.side, a.days, a.seed)
        rows.append({a.param: v, **summary})
        ts.append(pd.DataFrame(series).assign(**{a.param: v}))
        print(f"  done {a.param}={v:g}", file=sys.stderr, flush=True)
    df = pd.DataFrame(rows).set_index(a.param)
    print(df.to_string(float_format=lambda x: f"{x:.4g}"))
    if a.csv:
        df.to_csv(a.csv)
    if a.timeseries:
        pd.concat(ts, ignore_index=True)[[a.param, "day"] + [c for c in ts[0].columns if c not in ("day", a.param)]].to_csv(
            a.timeseries, index=False)


if __name__ == "__main__":
    main()
