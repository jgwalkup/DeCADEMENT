"""
Scan of the ongoing NH4 and PO4 inputs (Input_NH4 as N, Input_PO4 as P per grid cell per day) in the minimal model, to find the input at which the
N limitation of the bacterium is relieved. Km slopes 2.5/2.5, flat 15 C, 365 days.

Stage nh4: Input_NH4 on a log scale (and 0), Input_PO4 = 0.09 x Input_NH4, for the C 100 and C 200 substrates.
Stage po4: Input_PO4 = multiplier x (0.09 x Input_NH4) at a chosen Input_NH4 per substrate, to check for P limitation.
Relieved = C released by cells short of N below 2 % of the C uptake (year total) and the daily mean biomass N:C above the trigger on all days 31-365.

Usage:
    python nh4_scan.py nh4 ../analysis_output/nh4_scan_minimal [--seeds 0 1 2]
    python nh4_scan.py po4 ../analysis_output/nh4_scan_minimal --nh4 ref=0.003 crich=0.008 [--mult 0 0.25 0.5 1 2 4]
Writes OUT/<stage>_metrics.csv, OUT/<stage>_daily.csv and OUT/<stage>.png and prints the table.
"""
import argparse
import os
import sys
from multiprocessing import Pool

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from n_tracker import run  # noqa: E402

plt.rcParams.update({"font.size": 9.5, "axes.titlesize": 10, "axes.labelsize": 9.5, "legend.fontsize": 8.5, "figure.dpi": 150})
SUBS = {"ref": ("input_minimal_e1500_km25", 100.0, "C 100 (C:N 20)"), "crich": ("input_minimal_e1500_crich_km25", 200.0, "C 200 (C:N 40)")}
CLIMATE = os.path.join(HERE, "..", "input_minimal_e1500", "climate.csv")
P_PER_N = 0.09
TRIGGER = 0.1382  # N:C at which N / (C+N+P) = 0.12 for the bacterium (P:C 0.0135)
NH4_VALUES = [0, 0.0005, 0.0008, 0.00125, 0.002, 0.003, 0.005, 0.008, 0.0125, 0.02]


def job_run(job):
    sub, nh4, po4, seed = job
    df = run((os.path.abspath(os.path.join(HERE, "..", SUBS[sub][0])), os.path.abspath(CLIMATE), 1, seed, 20, [0, 0], [0, 0], 1.0, True,
              {"Input_NH4": float(nh4), "Input_PO4": float(po4)}))
    df.insert(0, "po4", po4)
    df.insert(0, "nh4", nh4)
    df.insert(0, "sub", sub)
    return df


def metrics(d, init_C):
    m = d.groupby("doy").mean(numeric_only=True)
    up = m.uptake_C_bac.sum()
    nc = m.mic_N_bac / m.mic_C_bac
    late = m.loc[31:]
    cells = d[d.doy >= 31]
    leach_n, in_n = d.groupby("seed").leach_N.sum().mean(), d.groupby("seed").input_N.sum().mean()
    leach_p, in_p = d.groupby("seed").leach_P.sum().mean(), d.groupby("seed").input_P.sum().mean()
    out = dict(
        overflow_pct=100 * m.overflow_C_bac.sum() / up, overflow_N_pct=100 * m.overflow_C_limN_bac.sum() / up, overflow_P_pct=100 * m.overflow_C_limP_bac.sum() / up,
        overflow_late_pct=100 * m.loc[201:].overflow_C_bac.sum() / m.loc[201:].uptake_C_bac.sum(),
        minNC_d31_365=(late.mic_N_bac / late.mic_C_bac).min(), NC_d365=nc.loc[365], pct_celldays_below_N_trigger=100 * cells.cells_belowN_bac.sum() / cells.cells_bac.sum(),
        pct_celldays_below_P_trigger=100 * cells.cells_belowP_bac.sum() / cells.cells_bac.sum(),
        biomass_C_d365=m.mic_C_bac.loc[365], substrate_degraded_pct=100 * (init_C - m.sub_C_litter.loc[365]) / init_C,
        N_uptake=m.uptake_N_bac.sum(), P_uptake=m.uptake_P_bac.sum(), CUE_model=(m.CUE * m.uptake_C_bac).sum() / up,
        input_N=in_n, leached_N=leach_n, leach_frac_N=leach_n / in_n if in_n > 0 else np.nan,
        input_P=in_p, leached_P=leach_p, leach_frac_P=leach_p / in_p if in_p > 0 else np.nan,
        NH4_pool_mean=m.mon_N_nh4.mean(), PO4_pool_mean=m.mon_P_po4.mean(), max_resid_N=d.resid_N.abs().max(), max_resid_P=d.resid_P.abs().max())
    out["effective_CUE"] = out["CUE_model"] - out["overflow_pct"] / 100
    out["relieved"] = bool(out["overflow_N_pct"] < 2 and out["minNC_d31_365"] >= TRIGGER)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("stage", choices=["nh4", "po4"])
    ap.add_argument("out_dir")
    ap.add_argument("--seeds", type=int, nargs="+", default=[0, 1, 2])
    ap.add_argument("--nh4", nargs="*", default=[], metavar="SUB=VALUE", help="po4 stage: Input_NH4 per substrate, e.g. ref=0.003 crich=0.008")
    ap.add_argument("--mult", type=float, nargs="+", default=[0, 0.25, 0.5, 1, 2, 4])
    ap.add_argument("--values", nargs="*", default=[], metavar="SUB=v1,v2,...", help="nh4 stage: Input_NH4 values per substrate instead of the default list")
    ap.add_argument("--tag", default="", help="suffix for the output file names")
    ap.add_argument("--procs", type=int, default=10)
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    if a.stage == "nh4":
        custom = {k: [float(x) for x in v.split(",")] for k, v in (x.split("=") for x in a.values)}
        conds = [(sub, v, P_PER_N * v) for sub in SUBS for v in custom.get(sub, [] if custom else NH4_VALUES)]
    else:
        chosen = {k: float(v) for k, v in (x.split("=") for x in a.nh4)}
        conds = [(sub, chosen[sub], P_PER_N * chosen[sub] * mu) for sub in SUBS for mu in a.mult]
    jobs = [(s, n, p, seed) for s, n, p in conds for seed in a.seeds]
    with Pool(min(len(jobs), a.procs)) as pool:
        parts = pool.map(job_run, jobs, chunksize=1)
    daily = pd.concat(parts, ignore_index=True)
    daily.to_csv(os.path.join(a.out_dir, f"{a.stage}{a.tag}_daily.csv"), index=False)
    rows = []
    for (sub, nh4, po4), d in daily.groupby(["sub", "nh4", "po4"]):
        r = metrics(d, SUBS[sub][1])
        r.update(sub=sub, nh4=nh4, po4=po4, po4_mult=po4 / (P_PER_N * nh4) if nh4 > 0 else np.nan)
        rows.append(r)
    met = pd.DataFrame(rows)
    # leaching attributable to the input = leaching minus the leaching of the zero-input run of the same substrate (P released by cells also leaches)
    bpath = os.path.join(a.out_dir, "nh4_metrics.csv")
    ref_tab = met[met.nh4 == 0] if (met.nh4 == 0).any() else (pd.read_csv(bpath).query("nh4 == 0") if os.path.exists(bpath) else met.iloc[0:0])
    base = ref_tab.set_index("sub")[["leached_N", "leached_P"]]
    for k in ("N", "P"):
        met[f"leach_attrib_pct_{k}"] = [100 * (r[f"leached_{k}"] - base.loc[r["sub"], f"leached_{k}"]) / r[f"input_{k}"]
                                         if (r["sub"] in base.index and r[f"input_{k}"] > 0) else np.nan for _, r in met.iterrows()]
    met = met[["sub", "nh4", "po4", "po4_mult"] + [c for c in met.columns if c not in ("sub", "nh4", "po4", "po4_mult")]].sort_values(["sub", "nh4", "po4"])
    met.to_csv(os.path.join(a.out_dir, f"{a.stage}{a.tag}_metrics.csv"), index=False)
    pd.set_option("display.width", 260)
    pd.set_option("display.max_columns", 40)
    show = ["sub", "nh4", "po4", "overflow_N_pct", "overflow_P_pct", "minNC_d31_365", "NC_d365", "pct_celldays_below_N_trigger", "biomass_C_d365",
            "substrate_degraded_pct", "effective_CUE", "leach_attrib_pct_N", "leach_attrib_pct_P", "relieved"]
    print(met[show].round(4).to_string(index=False))
    # figure
    xkey = "nh4" if a.stage == "nh4" else "po4_mult"
    fig, axes = plt.subplots(2, 4, figsize=(19, 8))
    for i, sub in enumerate(SUBS):
        t = met[met["sub"] == sub].sort_values(xkey)
        x = t[xkey].replace(0, np.nan) if a.stage == "nh4" else t[xkey]
        xs = np.where(t[xkey] == 0, t[xkey][t[xkey] > 0].min() / 3 if (t[xkey] > 0).any() else 1e-4, t[xkey]) if a.stage == "nh4" else t[xkey].to_numpy()
        ax = axes[i][0]
        ax.plot(xs, t.overflow_N_pct, "-o", color="#d6172a", ms=4, label="released, N-limited")
        ax.plot(xs, t.overflow_P_pct, "-s", color="#9467bd", ms=4, label="released, P-limited")
        ax.axhline(2, color="gray", ls=":")
        ax.set_ylabel("% of C taken up")
        ax.set_title(f"{SUBS[sub][2]}: C released by cells short of N or P")
        ax = axes[i][1]
        ax.plot(xs, t.minNC_d31_365, "-o", color="#1f77b4", ms=4, label="lowest daily N:C (days 31-365)")
        ax.plot(xs, t.NC_d365, "-s", color="#2ca02c", ms=4, label="N:C on day 365")
        ax.axhline(TRIGGER, color="gray", ls=":", label="trigger 0.138")
        ax.set_ylabel("Biomass N:C")
        ax.set_title("Biomass N:C")
        ax = axes[i][2]
        ax.plot(xs, t.biomass_C_d365, "-o", color="k", ms=4, label="biomass C, day 365")
        ax2 = ax.twinx()
        ax2.plot(xs, t.substrate_degraded_pct, "-s", color="#8c6d31", ms=4)
        ax2.set_ylabel("Substrate C degraded (%)", color="#8c6d31")
        ax.set_ylabel("Biomass C")
        ax.set_title("Growth")
        ax = axes[i][3]
        ax.plot(xs, t.leach_attrib_pct_N, "-o", color="#1f77b4", ms=4, label="N leached because of the input, % of the N input")
        ax.plot(xs, t.leach_attrib_pct_P, "-s", color="#9467bd", ms=4, label="P leached because of the input, % of the P input")
        ax.set_ylabel("% of the input")
        ax.set_ylim(0, 105)
        ax.set_title("Fate of the input")
        for a_ in axes[i]:
            a_.spines[["top"]].set_visible(False)
            if a.stage == "nh4":
                a_.set_xscale("log")
            a_.set_xlabel("Input_NH4 (N per cell per day; 0 plotted at the left)" if a.stage == "nh4" else "Input_PO4 as a multiple of 0.09 x Input_NH4")
    for a_ in axes[0]:
        a_.legend(frameon=False, fontsize=7.5)
    fig.suptitle(f"Minimal bacterium, Km slopes 2.5/2.5, flat 15 C: NH4 and PO4 input scan ({a.stage} stage, {len(a.seeds)} seeds)", fontsize=11, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(os.path.join(a.out_dir, f"{a.stage}{a.tag}.png"), dpi=150, bbox_inches="tight")


if __name__ == "__main__":
    main()
