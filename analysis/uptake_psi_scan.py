"""
Test of the moisture sensitivity of the uptake (Uptake_psi_slope: uptake factor = exp(slope x (Psi - wp_fc)) below wp_fc; default 0.5) in the minimal model under the
real climate (year 1 of input/climate.csv, 365 days), for the bacterium and the fungus on the C 100 and C 200 substrates, Km slopes 2.5/2.5, Input_NH4 0.02 and
Input_PO4 0.002 (as set in the input folders). The uptake factor at Psi = -5 MPa is exp(-4 x slope).

Usage: python uptake_psi_scan.py ../analysis_output/uptake_psi_scan [--slopes 0.5 0.35 0.25 0.1 0] [--seeds 0 1 2]
Writes OUT/e_daily_all.csv, OUT/metrics.csv, OUT/uptake_psi_scan.png and prints the table.
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
from e_tracker import run  # noqa: E402

plt.rcParams.update({"font.size": 9.5, "axes.titlesize": 10, "axes.labelsize": 9.5, "legend.fontsize": 8.5, "figure.dpi": 150})
CLIM = os.path.join(HERE, "..", "input", "climate.csv")
CASES = {("bacterium", "C 100"): ("input_minimal_e1500_km25", "bac", 100.0), ("bacterium", "C 200"): ("input_minimal_e1500_crich_km25", "bac", 200.0),
         ("fungus", "C 100"): ("input_minimal_e1500_km25_fb1", "fun", 100.0), ("fungus", "C 200"): ("input_minimal_e1500_crich_km25_fb1", "fun", 200.0)}
EMAX = 1400.0


def job(j):
    case, slope, seed = j
    inp = os.path.abspath(os.path.join(HERE, "..", CASES[case][0]))
    df = run((inp, os.path.abspath(CLIM), 1, seed, 20, 1.0, False, {"Uptake_psi_slope": float(slope)}))
    df.insert(0, "slope", slope)
    df.insert(0, "substrate", case[1])
    df.insert(0, "taxon", case[0])
    return df


def metrics(d, g, C0):
    m = d.groupby("doy").mean(numeric_only=True)
    ratio = m[f"E_uptake_{g}"] / m[f"E_cost_{g}"]
    up = m[f"uptake_C_{g}"].sum()
    dry = m.loc[201:293]
    ebio = m[f"E_uptake_{g}"] / m[f"biomass_C_{g}"]
    cell_days = d[f"cells_{g}"].sum()
    return dict(biomass_C_d365=m.loc[365, f"biomass_C_{g}"], degraded_pct=100 * (C0 - m.loc[365, "sub_C_litter"]) / C0,
                pct_celldays_Eshortfall=100 * d[f"cells_Eclip_{g}"].sum() / cell_days, days_Eratio_lt1=int((ratio < 1).sum()), min_Eratio=ratio.min(),
                min_EC_pct_Emax=100 * m[f"EC_{g}"].min() / EMAX, toE_pct_dry=100 * dry[f"toE_C_{g}"].sum() / dry[f"uptake_C_{g}"].sum(),
                retained_pct_dry=100 * (dry[f"uptake_BC_{g}"].sum() - dry[f"enzyme_C_{g}"].sum() - dry[f"overflow_C_{g}"].sum()) / dry[f"uptake_C_{g}"].sum(),
                overflow_pct=100 * m[f"overflow_C_{g}"].sum() / up, CUE_model=(m.CUE * m[f"uptake_C_{g}"]).sum() / up,
                spike_ratio=ebio.loc[285:315].max() / ebio.loc[260:280].mean(), uptake_C_total=up, E_resid=d.resid_E.abs().max())


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("out_dir")
    ap.add_argument("--slopes", type=float, nargs="+", default=[0.5, 0.35, 0.25, 0.1, 0.0])
    ap.add_argument("--seeds", type=int, nargs="+", default=[0, 1, 2])
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    jobs = [(c, s, sd) for c in CASES for s in a.slopes for sd in a.seeds]
    with Pool(min(len(jobs), 10)) as p:
        parts = p.map(job, jobs, chunksize=1)
    daily = pd.concat(parts, ignore_index=True)
    daily.to_csv(os.path.join(a.out_dir, "e_daily_all.csv"), index=False)
    rows = []
    for (t, sb, sl), d in daily.groupby(["taxon", "substrate", "slope"]):
        g, C0 = CASES[(t, sb)][1], CASES[(t, sb)][2]
        r = metrics(d, g, C0)
        r.update(taxon=t, substrate=sb, slope=sl, uptake_factor_at_psi_minus5=float(np.exp(-4 * sl)))
        rows.append(r)
    met = pd.DataFrame(rows)
    met = met[["taxon", "substrate", "slope", "uptake_factor_at_psi_minus5"] + [c for c in met.columns if c not in ("taxon", "substrate", "slope", "uptake_factor_at_psi_minus5")]]
    met = met.sort_values(["taxon", "substrate", "slope"], ascending=[True, True, False])
    met.to_csv(os.path.join(a.out_dir, "metrics.csv"), index=False)
    pd.set_option("display.width", 250)
    pd.set_option("display.max_columns", 30)
    print(met.drop(columns=["uptake_C_total", "E_resid"]).round(3).to_string(index=False))
    # figure: time series for the C 100 substrate (both taxa) and summary
    cm = dict(zip(sorted(a.slopes, reverse=True), plt.cm.viridis(np.linspace(0.05, 0.9, len(a.slopes)))))
    fig, axes = plt.subplots(2, 4, figsize=(19, 8), sharex="col")
    for i, (t, g) in enumerate((("bacterium", "bac"), ("fungus", "fun"))):
        dd = daily[(daily.taxon == t) & (daily.substrate == "C 100")]
        for sl in sorted(a.slopes, reverse=True):
            m = dd[dd.slope == sl].groupby("doy").mean(numeric_only=True)
            lab = f"slope {sl:g} (factor {np.exp(-4 * sl):.2f} at -5 MPa)"
            axes[i][0].plot(m.index, (m[f"E_uptake_{g}"] / m[f"E_cost_{g}"]).rolling(5, center=True, min_periods=1).mean(), color=cm[sl], label=lab)
            axes[i][1].plot(m.index, 100 * m[f"EC_{g}"] / EMAX, color=cm[sl])
            axes[i][2].plot(m.index, 100 * (m[f"toE_C_{g}"] / m[f"uptake_C_{g}"]).rolling(5, center=True, min_periods=1).mean(), color=cm[sl])
        axes[i][0].axhline(1, color="gray", ls=":")
        axes[i][0].set_ylim(0, 4)
        axes[i][0].set_title(f"{t}, C 100: energy produced / spent")
        axes[i][1].set_title("E:C as % of E_max")
        axes[i][2].set_title("% of C taken up routed to energy")
        axes[i][2].set_ylim(40, 100)
        ax = axes[i][3]
        for sb, ls in (("C 100", "-"), ("C 200", "--")):
            x = met[(met.taxon == t) & (met.substrate == sb)].sort_values("slope")
            ax.plot(x.slope, x.pct_celldays_Eshortfall, ls, color="#d6172a", marker="o", ms=4, label=f"{sb}: % of cell-days with an energy shortfall")
        ax.set_ylabel("% of cell-days with an energy shortfall", color="#d6172a")
        ax2 = ax.twinx()
        for sb, ls in (("C 100", "-"), ("C 200", "--")):
            x = met[(met.taxon == t) & (met.substrate == sb)].sort_values("slope")
            ax2.plot(x.slope, x.degraded_pct, ls, color="#8c6d31", marker="s", ms=4, label=f"{sb}: substrate degraded")
        ax2.set_ylabel("Substrate C degraded (%)", color="#8c6d31")
        ax.set_xlabel("Uptake_psi_slope (0.5 = current)")
        ax.set_title(f"{t}: outcome vs slope (solid C 100, dashed C 200)")
        for a_ in axes[i][:3]:
            a_.set_xlabel("Day")
        for a_ in axes[i]:
            a_.spines[["top"]].set_visible(False)
    axes[0][0].legend(frameon=False, fontsize=7.5, loc="upper right")
    fig.suptitle("Lowering the moisture sensitivity of the uptake: minimal model, real climate year 1 (3 seeds)", fontsize=11, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(os.path.join(a.out_dir, "uptake_psi_scan.png"), dpi=150, bbox_inches="tight")


if __name__ == "__main__":
    main()
