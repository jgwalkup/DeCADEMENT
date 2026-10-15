"""
Compare real-climate pulse runs (climate_run.py output) across warming offsets, HSP vs no-HSP control.

Usage:
    python climate_warming_summary.py ../analysis_output/climate_warming \\
        0=../analysis_output/climate_run_e1500_wetpsi 2=../analysis_output/climate_run_e1500_wetpsi_plus2 \\
        4=../analysis_output/climate_run_e1500_wetpsi_plus4 6=../analysis_output/climate_run_e1500_wetpsi_plus6
Writes <out_dir>/warming_summary.csv, warming_summary.png (end-of-pulse outcomes vs warming) and
warming_timeseries.png (within-year trajectories, one column per offset).
"""
import argparse
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from scenario_utils import scenario_color_map  # noqa: E402

plt.rcParams.update({
    "font.size": 10, "axes.titlesize": 11, "axes.labelsize": 10,
    "legend.fontsize": 9, "figure.dpi": 150,
})
BLUE = "#1f77b4"


def load(runs):
    parts = []
    for off, d in runs.items():
        t = pd.read_csv(os.path.join(d, "timeseries.csv"))
        t["offset"] = off
        parts.append(t)
    return pd.concat(parts, ignore_index=True)


def per_pulse(ts):
    last = ts[ts["doy"] == 365].set_index(["offset", "seed", "pulse", "hsp"])
    cue = ts[ts["doy"] > 335].groupby(["offset", "seed", "pulse", "hsp"])["CUE"].mean()
    out = last[["biomass_C", "energy_to_C", "alive_cells", "substrate_remaining_pct"]].join(cue)
    out["deg"] = 100 - out["substrate_remaining_pct"]
    return out.reset_index()


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("out_dir")
    ap.add_argument("runs", nargs="+", help="offset=dir pairs")
    a = ap.parse_args()
    runs = {float(r.split("=")[0]): r.split("=")[1] for r in a.runs}
    os.makedirs(a.out_dir, exist_ok=True)
    ts = load(runs)
    pp = per_pulse(ts)
    offs = sorted(runs)

    # paired HSP vs control per (offset, seed, pulse)
    piv = pp.pivot_table(index=["offset", "seed", "pulse"], columns="hsp",
                         values=["biomass_C", "energy_to_C", "CUE", "deg", "alive_cells"])
    rel = 100 * (piv[("biomass_C", True)] - piv[("biomass_C", False)]) / piv[("biomass_C", False)]
    dabs = piv[("biomass_C", True)] - piv[("biomass_C", False)]
    rows = []
    for o in offs:
        r = {"offset": o}
        for col in ("biomass_C", "energy_to_C", "CUE", "deg", "alive_cells"):
            for hsp, lab in ((True, "HSP"), (False, "ctrl")):
                v = piv.loc[o][(col, hsp)]
                r[f"{col}_{lab}"] = v.mean()
                r[f"{col}_{lab}_sd"] = v.std()
        r["rel_biomass_mean"] = rel.loc[o].mean()
        r["rel_biomass_sd"] = rel.loc[o].std()
        r["abs_biomass_diff_mean"] = dabs.loc[o].mean()
        rows.append(r)
    summ = pd.DataFrame(rows).set_index("offset")
    summ.to_csv(os.path.join(a.out_dir, "warming_summary.csv"))
    print(summ[["biomass_C_HSP", "biomass_C_ctrl", "rel_biomass_mean", "rel_biomass_sd", "abs_biomass_diff_mean",
                "CUE_HSP", "CUE_ctrl", "energy_to_C_HSP", "energy_to_C_ctrl", "deg_HSP", "deg_ctrl",
                "alive_cells_HSP", "alive_cells_ctrl"]].round(2).T.to_string())

    # ---- summary figure ----
    fig, axes = plt.subplots(2, 3, figsize=(14, 8))
    fig.suptitle("Warming scenarios (real daily temperature + offset, Ψ = -0.5 MPa): end-of-pulse outcomes, "
                 "mean ± SD over 5 seeds × 6 pulses", fontsize=11, fontweight="bold")
    panels = [("biomass_C", "Biomass C (day 365)", "C per grid cell"), ("CUE", "CUE (last 30 days)", "CUE"),
              ("energy_to_C", "Energy : C (day 365)", "living biomass"), ("deg", "Substrate degraded", "% of initial C"),
              ("alive_cells", "Cells alive (day 365)", "of 400")]
    for ax, (col, title, yl) in zip(axes.flat, panels):
        for lab, c, ls in (("HSP", BLUE, "-"), ("ctrl", "gray", "--")):
            ax.errorbar(offs, [summ.loc[o, f"{col}_{lab}"] for o in offs], yerr=[summ.loc[o, f"{col}_{lab}_sd"] for o in offs],
                        color=c, ls=ls, marker="o", ms=5, capsize=3, lw=1.6, label="HSP" if lab == "HSP" else "no HSP (control)")
        ax.set_title(title, fontsize=10)
        ax.set_ylabel(yl)
        if col == "CUE":
            ax.set_ylim(0, 1)
        ax.spines[["top", "right"]].set_visible(False)
    ax = axes.flat[5]
    ax.errorbar(offs, summ["rel_biomass_mean"], yerr=summ["rel_biomass_sd"], color=BLUE, marker="o", ms=5, capsize=3,
                lw=1.6)
    ax.axhline(0, color="gray", lw=0.8, ls="--")
    ax.set_title("Biomass change with HSP vs control", fontsize=10)
    ax.set_ylabel("% (paired by seed and pulse)")
    ax.spines[["top", "right"]].set_visible(False)
    for ax in axes[1]:
        ax.set_xlabel("Warming offset (°C)")
    for ax in axes[0]:
        ax.set_xlabel("Warming offset (°C)")
    axes[0][0].legend()
    fig.tight_layout()
    fig.savefig(os.path.join(a.out_dir, "warming_summary.png"), dpi=150, bbox_inches="tight")
    plt.close(fig)

    # ---- within-year trajectories ----
    m = ts.groupby(["offset", "hsp", "doy"], as_index=False)[["CUE", "energy_to_C", "biomass_C", "Temp"]].mean()
    rows_ = [("Temp", "Temperature (°C)", None), ("CUE", "CUE", (0, 1)), ("energy_to_C", "Energy : C", (0, None)),
             ("biomass_C", "Biomass C", (0, None))]
    fig, axes = plt.subplots(len(rows_), len(offs), figsize=(3.6 * len(offs), 10), sharex=True, squeeze=False)
    fig.suptitle("Within-year trajectories by warming offset: HSP (solid) vs no-HSP control (dashed), mean over seeds and pulses",
                 fontsize=11, fontweight="bold")
    for j, o in enumerate(offs):
        for i, (col, lab, yl) in enumerate(rows_):
            ax = axes[i][j]
            for hsp, ls, c in ((False, "--", "gray"), (True, "-", BLUE)):
                d = m[(m["offset"] == o) & (m["hsp"] == hsp)].sort_values("doy")
                y = d[col].rolling(7, center=True, min_periods=1).mean() if col in ("CUE", "energy_to_C") else d[col]
                if col == "Temp":
                    if hsp:
                        ax.plot(d["doy"], y, color="#d6172a", lw=1.2)
                else:
                    ax.plot(d["doy"], y, color=c, ls=ls, lw=1.5, label="HSP" if hsp else "no HSP (control)")
            if yl:
                ax.set_ylim(*yl)
            if i == 0:
                ax.set_title(f"+{o:g}°C")
            if j == 0:
                ax.set_ylabel(lab)
            if i == len(rows_) - 1:
                ax.set_xlabel("Day of pulse")
            ax.spines[["top", "right"]].set_visible(False)
    h, l = axes[1][0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncols=2, framealpha=0.85, bbox_to_anchor=(0.5, -0.02))
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(os.path.join(a.out_dir, "warming_timeseries.png"), dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("wrote", a.out_dir)


if __name__ == "__main__":
    main()
