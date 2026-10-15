"""
Compare MMRT and Arrhenius kinetics in the full model, with the outputs of full_gradient_resupply.py (HSP vs no-HSP control).

Usage:
    python compare_kinetics_resupply.py MMRT_DIR ARRHENIUS_DIR OUT_DIR
Only the seeds and temperatures present in the Arrhenius run are used from the MMRT run. Writes in OUT_DIR:
  kinetics_timeseries.png   CUE, biomass C and remaining substrate over time (rows) at each temperature (columns)
  kinetics_summary.png      end-of-pulse outcomes vs temperature, MMRT vs Arrhenius, HSP (solid) vs control (dashed)
and prints the numbers.
"""
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

plt.rcParams.update({"font.size": 10, "axes.titlesize": 11, "axes.labelsize": 10, "legend.fontsize": 9, "figure.dpi": 150})
COL = {"MMRT": "#1f77b4", "Arrhenius": "#ff7f0e"}
ROWS = [("CUE", "CUE", (0, 0.7)), ("biomass_C", "Biomass C per grid cell", (0, None)),
        ("substrate_remaining_pct", "Substrate C remaining (% of initial)", (0, 100))]


def load(mmrt, arr):
    a = pd.read_csv(os.path.join(arr, "timeseries.csv"))
    m = pd.read_csv(os.path.join(mmrt, "timeseries.csv"))
    m = m[m["seed"].isin(a["seed"].unique()) & m["Temp"].isin(a["Temp"].unique())]
    a["kin"], m["kin"] = "Arrhenius", "MMRT"
    return pd.concat([m, a], ignore_index=True), sorted(a["Temp"].unique()), a["seed"].nunique()


def end_values(ts):
    e = ts[ts["doy"] == 365].groupby(["kin", "hsp", "pulse", "Temp"])[
        ["biomass_C", "energy_to_C", "alive_rows", "substrate_remaining_pct"]].mean()
    c = ts[ts["doy"] > 335].groupby(["kin", "hsp", "pulse", "Temp"])["CUE"].mean()
    out = e.join(c).reset_index()
    out["degraded"] = 100 - out["substrate_remaining_pct"]
    return out


def timeseries_plot(ts, temps, n, out):
    fig, axes = plt.subplots(3, len(temps), figsize=(3.2 * len(temps), 8.8), sharex=True, squeeze=False)
    fig.suptitle(f"Full model, 365-day resupply: MMRT (blue) vs Arrhenius (orange), HSP (solid) vs no-HSP control (dashed), "
                 f"mean of {n} communities, 7-day mean", fontsize=11, fontweight="bold")
    for i, (col, lab, yl) in enumerate(ROWS):
        for j, t in enumerate(temps):
            ax = axes[i][j]
            for kin in COL:
                for hsp, ls in ((False, "--"), (True, "-")):
                    d = ts[(ts["Temp"] == t) & (ts["hsp"] == hsp) & (ts["kin"] == kin)].pivot(index="day", columns="seed", values=col)
                    ax.plot(d.index, d.mean(axis=1).rolling(7, center=True, min_periods=1).mean(), color=COL[kin], ls=ls, lw=1.4,
                            label=f"{kin} {'HSP' if hsp else 'control'}")
            for x in (365, 730):
                ax.axvline(x, color="lightgray", lw=0.8)
            ax.set_ylim(*yl)
            if i == 0:
                ax.set_title(f"{t:g}°C")
            if j == 0:
                ax.set_ylabel(lab)
            if i == len(ROWS) - 1:
                ax.set_xlabel("Day")
            ax.spines[["top", "right"]].set_visible(False)
    h, l = axes[0][0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncols=4, framealpha=0.85, bbox_to_anchor=(0.5, -0.01))
    fig.tight_layout()
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)


def summary_plot(e, n, out):
    fig, axes = plt.subplots(2, 4, figsize=(16, 7.5), sharex=True)
    fig.suptitle(f"MMRT vs Arrhenius, full model with resupply, day 365 of pulse 1 (top) and mean of pulses 2-3 (bottom); "
                 f"{n} communities", fontsize=11, fontweight="bold")
    for r, (name, pulses) in enumerate((("pulse 1", [1]), ("pulses 2-3", [2, 3]))):
        g = e[e["pulse"].isin(pulses)].groupby(["kin", "hsp", "Temp"], as_index=False).mean(numeric_only=True)
        for ax, (col, title) in zip(axes[r], (("biomass_C", "Biomass C"), ("rel", "Biomass change, HSP vs control (%)"),
                                              ("degraded", "Substrate degraded (% of initial C)"), ("CUE", "CUE (last 30 days)"))):
            for kin in COL:
                if col == "rel":
                    h = g[(g["kin"] == kin) & g["hsp"]].set_index("Temp")["biomass_C"]
                    c = g[(g["kin"] == kin) & ~g["hsp"]].set_index("Temp")["biomass_C"]
                    ax.plot(h.index, 100 * (h - c) / c, color=COL[kin], marker="o", ms=4, lw=1.6, label=kin)
                    continue
                for hsp, ls in ((True, "-"), (False, "--")):
                    d = g[(g["kin"] == kin) & (g["hsp"] == hsp)].sort_values("Temp")
                    ax.plot(d["Temp"], d[col], color=COL[kin], ls=ls, marker="o", ms=4, lw=1.5,
                            label=f"{kin} {'HSP' if hsp else 'control'}" if r == 0 and col == "biomass_C" else None)
            if col == "rel":
                ax.axhline(0, color="gray", lw=0.8, ls="--")
            ax.set_title(f"{title}, {name}", fontsize=10)
            ax.spines[["top", "right"]].set_visible(False)
            if r == 1:
                ax.set_xlabel("Temperature (°C)")
    h, l = axes[0][0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncols=4, framealpha=0.85, bbox_to_anchor=(0.5, -0.02))
    fig.tight_layout()
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)


def main(mmrt, arr, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    ts, temps, n = load(mmrt, arr)
    timeseries_plot(ts, temps, n, os.path.join(out_dir, "kinetics_timeseries.png"))
    e = end_values(ts)
    summary_plot(e, n, os.path.join(out_dir, "kinetics_summary.png"))
    pd.set_option("display.width", 220)
    for p in (1, 2, 3):
        t = e[e["pulse"] == p].pivot_table(index="Temp", columns=["kin", "hsp"], values=["biomass_C", "degraded", "CUE", "energy_to_C"])
        o = pd.DataFrame(index=t.index)
        for kin in COL:
            o[f"{kin[:3]}_bio_%"] = 100 * (t["biomass_C"][(kin, True)] - t["biomass_C"][(kin, False)]) / t["biomass_C"][(kin, False)]
            o[f"{kin[:3]}_bioH"] = t["biomass_C"][(kin, True)]
            o[f"{kin[:3]}_bioC"] = t["biomass_C"][(kin, False)]
            o[f"{kin[:3]}_degH"] = t["degraded"][(kin, True)]
            o[f"{kin[:3]}_degC"] = t["degraded"][(kin, False)]
            o[f"{kin[:3]}_CUEH"] = t["CUE"][(kin, True)]
            o[f"{kin[:3]}_CUEC"] = t["CUE"][(kin, False)]
            o[f"{kin[:3]}_ECC"] = t["energy_to_C"][(kin, False)]
        print(f"\npulse {p}, day 365 (MMR = MMRT, Arr = Arrhenius; bio_% = HSP vs control; H = HSP, C = control)")
        print(o.round(2).to_string())


if __name__ == "__main__":
    main(*sys.argv[1:4])
