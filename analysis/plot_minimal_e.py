"""
Energy balance and C allocation of the minimal model for two substrate stoichiometries, from e_tracker.py output (same layout as plot_minimal_n.py).
Usage: python plot_minimal_e.py ../analysis_output/minimal_e_km25_fb1 [bac|fun]   (expects ref/e_daily.csv and crich/e_daily.csv)
Writes minimal_e.png and prints a table.
"""
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

plt.rcParams.update({"font.size": 9.5, "axes.titlesize": 10, "axes.labelsize": 9.5, "legend.fontsize": 8.5, "figure.dpi": 150})


def _pop(flag):
    if flag in sys.argv:
        i = sys.argv.index(flag)
        v = sys.argv[i + 1]
        del sys.argv[i:i + 2]
        return v
    return None


CLIM = _pop("--climate")      # optional climate csv: adds a strip with temperature and water potential on top
LABEL = _pop("--label") or "flat 15 C"


def make_axes(figsize, nrows=2):
    figsize = (figsize[0], figsize[1] * nrows / 2)
    if not CLIM:
        return plt.subplots(nrows, 4, figsize=figsize, sharex=True, squeeze=False)
    fig = plt.figure(figsize=(figsize[0], figsize[1] + 1.6))
    gs = fig.add_gridspec(nrows + 1, 4, height_ratios=[0.38] + [1] * nrows)
    top = fig.add_subplot(gs[0, :])
    c = pd.read_csv(CLIM, index_col=0).iloc[:365]
    top.plot(np.arange(1, len(c) + 1), c["Temp"], color="#d6172a", lw=1.1)
    top.set_ylabel("Temperature (C)")
    t2 = top.twinx()
    t2.plot(np.arange(1, len(c) + 1), c["Psi"], color="#2166ac", lw=1.0)
    t2.set_ylabel("Psi (MPa)")
    top.set_xlim(0, 366)
    top.spines[["top"]].set_visible(False)
    axes = np.array([[fig.add_subplot(gs[r + 1, k]) for k in range(4)] for r in range(nrows)])
    for a in axes.flat:
        a.set_xlim(0, 366)
    return fig, axes
root = sys.argv[1]
G = sys.argv[2] if len(sys.argv) > 2 else "bac"
NAME = {"bac": "bacterium", "fun": "fungus"}[G]
EMAX = 1400.0
_LAB = {"ref": "IdealSub C 100, N 5, P 0.5 (C:N 20)", "crich": "IdealSub C 200, N 5, P 0.5 (C:N 40)"}
runs = {_LAB[k]: pd.read_csv(os.path.join(root, k, "e_daily.csv")) for k in ("ref", "crich") if os.path.exists(os.path.join(root, k, "e_daily.csv"))}
COST = [("maint", "maintenance", "#7f7f7f"), ("transport", "transport (active + transporters)", "#17becf"), ("enzyme", "enzymes", "#bcbd22"),
        ("hsp", "heat shock proteins", "#d6172a"), ("osmo", "osmolytes", "#9467bd")]


def sm(x, w=7):
    return x.rolling(w, center=True, min_periods=1).mean()


def shares(m):
    up = m[f"uptake_C_{G}"].replace(0, np.nan)
    toE = m[f"toE_C_{G}"] / up
    enz = m[f"enzyme_C_{G}"] / up
    ov = m[f"overflow_C_{G}"] / up
    kept = (m[f"uptake_BC_{G}"] - m[f"enzyme_C_{G}"] - m[f"overflow_C_{G}"]) / up
    return toE, enz, ov, kept


fig, axes = make_axes((19, 8), len(runs))
for i, (lab, d) in enumerate(runs.items()):
    m = d.groupby("doy").mean(numeric_only=True)
    ax = axes[i][0]
    ax.plot(m.index, m[f"EC_{G}"], color="k", label="E:C of the biomass")
    ax.axhline(EMAX, color="k", ls=":", lw=1)
    ax.set_ylim(0, 1500)
    ax.set_ylabel("E : C (dotted = E_max 1400)")
    ax2 = ax.twinx()
    ax2.plot(m.index, m[f"E_demand_{G}"], color="#1f77b4", label="fraction of uptake C routed to energy")
    ax2.set_ylim(0, 1)
    ax2.set_ylabel("Fraction of uptake C routed to energy", color="#1f77b4")
    ax.set_title(f"{lab}\nEnergy status of the cells")
    ax = axes[i][1]
    bio = m[f"biomass_C_{G}"].replace(0, np.nan)
    ax.stackplot(m.index, *[sm(m[f"E_{k}_{G}"] / bio) for k, _, _ in COST], colors=[c for _, _, c in COST], labels=[l for _, l, _ in COST])
    ax.plot(m.index, sm(m[f"E_uptake_{G}"] / bio), color="k", lw=1.3, label="energy produced from uptake")
    ax.set_ylabel("Energy per unit biomass C per day")
    ax.set_title("Energy produced and spent")
    ax = axes[i][2]
    toE, enz, ov, kept = shares(m)
    ax.stackplot(m.index, 100 * sm(toE), 100 * sm(enz), 100 * sm(ov), 100 * sm(kept.clip(lower=0)),
                 colors=["#1f77b4", "#bcbd22", "#d6172a", "#2ca02c"],
                 labels=["to energy (respired / fermented)", "into enzymes", "released by cells short of N or P", "retained as biomass"])
    ax.set_ylim(0, 100)
    ax.set_ylabel("% of C taken up")
    ax.set_title("Fate of the C taken up")
    ax = axes[i][3]
    ratio = m[f"E_uptake_{G}"] / m[f"E_cost_{G}"].replace(0, np.nan)
    ax.plot(m.index, sm(ratio), color="k", label="energy produced / energy spent")
    ax.axhline(1, color="gray", ls=":", lw=1)
    ax.set_ylabel("Energy produced / spent")
    ax3 = ax.twinx()
    ax3.plot(m.index, 100 * sm(ov), color="#d6172a", label="C released by cells short of N or P, % of uptake")
    ax3.set_ylim(0, None)
    ax3.set_ylabel("% of C taken up", color="#d6172a")
    nclip = int(d[[f"cells_Eclip_{G}"]].sum().iloc[0])
    ax.set_title(f"Energy balance and the C released ({nclip} cell-days with an energy shortfall)")
    for a in axes[i]:
        a.spines[["top"]].set_visible(False)
    if i == len(runs) - 1:
        for a in axes[i]:
            a.set_xlabel("Day")
axes[0][0].legend(frameon=False, loc="lower right", fontsize=7.5)
axes[0][1].legend(frameon=False, fontsize=7, loc="upper left")
axes[0][2].legend(frameon=False, fontsize=7, loc="lower left")
axes[0][3].legend(frameon=False, fontsize=7.5, loc="upper left")
_note = f", {sys.argv[3]}" if len(sys.argv) > 3 else ""
fig.suptitle(f"Minimal model ({NAME}), {LABEL}{_note}: energy balance and C allocation with a C:N 20 and a C:N 40 substrate (3 seeds)", fontsize=11, fontweight="bold")
fig.tight_layout(rect=(0, 0, 1, 0.96))
fig.savefig(os.path.join(root, "minimal_e.png"), dpi=150, bbox_inches="tight")

pd.set_option("display.width", 230)
for lab, d in runs.items():
    m = d.groupby("doy").mean(numeric_only=True)
    print(f"\n=== {lab} ({NAME}); max |daily energy residual| {d.resid_E.abs().max():.4f}; cell-days with energy shortfall {int(d[f'cells_Eclip_{G}'].sum())}; energy lost to the E cap {d[f'E_cap_{G}'].sum():.4f}")
    rows = {}
    for lo, hi in ((1, 30), (31, 100), (101, 200), (201, 365)):
        x = m.loc[lo:hi]
        up_c = x[f"uptake_C_{G}"].sum()
        rows[f"days {lo}-{hi}"] = {
            "E:C (mean)": x[f"EC_{G}"].mean(), "E:C as % of E_max": 100 * x[f"EC_{G}"].mean() / EMAX, "C routed to energy, % of uptake": 100 * x[f"toE_C_{G}"].sum() / up_c,
            "C into enzymes, %": 100 * x[f"enzyme_C_{G}"].sum() / up_c, "C released (N/P short), %": 100 * x[f"overflow_C_{G}"].sum() / up_c,
            "C retained as biomass, %": 100 * (x[f"uptake_BC_{G}"].sum() - x[f"enzyme_C_{G}"].sum() - x[f"overflow_C_{G}"].sum()) / up_c,
            "maintenance, % of energy spent": 100 * x[f"E_maint_{G}"].sum() / x[f"E_cost_{G}"].sum(), "transport, %": 100 * x[f"E_transport_{G}"].sum() / x[f"E_cost_{G}"].sum(),
            "HSP, %": 100 * x[f"E_hsp_{G}"].sum() / x[f"E_cost_{G}"].sum(), "enzymes+osmolytes, %": 100 * (x[f"E_enzyme_{G}"].sum() + x[f"E_osmo_{G}"].sum()) / x[f"E_cost_{G}"].sum(),
            "energy produced / spent": x[f"E_uptake_{G}"].sum() / x[f"E_cost_{G}"].sum(), "model CUE": x["CUE"].mean()}
    print(pd.DataFrame(rows).round(2).to_string())
