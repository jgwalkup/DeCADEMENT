"""
N pools and fluxes of the minimal (single-taxon, single-substrate) model for two substrate stoichiometries, from n_tracker.py output.
Usage: python plot_minimal_n.py ../analysis_output/minimal_n [bac|fun]   (expects ref/n_daily.csv and crich/n_daily.csv; group of the single taxon)
Writes minimal_n.png and prints a table.
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
_LAB = {"ref": "IdealSub C 100, N 5, P 0.5 (C:N 20)", "crich": "IdealSub C 200, N 5, P 0.5 (C:N 40)"}
runs = {_LAB[k]: pd.read_csv(os.path.join(root, k, "n_daily.csv")) for k in ("ref", "crich") if os.path.exists(os.path.join(root, k, "n_daily.csv"))}
TARGET = 0.16 / 0.825 if G == "bac" else 0.09 / 0.9
MINQ = 0.1633 if G == "bac" else 0.0617
POOLS = [("sub_N_litter", "substrate", "#8c6d31"), ("sub_N_dead", "dead microbes and enzymes", "#bcbd22"), ("mon_N_org", "organic monomer", "#17becf"),
         ("mon_N_nh4", "NH4", "#9edae5"), ("enz_N", "enzymes", "#7f7f7f"), (f"mic_N_{G}", NAME, "#1f77b4" if G == "bac" else "#d6172a")]


def sm(x, w=7):
    return x.rolling(w, center=True, min_periods=1).mean()


fig, axes = make_axes((19, 8), len(runs))
for i, (lab, d) in enumerate(runs.items()):
    m = d.groupby("doy").mean(numeric_only=True)
    ax = axes[i][0]
    ax.stackplot(m.index, *[m[c] for c, _, _ in POOLS], colors=[c for _, _, c in POOLS], labels=[l for _, l, _ in POOLS])
    ax.set_title(f"{lab}\nN pools (per cell)")
    ax = axes[i][1]
    nc = m[f"mic_N_{G}"] / m[f"mic_C_{G}"]
    ax.plot(m.index, nc, color="#1f77b4", label="N:C of the biomass")
    ax.axhline(TARGET, color="#1f77b4", ls=":", lw=1)
    ax.axhline(MINQ, color="gray", ls="--", lw=0.8)
    ax.set_ylim(0.04, 0.22)
    ax.set_ylabel("N : C (dotted target, dashed minimum)")
    ax2 = ax.twinx()
    ax2.plot(m.index, m[f"mic_C_{G}"], color="k", lw=1.2, label="biomass C")
    ax2.set_ylabel("Biomass C per cell")
    ax.set_title("Biomass N:C and biomass C")
    ax = axes[i][2]
    ax.plot(m.index, sm(m[f"uptake_N_{G}"]), color="#1f77b4", label="N uptake")
    ax.plot(m.index, sm(m["decay_N"]), color="k", ls=":", label="N released by decay")
    ax.plot(m.index, sm(m[f"release_N_{G}"]), color="#2ca02c", label="N released to NH4 by living cells")
    ax.plot(m.index, sm(m[f"death_N_{G}"]), color="#8c564b", label="N in cells killed")
    if "input_N" in m and m["input_N"].max() > 0:
        ax.plot(m.index, m["input_N"], color="#ff7f0e", ls="--", label="N input (Input_NH4)")
        ax.plot(m.index, sm(m["leach_N"]), color="#ff7f0e", label="N leached")
    ax.set_title("N fluxes per day")
    ax = axes[i][3]
    share = sm(m[f"overflow_C_{G}"]) / sm(m[f"uptake_C_{G}"]).replace(0, np.nan) * 100
    ax.plot(m.index, share, color="#d6172a", label="C released by cells short of N, % of C uptake")
    ax.set_ylabel("% of C uptake")
    ax.set_ylim(0, None)
    ax3 = ax.twinx()
    ax3.plot(m.index, sm(m["CUE"]), color="#7f7f7f", label="CUE")
    ax3.set_ylim(0, 1)
    ax3.set_ylabel("CUE")
    ax.set_title("C lost when N is short, and CUE")
    for a in axes[i]:
        a.spines[["top"]].set_visible(False)
    if i == len(runs) - 1:
        for a in axes[i]:
            a.set_xlabel("Day")
axes[0][0].legend(frameon=False, fontsize=7.5)
axes[0][2].legend(frameon=False, fontsize=7.5)
axes[0][3].legend(frameon=False, fontsize=7.5, loc="upper left")
_in = next(iter(runs.values()))
_note = ""
if "input_N" in _in and _in["input_N"].max() > 0:
    _note = f", Input_NH4 {_in['input_N'].max():g} and Input_PO4 {_in['input_P'].max():g}"
fig.suptitle(f"Minimal model ({NAME}), {LABEL}{_note}: N pools and fluxes with a C:N 20 and a C:N 40 substrate (3 seeds)", fontsize=11, fontweight="bold")
fig.tight_layout(rect=(0, 0, 1, 0.96))
fig.savefig(os.path.join(root, "minimal_n.png"), dpi=150, bbox_inches="tight")

pd.set_option("display.width", 220)
for lab, d in runs.items():
    m = d.groupby("doy").mean(numeric_only=True)
    print(f"\n=== {lab}  (max |daily N residual| {d.resid_N.abs().max():.1e}, N leached over the year {m.leach_N.sum():.4f})")
    rows = {}
    for day in (30, 100, 200, 300, 365):
        r = m.loc[day]
        rows[f"day {day}"] = dict(biomass_C=r[f"mic_C_{G}"], N_to_C=r[f"mic_N_{G}"] / r[f"mic_C_{G}"], substrate_C_left_pct=100 * r.sub_C_litter / m.sub_C_litter.iloc[0],
                                  substrate_N=r.sub_N_litter, dead_N=r.sub_N_dead, monomer_N=r.mon_N_org + r.mon_N_nh4, bacteria_N=r[f"mic_N_{G}"],
                                  energy_to_C=r.energy_to_C, CUE=r.CUE)
    print(pd.DataFrame(rows).round(3).to_string())
    up_c, ov_c, up_n = m[f"uptake_C_{G}"].sum(), m[f"overflow_C_{G}"].sum(), m[f"uptake_N_{G}"].sum()
    onset = m.index[(m[f"overflow_C_{G}"] / m[f"uptake_C_{G}"].replace(0, np.nan) > 0.02)].min()
    below = m.index[(m[f"mic_N_{G}"] / m[f"mic_C_{G}"] < 0.97 * TARGET)].min()
    print(f"year totals per cell: C uptake {up_c:.2f}, C released as overflow {ov_c:.2f} ({100*ov_c/up_c:.1f} % of uptake), N uptake {up_n:.3f}, N in cells killed {m[f"death_N_{G}"].sum():.3f}, "
          f"N released by living cells {m[f"release_N_{G}"].sum():.3f}; overflow > 2 % of uptake from day {onset}; biomass N:C below 97 % of target from day {below}")
