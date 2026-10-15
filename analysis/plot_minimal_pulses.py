"""
Pulse-resolved figure and table for multi-pulse minimal-model runs (n_tracker.py and e_tracker.py output, one line per pulse, mean of seeds).
Usage: python plot_minimal_pulses.py ../analysis_output/minimal_clim6_psi25 bac|fun [--climate CSV] [--title TEXT]
Expects n_<group>/n_daily.csv and e_<group>/e_daily.csv in the folder. Writes pulses_<group>.png and prints the per-pulse table.
"""
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

plt.rcParams.update({"font.size": 9.5, "axes.titlesize": 10, "axes.labelsize": 9.5, "legend.fontsize": 8.5, "figure.dpi": 150})


def pop(flag):
    if flag in sys.argv:
        i = sys.argv.index(flag)
        v = sys.argv[i + 1]
        del sys.argv[i:i + 2]
        return v
    return None


CLIM, TITLE = pop("--climate"), pop("--title") or ""
root, G = sys.argv[1], sys.argv[2]
NAME = {"bac": "bacterium", "fun": "fungus"}[G]
n = pd.read_csv(os.path.join(root, f"n_{G}", "n_daily.csv"))
e = pd.read_csv(os.path.join(root, f"e_{G}", "e_daily.csv"))
EMAX = 1400.0
pulses = sorted(n.pulse.unique())
cols = dict(zip(pulses, plt.cm.viridis(np.linspace(0.05, 0.92, len(pulses)))))
mn = n.groupby(["pulse", "doy"]).mean(numeric_only=True)
me = e.groupby(["pulse", "doy"]).mean(numeric_only=True)


def sm(x, w=5):
    return x.rolling(w, center=True, min_periods=1).mean()


nr = 3 if CLIM else 2
fig = plt.figure(figsize=(19, 8.5 + (1.4 if CLIM else 0)))
gs = fig.add_gridspec(nr, 4, height_ratios=([0.4] if CLIM else []) + [1, 1])
if CLIM:
    top = fig.add_subplot(gs[0, :])
    c = pd.read_csv(CLIM, index_col=0).iloc[:365]
    top.plot(np.arange(1, 366), c["Temp"], color="#d6172a", lw=1.1)
    top.set_ylabel("Temperature (C)")
    t2 = top.twinx()
    t2.plot(np.arange(1, 366), c["Psi"], color="#2166ac", lw=1.0)
    t2.set_ylabel("Psi (MPa)")
    top.set_xlim(0, 366)
    top.spines[["top"]].set_visible(False)
ax = [[fig.add_subplot(gs[(1 if CLIM else 0) + r, k]) for k in range(4)] for r in range(2)]
panels = [(0, 0, "Biomass C per cell", lambda p: mn.loc[p][f"mic_C_{G}"]),
          (0, 1, "E:C as % of E_max", lambda p: 100 * me.loc[p][f"EC_{G}"] / EMAX),
          (0, 2, "Energy produced / spent", lambda p: sm(me.loc[p][f"E_uptake_{G}"] / me.loc[p][f"E_cost_{G}"])),
          (0, 3, "% of C taken up routed to energy", lambda p: 100 * sm(me.loc[p][f"toE_C_{G}"] / me.loc[p][f"uptake_C_{G}"])),
          (1, 0, "Substrate C remaining (% of initial)", lambda p: 100 * mn.loc[p]["sub_C_litter"] / mn.loc[p]["sub_C_litter"].iloc[0]),
          (1, 1, "Biomass N:C", lambda p: mn.loc[p][f"mic_N_{G}"] / mn.loc[p][f"mic_C_{G}"]),
          (1, 2, "C released by cells short of N or P (% of uptake)", lambda p: 100 * sm(mn.loc[p][f"overflow_C_{G}"] / mn.loc[p][f"uptake_C_{G}"])),
          (1, 3, "CUE (model)", lambda p: sm(mn.loc[p]["CUE"]))]
for r, k, title, fn in panels:
    a = ax[r][k]
    for p in pulses:
        a.plot(mn.loc[p].index, fn(p), color=cols[p], lw=1.2, label=f"pulse {int(p)}")
    a.set_title(title, fontsize=9.5)
    a.spines[["top", "right"]].set_visible(False)
    a.set_xlim(0, 366)
    if r == 1:
        a.set_xlabel("Day of pulse")
ax[0][2].axhline(1, color="gray", ls=":")
ax[0][1].axhline(100, color="k", ls=":", lw=0.8)
ax[0][0].legend(frameon=False, fontsize=7.5, ncols=2)
fig.suptitle(f"Minimal model ({NAME}), C 200 substrate, {len(pulses)} pulses of the real climate{(', ' + TITLE) if TITLE else ''} (mean of {n.seed.nunique()} seeds)", fontsize=11, fontweight="bold")
fig.tight_layout(rect=(0, 0, 1, 0.96))
fig.savefig(os.path.join(root, f"pulses_{G}.png"), dpi=150, bbox_inches="tight")

rows = {}
for p in pulses:
    x, y = mn.loc[p], me.loc[p]
    up = x[f"uptake_C_{G}"].sum()
    ratio = y[f"E_uptake_{G}"] / y[f"E_cost_{G}"]
    short = e[(e.pulse == p)][f"cells_Eclip_{G}"].sum() / n[n.pulse == p][f"cells_{G}"].sum()
    rows[f"pulse {int(p)}"] = dict(biomass_day1=x[f"mic_C_{G}"].iloc[0], biomass_d365=x[f"mic_C_{G}"].iloc[-1], cells_d365=x[f"cells_{G}"].iloc[-1],
                                   degraded_pct=100 * (1 - x["sub_C_litter"].iloc[-1] / x["sub_C_litter"].iloc[0]), N_to_C_d365=x[f"mic_N_{G}"].iloc[-1] / x[f"mic_C_{G}"].iloc[-1],
                                   overflow_pct=100 * x[f"overflow_C_{G}"].sum() / up, pct_celldays_Eshortfall=100 * short, days_Eratio_lt1=int((ratio < 1).sum()),
                                   min_Eratio=ratio.min(), min_EC_pct_Emax=100 * y[f"EC_{G}"].min() / EMAX,
                                   toE_pct_d201_293=100 * y.loc[201:293][f"toE_C_{G}"].sum() / y.loc[201:293][f"uptake_C_{G}"].sum(), CUE_model=(x["CUE"] * x[f"uptake_C_{G}"]).sum() / up,
                                   N_leached=x["leach_N"].sum())
pd.set_option("display.width", 220)
print(pd.DataFrame(rows).round(3).to_string())
