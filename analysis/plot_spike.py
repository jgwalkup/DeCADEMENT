"""
Figure for diagnose_spike.py output. Usage: python plot_spike.py ../analysis_output/spike_diagnosis/spike.csv
Writes spike_diagnosis.png next to the csv.
"""
import os
import sys

import matplotlib.pyplot as plt
import pandas as pd

plt.rcParams.update({"font.size": 9.5, "axes.titlesize": 10, "axes.labelsize": 9.5, "legend.fontsize": 8.5, "figure.dpi": 150})
d = pd.read_csv(sys.argv[1])
m = d.groupby(["kind", "doy"]).mean(numeric_only=True)
r = m.loc["real"].loc[200:340]
fig, ax = plt.subplots(2, 2, figsize=(13, 8.2), sharex=True)
a = ax[0][0]
a.plot(r.index, r["Psi"], color="#2166ac", label="Psi (MPa)")
a.set_ylabel("Psi (MPa)", color="#2166ac")
a2 = a.twinx()
a2.plot(r.index, r["f_uptake"], color="k", label="moisture factor on uptake")
a2.plot(r.index, r["f_enzyme"], color="#8c564b", ls="--", label="moisture factor on enzyme Vmax")
a2.plot(r.index, r["Temp"] / 30, color="#d6172a", lw=0.8, alpha=0.7, label="temperature / 30 C")
a2.set_ylim(0, 1.05)
a2.legend(frameon=False, loc="lower right", fontsize=7.5)
a.set_title("Rewetting: Psi rises from -5 to -0.5 MPa on days 294-296")
a = ax[0][1]
a.plot(r.index, r["mon_pool_before"], color="#17becf", label="monomer C pool seen by uptake")
a.plot(r.index, 10 * r["decay_C"], color="#8c6d31", label="C released by decay per day (x10)")
a.set_ylabel("C per cell")
a.legend(frameon=False)
a.set_title("Monomers pile up during the dry season and are drawn down after rewetting")
a = ax[1][0]
a.plot(r.index, r["uptake_C_per_biomass"], color="k", label="C taken up per unit biomass")
a.set_ylabel("C uptake per biomass per day")
a3 = a.twinx()
a3.plot(r.index, 100 * r["frac_pool_taken"], color="#1f77b4", label="share of the monomer pool taken per day (%)")
a3.set_ylabel("% of the pool taken per day", color="#1f77b4")
a.set_title("Uptake jumps about fivefold")
a = ax[1][1]
cl = {"real": ("#000000", "real climate"), "psi_flat": ("#2ca02c", "Psi -0.5 all year (real T)"), "no_rewet": ("#d6172a", "no rewetting (Psi held at its day-285 value)"),
      "t_flat_late": ("#ff7f0e", "real Psi, T 15 C from day 285")}
for k, (c, lab) in cl.items():
    x = m.loc[k].loc[200:340]
    a.plot(x.index, x["E_per_biomass"], color=c, label=lab)
a.set_ylabel("Energy produced per unit biomass per day")
a.legend(frameon=False, fontsize=7.5)
a.set_title("The spike needs the rewetting, not the warm days")
for row in ax:
    for a in row:
        a.set_xlabel("Day")
        a.axvspan(293.5, 296.5, color="gray", alpha=0.12, lw=0)
        a.spines[["top"]].set_visible(False)
fig.suptitle("Minimal bacterium (C 100 substrate), real climate: what drives the burst of energy production at day 295-300 (mean of 3 seeds)", fontsize=11, fontweight="bold")
fig.tight_layout(rect=(0, 0, 1, 0.96))
out = os.path.join(os.path.dirname(os.path.abspath(sys.argv[1])), "spike_diagnosis.png")
fig.savefig(out, dpi=150, bbox_inches="tight")
print("wrote", out)
