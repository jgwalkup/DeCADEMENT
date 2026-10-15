"""Figure for full_climate_run.py output: climate on top; below, one column per setting and one line per pulse (mean of seeds).
Usage: python plot_full_climate.py ../analysis_output/full_climate_osmo"""
import os
import sys

import matplotlib.pyplot as plt
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from scenario_utils import scenario_color_map  # noqa: E402

plt.rcParams.update({"font.size": 10, "axes.titlesize": 11, "axes.labelsize": 10, "legend.fontsize": 9, "figure.dpi": 150})
d = sys.argv[1]
s = pd.read_csv(os.path.join(d, "timeseries.csv"))
names = {"input_e1500_enz2": "No Osmo_tol_ref (tolerance 0)", "input_e1500_enz2_osmo": "Osmo_tol_ref 0.055 (tolerance 1)"}
pulses = sorted(s["pulse"].unique())
col = scenario_color_map([str(p) for p in pulses], palette="discrete")
m = s.groupby(["input", "pulse", "doy"], as_index=False)[["alive_rows", "taxa_alive", "biomass_C", "energy_to_C", "CUE", "Temp", "Psi"]].mean()
rows = [("alive_rows", "Cell-taxa alive"), ("taxa_alive", "Taxa present (of ~93)"), ("biomass_C", "Biomass C per grid cell"), ("energy_to_C", "Energy : C")]
fig = plt.figure(figsize=(12, 13))
gs = fig.add_gridspec(len(rows) + 1, 2, height_ratios=[0.8] + [1] * len(rows))
c = m[(m["input"] == "input_e1500_enz2") & (m["pulse"] == 1)].sort_values("doy")
ax = fig.add_subplot(gs[0, :])
ax.plot(c["doy"], c["Temp"], color="#d6172a", lw=1.2)
ax.set_ylabel("Temperature (°C)")
ax2 = ax.twinx()
ax2.plot(c["doy"], c["Psi"], color="#2166ac", lw=1.0, alpha=0.8)
ax2.set_ylabel("Ψ (MPa)")
ax.set_title("Real climate (one year, repeated each pulse; day 1 = 15 Dec)")
for a in (ax, ax2):
    a.spines[["top"]].set_visible(False)
for j, inp in enumerate(names):
    for i, (colname, lab) in enumerate(rows):
        a = fig.add_subplot(gs[i + 1, j])
        for p in pulses:
            q = m[(m["input"] == inp) & (m["pulse"] == p)].sort_values("doy")
            y = q[colname].rolling(7, center=True, min_periods=1).mean() if colname == "energy_to_C" else q[colname]
            a.plot(q["doy"], y, color=col[str(p)], lw=1.2, label=f"pulse {p}")
        if i == 0:
            a.set_title(names[inp])
        if j == 0:
            a.set_ylabel(lab)
        if i == len(rows) - 1:
            a.set_xlabel("Day of pulse")
        a.set_ylim(bottom=0)
        a.spines[["top", "right"]].set_visible(False)
h, l = fig.axes[2].get_legend_handles_labels()
fig.legend(h, l, loc="lower center", ncols=len(pulses), title="Pulse (mean of 3 communities)", framealpha=0.85, bbox_to_anchor=(0.5, -0.01))
fig.suptitle("Full model, real climate: drought tolerance off vs on (20x20, input_e1500_enz2)", fontsize=11, fontweight="bold", y=1.0)
fig.tight_layout(rect=(0, 0.02, 1, 0.985))
fig.savefig(os.path.join(d, "full_climate_osmo.png"), dpi=150, bbox_inches="tight")
print("wrote", os.path.join(d, "full_climate_osmo.png"))
