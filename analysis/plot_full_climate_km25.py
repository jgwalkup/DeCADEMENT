"""
Real-climate figure for full_climate_run.py output of input_e1500_enz2_km25 (Osmo_tol_lower 0, Osmo_tol_ref 0), next to the earlier
run at the old slopes (7/5) with Osmo_tol_ref 0.055 (tolerance 1 for all taxa; analysis_output/full_climate_osmo).
Usage: python plot_full_climate_km25.py ../analysis_output/full_climate_km25_osmo0 ../analysis_output/full_climate_osmo
Writes <first dir>/full_climate_km25_osmo0.png. Columns: slopes 7/5 (previous) vs 2.5/2.5; one line per pulse (mean over seeds).
"""
import os
import sys

import matplotlib.pyplot as plt
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from scenario_utils import scenario_color_map  # noqa: E402

plt.rcParams.update({"font.size": 10, "axes.titlesize": 11, "axes.labelsize": 10, "legend.fontsize": 9, "figure.dpi": 150})
new_dir, old_dir = sys.argv[1], sys.argv[2]
new = pd.read_csv(os.path.join(new_dir, "timeseries.csv"))
old = pd.read_csv(os.path.join(old_dir, "timeseries.csv"))
old = old[(old["input"] == "input_e1500_enz2_osmo") & (old["seed"] <= 2)]  # same three communities as seeds 0-2 of the new run
sets = [("Previous: slopes 7 / 5, Osmo_tol_ref 0.055 (tolerance 1), 3 communities", old),
        (f"Slopes 2.5 / 2.5, Osmo_tol_lower 0 / Osmo_tol_ref 0 (tolerance 1), {new['seed'].nunique()} communities", new)]
pulses = sorted(new["pulse"].unique())
col = scenario_color_map([str(p) for p in pulses], palette="discrete")
rows = [("alive_rows", "Cell-taxa alive"), ("taxa_alive", "Taxa present"), ("biomass_C", "Biomass C per grid cell"),
        ("substrate_remaining_pct", "Substrate C remaining (% of initial)"), ("energy_to_C", "Energy : C")]
fig = plt.figure(figsize=(12, 15))
gs = fig.add_gridspec(len(rows) + 1, 2, height_ratios=[0.8] + [1] * len(rows))
c = new[new["pulse"] == 1].groupby("doy")[["Temp", "Psi"]].mean().reset_index()
ax = fig.add_subplot(gs[0, :])
ax.plot(c["doy"], c["Temp"], color="#d6172a", lw=1.2)
ax.set_ylabel("Temperature (°C)")
ax2 = ax.twinx()
ax2.plot(c["doy"], c["Psi"], color="#2166ac", lw=1.0, alpha=0.8)
ax2.set_ylabel("Ψ (MPa)")
ax.set_title("Real climate (one year, repeated each pulse; day 1 = 15 Dec)")
for a in (ax, ax2):
    a.spines[["top"]].set_visible(False)
for j, (title, d) in enumerate(sets):
    m = d.groupby(["pulse", "doy"], as_index=False)[[r[0] for r in rows]].mean()
    for i, (colname, lab) in enumerate(rows):
        a = fig.add_subplot(gs[i + 1, j])
        for p in pulses:
            q = m[m["pulse"] == p].sort_values("doy")
            y = q[colname].rolling(7, center=True, min_periods=1).mean() if colname == "energy_to_C" else q[colname]
            a.plot(q["doy"], y, color=col[str(p)], lw=1.2, label=f"pulse {p}")
        if i == 0:
            a.set_title(title, fontsize=9)
        if j == 0:
            a.set_ylabel(lab)
        if i == len(rows) - 1:
            a.set_xlabel("Day of pulse")
        a.set_ylim(bottom=0)
        a.spines[["top", "right"]].set_visible(False)
h, l = fig.axes[2].get_legend_handles_labels()
fig.legend(h, l, loc="lower center", ncols=len(pulses), framealpha=0.85, bbox_to_anchor=(0.5, -0.01))
fig.suptitle("Full model, real climate, drought tolerance 1 for all taxa: old vs new Vmax-Km slopes (20x20)", fontsize=11,
             fontweight="bold", y=1.0)
fig.tight_layout(rect=(0, 0.02, 1, 0.985))
out = os.path.join(new_dir, "full_climate_km25_osmo0.png")
fig.savefig(out, dpi=150, bbox_inches="tight")
print("wrote", out)
