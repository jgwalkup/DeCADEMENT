"""Daily CUE over the year at each constant temperature, for several HSP settings and the no-HSP control.

    python analysis/plot_hsp_cue_timeseries.py OUT_PNG CONTROL_DIR LABEL=DIR [LABEL=DIR ...]

Each directory is a kinetics_gradient.py run (MMRT, flat wet climate, one pulse). One panel per temperature; lines are the
settings in the order given (light to dark) and the control (grey, dashed). CUE is the model's daily CUE, mean of seeds,
7-day centred running mean, shown only while more than 1 % of the substrate is left (after that uptake is near zero and
daily CUE is 0). A tick on the x axis marks the day 95 % of the substrate is degraded. Writes OUT_PNG and OUT_PNG with .csv
(daily values before smoothing).
"""
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

RAMP = ["#6da7ec", "#2a78d6", "#0d366b", "#06203f"]   # ordinal blue steps (light -> dark)
CONTROL = "#8a8a85"


def daily(d):
    t = pd.read_csv(f"{d}/timeseries.csv.gz")
    return t.groupby(["Temp_set", "doy"])[["CUE", "substrate_remaining_pct", "C_uptake"]].mean().reset_index()


def main(out, ctl_dir, specs):
    sets = [("Control (no HSP)", ctl_dir)] + [tuple(s.split("=", 1)) for s in specs]
    data = {lab: daily(d) for lab, d in sets}
    temps = sorted(set.intersection(*(set(v.Temp_set) for v in data.values())))
    pd.concat([v.assign(setting=k) for k, v in data.items()]).round(4).to_csv(out.rsplit(".", 1)[0] + ".csv", index=False)
    n = len(temps); nc = 3; nr = int(np.ceil(n / nc))
    fig, ax = plt.subplots(nr, nc, figsize=(14, 3.3 * nr), sharex=True, sharey=True)
    ax = np.atleast_1d(ax).ravel()
    for i, T in enumerate(temps):
        a = ax[i]
        for j, (lab, _) in enumerate(sets):
            g = data[lab][data[lab].Temp_set == T].sort_values("doy")
            cue = g.CUE.rolling(7, center=True, min_periods=1).mean().where(g.substrate_remaining_pct > 1)
            col, ls, lw = (CONTROL, "--", 1.4) if j == 0 else (RAMP[j - 1], "-", 1.8)
            a.plot(g.doy, cue, color=col, ls=ls, lw=lw, label=lab)
            deg = 100 - g.substrate_remaining_pct.to_numpy()
            if (deg >= 95).any():
                a.plot(g.doy.to_numpy()[np.argmax(deg >= 95)], 0.02, marker="|", ms=10, mew=2, color=col, clip_on=False)
        a.set_title(f"{T:g} C", fontsize=10, loc="left")
        a.set_ylim(0, 0.5)
        a.set_xlim(0, 365)
        a.grid(alpha=0.3, lw=0.6)
        a.spines[["top", "right"]].set_visible(False)
        if i % nc == 0:
            a.set_ylabel("Daily CUE (7-day mean)")
        if i >= n - nc:
            a.set_xlabel("Day of year")
    for a in ax[n:]:
        a.set_visible(False)
    h, l = ax[0].get_legend_handles_labels()
    fig.legend(h, l, loc="upper center", ncol=len(l), fontsize=9, frameon=False, bbox_to_anchor=(0.5, 0.955))
    fig.suptitle("Daily CUE by HSP energy cost across constant temperatures (bacterial minimal model, final set, MMRT, wet Psi, "
                 "one pulse, R_endo on; mean of 2 seeds)\nLines end when < 1 % of the substrate is left; ticks on the x axis "
                 "mark the day 95 % is degraded", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    fig.savefig(out, dpi=130)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3:])
