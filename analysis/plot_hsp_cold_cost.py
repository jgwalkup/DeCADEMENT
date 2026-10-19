"""HSP cold cost: HSP vs control (HSP production off) across a constant-temperature gradient, for several HSP settings.

    python analysis/plot_hsp_cold_cost.py OUT_PNG LABEL=HSP_DIR:CONTROL_DIR [LABEL=HSP_DIR:CONTROL_DIR ...]

Each directory is a kinetics_gradient.py run (MMRT, flat wet climate, one pulse). Panels: mean biomass C over the year,
HSP vs its control (%); day 50 % degraded for HSP runs (solid) and controls (dashed). Mean of seeds. Writes OUT_PNG and
OUT_PNG with .csv.
"""
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402


def summ(d):
    t = pd.read_csv(f"{d}/timeseries.csv.gz")
    out = []
    for (T, s), g in t.groupby(["Temp_set", "seed"]):
        g = g.sort_values("doy")
        dg = 100 - g.substrate_remaining_pct.to_numpy(); dd = g.doy.to_numpy()
        out.append(dict(T=T, biomass_mean=g.biomass_C.mean(), degraded_365=dg[-1],
                        day_50pct=dd[np.argmax(dg >= 50)] if (dg >= 50).any() else np.nan))
    return pd.DataFrame(out).groupby("T").mean()


def main(out, specs):
    rows, cols = [], plt.rcParams["axes.prop_cycle"].by_key()["color"]
    fig, ax = plt.subplots(1, 2, figsize=(13, 5))
    for i, spec in enumerate(specs):
        lab, dirs = spec.split("=", 1)
        h, c = (summ(x) for x in dirs.split(":"))
        T = h.index.intersection(c.index)
        diff = 100 * (h.loc[T, "biomass_mean"] / c.loc[T, "biomass_mean"] - 1)
        rows.append(pd.DataFrame({"setting": lab, "T": T, "biomass_vs_control_pct": diff.values,
                                  "day_50pct_hsp": h.loc[T, "day_50pct"].values,
                                  "day_50pct_control": c.loc[T, "day_50pct"].values,
                                  "degraded_365_hsp": h.loc[T, "degraded_365"].values,
                                  "degraded_365_control": c.loc[T, "degraded_365"].values}))
        ax[0].plot(T, diff, marker="o", ms=4, lw=1.7, color=cols[i % len(cols)], label=lab)
        ax[1].plot(T, h.loc[T, "day_50pct"], marker="o", ms=4, lw=1.7, color=cols[i % len(cols)], label=lab)
        ax[1].plot(T, c.loc[T, "day_50pct"], ls="--", lw=1.1, color=cols[i % len(cols)])
    pd.concat(rows).round(2).to_csv(out.rsplit(".", 1)[0] + ".csv", index=False)
    ax[0].axhline(0, color="k", lw=0.8)
    ax[0].set_ylim(-35, 60)
    ax[0].axvspan(5.3, 30.7, color="#3b6fb6", alpha=0.07, lw=0)
    ax[0].text(5.6, -32, "ambient range", fontsize=8, color="#3b6fb6")
    ax[0].set_ylabel("Mean biomass C, HSP vs control (%)\n(clipped at +60 %; controls collapse at 30 C)")
    ax[1].set_ylabel("Day 50 % degraded (solid HSP, dashed control)")
    for a in ax:
        a.set_xlabel("Constant temperature (C)")
        a.grid(alpha=0.3, lw=0.6)
        a.spines[["top", "right"]].set_visible(False)
    ax[0].legend(fontsize=8, frameon=False, loc="upper left")
    fig.suptitle("HSP cold cost across a constant-temperature gradient (bacterial minimal model, final set, MMRT, wet Psi, "
                 "one pulse, R_endo on; mean of 2 seeds)", fontsize=10.5)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(out, dpi=130)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2:])
