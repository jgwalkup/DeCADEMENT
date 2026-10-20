"""Substrate-specific decomposition and uptake in full-model runs (full_single_taxon.py output).

    python analysis/plot_full_substrates.py OUT_PNG TITLE LABEL=RUN_DIR:RUN_LABEL [...]
        [--climate ambient_temp] [--kinetics MMRT]

One row per setting (LABEL is the row title; RUN_LABEL the label in RUN_DIR/timeseries.csv.gz). Columns: each litter
substrate's C remaining (% of its initial C); litter C remaining per grid cell, stacked by substrate; C taken up per grid
cell per day by monomer (7-day mean, stacked; necromass monomers included). Mean of seeds. Writes OUT_PNG and a CSV of the
plotted means next to it.
"""
import argparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

COLORS = {"Cellulose": "#1c5cab", "Hemicellulose": "#6fa8e8", "Starch": "#b7d4f5", "Lignin": "#5c4033",
          "Chitin": "#e08a1e", "Protein1": "#2e8b57", "Protein2": "#5fb37f", "Protein3": "#9ed4b0",
          "OrgP1": "#9b59b6", "OrgP2": "#c9a3db", "DeadMic": "#8a8a85", "DeadEnz": "#c4c4bf"}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("out")
    ap.add_argument("title")
    ap.add_argument("specs", nargs="+")
    ap.add_argument("--climate", default="ambient_temp")
    ap.add_argument("--kinetics", default="MMRT")
    a = ap.parse_args()
    rows, cache, out = [], {}, []
    for spec in a.specs:
        lab, rest = spec.split("=", 1)
        d, rl = rest.rsplit(":", 1)
        if d not in cache:
            cache[d] = pd.read_csv(f"{d}/timeseries.csv.gz")
        t = cache[d]
        g = t[(t.label == rl) & (t.climate == a.climate) & (t.kinetics == a.kinetics)].groupby("doy").mean(
            numeric_only=True)
        rows.append((lab, g))
        out.append(g.assign(setting=lab))
    pd.concat(out).round(5).to_csv(a.out.rsplit(".", 1)[0] + ".csv")
    subs = [c[2:] for c in rows[0][1].columns if c.startswith("S_") and c[2:] not in ("DeadMic", "DeadEnz")]
    mons = [c for c in rows[0][1].columns if c.startswith("U_")]
    fig, ax = plt.subplots(len(rows), 3, figsize=(17, 3.6 * len(rows) + 0.8), squeeze=False, sharex=True)
    for i, (lab, g) in enumerate(rows):
        x = g.index
        for s in subs:
            ax[i, 0].plot(x, 100 * g[f"S_{s}"] / g[f"S_{s}"].iloc[0], color=COLORS.get(s), lw=1.7, label=s)
        ax[i, 1].stackplot(x, *[g[f"S_{s}"] for s in subs], colors=[COLORS.get(s) for s in subs], labels=subs)
        up = [g[m].rolling(7, center=True, min_periods=1).mean() for m in mons]
        names = [m.replace("U_Mon_", "") for m in mons]
        ax[i, 2].stackplot(x, *up, colors=[COLORS.get(n) for n in names], labels=names)
        ax[i, 0].set_ylabel(f"{lab}\n\nSubstrate C remaining (%)")
        ax[i, 1].set_ylabel("Litter C per grid cell")
        ax[i, 2].set_ylabel("C uptake per grid cell per day")
        for a_ in ax[i]:
            a_.grid(alpha=0.3, lw=0.6)
            a_.spines[["top", "right"]].set_visible(False)
            a_.set_xlim(0, 365)
    for a_ in ax[-1]:
        a_.set_xlabel("Day of year")
    ax[0, 0].legend(fontsize=7.5, frameon=False, ncol=2, loc="lower left")
    ax[0, 2].legend(fontsize=7.5, frameon=False, ncol=2, loc="upper right")
    fig.suptitle(a.title, fontsize=10.5)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(a.out, dpi=120)


if __name__ == "__main__":
    main()
