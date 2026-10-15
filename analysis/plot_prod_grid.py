"""
Heatmaps over a two-parameter grid from sweep_hsp_effect.py seed replicates (seed*.csv):
cold penalty (mean change in biomass vs no-HSP control at the cold temperatures), warm-side changes,
and CUE, pooled as seed means. Styled like plot_temperature.py heatmaps (YlOrBr, annotated cells).

Usage:
    python plot_prod_grid.py ../analysis_output/hsp_prod_e1500 --rows HSP_Induci_Prod --cols HSP_Consti_Prod \\
        --baseline 0.3 0.0003
Writes <dir>/prod_grid.png and <dir>/prod_grid_summary.csv.
"""
import argparse
import glob
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

plt.rcParams.update({
    "font.size": 10, "axes.titlesize": 11, "axes.labelsize": 10,
    "legend.fontsize": 9, "figure.dpi": 150,
})


def vfmt(v):
    return f"{v:.3g}"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("sweep_dir")
    ap.add_argument("--rows", required=True)
    ap.add_argument("--cols", required=True)
    ap.add_argument("--baseline", type=float, nargs=2, default=None, help="row value, column value to mark with *")
    ap.add_argument("--cold", type=float, nargs="+", default=[4, 8])
    a = ap.parse_args()

    s = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(os.path.join(a.sweep_dir, "seed*.csv")))], ignore_index=True)
    n_seeds = len(glob.glob(os.path.join(a.sweep_dir, "seed*.csv")))
    key = [a.rows, a.cols]
    m = s.groupby(key + ["Temp"]).agg(rel=("rel_biomC", "mean"), CUE=("CUE", "mean"), ctrl=("ref_CUE", "mean"),
                                       EC=("EC_end", "mean"), HSP=("HSP_prop", "mean")).reset_index()

    def at(t, col):
        return m[m["Temp"] == t].pivot(index=a.rows, columns=a.cols, values=col)

    cold = sum(at(t, "rel") for t in a.cold) / len(a.cold)
    cue_mean = m[m["Temp"] <= 32].groupby(key)["CUE"].mean().unstack(a.cols)
    panels = [
        (cold, f"Cold penalty (mean of {', '.join(f'{t:g}' for t in a.cold)}°C)", "% change in biomass vs no HSP"),
        (at(16, "rel"), "Biomass change at 16°C", "% vs no HSP"),
        (at(28, "rel"), "Biomass change at 28°C", "% vs no HSP"),
        (at(32, "rel"), "Biomass change at 32°C", "% vs no HSP"),
        (cue_mean, "Mean CUE, 4-32°C (HSP run)", "CUE"),
        (at(4, "HSP"), "HSP proportion at 4°C", "HSP mass / total biomass"),
    ]
    rows, cols = list(cold.index), list(cold.columns)
    fig, axes = plt.subplots(2, 3, figsize=(14, 8.5))
    fig.suptitle(f"HSP production rates: cold, warm and CUE responses (mean of {n_seeds} seeds)",
                 fontsize=11, fontweight="bold")
    for ax, (df, title, cb) in zip(axes.flat, panels):
        z = df.loc[rows, cols].values
        im = ax.imshow(z, aspect="auto", cmap="YlOrBr", origin="lower")
        ax.set_xticks(range(len(cols)), [vfmt(c) for c in cols])
        ax.set_yticks(range(len(rows)), [vfmt(r) for r in rows])
        ax.set_xlabel(a.cols)
        ax.set_ylabel(a.rows)
        ax.set_title(title, fontsize=10)
        vmax = np.nanmax(z)
        for i in range(len(rows)):
            for j in range(len(cols)):
                v = z[i, j]
                txt = f"{v:.3f}" if title.startswith(("Mean CUE", "HSP")) else f"{v:.0f}"
                if a.baseline and rows[i] == a.baseline[0] and cols[j] == a.baseline[1]:
                    txt += "*"
                ax.text(j, i, txt, ha="center", va="center", fontsize=8,
                        color="white" if v > 0.6 * vmax else "black")
        fig.colorbar(im, ax=ax, label=cb, fraction=0.046, pad=0.04)
    fig.text(0.5, -0.01, "* = current setting. 20x20 grid, 365 days; control = no HSP production.",
             fontsize=8, color="dimgray", ha="center")
    fig.tight_layout()
    fig.savefig(os.path.join(a.sweep_dir, "prod_grid.png"), dpi=150, bbox_inches="tight")
    out = m.pivot_table(index=key, columns="Temp", values=["rel", "CUE"])
    out.columns = [f"{c}_{t:g}" for c, t in out.columns]
    out.to_csv(os.path.join(a.sweep_dir, "prod_grid_summary.csv"))
    print("wrote prod_grid.png and prod_grid_summary.csv")


if __name__ == "__main__":
    main()
