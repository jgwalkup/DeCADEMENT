"""Two-parameter grid sweep of the minimal model (scenario climates, MMRT and Arrhenius), with heat-map figures.

    python analysis/grid2d_sweep.py run OUT_DIR --x NAME=v1,v2,... --y NAME=v1,v2,... [--procs 4] [--input DIR]
        [--climates CSV ...] [--kinetics 1 0] [--seeds 0 1] [--side 20] [--days 365]
    python analysis/grid2d_sweep.py plot OUT_DIR

run: every x value x y value x climate x kinetics (runtime Kinetics_Model: 1 MMRT, 0 Arrhenius) x seed, one 365-day pulse,
    side x side grid, R_endo 1, input folder otherwise unchanged. NAME may join several parameters.csv rows with "+"
    (all set to the same value) or name a monomer_energy.csv cell as monomer_energy:<row>:<column>. Appends to
    OUT_DIR/runs.csv and skips runs already there (restartable). Outputs per run as in morris_minimal.py.
plot: seed means; OUT_DIR/grid2d.png (heat maps of % degraded, day 50 %, CUE, endogenous-respiration share and substrate
    left, one column per climate x kinetics) and OUT_DIR/grid2d_warming.png (+5 C and +10 C minus ambient); tables in
    OUT_DIR/grid2d_means.csv. The input's baseline values are outlined.
"""
import argparse
import itertools
import os
import sys
from multiprocessing import Pool

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from morris_minimal import DEFAULT_CLIMATES, KIN, metrics, _name  # noqa: E402

KEYS = ["x", "y", "climate", "kinetics", "seed"]
CLIM = [("ambient_temp", "Ambient"), ("moderate_warm", "+5 °C"), ("extreme_warm", "+10 °C")]
ROWS = [("degraded_365", "% degraded, day 365", "Blues"), ("day_50pct", "Day 50 % degraded\n(366: not reached)", "Blues_r"),
        ("CUE_annual", "Annual CUE", "Blues"), ("endo_share_pct", "Endogenous respiration\n(% of respiration)", "Oranges"),
        ("substrate_end", "Substrate left, day 365 (%)", "Oranges")]


def _parse(spec):
    name, vals = spec.split("=")
    return name, [float(v) for v in vals.split(",")]


def _job(args):
    from sweep_minimal import run_one
    xn, xv, yn, yv, src, clim, kin, seed, side, days = args
    changes = {"R_endo": 1.0}
    for n, v in ((xn, xv), (yn, yv)):
        for row in n.split("+"):
            changes[row] = v
    _, series = run_one(src, changes, side, days, seed, climate=clim, runtime={"Kinetics_Model": kin})
    return dict(x_name=xn, x=xv, y_name=yn, y=yv, climate=_name(clim), kinetics=KIN[kin], seed=seed,
                **metrics(series, side))


def run(a):
    (xn, xs), (yn, ys) = _parse(a.x), _parse(a.y)
    src = os.path.abspath(a.input)
    os.makedirs(a.out_dir, exist_ok=True)
    out = os.path.join(a.out_dir, "runs.csv")
    done = set()
    if os.path.exists(out):
        prev = pd.read_csv(out)
        done = set(zip(prev.x.round(10), prev.y.round(10), prev.climate, prev.kinetics, prev.seed))
    jobs = [(xn, xv, yn, yv, src, os.path.abspath(c), k, s, a.side, a.days)
            for xv, yv, c, k, s in itertools.product(xs, ys, a.climates, a.kinetics, a.seeds)
            if (round(xv, 10), round(yv, 10), _name(c), KIN[k], s) not in done]
    print(f"{len(jobs)} runs to do -> {out}", flush=True)
    with Pool(a.procs) as pool:
        for n, r in enumerate(pool.imap_unordered(_job, jobs), 1):
            pd.DataFrame([r]).to_csv(out, mode="a", header=not os.path.exists(out), index=False)
            if n % 10 == 0 or n == len(jobs):
                print(f"  {n}/{len(jobs)}", flush=True)


def _baseline(src, name):
    row = name.split("+")[0]
    if row.startswith("monomer_energy:"):
        _, r, c = row.split(":")
        return float(pd.read_csv(os.path.join(src, "monomer_energy.csv"), index_col=0).loc[r, c])
    p = pd.read_csv(os.path.join(src, "parameters.csv"), header=None, index_col=0)[1]
    return float(p[row])


def _heat(ax, m, cmap, fmt, vmin=None, vmax=None, base=None, center=None):
    import matplotlib.colors as mcolors
    norm = mcolors.TwoSlopeNorm(vcenter=0, vmin=min(-1e-9, vmin), vmax=max(1e-9, vmax)) if center is not None else None
    im = ax.imshow(m.to_numpy(), origin="lower", aspect="auto", cmap=cmap, norm=norm,
                   vmin=None if norm else vmin, vmax=None if norm else vmax)
    for (i, j), v in np.ndenumerate(m.to_numpy()):
        if np.isfinite(v):
            ax.text(j, i, fmt.format(v), ha="center", va="center", fontsize=6.5,
                    color="white" if (abs(im.norm(v) - 0.5) > 0.3 if center is not None else im.norm(v) > 0.6) else "black")
    ax.set_xticks(range(m.shape[1]), [f"{c:g}" for c in m.columns], fontsize=7)
    ax.set_yticks(range(m.shape[0]), [f"{r:g}" for r in m.index], fontsize=7)
    if base is not None and base[0] in m.columns and base[1] in m.index:
        j, i = list(m.columns).index(base[0]), list(m.index).index(base[1])
        ax.add_patch(__import__("matplotlib.patches", fromlist=["Rectangle"]).Rectangle(
            (j - 0.5, i - 0.5), 1, 1, fill=False, ec="#d62728", lw=1.6))
    return im


def plot(a):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    r = pd.read_csv(os.path.join(a.out_dir, "runs.csv"))
    xn, yn = r.x_name.iloc[0], r.y_name.iloc[0]
    src = os.path.join(HERE, "..", "input_minimal_e1500_crich_km25_in")
    try:
        base = (_baseline(src, xn), _baseline(src, yn))
    except Exception:
        base = None
    outs = [o for o, _, _ in ROWS]
    m = r.groupby(["climate", "kinetics", "x", "y"])[outs].mean().reset_index()
    m.to_csv(os.path.join(a.out_dir, "grid2d_means.csv"), index=False)
    cols = [(c, cn, k) for c, cn in CLIM for k in ("MMRT", "Arrhenius") if ((m.climate == c) & (m.kinetics == k)).any()]

    fig, axs = plt.subplots(len(ROWS), len(cols), figsize=(3.1 * len(cols) + 1, 2.7 * len(ROWS)), squeeze=False)
    for i, (o, lab, cmap) in enumerate(ROWS):
        vmin, vmax = m[o].min(), m[o].max()
        for j, (c, cn, k) in enumerate(cols):
            g = m[(m.climate == c) & (m.kinetics == k)].pivot(index="y", columns="x", values=o)
            fmt = "{:.2f}" if o == "CUE_annual" else "{:.0f}"
            im = _heat(axs[i, j], g, cmap, fmt, vmin, vmax, base)
            if i == 0:
                axs[i, j].set_title(f"{cn}, {k}", fontsize=9)
            if j == 0:
                axs[i, j].set_ylabel(f"{lab}\n\n{yn.split('+')[0]}", fontsize=8)
            if i == len(ROWS) - 1:
                axs[i, j].set_xlabel(xn.split("+")[0], fontsize=8)
        fig.colorbar(im, ax=axs[i, :].tolist(), fraction=0.015, pad=0.01)
    fig.suptitle(f"{xn.split('+')[0]} x {yn.split('+')[0]} (bacterial minimal model, one pulse, R_endo on; mean of seeds; "
                 "red outline: current baseline)", fontsize=10)
    fig.savefig(os.path.join(a.out_dir, "grid2d.png"), dpi=120, bbox_inches="tight")

    # warming responses
    w = m.pivot_table(index=["kinetics", "x", "y"], columns="climate", values=outs)
    wr = [("degraded_365", "Change in % degraded"), ("day_50pct", "Change in day 50 %"),
          ("CUE_annual", "Change in CUE"), ("endo_share_pct", "Change in endogenous-respiration share")]
    wcols = [(warm, wn, k) for warm, wn in (("moderate_warm", "+5 °C"), ("extreme_warm", "+10 °C"))
             for k in ("MMRT", "Arrhenius") if (o := ("degraded_365", warm)) in w.columns]
    if wcols:
        fig, axs = plt.subplots(len(wr), len(wcols), figsize=(3.1 * len(wcols) + 1, 2.7 * len(wr)), squeeze=False)
        for i, (o, lab) in enumerate(wr):
            d = {(warm, k): (w[(o, warm)] - w[(o, "ambient_temp")]).xs(k, level="kinetics") for warm, _, k in wcols}
            lim = max(abs(pd.concat(d.values()).min()), abs(pd.concat(d.values()).max()))
            for j, (warm, wn, k) in enumerate(wcols):
                g = d[(warm, k)].unstack("x")
                fmt = "{:+.3f}" if o == "CUE_annual" else "{:+.0f}"
                im = _heat(axs[i, j], g, "RdBu_r" if o != "degraded_365" else "RdBu", fmt, -lim, lim, base, center=0)
                if i == 0:
                    axs[i, j].set_title(f"{wn} minus ambient, {k}", fontsize=9)
                if j == 0:
                    axs[i, j].set_ylabel(f"{lab}\n\n{yn.split('+')[0]}", fontsize=8)
                if i == len(wr) - 1:
                    axs[i, j].set_xlabel(xn.split("+")[0], fontsize=8)
            fig.colorbar(im, ax=axs[i, :].tolist(), fraction=0.015, pad=0.01)
        fig.suptitle(f"Warming responses: {xn.split('+')[0]} x {yn.split('+')[0]} (red outline: current baseline)",
                     fontsize=10)
        fig.savefig(os.path.join(a.out_dir, "grid2d_warming.png"), dpi=120, bbox_inches="tight")
    print(f"{len(r)} runs -> grid2d.png, grid2d_warming.png, grid2d_means.csv")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run")
    r.add_argument("out_dir")
    r.add_argument("--x", required=True, metavar="NAME=v1,v2,...")
    r.add_argument("--y", required=True, metavar="NAME=v1,v2,...")
    r.add_argument("--procs", type=int, default=4)
    r.add_argument("--input", default=os.path.join(HERE, "..", "input_minimal_e1500_crich_km25_in"))
    r.add_argument("--climates", nargs="+", default=DEFAULT_CLIMATES)
    r.add_argument("--kinetics", type=int, nargs="+", default=[1, 0])
    r.add_argument("--seeds", type=int, nargs="+", default=[0, 1])
    r.add_argument("--side", type=int, default=20)
    r.add_argument("--days", type=int, default=365)
    p = sub.add_parser("plot")
    p.add_argument("out_dir")
    a = ap.parse_args()
    {"run": run, "plot": plot}[a.cmd](a)


if __name__ == "__main__":
    main()
