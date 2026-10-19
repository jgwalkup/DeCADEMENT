"""Morris elementary-effects screening of the minimal model (both kinetics, scenario climates).

Three steps, so the run step can be spread over HPC array tasks:

    python analysis/morris_minimal.py design OUT_DIR [--factors analysis/morris_factors_minimal.csv] [--r 10]
        [--levels 4] [--pool 0] [--seed 1]
    python analysis/morris_minimal.py run OUT_DIR [--task-id 0 --n-tasks 1] [--procs 2] [--input DIR]
        [--climates CSV ...] [--kinetics 1 0] [--seeds 0 1] [--side 20] [--days 365]
    python analysis/morris_minimal.py analyze OUT_DIR [--boot 1000]

design: r trajectories of k+1 points each on a p-level grid in [0, 1]^k (Morris 1991; step delta = p / (2 (p - 1))),
    each changing one factor at a time; with --pool N > r, N candidate trajectories are drawn and r of them chosen to
    maximise their spread (Campolongo et al. 2007, greedy). Unit values map to the factor ranges in the factor table
    (scale log: geometric between low and high; linear: arithmetic). Writes OUT_DIR/design.csv and a copy of the table.
    A factor's "rows" may join several parameters.csv rows with "+" (all set to the same value) or name a
    monomer_energy.csv cell as monomer_energy:<row>:<column>.
run: every design point x climate x kinetics (runtime Kinetics_Model: 1 MMRT, 0 Arrhenius) x seed, one 365-day pulse,
    side x side grid, R_endo 1, input folder otherwise unchanged. Jobs are dealt to tasks round-robin (job j goes to task
    j mod n_tasks); each task writes OUT_DIR/runs/task_<id>.csv and skips jobs already in it (restartable). On SLURM, use
    --task-id $SLURM_ARRAY_TASK_ID --n-tasks <array size> --procs $SLURM_CPUS_PER_TASK.
analyze: seeds are averaged per point (all points use the same seeds), elementary effects EE = dy / delta along each
    trajectory, then per output, climate and kinetics: mu, mu* (mean |EE|), sigma (SD of EE) and a bootstrap 95 % interval
    of mu*; also for the warming responses (+5 C minus ambient, +10 C minus ambient, same point and kinetics).
    mu* and sigma are also given relative to the SD of the output over all design points (comparable across outputs).
    Writes OUT_DIR/morris_indices.csv, OUT_DIR/outputs_by_point.csv and OUT_DIR/morris_mu_star.png.

Outputs per run: % degraded at day 365; day of 50 % / 95 % degraded (366 if not reached); uptake-weighted annual CUE;
whole-year endogenous respiration C per grid cell and its share of respiration; substrate left (%) on the first day
endogenous respiration exceeds 0.1 % of biomass C (0 if it never does); substrate left at day 365; maximum and mean
biomass C; drought kills per grid cell.
"""
import argparse
import glob
import itertools
import os
import sys
from multiprocessing import Pool

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

CLIM_DIR = os.path.join(HERE, "..", "input", "climate_data")
DEFAULT_CLIMATES = [os.path.join(CLIM_DIR, f) for f in
                    ("climate_ambient_temp.csv", "climate_moderate_warm.csv", "climate_extreme_warm.csv")]
KIN = {1: "MMRT", 0: "Arrhenius"}
OUTPUTS = ["degraded_365", "day_50pct", "day_95pct", "CUE_annual", "endo_C_year", "endo_share_pct",
           "substrate_at_endo_onset", "substrate_end", "biomass_max", "biomass_mean", "drought_kills_per_cell"]
KEYS = ["point", "climate", "kinetics", "seed"]


def _name(path):
    return os.path.splitext(os.path.basename(path))[0].replace("climate_", "")


# ── design ───────────────────────────────────────────────────────────────────────────────────────────────────────────
def _trajectory(k, p, rng):
    """One Morris trajectory: (k+1, k) unit-space points and the factor changed at each step (Morris 1991)."""
    delta = p / (2 * (p - 1))
    grid = np.arange(p) / (p - 1)
    base = rng.choice(grid[grid <= 1 - delta + 1e-12], size=k)        # start so that +delta stays inside [0, 1]
    B = np.tril(np.ones((k + 1, k)), -1)                                # strictly lower triangular
    D = np.diag(rng.choice([-1, 1], size=k))
    P = np.eye(k)[rng.permutation(k)]
    J = np.ones((k + 1, k))
    X = (J * base + (delta / 2) * ((2 * B - J) @ D + J)) @ P
    X = np.clip(np.round(X * (p - 1)) / (p - 1), 0, 1)
    changed = [int(np.flatnonzero(np.abs(X[i + 1] - X[i]) > 1e-9)[0]) for i in range(k)]
    return X, changed


def _distances(trajs):
    """Pairwise trajectory distances (sum of point-to-point Euclidean distances; Campolongo et al. 2007)."""
    n = len(trajs)
    D = np.zeros((n, n))
    for a, b in itertools.combinations(range(n), 2):
        D[a, b] = D[b, a] = np.sqrt(((trajs[a][0][:, None, :] - trajs[b][0][None, :, :]) ** 2).sum(axis=2)).sum()
    return D


def to_value(u, f):
    lo, hi = float(f["low"]), float(f["high"])
    if f["scale"] == "log":
        return float(np.exp(np.log(lo) + u * (np.log(hi) - np.log(lo))))
    return float(lo + u * (hi - lo))


def design(a):
    fac = pd.read_csv(a.factors)
    k, p = len(fac), a.levels
    rng = np.random.default_rng(a.seed)
    n_cand = max(a.pool, a.r)
    cand = [_trajectory(k, p, rng) for _ in range(n_cand)]
    if n_cand > a.r:      # greedy: repeatedly drop the trajectory whose removal keeps the largest spread
        D2 = _distances(cand) ** 2
        idx = list(range(n_cand))
        while len(idx) > a.r:
            sub = D2[np.ix_(idx, idx)]
            idx.pop(int(np.argmin(sub.sum(axis=1))))   # removing j loses its row sum of squared distances
        cand = [cand[i] for i in idx]
    rows = []
    for t, (X, changed) in enumerate(cand):
        for s in range(k + 1):
            row = dict(point=t * (k + 1) + s, traj=t, step=s, changed=fac.factor[changed[s - 1]] if s else "")
            for j, f in fac.iterrows():
                row[f"u_{f.factor}"] = X[s, j]
                row[f.factor] = to_value(X[s, j], f)
            rows.append(row)
    os.makedirs(a.out_dir, exist_ok=True)
    pd.DataFrame(rows).to_csv(os.path.join(a.out_dir, "design.csv"), index=False)
    fac.to_csv(os.path.join(a.out_dir, "factors.csv"), index=False)
    with open(os.path.join(a.out_dir, "design_meta.txt"), "w") as fh:
        fh.write(f"k={k} r={a.r} levels={p} delta={p / (2 * (p - 1))} pool={a.pool} seed={a.seed}\n")
    print(f"{len(rows)} points ({a.r} trajectories x {k + 1}) -> {a.out_dir}/design.csv")


# ── run ──────────────────────────────────────────────────────────────────────────────────────────────────────────────
def metrics(series, side):
    s = pd.DataFrame(series)
    sub, bio = s.substrate_remaining_pct.to_numpy(), s.biomass_C.to_numpy()
    deg, doy = 100 - sub, np.arange(1, len(s) + 1)
    w = s.C_uptake.clip(lower=0).to_numpy()
    endo = s.R_endo_C.to_numpy()
    onset = endo > 0.001 * bio
    resp = float(s.respiration.sum())
    return dict(
        degraded_365=float(deg[-1]),
        day_50pct=float(doy[np.argmax(deg >= 50)]) if (deg >= 50).any() else 366.0,
        day_95pct=float(doy[np.argmax(deg >= 95)]) if (deg >= 95).any() else 366.0,
        CUE_annual=float(np.average(s.CUE.fillna(0), weights=w)) if w.sum() > 0 else 0.0,
        endo_C_year=float(endo.sum()),
        endo_share_pct=100 * float(endo.sum()) / resp if resp > 0 else 0.0,
        substrate_at_endo_onset=float(sub[np.argmax(onset)]) if onset.any() else 0.0,
        substrate_end=float(sub[-1]),
        biomass_max=float(bio.max()),
        biomass_mean=float(bio.mean()),
        drought_kills_per_cell=float(s.killed.sum()) / side ** 2,
    )


def _job(args):
    from sweep_minimal import run_one
    point, changes, src, clim, kin, seed, side, days = args
    _, series = run_one(src, {"R_endo": 1.0, **changes}, side, days, seed, climate=clim,
                        runtime={"Kinetics_Model": kin})
    return dict(point=point, climate=_name(clim), kinetics=KIN[kin], seed=seed, **metrics(series, side))


def run(a):
    des = pd.read_csv(os.path.join(a.out_dir, "design.csv"))
    fac = pd.read_csv(os.path.join(a.out_dir, "factors.csv"))
    src = os.path.abspath(a.input)
    jobs = []
    for _, d in des.iterrows():
        changes = {}
        for _, f in fac.iterrows():
            for row in f["rows"].split("+"):
                changes[row] = d[f.factor]
        for clim, kin, seed in itertools.product(a.climates, a.kinetics, a.seeds):
            jobs.append((int(d.point), changes, src, os.path.abspath(clim), kin, seed, a.side, a.days))
    mine = [j for i, j in enumerate(jobs) if i % a.n_tasks == a.task_id]
    os.makedirs(os.path.join(a.out_dir, "runs"), exist_ok=True)
    out = os.path.join(a.out_dir, "runs", f"task_{a.task_id:04d}.csv")
    done = set()
    if os.path.exists(out):
        prev = pd.read_csv(out)
        done = set(map(tuple, prev[KEYS].astype(str).to_numpy()))
    todo = [j for j in mine if (str(j[0]), _name(j[3]), KIN[j[4]], str(j[5])) not in done]
    print(f"task {a.task_id}/{a.n_tasks}: {len(mine)} jobs, {len(todo)} to run", flush=True)
    with Pool(a.procs) as pool:
        for n, r in enumerate(pool.imap_unordered(_job, todo), 1):
            pd.DataFrame([r]).to_csv(out, mode="a", header=not os.path.exists(out), index=False)
            if n % 10 == 0 or n == len(todo):
                print(f"  {n}/{len(todo)}", flush=True)


# ── analyze ──────────────────────────────────────────────────────────────────────────────────────────────────────────
def _indices(ee, boot, rng):
    ee = np.asarray(ee, dtype=float)
    ee = ee[np.isfinite(ee)]
    if ee.size == 0:
        return dict(mu=np.nan, mu_star=np.nan, sigma=np.nan, mu_star_lo=np.nan, mu_star_hi=np.nan, n=0)
    bs = np.abs(ee[rng.integers(0, ee.size, size=(boot, ee.size))]).mean(axis=1)
    return dict(mu=ee.mean(), mu_star=np.abs(ee).mean(), sigma=ee.std(ddof=1) if ee.size > 1 else 0.0,
                mu_star_lo=np.percentile(bs, 2.5), mu_star_hi=np.percentile(bs, 97.5), n=ee.size)


def analyze(a):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    des = pd.read_csv(os.path.join(a.out_dir, "design.csv"))
    fac = pd.read_csv(os.path.join(a.out_dir, "factors.csv"))
    k = len(fac)
    meta = open(os.path.join(a.out_dir, "design_meta.txt")).read()
    p = int(meta.split("levels=")[1].split()[0])
    delta = p / (2 * (p - 1))
    runs = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(os.path.join(a.out_dir, "runs", "task_*.csv")))])
    runs = runs.drop_duplicates(KEYS)
    pts = runs.groupby(["point", "climate", "kinetics"])[OUTPUTS].mean().reset_index()
    # warming responses as extra pseudo-climates
    wide = pts.set_index(["point", "kinetics", "climate"])[OUTPUTS].unstack("climate")
    extra = []
    for warm, lab in (("moderate_warm", "+5 C minus ambient"), ("extreme_warm", "+10 C minus ambient")):
        if ("degraded_365", warm) in wide.columns and ("degraded_365", "ambient_temp") in wide.columns:
            dif = pd.DataFrame({o: wide[(o, warm)] - wide[(o, "ambient_temp")] for o in OUTPUTS}).reset_index()
            extra.append(dif.assign(climate=lab))
    pts = pd.concat([pts] + extra, ignore_index=True)
    pts = pts.merge(des[["point", "traj", "step", "changed"] + [f"u_{f}" for f in fac.factor]], on="point")
    pts.to_csv(os.path.join(a.out_dir, "outputs_by_point.csv"), index=False)

    rng = np.random.default_rng(0)
    rows = []
    for (cl, kn), g in pts.groupby(["climate", "kinetics"]):
        g = g.set_index("point")
        ee = {f: {o: [] for o in OUTPUTS} for f in fac.factor}
        for t, gt in g.groupby("traj"):
            gt = gt.sort_values("step")
            if len(gt) != k + 1:
                continue
            for i in range(1, k + 1):
                f = gt.changed.iloc[i]
                du = gt[f"u_{f}"].iloc[i] - gt[f"u_{f}"].iloc[i - 1]
                for o in OUTPUTS:
                    ee[f][o].append((gt[o].iloc[i] - gt[o].iloc[i - 1]) / np.sign(du) / delta)
        for o in OUTPUTS:
            sd = g[o].std()
            for f in fac.factor:
                r = _indices(ee[f][o], a.boot, rng)
                rows.append(dict(climate=cl, kinetics=kn, output=o, factor=f, group=fac.set_index("factor").group[f],
                                 **r, mu_star_rel=r["mu_star"] / sd if sd > 0 else np.nan,
                                 sigma_rel=r["sigma"] / sd if sd > 0 else np.nan))
    ind = pd.DataFrame(rows)
    ind.to_csv(os.path.join(a.out_dir, "morris_indices.csv"), index=False)

    # figure: relative mu* per factor (rows) for each output (columns); one panel per climate x kinetics
    sel = ["degraded_365", "day_50pct", "CUE_annual", "endo_share_pct", "substrate_end", "biomass_max",
           "drought_kills_per_cell"]
    panels = [(c, kn) for c in ["ambient_temp", "moderate_warm", "extreme_warm", "+5 C minus ambient",
                                "+10 C minus ambient"] for kn in ["MMRT", "Arrhenius"]
              if ((ind.climate == c) & (ind.kinetics == kn)).any()]
    order = (ind[ind.output.isin(sel)].groupby("factor").mu_star_rel.mean().sort_values(ascending=False).index)
    fig, axs = plt.subplots(int(np.ceil(len(panels) / 2)), 2, figsize=(15, 4.6 * np.ceil(len(panels) / 2)), squeeze=False)
    for ax, (c, kn) in zip(axs.ravel(), panels):
        m = (ind[(ind.climate == c) & (ind.kinetics == kn) & ind.output.isin(sel)]
             .pivot(index="factor", columns="output", values="mu_star_rel").reindex(index=order, columns=sel))
        im = ax.imshow(m.to_numpy(), aspect="auto", cmap="Blues", vmin=0, vmax=np.nanpercentile(m.to_numpy(), 98))
        ax.set_yticks(range(len(order)), order, fontsize=7)
        ax.set_xticks(range(len(sel)), sel, rotation=35, ha="right", fontsize=7)
        ax.set_title(f"{c}, {kn}", fontsize=9)
        fig.colorbar(im, ax=ax, fraction=0.03, label="mu* / SD of output")
    for ax in axs.ravel()[len(panels):]:
        ax.set_visible(False)
    fig.suptitle("Morris screening, minimal model: mu* relative to the output's SD over all design points "
                 "(rows sorted by mean relative mu*)", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(os.path.join(a.out_dir, "morris_mu_star.png"), dpi=120)
    print(f"{len(runs)} runs, {pts.point.nunique()} points -> morris_indices.csv")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    d = sub.add_parser("design")
    d.add_argument("out_dir")
    d.add_argument("--factors", default=os.path.join(HERE, "morris_factors_minimal.csv"))
    d.add_argument("--r", type=int, default=10)
    d.add_argument("--levels", type=int, default=4)
    d.add_argument("--pool", type=int, default=0)
    d.add_argument("--seed", type=int, default=1)
    r = sub.add_parser("run")
    r.add_argument("out_dir")
    r.add_argument("--task-id", type=int, default=int(os.environ.get("SLURM_ARRAY_TASK_ID", 0)))
    r.add_argument("--n-tasks", type=int, default=1)
    r.add_argument("--procs", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", 2)))
    r.add_argument("--input", default=os.path.join(HERE, "..", "input_minimal_e1500_crich_km25_in"))
    r.add_argument("--climates", nargs="+", default=DEFAULT_CLIMATES)
    r.add_argument("--kinetics", type=int, nargs="+", default=[1, 0])
    r.add_argument("--seeds", type=int, nargs="+", default=[0, 1])
    r.add_argument("--side", type=int, default=20)
    r.add_argument("--days", type=int, default=365)
    z = sub.add_parser("analyze")
    z.add_argument("out_dir")
    z.add_argument("--boot", type=int, default=1000)
    a = ap.parse_args()
    {"design": design, "run": run, "analyze": analyze}[a.cmd](a)


if __name__ == "__main__":
    main()
