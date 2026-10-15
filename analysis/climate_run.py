"""
Multi-pulse runs of the minimal scenario under a real climate file, with and without HSP production.

Mirrors src/dementpy.py's pulse loop: each pulse is `end_time` days; at the end of a pulse the real
Grid.repopulation() resets substrate, monomers and enzymes and re-seeds the microbial community
(a stub stands in for the Output class, which repopulation() only uses for per-taxon cell counts).

Usage:
    python climate_run.py ../input_minimal ../input/climate.csv ../analysis_output/climate_run_e1500 \\
        [--pulses 6] [--seeds 0 1 2 3 4] [--e-resp 1500] [--plot-only]

Parameter set (E_resp 1500): intercept 45k, slope -90k, HSP cost 300, Maint_E_Rate_Const 2.5, inducible-dominated HSP.
Writes <out_dir>/timeseries.csv and climate_run_timeseries.png, climate_run_pulses.png.
"""
import argparse
import os
import shutil
import sys
import tempfile
from multiprocessing import Pool

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "src"))
sys.path.insert(0, HERE)
from initialization import initialize_data  # noqa: E402
from grid import Grid  # noqa: E402
from sweep_minimal import set_params, set_runtime  # noqa: E402
from scenario_utils import scenario_color_map  # noqa: E402

plt.rcParams.update({
    "font.size": 10, "axes.titlesize": 11, "axes.labelsize": 10,
    "legend.fontsize": 9, "figure.dpi": 150,
})

PARAMS = {"HSP_Induci_Prod": 0.3, "HSP_Consti_Prod": 0.0003, "Maint_E_Rate_Const": 2.5,
          "Maint_E_HSP_Ea_Intercept": 45000, "Maint_E_HSP_Ea_Slope": -90000, "HSP_prod_E_cost": 300}
OFF = {"HSP_Induci_Prod": 0.0, "HSP_Consti_Prod": 0.0}
SIDE = 20
CYCLE = 365


class _StubOutput:
    """Minimal stand-in for Output: repopulation() reads Taxon_count_repop and calls record_pulse_init."""

    def __init__(self, n_taxa, total_days):
        self.Taxon_count_repop = pd.DataFrame(np.zeros((n_taxa, total_days + 2)))

    def record_pulse_init(self, grid, next_pulse):
        pass


def run_pulses(args):
    src, climate, changes, e_resp, pulses, seed, hsp, SIDE = args
    tmp = tempfile.mkdtemp()
    work = os.path.join(tmp, "in")
    shutil.copytree(src, work)
    shutil.copy(climate, os.path.join(work, "climate.csv"))
    set_params(os.path.join(work, "parameters.csv"), changes)
    set_runtime(os.path.join(work, "runtime.txt"),
                dict(x=SIDE, y=SIDE, gridsize=SIDE * SIDE, pulse=pulses, end_time=CYCLE))
    p = os.path.join(work, "monomer_energy.csv")
    out = []
    for l in open(p).read().split("\n"):
        c = l.split(",")
        out.append(f"{c[0]},{e_resp},{c[2]},{c[3]}" if c[0] == "Mon_IdealSub" else l)
    open(p, "w").write("\n".join(out))
    clim = pd.read_csv(os.path.join(work, "climate.csv"), index_col=0)
    cwd = os.getcwd()
    os.chdir(work)
    rows = []
    try:
        np.random.seed(seed)
        rt = pd.read_csv("runtime.txt", header=None, index_col=0, sep="\t").drop("output_format").astype(float)
        data = initialize_data(rt)
        eco = Grid(rt, data)
        G = SIDE * SIDE
        sub = [s for s in data["Substrate_index"] if s not in ("DeadMic", "DeadEnz")][0]
        sub0 = float(np.asarray(eco.Substrates.loc[sub, "C"]).sum())
        stub = _StubOutput(1, pulses * CYCLE)
        ps = lambda x: float(np.asarray(x).sum())
        for pl in range(pulses):
            for i in range(pl * CYCLE, (pl + 1) * CYCLE):
                eco.degradation(pl, i)
                eco.uptake(pl, i)
                eco.metabolism(i)
                eco.mortality(i)
                eco.reproduction(i)
                C = np.asarray(eco.Microbes["C"])
                E = np.asarray(eco.Microbes["E"])
                alive = C > 0
                n = int(alive.sum())
                stub.Taxon_count_repop.iloc[0, i + 1] = n
                nn = max(n, 1)
                rows.append(dict(
                    hsp=hsp, seed=seed, pulse=pl + 1, day=i + 1, doy=i - pl * CYCLE + 1,
                    Temp=clim["Temp"].iloc[i], Psi=clim["Psi"].iloc[i], FM=clim["FM"].iloc[i],
                    CUE=float(eco.CUE_system),
                    energy_to_C=E[alive].sum() / C[alive].sum() if n else np.nan,
                    biomass_C=C.sum() / G, alive_cells=n,
                    substrate_remaining_pct=100 * ps(eco.Substrates.loc[sub, "C"]) / sub0,
                    net_energy=ps(eco.Microbe_E_Gain) / nn,
                    HSP_prop=float(np.asarray(eco.HSP_prop)[alive].mean()) if n else np.nan,
                    maint_Ea=float(np.asarray(eco.maint_E_cost_Ea)[alive].mean()) if n else np.nan,
                    osmo_E_cost=ps(eco.Taxon_Osmo_E_Cost_total) / nn,
                    E_cost_total=ps(eco.Microbe_E_Cost) / nn,
                ))
                if i == (pl + 1) * CYCLE - 1 and pl < pulses - 1:
                    eco.repopulation(stub, i, True)
    finally:
        os.chdir(cwd)
        shutil.rmtree(tmp, ignore_errors=True)
    return pd.DataFrame(rows)


def plots(ts, out_dir, pulses, n_seeds):
    pcol = scenario_color_map([str(p) for p in range(1, pulses + 1)], palette="discrete")
    one = ts[(ts["pulse"] == 1) & (ts["seed"] == ts["seed"].min()) & ts["hsp"]]
    m = ts.groupby(["hsp", "pulse", "doy"], as_index=False)[
        ["CUE", "energy_to_C", "biomass_C", "alive_cells", "substrate_remaining_pct", "HSP_prop"]].mean()
    rows = [("CUE", "CUE", (0, 1)), ("energy_to_C", "Energy : C", (0, None)), ("biomass_C", "Biomass C", (0, None)),
            ("alive_cells", "Cells alive", (0, 1.03 * ts["alive_cells"].max())), ("substrate_remaining_pct", "Substrate remaining (%)", (0, 100))]
    fig = plt.figure(figsize=(13, 15))
    gs = fig.add_gridspec(len(rows) + 1, 2, height_ratios=[0.9] + [1] * len(rows))
    axc = fig.add_subplot(gs[0, :])
    axc.plot(one["doy"], one["Temp"], color="#d6172a", lw=1.3, label="Temperature (°C)")
    axc.set_ylabel("Temperature (°C)")
    ax2 = axc.twinx()
    ax2.plot(one["doy"], one["Psi"], color="#2166ac", lw=1.0, alpha=0.8, label="Water potential Ψ")
    ax2.set_ylabel("Ψ (MPa)")
    axc.set_title("Climate (one year, repeated each pulse; day 1 = 15 Dec)")
    axc.spines[["top"]].set_visible(False)
    ax2.spines[["top"]].set_visible(False)
    for i, (col, lab, yl) in enumerate(rows):
        for j, hsp in enumerate((True, False)):
            ax = fig.add_subplot(gs[i + 1, j])
            for pl in range(1, pulses + 1):
                d = m[(m["hsp"] == hsp) & (m["pulse"] == pl)].sort_values("doy")
                y = d[col].rolling(7, center=True, min_periods=1).mean() if col in ("CUE", "energy_to_C") else d[col]
                ax.plot(d["doy"], y, color=pcol[str(pl)], lw=1.2, label=f"pulse {pl}")
            ax.set_ylim(*yl)
            if i == 0:
                ax.set_title("With HSP production" if hsp else "No HSP production (control)")
            if j == 0:
                ax.set_ylabel(lab)
            if i == len(rows) - 1:
                ax.set_xlabel("Day of pulse")
            ax.spines[["top", "right"]].set_visible(False)
    h, l = fig.axes[2].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncols=pulses, title="Pulse (mean of seeds)", framealpha=0.85,
               bbox_to_anchor=(0.5, -0.035))
    fig.suptitle(f"Real-climate runs, E_resp 1500 set (mean of {n_seeds} seeds)", fontsize=11, fontweight="bold", y=1.0)
    fig.tight_layout(rect=(0, 0, 1, 0.985))
    fig.savefig(os.path.join(out_dir, "climate_run_timeseries.png"), dpi=150, bbox_inches="tight")
    plt.close(fig)

    # per-pulse end-of-pulse summary (seed mean +/- SD)
    last = ts[ts["doy"] == CYCLE]
    cue = ts[ts["doy"] > CYCLE - 30].groupby(["hsp", "seed", "pulse"])["CUE"].mean().reset_index()
    panels = [("biomass_C", "Biomass C (day 365 of pulse)", last), ("CUE", "CUE (last 30 days)", cue),
              ("alive_cells", "Cells alive", last), ("substrate_remaining_pct", "Substrate remaining (%)", last)]
    fig, axes = plt.subplots(1, 4, figsize=(16, 4.2))
    for ax, (col, title, d) in zip(axes, panels):
        for hsp, c, ls in ((True, "#1f77b4", "-"), (False, "gray", "--")):
            g = d[d["hsp"] == hsp].groupby("pulse")[col]
            ax.errorbar(g.mean().index, g.mean().values, yerr=g.std().values, color=c, ls=ls, marker="o", ms=4,
                        capsize=2, lw=1.5, label="HSP" if hsp else "no HSP (control)")
        ax.set_title(title, fontsize=10)
        ax.set_xlabel("Pulse")
        if col == "CUE":
            ax.set_ylim(0, 1)
        ax.spines[["top", "right"]].set_visible(False)
    axes[0].legend()
    fig.suptitle("End-of-pulse outcomes across pulses (mean ± SD over seeds)", fontsize=11, fontweight="bold")
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "climate_run_pulses.png"), dpi=150, bbox_inches="tight")
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("input_folder")
    ap.add_argument("climate_file")
    ap.add_argument("out_dir")
    ap.add_argument("--pulses", type=int, default=6)
    ap.add_argument("--seeds", type=int, nargs="+", default=[0, 1, 2, 3, 4])
    ap.add_argument("--e-resp", type=float, default=1500)
    ap.add_argument("--side", type=int, default=20, help="grid side (side x side cells)")
    ap.add_argument("--set", action="append", default=[], metavar="NAME=VALUE",
                    help="override a parameter of the HSP/control parameter set (repeatable)")
    ap.add_argument("--plot-only", action="store_true")
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    csv = os.path.join(a.out_dir, "timeseries.csv")
    if not a.plot_only:
        src, clim = os.path.abspath(a.input_folder), os.path.abspath(a.climate_file)
        params = {**PARAMS, **{k: float(v) for k, v in (x.split("=") for x in a.set)}}
        jobs = []
        for s in a.seeds:
            jobs.append((src, clim, dict(params), a.e_resp, a.pulses, s, True, a.side))
            jobs.append((src, clim, {**params, **OFF}, a.e_resp, a.pulses, s, False, a.side))
        print(f"{len(jobs)} runs of {a.pulses} pulses", file=sys.stderr, flush=True)
        with Pool(min(len(jobs), os.cpu_count())) as pool:
            parts = pool.map(run_pulses, jobs, chunksize=1)
        pd.concat(parts, ignore_index=True).to_csv(csv, index=False)
    ts = pd.read_csv(csv)
    plots(ts, a.out_dir, int(ts["pulse"].max()), ts["seed"].nunique())
    print("wrote figures in", a.out_dir)


if __name__ == "__main__":
    main()
