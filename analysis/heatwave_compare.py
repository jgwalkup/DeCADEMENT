#!/usr/bin/env python3
"""
heatwave_compare.py
===================
Paired comparison of heatwave scenarios against the baseline, community by
community, for the branch year (pulse 6 in the heatwave runs).

Every community's scenarios start from the same spin-up state, so the right
statistic for a scenario effect is the paired difference per community
(scenario − baseline), not the difference of scenario means with the
across-community spread as its uncertainty.  For each scenario and each
metric this script reports the mean paired difference, its SD, SE and 95 %
confidence interval, the relative change, the share of communities whose
difference is negative, a paired t-test and a Wilcoxon signed-rank test, and
— as a reference for what a one-year change looks like without a heatwave —
the SD across communities of the baseline's change from the previous pulse.

Metrics are the per-pulse community scalars written by export_summary_data.py
(community_pulse_metrics.csv): net CUE, C uptake / retained / respired,
endogenous respiration C, rebalancing C, unpaid energy, biomass C (mean and
end of pulse), substrate C at the end, taxa alive and cells at the end.

Input (one of)
--------------
  --summary-dir DIR   export_summary_data.py output root holding
                      community_pulse_metrics.csv (preferred: all pulses)
  --base DIR          raw simulation output (<base>/<scenario>/<community>_s<seed>/);
                      the same metrics are computed for the branch year from the
                      scenario folders, with -j workers.  Only communities present
                      in every scenario are used.

Output
------
  <out-dir>/heatwave_compare.csv         one row per scenario × metric
  <out-dir>/heatwave_compare_pairs.csv   per-community paired differences (long)
  <out-dir>/heatwave_compare.png         histograms of the paired differences (--fig)

Usage
-----
  python analysis/heatwave_compare.py --summary-dir summary_20261026 --fig
  python analysis/heatwave_compare.py --base sim_heatwave_20261026 -j 16 --fig
"""

from __future__ import annotations

import argparse
import math
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

_HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(_HERE))
from branch_output import branch_info, initialization_dir          # noqa: E402
from export_summary_data import PULSE_METRICS, pulse_metrics        # noqa: E402
from scenario_utils import discover_scenarios, pick_runs, read_runtime  # noqa: E402

try:
    from scipy import stats as _stats
except ImportError:                                               # p-values become NaN
    _stats = None


# ── Raw mode: per-run metrics for the branch year ─────────────────────────────

def _read(run_dir: Path, name: str, index_col=0):
    p = run_dir / f"{name}.csv"
    if not p.exists():
        return None
    df = pd.read_csv(p, index_col=index_col)
    df.columns = [int(c) if str(c).lstrip("-").isdigit() else c for c in df.columns]
    return df


def _e_yield(run_dir: Path) -> float | None:
    try:
        return float(pd.read_csv(initialization_dir(run_dir) / "MonEnergy_Attrib.csv",
                                 index_col=0).loc["Mon_DeadMic", "E_resp"])
    except Exception:
        return None


def metrics_from_run(run_dir: Path) -> pd.DataFrame | None:
    """pulse_metrics() from the tables in one run folder, read as they are.  A
    split scenario folder holds the branch days only, so only the branch pulse(s)
    come out; an unsplit folder gives every pulse."""
    run_dir = Path(run_dir)
    rt = read_runtime(run_dir)
    tabs = {n: _read(run_dir, n) for n in
            ["C_Uptake_Taxon", "BC_Uptake_Taxon", "Enzyme_C_Taxon", "R_endo_C_Taxon",
             "R_stoich_E_Taxon", "E_unpaid_Taxon", "Taxon_count", "MicrobesSeries",
             "SubstratesSeries", "RespSeries"]}
    cnt = tabs["Taxon_count"]
    if cnt is None or cnt.shape[1] < 2:                 # older output: final count only
        rp = _read(run_dir, "Taxon_count_repop")
        if rp is None:
            return None
        info = branch_info(run_dir)
        start = (info["branch_day"] + 1) if info else 1
        rp.columns = range(start, start + rp.shape[1])
        cnt = rp
    resp = tabs["RespSeries"]
    resp = resp.iloc[0] if resp is not None else None
    return pulse_metrics(tabs["C_Uptake_Taxon"], tabs["BC_Uptake_Taxon"], tabs["Enzyme_C_Taxon"],
                         tabs["R_endo_C_Taxon"], tabs["R_stoich_E_Taxon"], tabs["E_unpaid_Taxon"],
                         cnt, tabs["MicrobesSeries"], tabs["SubstratesSeries"], resp,
                         _e_yield(run_dir), rt["pulse"], rt["end_time"])


def _one(args):
    label, rd = args
    try:
        m = metrics_from_run(rd)
    except Exception as e:
        print(f"  [{rd.name} / {label}] FAILED: {e}", flush=True)
        return None
    if m is None:
        return None
    m = m.rename_axis("pulse").reset_index()
    m.insert(0, "community", rd.name)
    m.insert(0, "scenario", label)
    return m


def metrics_from_base(base: Path, n_jobs: int) -> pd.DataFrame:
    scenarios = discover_scenarios(base)
    if not scenarios:
        sys.exit(f"ERROR: no scenario folders under {base}")
    runs = {label: {rd.name: rd for rd in pick_runs(sdir, min_cols=100)} for label, sdir in scenarios}
    common = set.intersection(*(set(r) for r in runs.values()))
    skipped = set.union(*(set(r) for r in runs.values())) - common
    if skipped:
        print(f"{len(skipped)} communities missing from some scenario, skipped: {sorted(skipped)[:10]}"
              f"{' …' if len(skipped) > 10 else ''}")
    jobs = [(label, runs[label][c]) for label in runs for c in sorted(common)]
    print(f"Computing branch-year metrics for {len(common)} communities × {len(runs)} scenarios "
          f"({n_jobs} workers)...", flush=True)
    if n_jobs > 1:
        with ProcessPoolExecutor(max_workers=n_jobs) as ex:
            parts = list(ex.map(_one, jobs, chunksize=4))
    else:
        parts = [_one(j) for j in jobs]
    parts = [p for p in parts if p is not None]
    if not parts:
        sys.exit("ERROR: no metrics computed")
    return pd.concat(parts, ignore_index=True)


# ── Paired statistics ─────────────────────────────────────────────────────────

def _t_ci(d: np.ndarray, alpha: float = 0.05):
    n = d.size
    if n < 2:
        return np.nan, np.nan, np.nan
    se = d.std(ddof=1) / math.sqrt(n)
    tq = _stats.t.ppf(1 - alpha / 2, n - 1) if _stats else 1.96
    return se, d.mean() - tq * se, d.mean() + tq * se


def paired_table(pm: pd.DataFrame, baseline: str, pulse: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Per-scenario × metric paired statistics and the per-community differences."""
    scen = [s for s in pm.scenario.unique() if s != baseline]
    if baseline not in set(pm.scenario):
        sys.exit(f"ERROR: baseline scenario {baseline!r} not in the table; have {sorted(pm.scenario.unique())}")
    cur = pm[pm.pulse == pulse].set_index(["scenario", "community"])
    base = cur.loc[baseline]

    # Reference: the baseline's own change from the previous pulse (no heatwave)
    pn = int(pulse.split()[-1])
    prev = f"Pulse {pn - 1}"
    ref = None
    if prev in set(pm.pulse):
        bprev = pm[(pm.pulse == prev) & (pm.scenario == baseline)].set_index("community")
        common = base.index.intersection(bprev.index)
        ref = (base.loc[common, PULSE_METRICS] - bprev.loc[common, PULSE_METRICS])

    rows, pairs = [], []
    for s in scen:
        other = cur.loc[s]
        common = base.index.intersection(other.index)
        diff = other.loc[common, PULSE_METRICS] - base.loc[common, PULSE_METRICS]
        long = diff.rename_axis("community").reset_index().melt(id_vars="community", var_name="metric", value_name="diff")
        long.insert(0, "scenario", s)
        pairs.append(long)
        for m in PULSE_METRICS:
            d = diff[m].to_numpy(float)
            d = d[np.isfinite(d)]
            b = base.loc[common, m].to_numpy(float)
            se, lo, hi = _t_ci(d)
            if _stats and d.size > 1:
                t_p = float(_stats.ttest_1samp(d, 0.0).pvalue)
                try:
                    w_p = float(_stats.wilcoxon(d).pvalue) if np.any(d != 0) else 1.0
                except ValueError:
                    w_p = np.nan
            else:
                t_p = w_p = np.nan
            bmean = float(np.nanmean(b))
            rows.append(dict(
                scenario=s, metric=m, n=int(d.size),
                baseline_mean=bmean, scenario_mean=float(np.nanmean(other.loc[common, m])),
                diff_mean=float(d.mean()) if d.size else np.nan,
                diff_sd=float(d.std(ddof=1)) if d.size > 1 else np.nan,
                diff_se=se, ci95_lo=lo, ci95_hi=hi,
                rel_change_pct=(100.0 * d.mean() / bmean) if (d.size and bmean) else np.nan,
                frac_negative=float((d < 0).mean()) if d.size else np.nan,
                t_p=t_p, wilcoxon_p=w_p,
                baseline_prev_pulse_change_sd=(float(ref[m].std(ddof=1)) if ref is not None else np.nan),
                baseline_prev_pulse_change_mean=(float(ref[m].mean()) if ref is not None else np.nan),
            ))
    return pd.DataFrame(rows), pd.concat(pairs, ignore_index=True)


def spinup_check(pm: pd.DataFrame, baseline: str, pulse: str) -> None:
    """Pulses before the branch must be identical across scenarios (shared spin-up)."""
    pn = int(pulse.split()[-1])
    pre = pm[pm.pulse.map(lambda p: int(p.split()[-1]) < pn)]
    if pre.empty:
        return
    base = pre[pre.scenario == baseline].set_index(["community", "pulse"])[PULSE_METRICS]
    worst = 0.0
    for s in pre.scenario.unique():
        if s == baseline:
            continue
        o = pre[pre.scenario == s].set_index(["community", "pulse"])[PULSE_METRICS]
        common = base.index.intersection(o.index)
        worst = max(worst, float(np.nanmax(np.abs(o.loc[common].to_numpy(float) - base.loc[common].to_numpy(float)))))
    print(f"Spin-up check: max |scenario − baseline| over pulses before {pulse}: {worst:.3g} "
          f"({'identical, as expected' if worst == 0 else 'NOT identical — runs are not paired'})\n")


def _fmt(df: pd.DataFrame) -> str:
    cols = ["scenario", "metric", "n", "baseline_mean", "diff_mean", "diff_se", "ci95_lo", "ci95_hi",
            "rel_change_pct", "frac_negative", "t_p", "wilcoxon_p", "baseline_prev_pulse_change_sd"]
    out = df[cols].copy()
    for c in cols[3:]:
        out[c] = out[c].map(lambda v: f"{v:.4g}" if pd.notna(v) else "—")
    return out.to_string(index=False)


def make_figure(pairs: pd.DataFrame, table: pd.DataFrame, path: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    metrics = ["net_cue", "c_resp", "c_endo_resp", "biomass_C_end"]
    scen = list(pairs.scenario.unique())
    fig, axes = plt.subplots(1, len(metrics), figsize=(4.2 * len(metrics), 3.6))
    for ax, m in zip(axes, metrics):
        sub = pairs[pairs.metric == m]
        vals = sub["diff"].to_numpy(float)
        vals = vals[np.isfinite(vals)]
        if vals.size == 0:
            ax.set_visible(False)
            continue
        bins = np.linspace(np.nanmin(vals), np.nanmax(vals), 41)
        for s in scen:
            d = sub[sub.scenario == s]["diff"].to_numpy(float)
            ax.hist(d[np.isfinite(d)], bins=bins, alpha=0.55, label=s.replace("climate_", ""))
            r = table[(table.scenario == s) & (table.metric == m)].iloc[0]
            ax.axvline(r.diff_mean, color=ax.patches[-1].get_facecolor(), lw=1.5)
        ax.axvline(0, color="k", lw=0.8, ls="--")
        ax.set_title(m)
        ax.set_xlabel("scenario − baseline (per community)")
    axes[0].set_ylabel("communities")
    axes[0].legend(frameon=False, fontsize=8)
    fig.suptitle("Paired branch-year differences vs. baseline", y=1.02)
    fig.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    print(f"Figure: {path}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--summary-dir", metavar="DIR", help="export_summary_data.py output root")
    src.add_argument("--base", metavar="DIR", help="raw simulation output root")
    ap.add_argument("--baseline", default="climate_ambient_temp")
    ap.add_argument("--pulse", type=int, default=None,
                    help="pulse to compare (default: the last pulse present)")
    ap.add_argument("--out-dir", metavar="DIR", default=None,
                    help="where to write results (default: --summary-dir, or <base>/heatwave_compare)")
    ap.add_argument("--jobs", "-j", type=int, default=1, help="workers for --base mode (0 = all cores)")
    ap.add_argument("--fig", action="store_true", help="also write heatwave_compare.png")
    a = ap.parse_args()

    if a.summary_dir:
        root = Path(a.summary_dir)
        f = root / "community_pulse_metrics.csv"
        if not f.exists():
            sys.exit(f"ERROR: {f} not found — re-run export_summary_data.py (current version writes it), "
                     f"or use --base to compute from the raw runs.")
        pm = pd.read_csv(f)
        out = Path(a.out_dir) if a.out_dir else root
    else:
        import os
        n_jobs = a.jobs if a.jobs > 0 else (int(os.environ.get("SLURM_CPUS_PER_TASK", 0)) or os.cpu_count() or 1)
        pm = metrics_from_base(Path(a.base), n_jobs)
        out = Path(a.out_dir) if a.out_dir else Path(a.base) / "heatwave_compare"
    out.mkdir(parents=True, exist_ok=True)

    pulses = sorted(pm.pulse.unique(), key=lambda p: int(p.split()[-1]))
    pulse = f"Pulse {a.pulse}" if a.pulse else pulses[-1]
    if pulse not in pulses:
        sys.exit(f"ERROR: {pulse} not in the table; have {pulses}")
    n_comm = pm[pm.pulse == pulse].groupby("scenario").community.nunique()
    print(f"Comparing {pulse} against {a.baseline}; communities per scenario: "
          + ", ".join(f"{s}={n}" for s, n in n_comm.items()) + "\n")
    spinup_check(pm, a.baseline, pulse)

    table, pairs = paired_table(pm, a.baseline, pulse)
    table.to_csv(out / "heatwave_compare.csv", index=False)
    pairs.to_csv(out / "heatwave_compare_pairs.csv", index=False)
    print(_fmt(table))
    print(f"\nWrote {out / 'heatwave_compare.csv'} and heatwave_compare_pairs.csv")
    if _stats is None:
        print("(scipy not available: p-values omitted, CI uses z=1.96)")
    if a.fig:
        make_figure(pairs, table, out / "heatwave_compare.png")


if __name__ == "__main__":
    main()
