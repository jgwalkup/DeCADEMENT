#!/usr/bin/env python3
"""Run the REAL compute_run on each community and report any per-taxon frame
whose row index is non-unique — i.e. exactly what would crash the fold.

This catches duplicates in DERIVED frames (energy_C, CU-BC, hsp_*, fractions,
etc.), not just the raw output files diag_dupes.py inspects.

Usage:
  python diag_fold_dupes.py --scan <temp_gradient_dir_or_scenario_dir>
  python diag_fold_dupes.py <one_community_dir>
"""
import sys
from pathlib import Path

# import compute_run from the sibling export script
sys.path.insert(0, str(Path(__file__).resolve().parent))
import export_summary_data as E
import pandas as pd


def check_run(run):
    """Return list of (key, dup_labels) for frames with non-unique row index."""
    d = E.compute_run(run)
    if d is None:
        return "skipped"           # compute_run itself rejected/guarded it
    bad = []
    for key, val in d.items():
        if isinstance(val, pd.DataFrame) and not val.index.is_unique:
            dups = val.index[val.index.duplicated(keep=False)].unique().tolist()
            bad.append((key, dups))
    return bad


def is_run_dir(p):
    return (p / "MicrobesSeries.csv").exists() or (p / "Initialization").is_dir()


if len(sys.argv) >= 3 and sys.argv[1] == "--scan":
    root = Path(sys.argv[2])
    runs = sorted({p.parent for p in root.rglob("Initialization") if p.is_dir()})
    if not runs:
        runs = sorted(p for p in root.rglob("community_*") if is_run_dir(p))
    print(f"Running compute_run on {len(runs)} communities under {root}\n")
    n_bad = n_skip = 0
    for r in runs:
        res = check_run(r)
        if res == "skipped":
            n_skip += 1
            print(f"  [skip]  {r.parent.name}/{r.name}")
        elif res:
            n_bad += 1
            for key, dups in res:
                print(f"  [DUP]   {r.parent.name}/{r.name}  frame '{key}': "
                      f"{dups[:10]}{' …' if len(dups) > 10 else ''}")
    print(f"\nSummary: {n_bad} communities with duplicate frames, "
          f"{n_skip} skipped, {len(runs)-n_bad-n_skip} clean / {len(runs)} total.")
else:
    run = Path(sys.argv[1])
    res = check_run(run)
    if res == "skipped":
        print(f"{run.name}: compute_run returned None (guarded/skipped).")
    elif not res:
        print(f"{run.name}: all frames have unique indices (clean).")
    else:
        print(f"{run.name}: duplicate-index frames:")
        for key, dups in res:
            print(f"  frame '{key}': {dups[:10]}{' …' if len(dups) > 10 else ''}")
