#!/usr/bin/env python3
"""Scan community outputs for duplicate taxon labels.

Usage:
  # one community
  python diag_dupes.py <community_dir>
  # every community under a scenario dir or a whole temp_gradient dir
  python diag_dupes.py --scan <dir>
"""
import sys, pandas as pd
from pathlib import Path

def dup_labels(path):
    if not path.exists():
        return None
    idx = pd.read_csv(path, index_col=0).index
    return idx[idx.duplicated(keep=False)].unique().tolist()

def check_one(run, verbose=True):
    """Return dict of file->dup-labels for the files that have duplicates."""
    findings = {}
    for rel in ["MicrobesSeries.csv", "CUE_TaxonSeries.csv",
                "Initialization/EnzGenes.csv", "Initialization/Microbes_pp.csv"]:
        d = dup_labels(run / rel)
        if d:
            findings[rel] = d
    if verbose:
        print(f"Inspecting: {run}")
        for rel in ["MicrobesSeries.csv", "CUE_TaxonSeries.csv",
                    "Initialization/EnzGenes.csv", "Initialization/Microbes_pp.csv"]:
            p = run / rel
            if not p.exists():
                print(f"  [missing] {rel}"); continue
            idx = pd.read_csv(p, index_col=0).index
            d = idx[idx.duplicated(keep=False)].unique().tolist()
            flag = f"  <-- DUP {d[:10]}" if d else ""
            print(f"  {rel:34s} rows={len(idx):4d} unique={idx.nunique():4d}{flag}")
    return findings

def is_run_dir(p):
    return (p / "MicrobesSeries.csv").exists() or (p / "Initialization").is_dir()

if len(sys.argv) >= 3 and sys.argv[1] == "--scan":
    root = Path(sys.argv[2])
    # find all run dirs (community_* with output) under root, any depth up to 2
    runs = sorted({p.parent if p.name == "Initialization" else p
                   for p in root.rglob("Initialization") if p.is_dir()})
    if not runs:
        runs = sorted(p for p in root.rglob("community_*") if is_run_dir(p))
    print(f"Scanning {len(runs)} community dirs under {root}\n")
    bad = []
    for r in runs:
        f = check_one(r, verbose=False)
        if f:
            bad.append((r, f))
    if not bad:
        print("No duplicate taxon labels found in any community.")
    else:
        print(f"{len(bad)} community(ies) with duplicate taxon labels:\n")
        for r, f in bad:
            scen = r.parent.name
            for rel, dups in f.items():
                print(f"  {scen}/{r.name}  {rel}: {dups[:10]}{' …' if len(dups)>10 else ''}")
    print(f"\nSummary: {len(bad)}/{len(runs)} communities affected.")
else:
    check_one(Path(sys.argv[1]))
