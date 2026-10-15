#!/usr/bin/env python3
"""
extract_community_taxa.py
=========================
Extract taxon IDs and membership counts for a range of communities
from a DEMENTpy community pool.

For each community in the specified range, reads Microbes_pp.csv to get
the taxon IDs present in that community. Produces two output files:

  <output-dir>/community_taxa.csv
      One row per (community, taxon) with community ID, taxon ID,
      fungi/bacteria classification (fb), and initial biomass (C, N, P, E).

  <output-dir>/taxon_summary.csv
      One row per unique taxon across all communities with:
        - taxon_id        : original pool ID (e.g. Tax47)
        - n_communities   : how many communities this taxon appears in
        - community_ids   : comma-separated list of community IDs

USAGE
-----
    python analysis/extract_community_taxa.py \\
        --community-dir community_pool_2026051301/communities \\
        --range         1 30 \\
        --output-dir    analysis_output/community_taxa

OPTIONS
-------
    --community-dir DIR    Directory containing community_NNNN/ subdirectories
    --range         START END
                           Inclusive range of community indices to process
                           (default: all communities found)
    --output-dir    DIR    Directory to write output CSV files (default: .)
    --fb-file       FILE   Path to fb.csv inside each community dir, used to
                           classify fungi (1) vs bacteria (0) (default: fb.csv)
    --quiet                Suppress progress output

EXAMPLES
--------
    # Extract taxa for communities 1-30
    python analysis/extract_community_taxa.py \\
        --community-dir community_pool_2026051301/communities \\
        --range 1 30 \\
        --output-dir analysis_output/community_taxa

    # All communities in the pool
    python analysis/extract_community_taxa.py \\
        --community-dir community_pool_2026051301/communities \\
        --output-dir analysis_output/community_taxa
"""

import argparse
import sys
from pathlib import Path
from collections import defaultdict

import numpy as np
import pandas as pd


# ── Community directory discovery ─────────────────────────────────────────────

def find_community_dirs(community_dir: Path, start: int, end: int) -> list[tuple[int, Path]]:
    """Return sorted list of (community_index, path) for directories in [start, end]."""
    dirs = sorted(community_dir.glob("community_*"))
    result = []
    for d in dirs:
        if not d.is_dir():
            continue
        # Extract numeric index from name (community_0001 → 1)
        try:
            idx = int(d.name.split("_")[-1])
        except ValueError:
            continue
        if start <= idx <= end:
            result.append((idx, d))
    return result


# ── Per-community extraction ───────────────────────────────────────────────────

def extract_community(comm_idx: int, comm_dir: Path) -> pd.DataFrame | None:
    """
    Read Microbes_pp.csv and fb.csv from a community directory.

    Returns a DataFrame with columns:
        community_id, taxon_id, fb (0=bacteria, 1=fungi), C, N, P, E
    Returns None if required files are missing.
    """
    mic_path = comm_dir / "Microbes_pp.csv"
    fb_path  = comm_dir / "fb.csv"

    if not mic_path.exists():
        print(f"  WARNING: Microbes_pp.csv not found in {comm_dir}", file=sys.stderr)
        return None

    try:
        mic = pd.read_csv(mic_path, index_col=0)
    except Exception as e:
        print(f"  WARNING: could not read {mic_path}: {e}", file=sys.stderr)
        return None

    # fb: fungi/bacteria index (0 = bacteria, 1 = fungi)
    fb_vals = None
    if fb_path.exists():
        try:
            fb_raw  = pd.read_csv(fb_path, header=None)
            fb_vals = fb_raw.iloc[:, 0].values
        except Exception:
            pass

    rows = []
    for i, taxon_id in enumerate(mic.index):
        row = {
            "community_id": comm_idx,
            "taxon_id":     str(taxon_id),
            "fb":           int(fb_vals[i]) if fb_vals is not None and i < len(fb_vals) else np.nan,
        }
        for col in ["C", "N", "P", "E"]:
            row[col] = float(mic.loc[taxon_id, col]) if col in mic.columns else np.nan
        rows.append(row)

    return pd.DataFrame(rows)


# ── Main ───────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Extract taxon IDs and counts for a range of DEMENTpy communities.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "--community-dir", required=True, metavar="DIR",
        help="Directory containing community_NNNN/ subdirectories.",
    )
    parser.add_argument(
        "--range", nargs=2, type=int, default=None, metavar=("START", "END"),
        help="Inclusive range of community indices to process. Default: all.",
    )
    parser.add_argument(
        "--output-dir", default=".", metavar="DIR",
        help="Directory to write output CSV files (default: current directory).",
    )
    parser.add_argument(
        "--quiet", action="store_true",
        help="Suppress progress output.",
    )
    args = parser.parse_args()

    community_dir = Path(args.community_dir)
    if not community_dir.exists():
        print(f"ERROR: --community-dir not found: {community_dir}", file=sys.stderr)
        sys.exit(1)

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # ── Determine range ────────────────────────────────────────────────────────
    if args.range:
        start, end = args.range
    else:
        # Auto-detect from directory names
        all_dirs = sorted(community_dir.glob("community_*"))
        indices  = []
        for d in all_dirs:
            try:
                indices.append(int(d.name.split("_")[-1]))
            except ValueError:
                pass
        if not indices:
            print(f"ERROR: no community directories found in {community_dir}", file=sys.stderr)
            sys.exit(1)
        start, end = min(indices), max(indices)

    comm_dirs = find_community_dirs(community_dir, start, end)
    if not comm_dirs:
        print(f"ERROR: no communities found in range {start}–{end}", file=sys.stderr)
        sys.exit(1)

    if not args.quiet:
        print(f"Community dir : {community_dir}")
        print(f"Range         : {start}–{end}  ({len(comm_dirs)} communities)")
        print(f"Output dir    : {output_dir}")
        print()

    # ── Extract ────────────────────────────────────────────────────────────────
    all_rows = []
    n_failed = 0

    for i, (comm_idx, comm_dir) in enumerate(comm_dirs):
        if not args.quiet:
            print(f"\r  Processing community {comm_idx:04d} ({i+1}/{len(comm_dirs)})...",
                  end="", flush=True)

        df = extract_community(comm_idx, comm_dir)
        if df is not None:
            all_rows.append(df)
        else:
            n_failed += 1

    if not args.quiet:
        print()  # newline after progress

    if not all_rows:
        print("ERROR: no data extracted — check community directory structure.", file=sys.stderr)
        sys.exit(1)

    # ── Community-taxon table ──────────────────────────────────────────────────
    comm_taxa = pd.concat(all_rows, ignore_index=True)
    comm_taxa["fb_label"] = comm_taxa["fb"].map({0: "bacteria", 1: "fungi"})

    out_ct = output_dir / "community_taxa.csv"
    comm_taxa.to_csv(out_ct, index=False)

    # ── Taxon summary ──────────────────────────────────────────────────────────
    # Count communities per taxon and list which ones
    taxon_communities = defaultdict(list)
    for _, row in comm_taxa.iterrows():
        taxon_communities[row["taxon_id"]].append(int(row["community_id"]))

    summary_rows = []
    for taxon_id, comm_list in taxon_communities.items():
        # Get fb label from first occurrence
        fb_val = comm_taxa.loc[
            comm_taxa["taxon_id"] == taxon_id, "fb"
        ].iloc[0]
        fb_label = {0: "bacteria", 1: "fungi"}.get(int(fb_val) if not np.isnan(fb_val) else -1, "unknown")

        # Mean initial biomass across communities (should be identical — same traits)
        mean_C = comm_taxa.loc[comm_taxa["taxon_id"] == taxon_id, "C"].mean()

        summary_rows.append({
            "taxon_id":       taxon_id,
            "fb":             int(fb_val) if not np.isnan(fb_val) else np.nan,
            "fb_label":       fb_label,
            "n_communities":  len(comm_list),
            "community_ids":  ",".join(str(c) for c in sorted(comm_list)),
            "mean_C_init":    round(mean_C, 6),
        })

    taxon_summary = (
        pd.DataFrame(summary_rows)
        .sort_values("n_communities", ascending=False)
        .reset_index(drop=True)
    )

    out_ts = output_dir / "taxon_summary.csv"
    taxon_summary.to_csv(out_ts, index=False)

    # ── Report ─────────────────────────────────────────────────────────────────
    n_comm     = comm_taxa["community_id"].nunique()
    n_taxa     = comm_taxa["taxon_id"].nunique()
    n_bacteria = taxon_summary[taxon_summary["fb_label"] == "bacteria"].shape[0]
    n_fungi    = taxon_summary[taxon_summary["fb_label"] == "fungi"].shape[0]
    mean_per_comm = comm_taxa.groupby("community_id")["taxon_id"].count().mean()
    mean_comm_per_taxon = taxon_summary["n_communities"].mean()

    print(f"\n{'─'*55}")
    print(f"  Communities processed : {n_comm}  ({n_failed} failed)")
    print(f"  Unique taxa found     : {n_taxa}")
    print(f"    Bacteria            : {n_bacteria}")
    print(f"    Fungi               : {n_fungi}")
    print(f"  Mean taxa/community   : {mean_per_comm:.1f}")
    print(f"  Mean communities/taxon: {mean_comm_per_taxon:.2f}")
    print(f"  Most common taxon     : {taxon_summary.iloc[0]['taxon_id']} "
          f"(in {taxon_summary.iloc[0]['n_communities']} communities)")
    print(f"  Singletons (1 comm)   : "
          f"{(taxon_summary['n_communities'] == 1).sum()}")
    print(f"{'─'*55}")
    print(f"\n  Saved:")
    print(f"    {out_ct}")
    print(f"    {out_ts}")


if __name__ == "__main__":
    main()
