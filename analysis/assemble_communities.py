"""
assemble_communities.py
=======================
Generate a fixed set of microbial communities by drawing from a pre-generated
taxon pool. Produces two sets of community directories:

  <output-dir>/single_taxon/taxon_<ID>/   -- one directory per taxon in the pool
                                             each containing a single-taxon community
  <output-dir>/communities/community_<ID>/ -- one directory per ensemble community
                                             each containing m taxa drawn without
                                             replacement from the pool

Each community directory contains one CSV per microbial trait matrix, in the
exact format expected by the modified initialize_data() in initialization.py.
All communities are generated from a single seed so the full set is fixed and
reproducible before any simulations are run.

Expected project layout
-----------------------
    Decoupled-Energy-DEMENTpy/   <- run from here
    ├── src/                     <- --src-dir
    ├── input/                   <- --input-dir
    └── analysis/
        └── assemble_communities.py

Usage
-----
    python analysis/assemble_communities.py \\
        --src-dir       src \\
        --input-dir     input \\
        --output-dir    community_pool \\
        --n-communities 1200 \\
        --community-size 100 \\
        --seed          99 \\
        --mode          both

Options
-------
  --src-dir        Directory containing DEMENTpy .py source files (default: src)
  --input-dir      Directory containing input CSVs and runtime.txt (default: input)
  --output-dir     Root output directory (default: community_pool)
  --runtime        runtime.txt filename inside --input-dir (default: runtime.txt)
  --n-communities  Number of ensemble communities to generate (default: 1000)
  --community-size Number of taxa per ensemble community (default: 100)
  --seed           Random seed for reproducibility (default: None)
  --mode           What to generate: 'single', 'ensemble', or 'both' (default: both)

Community directory contents (CSV files)
-----------------------------------------
Each community directory contains the following files, one per trait matrix.
These match exactly the keys expected by the modified initialize_data():

  fb.csv                  -- fungi/bacteria index (1D array)
  Microbes_pp.csv         -- microbial pool before grid placement (C,N,P,E per taxon)
  MinRatios.csv           -- minimum cell stoichiometric quotas
  EnzGenes.csv            -- enzyme gene presence/absence matrix
  OsmoGenes.csv           -- osmolyte gene presence/absence matrix
  HSPGenes.csv            -- HSP gene presence/absence matrix
  UptakeGenes.csv         -- transporter gene presence/absence matrix
  UptakeGenes_trait.csv   -- per-gene transporter cost (scalar per taxon)
  UptakeGenesCost.csv     -- per-gene transporter cost matrix
  UptakeProdCost.csv      -- alias of UptakeGenesCost (matches Data_Dictionary)
  EnzProdConsti_trait.csv -- constitutive enzyme production rate (scalar per taxon)
  EnzProdInduci_trait.csv -- inducible enzyme production rate (scalar per taxon)
  EnzProdConstit.csv      -- constitutive enzyme cost matrix
  EnzProdInduce.csv       -- inducible enzyme cost matrix
  HSPProdConsti_trait.csv -- constitutive HSP production rate (scalar per taxon)
  HSPProdInduci_trait.csv -- inducible HSP production rate (scalar per taxon)
  HSPProdConsti.csv       -- constitutive HSP cost matrix
  HSPProdInduci.csv       -- inducible HSP cost matrix
  OsmoProdConsti_trait.csv -- constitutive osmolyte production rate (scalar per taxon)
  OsmoProdInduci_trait.csv -- inducible osmolyte production rate (scalar per taxon)
  OsmoProdConsti.csv      -- constitutive osmolyte cost matrix
  OsmoProdInduci.csv      -- inducible osmolyte cost matrix
  TaxDroughtTol.csv       -- taxon-specific drought tolerance
  basal_death_prob.csv    -- basal death probability (1D array)
  death_rate.csv          -- death rate coefficient (1D array)
  (taxon IDs are preserved from the original pool — no membership file needed)
"""

import argparse
import os
import sys
import numpy as np
import pandas as pd
from pathlib import Path


# ---------------------------------------------------------------------------
# Runtime / parameter loading  (mirrors generate_trait_pool.py)
# ---------------------------------------------------------------------------

def load_runtime(path):
    rt = pd.read_csv(path, sep='\t', header=None, index_col=0)
    rt.columns = [1]
    numeric_mask = pd.to_numeric(rt[1], errors='coerce').notna()
    rt = rt.loc[numeric_mask].copy()
    rt[1] = pd.to_numeric(rt[1])
    return rt


def override_runtime(runtime, n_taxa):
    rt = runtime.copy()
    rt.loc['n_taxa',   1] = float(n_taxa)
    rt.loc['gridsize', 1] = 1.0
    return rt


# ---------------------------------------------------------------------------
# Full pool initialization
# Runs the complete Microbe initialization for all N taxa in the pool.
# Returns a dict of all trait arrays keyed by their Data_Dictionary names.
# ---------------------------------------------------------------------------

def initialize_pool(runtime, src_dir, input_dir):
    """
    Run the DEMENTpy Microbe initialization for the full taxon pool
    (gridsize=1, no spatial expansion). Returns a dict of all per-taxon
    trait DataFrames / Series / arrays needed to reconstruct any community.
    """
    sys.path.insert(0, src_dir)
    from substrate import Substrate
    from monomer   import Monomer
    from enzyme    import Enzyme
    from microbe   import Microbe

    def csv(name, **kwargs):
        return pd.read_csv(os.path.join(input_dir, name), **kwargs)

    parameters      = csv('parameters.csv',         header=None, index_col=0).astype('float32')
    substrates_init = csv('initial_substrates.csv', header=0,    index_col=0).astype('float32')
    sub_mon_input   = csv('sub_mon_inputs.csv',     header=0,    index_col=0).astype('float32')
    Ea_input        = csv('enzyme_ea.csv',           header=0,    index_col=0).astype('float32')
    ST0_input       = csv('enzyme_S.csv',            header=0,    index_col=0).astype('float32')
    mon_energy      = csv('monomer_energy.csv',      header=0,    index_col=0).astype('float32')

    # Substrate → Monomer → Enzyme (needed for UptakeGenes)
    Substrates = Substrate(runtime, parameters, substrates_init)
    substrates_req_enzyme        = Substrates.substrate_degradation_enzyme()
    substrates_produced_monomers = Substrates.substrate_produced_monomer()
    substrates_index             = Substrates.substrate_index

    Monomers = Monomer(runtime, parameters, mon_energy)
    monomers_initial_pool = Monomers.monomer_initialization(Substrates.Substrates_start)
    monomers_index        = monomers_initial_pool.index

    Enzymes = Enzyme(runtime, parameters, substrates_index, monomers_index)
    enzymes_attributes = Enzymes.enzyme_attributes()

    # Microbe — full trait initialization for N taxa
    Microbes = Microbe(runtime, parameters)

    # Fix: self.fb is drawn as int8 in Microbe.__init__, which overflows when
    # n_taxa > 127 and sum(fb) > 127 (used in np.tile in minimum_cell_quota).
    # Cast to int32 here before any method that sums fb is called.
    Microbes.fb = Microbes.fb.astype('int32')

    microbial_pool   = Microbes.microbial_pool_initialization()   # (Microbes_pp, fb_grid)
    microbial_min_ratios = Microbes.minimum_cell_quota()

    EnzGenes    = Microbes.microbe_enzyme_gene()
    OsmoGenes   = Microbes.microbe_osmolyte_gene()
    HSPGenes    = Microbes.microbe_heat_stress_gene()
    UptakeGenes = Microbes.microbe_uptake_gene(
        substrates_req_enzyme, EnzGenes, substrates_produced_monomers
    )

    enz_prod    = Microbes.microbe_enzproduction_rate(EnzGenes, enzymes_attributes)
    upt_cost    = Microbes.microbe_uptake_cost(UptakeGenes)
    hsp_prod    = Microbes.microbe_heat_stress_protein_production_rate(HSPGenes)
    osmo_prod   = Microbes.microbe_osmoproduction_rate(OsmoGenes)
    drought_tol = Microbes.microbe_drought_tol(osmo_prod[2], osmo_prod[3])
    mortality   = Microbes.microbe_mortality(microbial_pool[1])

    pool = {
        # Fungi/bacteria index and microbial pool
        'fb':                   microbial_pool[1],          # 1D array, length n_taxa
        'Microbes_pp':          microbial_pool[0],          # DataFrame (n_taxa, 4)
        'MinRatios':            microbial_min_ratios,       # DataFrame (n_taxa, 4)
        # Gene matrices
        'EnzGenes':             EnzGenes,
        'OsmoGenes':            OsmoGenes,
        'HSPGenes':             HSPGenes,
        'UptakeGenes':          UptakeGenes,
        # Transporter production
        'UptakeGenes_trait':    upt_cost[0],                # Series (n_taxa,)
        'UptakeGenesCost':      upt_cost[1],                # DataFrame (n_taxa, n_uptake)
        'UptakeProdCost':       upt_cost[1],                # alias
        # Enzyme production
        'EnzProdConsti_trait':  enz_prod[0],                # Series (n_taxa,)
        'EnzProdInduci_trait':  enz_prod[1],                # Series (n_taxa,)
        'EnzProdConstit':       enz_prod[2],                # DataFrame (n_taxa, n_enzymes)
        'EnzProdInduce':        enz_prod[3],                # DataFrame (n_taxa, n_enzymes)
        # HSP production
        'HSPProdConsti_trait':  hsp_prod[0],                # Series (n_taxa,)
        'HSPProdInduci_trait':  hsp_prod[1],                # Series (n_taxa,)
        'HSPProdConsti':        hsp_prod[2],                # DataFrame (n_taxa, n_HSP)
        'HSPProdInduci':        hsp_prod[3],                # DataFrame (n_taxa, n_HSP)
        # Osmolyte production
        'OsmoProdConsti_trait': osmo_prod[0],               # Series (n_taxa,)
        'OsmoProdInduci_trait': osmo_prod[1],               # Series (n_taxa,)
        'OsmoProdConsti':       osmo_prod[2],               # DataFrame (n_taxa, n_osmolytes)
        'OsmoProdInduci':       osmo_prod[3],               # DataFrame (n_taxa, n_osmolytes)
        # Drought tolerance and mortality
        'TaxDroughtTol':        drought_tol,                # Series (n_taxa,)
        'basal_death_prob':     mortality[0],               # 1D array
        'death_rate':           mortality[1],               # 1D array
    }

    # Store taxon index and req_enzyme for use in main()
    pool['_index']       = Microbes.index          # ['Tax1', 'Tax2', ...]
    pool['_req_enzyme']  = substrates_req_enzyme   # enzyme-substrate mapping

    return pool


# ---------------------------------------------------------------------------
# Community subsetting
# Slice all trait arrays in the pool to the selected taxon indices.
# Re-labels taxa Tax1…Taxm for compatibility with the model.
# ---------------------------------------------------------------------------

def subset_community(pool, selected_positions):
    """
    Subset the full pool to the taxa at `selected_positions` (integer positions,
    0-based). Retains original taxon IDs (e.g. Tax47, Tax312) so that results
    can be linked back to the full pool without an intermediate join.

    Returns a dict with the same keys as pool (minus '_index' and '_req_enzyme'), containing
    only the rows/elements corresponding to the selected taxa.
    """
    community = {}

    # Keys prefixed with '_' are metadata, not per-taxon arrays
    _skip = {'_index', '_req_enzyme'}

    for key, obj in pool.items():
        if key in _skip:
            continue

        if isinstance(obj, pd.DataFrame):
            community[key] = obj.iloc[selected_positions].copy()

        elif isinstance(obj, pd.Series):
            community[key] = obj.iloc[selected_positions].copy()

        elif isinstance(obj, np.ndarray):
            # 1D arrays: fb, basal_death_prob, death_rate
            community[key] = obj[selected_positions].copy()

        else:
            community[key] = obj   # scalars / unexpected types

    return community


# ---------------------------------------------------------------------------
# Saving a community
# ---------------------------------------------------------------------------

def save_community(community, out_dir):
    """Save all trait arrays in a community dict to CSVs in out_dir."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    for key, obj in community.items():
        path = out_dir / f"{key}.csv"

        if isinstance(obj, pd.DataFrame):
            obj.to_csv(path)
        elif isinstance(obj, pd.Series):
            obj.to_csv(path, header=True)
        elif isinstance(obj, np.ndarray):
            np.savetxt(path, obj, delimiter=',')
        else:
            # Unexpected type — skip with a note
            pass


# ---------------------------------------------------------------------------
# Progress display
# ---------------------------------------------------------------------------

def progress(current, total, prefix=''):
    pct   = 100 * current / total
    bar   = '#' * int(pct / 2) + '-' * (50 - int(pct / 2))
    print(f"\r{prefix}[{bar}] {current}/{total} ({pct:.0f}%)", end='', flush=True)
    if current == total:
        print()


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Assemble DeCA-DEMENT communities from a pre-initialized taxon pool.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument('--src-dir',        type=str, default='src',
                        help='Directory containing DEMENTpy .py source files')
    parser.add_argument('--input-dir',      type=str, default='input',
                        help='Directory containing input CSVs and runtime.txt')
    parser.add_argument('--output-dir',     type=str, default='community_pool',
                        help='Root output directory for community files')
    parser.add_argument('--runtime',        type=str, default='runtime.txt',
                        help='runtime.txt filename inside --input-dir')
    parser.add_argument('--n-taxa',         type=int, default=None,
                        help='Pool size N (overrides runtime.txt n_taxa)')
    parser.add_argument('--n-communities',  type=int, default=1000,
                        help='Number of ensemble communities to generate')
    parser.add_argument('--community-size', type=int, default=100,
                        help='Number of taxa per ensemble community (m)')
    parser.add_argument('--seed',           type=int, default=None,
                        help='Random seed for reproducibility')
    parser.add_argument('--mode',           type=str, default='both',
                        choices=['single', 'ensemble', 'both'],
                        help="What to generate: 'single' (one taxon per community), "
                             "'ensemble' (m-taxon communities), or 'both'")
    parser.add_argument('--super-community', action='store_true',
                        help='Also save a super-community containing all N taxa '
                             'in <output-dir>/super_community/. Useful for '
                             'trait-space analysis of the full pool. Not included '
                             'in ensemble simulations.')
    args = parser.parse_args()

    # ── Seed ────────────────────────────────────────────────────────────────
    if args.seed is not None:
        np.random.seed(args.seed)
        print(f"Random seed: {args.seed}")

    # ── Path validation ──────────────────────────────────────────────────────
    src_dir    = os.path.abspath(args.src_dir)
    input_dir  = os.path.abspath(args.input_dir)
    output_dir = Path(os.path.abspath(args.output_dir))

    for label, path in [('--src-dir', src_dir), ('--input-dir', input_dir)]:
        if not os.path.isdir(path):
            print(f"ERROR: {label} '{path}' does not exist.")
            sys.exit(1)

    runtime_path = os.path.join(input_dir, args.runtime)
    if not os.path.isfile(runtime_path):
        print(f"ERROR: runtime file '{runtime_path}' not found.")
        sys.exit(1)

    # ── Runtime ─────────────────────────────────────────────────────────────
    runtime = load_runtime(runtime_path)
    n_taxa  = args.n_taxa if args.n_taxa is not None else int(float(runtime.loc['n_taxa', 1]))

    if args.mode in ('ensemble', 'both'):
        if args.community_size >= n_taxa:
            print(f"ERROR: --community-size ({args.community_size}) must be "
                  f"less than pool size ({n_taxa}).")
            sys.exit(1)

    print(f"Pool size (N)      : {n_taxa}")
    if args.mode in ('ensemble', 'both'):
        print(f"Community size (m) : {args.community_size}")
        print(f"N communities      : {args.n_communities}")
    print(f"Mode               : {args.mode}")

    output_dir.mkdir(parents=True, exist_ok=True)

    # ── Initialize pool ──────────────────────────────────────────────────────
    print("\nInitializing taxon pool (running DEMENTpy trait sampling)...")
    runtime_pool = override_runtime(runtime, n_taxa)
    pool = initialize_pool(runtime_pool, src_dir, input_dir)
    print(f"Pool initialized: {n_taxa} taxa")

    # ── Save ReqEnz.csv to output dir root ──────────────────────────────────
    # ReqEnz is a property of the model configuration, not of any individual
    # community — save it once at the pool level so community_summary.py can
    # find it without needing --input-dir or --src-dir.
    req_enz_path = output_dir / 'ReqEnz.csv'
    pool['_req_enzyme'].to_csv(req_enz_path)
    print(f"ReqEnz.csv saved to: {req_enz_path}")

    # ── Single-taxon communities ─────────────────────────────────────────────
    if args.mode in ('single', 'both'):
        single_dir = output_dir / 'single_taxon'
        single_dir.mkdir(parents=True, exist_ok=True)
        print(f"\nGenerating {n_taxa} single-taxon communities → {single_dir}")

        for i in range(n_taxa):
            taxon_id  = pool['_index'][i]
            comm_dir  = single_dir / f"taxon_{taxon_id}"
            community = subset_community(pool, [i])
            save_community(community, comm_dir)
            progress(i + 1, n_taxa, prefix='  Single-taxon: ')

        print(f"  Single-taxon communities saved.")

    # ── Ensemble communities ─────────────────────────────────────────────────
    if args.mode in ('ensemble', 'both'):
        ensemble_dir = output_dir / 'communities'
        ensemble_dir.mkdir(parents=True, exist_ok=True)
        print(f"\nGenerating {args.n_communities} ensemble communities → {ensemble_dir}")

        # Pre-generate all community draws so they are fixed before any
        # simulation is run. Each row is an unordered set of taxon positions.
        all_draws = [
            np.random.choice(n_taxa, size=args.community_size, replace=False)
            for _ in range(args.n_communities)
        ]

        for k, positions in enumerate(all_draws):
            comm_id  = k + 1
            comm_dir = ensemble_dir / f"community_{comm_id:04d}"
            community = subset_community(pool, positions)
            save_community(community, comm_dir)
            progress(k + 1, args.n_communities, prefix='  Ensemble:     ')

        print(f"  Ensemble communities saved.")

    # ── Super-community (all N taxa) ────────────────────────────────────────
    if args.super_community:
        super_dir = output_dir / 'super_community'
        super_dir.mkdir(parents=True, exist_ok=True)
        print(f"\nGenerating super-community ({n_taxa} taxa) → {super_dir}")
        all_positions = list(range(n_taxa))
        super_comm = subset_community(pool, all_positions)
        save_community(super_comm, super_dir)
        print(f"  Super-community saved.")

    # ── Summary ──────────────────────────────────────────────────────────────
    print(f"\nDone. Output written to: {output_dir}")
    if args.mode in ('single', 'both'):
        print(f"  {n_taxa} single-taxon communities in: {output_dir / 'single_taxon'}")
    if args.mode in ('ensemble', 'both'):
        print(f"  {args.n_communities} ensemble communities in: {output_dir / 'communities'}")
    if args.super_community:
        print(f"  Super-community (N={n_taxa}) in: {output_dir / 'super_community'}")


if __name__ == '__main__':
    main()
