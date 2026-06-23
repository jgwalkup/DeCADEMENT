"""
generate_trait_pool.py
======================
Generate a pool of taxa by running the DEMENTpy initialization sequence up to
(but not including) spatial grid expansion. Imports and calls the actual
Substrate, Monomer, Enzyme, and Microbe classes directly, so trait assignments
are always identical to what the model itself would produce.

Expected project layout
-----------------------
    Decoupled-Energy-DEMENTpy/   <- run from here
    ├── src/                     <- --src-dir  (Python source files)
    ├── input/                   <- --input-dir (CSV input files + runtime.txt)
    └── analysis/
        └── generate_trait_pool.py

Usage
-----
    python analysis/generate_trait_pool.py \\
        --src-dir src \\
        --input-dir input \\
        --output-dir taxon_output \\
        --n-taxa 1000 \\
        --seed 42

All paths are resolved relative to the current working directory.

Options
-------
  --src-dir        Directory containing DEMENTpy .py source files (default: src)
  --input-dir      Directory containing input CSVs and runtime.txt (default: input)
  --output-dir     Directory for output files; created if absent (default: taxon_output)
  --n-taxa         Number of taxa to generate (overrides runtime.txt value)
  --runtime        runtime.txt filename, looked up inside --input-dir (default: runtime.txt)
  --seed           Random seed for reproducibility (default: None)
  --community-size Community size for coverage diagnostics (default: 100)

Output
------
  <output-dir>/trait_pool.csv   — one row per taxon, columns described below

Output columns
--------------
  taxon                : taxon label (Tax1 ... TaxN)

  -- Enzyme traits --
  n_enz_genes          : number of enzyme genes held by taxon
  enz_consti_rate      : per-gene constitutive enzyme production rate (LHS draw)
  enz_induci_rate      : per-gene inducible enzyme production rate (LHS draw)
  total_enz_consti     : n_enz_genes x enz_consti_rate
  total_enz_induci     : n_enz_genes x enz_induci_rate

  -- Transporter traits --
  n_uptake_genes       : number of monomer transporter genes held by taxon
  uptake_c_cost        : per-transporter C cost (LHS draw)
  total_uptake_cost    : n_uptake_genes x uptake_c_cost

  -- HSP traits --
  HSP_gene             : 0 or 1 (presence/absence of HSP gene)
  HSP_consti_rate      : constitutive HSP production rate (0 if HSP_gene=0)
  HSP_induci_rate      : inducible HSP production rate (0 if HSP_gene=0)

  -- Osmolyte traits --
  n_osmo_genes         : number of osmolyte genes held by taxon
  osmo_consti_rate     : per-gene constitutive osmolyte production rate (LHS draw)
  osmo_induci_rate     : per-gene inducible osmolyte production rate (LHS draw)
  total_osmo_consti    : n_osmo_genes x osmo_consti_rate
  total_osmo_induci    : n_osmo_genes x osmo_induci_rate
"""

import argparse
import sys
import os
import numpy as np
import pandas as pd


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def load_runtime(path):
    """
    Load runtime.txt as a single-column DataFrame indexed by parameter name.
    Non-numeric rows (e.g. output_format, kinetics_model) are dropped before
    the DataFrame is passed to DEMENTpy classes, mirroring what dementpy.py does.
    """
    rt = pd.read_csv(path, sep='\t', header=None, index_col=0)
    rt.columns = [1]
    numeric_mask = pd.to_numeric(rt[1], errors='coerce').notna()
    rt = rt.loc[numeric_mask].copy()
    rt[1] = pd.to_numeric(rt[1])
    return rt


def override_runtime(runtime, n_taxa):
    """Return a copy of runtime with n_taxa overridden and gridsize forced to 1."""
    rt = runtime.copy()
    rt.loc['n_taxa',   1] = float(n_taxa)
    rt.loc['gridsize', 1] = 1.0
    return rt


# ---------------------------------------------------------------------------
# Coverage diagnostics
# ---------------------------------------------------------------------------

def coverage_report(trait_pool, n_community=100, min_appearances=10):
    """
    Print trait-space coverage summary and community-sampling coverage table.

    Coverage is assessed using a Binomial model: each taxon has probability
    p = m/N of being drawn into any single community of size m from a pool
    of N. Across k communities, appearances per taxon ~ Binomial(k, p).

    The diagnostic reports:
      - E[appearances]: expected times each taxon appears across k communities
      - P(taxon < t):   probability any single taxon falls below the threshold
      - P(all >= t):    probability every taxon in the pool meets the threshold
                        = (1 - P(taxon < t))^N

    P(all >= t) is the most useful summary: it is the probability that your
    least-sampled taxon still meets the minimum appearance threshold.
    """
    from scipy.stats import binom as binom_dist

    n_taxa        = len(trait_pool)
    producers     = int(trait_pool['HSP_gene'].sum())
    non_producers = n_taxa - producers

    print("\n=== Trait Pool Summary ===")
    print(f"  Total taxa generated : {n_taxa}")
    print(f"  HSP producers        : {producers}  ({100*producers/n_taxa:.1f}%)")
    print(f"  HSP non-producers    : {non_producers}  ({100*non_producers/n_taxa:.1f}%)")

    print("\n  Acquisition axis (total_enz_induci):")
    for grp, label in [
        (trait_pool,                           'All taxa'),
        (trait_pool[trait_pool.HSP_gene == 0], 'HSP=0'),
        (trait_pool[trait_pool.HSP_gene == 1], 'HSP=1'),
    ]:
        if len(grp) == 0:
            continue
        q = grp['total_enz_induci'].quantile([0.25, 0.5, 0.75]).values
        print(f"    {label:18s}  Q1={q[0]:.5f}  median={q[1]:.5f}  Q3={q[2]:.5f}")

    prod = trait_pool[trait_pool.HSP_gene == 1]
    if len(prod) > 0:
        print("\n  HSP production intensity (producers only):")
        for col, label in [('HSP_consti_rate', 'Constitutive'),
                            ('HSP_induci_rate', 'Inducible')]:
            q = prod[col].quantile([0.25, 0.5, 0.75]).values
            print(f"    {label:18s}  Q1={q[0]:.5f}  median={q[1]:.5f}  Q3={q[2]:.5f}")

    # Binomial coverage model
    p = n_community / n_taxa   # P(taxon drawn into one community)

    print(f"\n=== Community Sampling Coverage ===")
    print(f"  Pool size (N)              : {n_taxa}")
    print(f"  Community size (m)         : {n_community}")
    print(f"  Min. appearances threshold : {min_appearances}")
    print(f"  P(drawn per community)     : {p:.4f}")
    print()
    print(f"  Each taxon appears Binomial(k, {p:.4f}) times across k communities.")
    print(f"  P(all >= t) is the probability every taxon meets the threshold —")
    print(f"  the key diagnostic for adequate individual-taxon coverage.")
    print()
    print(f"  {'k':>6}  {'E[appearances]':>16}  {'P(taxon < t)':>14}  {'P(all taxa >= t)':>18}")
    print(f"  {'':-<6}  {'':-<16}  {'':-<14}  {'':-<18}")
    for k in [10, 25, 50, 100, 250, 500, 1000]:
        p_under  = binom_dist.cdf(min_appearances - 1, k, p)
        p_all_ok = (1 - p_under) ** n_taxa
        e_appear = k * p
        flag = "  <-- 95% threshold" if p_all_ok >= 0.95 and (k == 10 or
               (1 - binom_dist.cdf(min_appearances - 1, k - 1, p)) ** n_taxa < 0.95) else ""
        print(f"  {k:>6}  {e_appear:>16.1f}  {p_under:>14.6f}  {p_all_ok:>18.6f}{flag}")
    print()


# ---------------------------------------------------------------------------
# Main generation logic
# ---------------------------------------------------------------------------

def generate_trait_pool(runtime, src_dir, input_dir):
    """
    Run the DEMENTpy initialization sequence for a single-cell grid and return
    a per-taxon trait summary DataFrame.

    Parameters
    ----------
    runtime   : DataFrame — runtime parameters (numeric rows only, gridsize=1)
    src_dir   : str — absolute path to directory containing DEMENTpy .py files
    input_dir : str — absolute path to directory containing input CSV files
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

    # Substrate
    Substrates = Substrate(runtime, parameters, substrates_init)
    substrates_req_enzyme        = Substrates.substrate_degradation_enzyme()
    substrates_produced_monomers = Substrates.substrate_produced_monomer()
    substrates_index             = Substrates.substrate_index

    # Monomer
    Monomers = Monomer(runtime, parameters, mon_energy)
    monomers_initial_pool = Monomers.monomer_initialization(Substrates.Substrates_start)
    monomers_index        = monomers_initial_pool.index

    # Enzyme
    Enzymes = Enzyme(runtime, parameters, substrates_index, monomers_index)
    enzymes_attributes = Enzymes.enzyme_attributes()

    # Microbe
    Microbes = Microbe(runtime, parameters)
    index    = Microbes.index

    EnzGenes    = Microbes.microbe_enzyme_gene()
    OsmoGenes   = Microbes.microbe_osmolyte_gene()
    HSPGenes    = Microbes.microbe_heat_stress_gene()
    UptakeGenes = Microbes.microbe_uptake_gene(
        substrates_req_enzyme, EnzGenes, substrates_produced_monomers
    )

    enz_consti_s, enz_induci_s, _, _ = \
        Microbes.microbe_enzproduction_rate(EnzGenes, enzymes_attributes)

    upt_cost_s, _ = \
        Microbes.microbe_uptake_cost(UptakeGenes)

    hsp_consti_s, hsp_induci_s, _, _ = \
        Microbes.microbe_heat_stress_protein_production_rate(HSPGenes)

    osmo_consti_s, osmo_induci_s, _, _ = \
        Microbes.microbe_osmoproduction_rate(OsmoGenes)

    # Aggregate gene counts
    n_enz_genes  = EnzGenes.sum(axis=1).values
    n_upt_genes  = UptakeGenes.sum(axis=1).values
    n_osmo_genes = OsmoGenes.sum(axis=1).values
    hsp_gene     = HSPGenes.sum(axis=1).values   # 0 or 1 when n_HSP=1

    trait_pool = pd.DataFrame({
        'taxon':             index,
        # Enzyme
        'n_enz_genes':       n_enz_genes.astype(int),
        'enz_consti_rate':   enz_consti_s.values,
        'enz_induci_rate':   enz_induci_s.values,
        'total_enz_consti':  (n_enz_genes * enz_consti_s.values),
        'total_enz_induci':  (n_enz_genes * enz_induci_s.values),
        # Transporter
        'n_uptake_genes':    n_upt_genes.astype(int),
        'uptake_c_cost':     upt_cost_s.values,
        'total_uptake_cost': (n_upt_genes * upt_cost_s.values),
        # HSP
        'HSP_gene':          hsp_gene.astype(int),
        'HSP_consti_rate':   hsp_consti_s.values,
        'HSP_induci_rate':   hsp_induci_s.values,
        # Osmolyte
        'n_osmo_genes':      n_osmo_genes.astype(int),
        'osmo_consti_rate':  osmo_consti_s.values,
        'osmo_induci_rate':  osmo_induci_s.values,
        'total_osmo_consti': (n_osmo_genes * osmo_consti_s.values),
        'total_osmo_induci': (n_osmo_genes * osmo_induci_s.values),
    })

    # Zero HSP rates for non-producers — gene absent means rates are never
    # expressed; mask here so downstream analysis is not misled by LHS draws
    trait_pool.loc[trait_pool.HSP_gene == 0,
                   ['HSP_consti_rate', 'HSP_induci_rate']] = 0.0

    return trait_pool


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Generate a fixed trait pool for DeCA-DEMENT guild analysis.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument('--src-dir',       type=str, default='src',
                        help='Directory containing DEMENTpy .py source files')
    parser.add_argument('--input-dir',     type=str, default='input',
                        help='Directory containing input CSVs and runtime.txt')
    parser.add_argument('--output-dir',    type=str, default='taxon_output',
                        help='Directory for output files (created if absent)')
    parser.add_argument('--runtime',       type=str, default='runtime.txt',
                        help='runtime.txt filename, looked up inside --input-dir')
    parser.add_argument('--n-taxa',        type=int, default=None,
                        help='Number of taxa to generate (overrides runtime.txt)')
    parser.add_argument('--seed',          type=int, default=None,
                        help='Random seed for reproducibility')
    parser.add_argument('--community-size',  type=int, default=100,
                        help='Community size for coverage diagnostics (default: 100)')
    parser.add_argument('--min-appearances', type=int, default=10,
                        help='Minimum times each taxon must appear across communities (default: 10)')
    args = parser.parse_args()

    if args.seed is not None:
        np.random.seed(args.seed)
        print(f"Random seed: {args.seed}")

    # Resolve all paths relative to cwd
    src_dir    = os.path.abspath(args.src_dir)
    input_dir  = os.path.abspath(args.input_dir)
    output_dir = os.path.abspath(args.output_dir)
    os.makedirs(output_dir, exist_ok=True)

    # Validate directories
    for label, path in [('--src-dir', src_dir), ('--input-dir', input_dir)]:
        if not os.path.isdir(path):
            print(f"ERROR: {label} '{path}' does not exist or is not a directory.")
            sys.exit(1)

    runtime_path = os.path.join(input_dir, args.runtime)
    if not os.path.isfile(runtime_path):
        print(f"ERROR: runtime file '{runtime_path}' not found.")
        sys.exit(1)

    runtime = load_runtime(runtime_path)
    n_taxa  = args.n_taxa if args.n_taxa is not None else int(float(runtime.loc['n_taxa', 1]))
    print(f"Generating trait pool: {n_taxa} taxa")

    runtime_pool = override_runtime(runtime, n_taxa)
    trait_pool   = generate_trait_pool(runtime_pool, src_dir, input_dir)

    out_path = os.path.join(output_dir, 'trait_pool.csv')
    trait_pool.to_csv(out_path, index=False)
    print(f"Trait pool saved to: {out_path}")

    coverage_report(trait_pool, n_community=args.community_size,
                    min_appearances=args.min_appearances)


if __name__ == '__main__':
    main()
