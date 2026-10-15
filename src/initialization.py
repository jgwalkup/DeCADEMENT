"""
This module, with only one function "initialized_data(), initializes data related to
substrate, monomer, enzyme, and microbe, as well as their distribution on the spatial grid,
preceding the actual decompostion-related computations.
"""

import pandas as pd
import numpy as np
import xarray as xr

# new
import warnings
import numbers
from pathlib import Path
from typing import Union  #  Use Union for Python < 3.10 compatibility
# to here

from substrate import Substrate
from monomer import Monomer
from enzyme import Enzyme
from microbe import Microbe
from utility import expand


def _optional(parameters, name, default):
    """Value of an optional parameters.csv row; prints a notice when the row is absent and the default is used."""
    if name in parameters.index:
        return float(parameters.loc[name, 1])
    print(f"  parameters.csv has no '{name}' row -- using the default {default}")
    return float(default)


def initialize_data(runtime_parameters, climate_file="climate.csv", community_dir=None):
    """
    Parameters:
        runtime_parameters: user-specified parameters setting up the system;
                            all other paras loaded by reading the parameters.csv
        climate_file:       path to climate forcing CSV (default: 'climate.csv')
        community_dir:      optional path to a pre-assembled community directory
                            produced by assemble_communities.py. When provided,
                            microbial trait arrays are loaded from this directory
                            instead of being drawn randomly, and n_taxa in
                            runtime_parameters is updated to match the community.
                            When None (default), behaviour is identical to the
                            original model: taxa are drawn randomly as before.
    Return:
        Data_Dictionary: a dictionary of all variables that feeds the grid.py module
    """

    # Load all input files
    parameters = pd.read_csv("parameters.csv", header=None, index_col=0).astype(
        "float32"
    )  # parameters
    substrates_init = pd.read_csv(
        "initial_substrates.csv", header=0, index_col=0
    ).astype("float32")  # initial substrates
    sub_mon_input = pd.read_csv("sub_mon_inputs.csv", header=0, index_col=0).astype(
        "float32"
    )  # inputs of substrates and monomers
    Ea_input = pd.read_csv("enzyme_ea.csv", header=0, index_col=0).astype(
        "float32"
    )  # enzyme activation energy
    mon_energy = pd.read_csv("monomer_energy.csv", header=0, index_col=0).astype(
        "float32"
    )  # enzyme activation energy
    # HT0_input        = pd.read_csv("enzyme_H.csv",          header=0,    index_col=0).astype('float32')   # enzyme enthalpy
    ST0_input = pd.read_csv("enzyme_S.csv", header=0, index_col=0).astype(
        "float32"
    )  # enzyme enthalpy
    climate = pd.read_csv(climate_file, header=0, index_col=0)  # climate forcings

    # daily temperature and water potential
    daily_temp = climate["Temp"].astype("float32")  # temperaure series
    daily_psi = climate["Psi"].astype("float32")  # water potential series

    # ...an instance of Substrate class
    Substrates = Substrate(runtime_parameters, parameters, substrates_init)
    # ...substrate initial pool size
    substrates_initial_pool = Substrates.Substrates_start
    # ...substrate input rate
    substrates_input_rate = Substrates.substrate_input(sub_mon_input)
    # ...substrates-produced monomers
    substrates_produced_monomers = Substrates.substrate_produced_monomer()
    # ...substrates degradation required enzymes
    substrates_req_enzyme = Substrates.substrate_degradation_enzyme()
    # ...substrates index
    substrates_index = Substrates.substrate_index

    # ...an instance of Monomer class
    Monomers = Monomer(runtime_parameters, parameters, mon_energy)
    # ...monomers initial pool size
    monomers_initial_pool = Monomers.monomer_initialization(substrates_initial_pool)
    # ...Monomer energy attributes
    monomer_energy_characteristics = Monomers.Monomer_energy
    # ...initial monomer ratios
    monomer_ratio_initial = Monomers.monomer_ratios(monomers_initial_pool)
    # ...monomers input rate
    monomers_input_rate = Monomers.monomer_input_rate(sub_mon_input)
    # ...monomers uptake required enzymes
    monomers_uptake_reqenzyme = Monomers.monomer_uptake_reqenzyme()
    # ...monomer index
    monomers_index = monomers_initial_pool.index

    # ...an instance of Enzyme class
    Enzymes = Enzyme(runtime_parameters, parameters, substrates_index, monomers_index)
    # ...enzyme initial pool size:0
    enzymes_initial_pool = Enzymes.enzyme_pool_initialization()
    # ...enzyme attributes
    enzymes_attributes = Enzymes.enzyme_attributes()
    # ...enzymes of substrate degradation Ea
    enzymes_Ea = Enzymes.enzyme_Ea(Ea_input)
    # ...enzymes of substrate degradation HT0 for calculating Vmax
    # enzymes_Vmax_HT0 = Enzymes.enzyme_HT0(HT0_input).T
    enzymes_Vmax_HT0 = Enzymes.enzyme_HT0().T
    # ...enzymes of substrate degradation ST0 for calculating Vmax
    # enzymes_Vmax_ST0 = Enzymes.enzyme_ST0(enzymes_Vmax_HT0)
    enzymes_Vmax_ST0 = Enzymes.enzyme_ST0(ST0_input).T
    # ...enzymes of substrate degradation Cp for calculating Vmax
    enzymes_Vmax_Cp = Enzymes.enzyme_Cp().T
    # ...enzymes of substrate degradation HT0 for calculating Km
    # enzymes_Km_HT0 = Enzymes.enzyme_HT0(HT0_input)
    enzymes_Km_HT0 = Enzymes.enzyme_HT0()
    # ...enzymes of substrate degradation ST0 for calculating Km
    # enzymes_Km_ST0 = Enzymes.enzyme_ST0(enzymes_Km_HT0)
    enzymes_Km_ST0 = Enzymes.enzyme_ST0(ST0_input)
    # ...enzymes of substrate degradation Cp for calculating Km
    enzymes_Km_Cp = Enzymes.enzyme_Cp()
    # ...monomers uptake enzyme Ea
    enzymes_uptake_Ea = Enzymes.enzyme_uptake_Ea()
    # ...enzymes of substrate degradation Vmax
    enzymes_Vmax, enzymes_Vmax_T = Enzymes.enzyme_Vmax(substrates_req_enzyme)
    # ...enzymes of substrate degradation substrate specificity-efficiency tradeoff
    enzymes_sub_spec, enzymes_sub_spec_T = Enzymes.enzyme_Vmax_sub_spec(
        substrates_req_enzyme
    )
    # ...monomers uptake enzyme Vmax
    enzymes_uptake_Vmax = Enzymes.enzyme_uptake_Vmax(monomers_uptake_reqenzyme)
    # ...enzymes of substrate degradation Km
    enzymes_Km = Enzymes.enzyme_Km(enzymes_Vmax)
    # ...monomers uptake enzyme Km
    enzymes_uptake_Km = Enzymes.enzyme_uptake_Km(enzymes_uptake_Vmax)
    # ...monomers of substrate degradation Cp
    enzymes_uptake_Cp = Enzymes.enzyme_uptake_Cp()
    # ...monomers of substrate degradation ST0
    enzymes_uptake_HT0 = Enzymes.enzyme_uptake_HT0()
    # ...monomers of substrate degradation ST0
    enzymes_uptake_ST0 = Enzymes.enzyme_uptake_ST0()
    # ...monomers of substrate degradation substrate specificity-efficiency tradeoff
    enzymes_uptake_sub_spec = Enzymes.uptake_sub_spec(monomers_uptake_reqenzyme)

    # ── Microbial community ──────────────────────────────────────────────────
    # When community_dir is provided, load pre-assembled trait arrays from that
    # directory (produced by assemble_communities.py) and override n_taxa in
    # runtime_parameters to match the loaded community size. Grid placement
    # (microbial_community_initialization) still runs as normal.
    # When community_dir is None, behaviour is identical to the original model.

    if community_dir is not None:
        # ── Load pre-assembled community ─────────────────────────────────────
        cdir = Path(community_dir)
        if not cdir.is_dir():
            raise FileNotFoundError(
                f"community_dir '{community_dir}' does not exist or is not a directory."
            )

        def load_csv(name, **kwargs):
            p = cdir / f"{name}.csv"
            if not p.exists():
                raise FileNotFoundError(
                    f"Expected community file '{p}' not found. "
                    f"Was this community generated by assemble_communities.py?"
                )
            return pd.read_csv(p, **kwargs)

        # Infer n_taxa and taxon index from EnzGenes (always present)
        _EnzGenes = load_csv("EnzGenes", header=0, index_col=0)
        n_taxa_community = len(_EnzGenes)
        taxon_index = list(_EnzGenes.index)

        # Override n_taxa in runtime so downstream code uses the right count
        runtime_parameters = runtime_parameters.copy()
        runtime_parameters.loc["n_taxa", 1] = float(n_taxa_community)
        print(f"Loading pre-assembled community from: {cdir}")
        print(f"  Community size: {n_taxa_community} taxa")

        # ...an instance of Microbe class (needed for grid placement methods)
        Microbes = Microbe(runtime_parameters, parameters)

        # Load per-taxon arrays from CSV files
        _fb_taxon = np.loadtxt(cdir / "fb.csv", delimiter=",").astype("int8")
        _microbes_pp = load_csv("Microbes_pp", header=0, index_col=0).astype("float32")

        # Grid-expand fb and Microbes_pp — microbial_pool_initialization()
        # normally does this via np.tile before passing to
        # microbial_community_initialization(). We must replicate that here.
        gridsize = int(runtime_parameters.loc["gridsize", 1])
        _fb_grid = np.tile(_fb_taxon, gridsize)  # (n_taxa,) -> (n_taxa*gridsize,)

        # Expand Microbes_pp: tile rows to fill the full spatial grid,
        # preserving the per-taxon index repeated gridsize times
        _microbes_pp_grid = pd.DataFrame(
            np.tile(_microbes_pp.values, (gridsize, 1)),
            index=list(_microbes_pp.index) * gridsize,
            columns=_microbes_pp.columns,
            dtype="float32",
        )
        microbial_pool = (_microbes_pp_grid, _fb_grid)

        # ...Microbial community initialization (grid placement — unchanged)
        microbial_community = Microbes.microbial_community_initialization(
            microbial_pool[0], microbial_pool[1]
        )

        # ...Microbial minimum ratios
        microbial_min_ratios = load_csv("MinRatios", header=0, index_col=0).astype(
            "float32"
        )

        # ...Microbial enzyme genes
        microbial_enzyme_gene = _EnzGenes.astype("int8")

        # ...Microbial osmolyte genes
        microbial_osmolyte_gene = load_csv("OsmoGenes", header=0, index_col=0).astype(
            "int8"
        )

        # ...Microbial heat stress genes
        microbial_heat_stress_gene = load_csv("HSPGenes", header=0, index_col=0).astype(
            "int8"
        )

        # ...Microbial uptake genes
        microbial_uptake_gene = load_csv("UptakeGenes", header=0, index_col=0).astype(
            "int8"
        )

        # ...Microbial uptake cost — tuple(series, dataframe) to match Microbe method output
        _upt_trait = (
            load_csv("UptakeGenes_trait", header=0, index_col=0)
            .squeeze()
            .astype("float32")
        )
        _upt_cost = load_csv("UptakeGenesCost", header=0, index_col=0).astype("float32")
        microbial_uptake_cost = (_upt_trait, _upt_cost)

        # ...Microbial enzyme production rate — tuple(consti_s, induci_s, consti_df, induci_df)
        _enz_consti_s = (
            load_csv("EnzProdConsti_trait", header=0, index_col=0)
            .squeeze()
            .astype("float32")
        )
        _enz_induci_s = (
            load_csv("EnzProdInduci_trait", header=0, index_col=0)
            .squeeze()
            .astype("float32")
        )
        _enz_consti_df = load_csv("EnzProdConstit", header=0, index_col=0).astype(
            "float32"
        )
        _enz_induci_df = load_csv("EnzProdInduce", header=0, index_col=0).astype(
            "float32"
        )
        microbial_enzyme_prod_rate = (
            _enz_consti_s,
            _enz_induci_s,
            _enz_consti_df,
            _enz_induci_df,
        )

        # ...Microbial heat stress protein initial concentration
        microbial_heat_stress_protein_init = Microbes.HSP_initialization(
            microbial_community[0], microbial_heat_stress_gene
        )

        # ...Microbial heat stress protein production rate — tuple(consti_s, induci_s, consti_df, induci_df)
        _hsp_consti_s = (
            load_csv("HSPProdConsti_trait", header=0, index_col=0)
            .squeeze()
            .astype("float32")
        )
        _hsp_induci_s = (
            load_csv("HSPProdInduci_trait", header=0, index_col=0)
            .squeeze()
            .astype("float32")
        )
        _hsp_consti_df = load_csv("HSPProdConsti", header=0, index_col=0).astype(
            "float32"
        )
        _hsp_induci_df = load_csv("HSPProdInduce", header=0, index_col=0).astype(
            "float32"
        )
        microbial_heat_stress_protein_prod_rate = (
            _hsp_consti_s,
            _hsp_induci_s,
            _hsp_consti_df,
            _hsp_induci_df,
        )

        # ...Microbial osmolyte production rate — tuple(consti_s, induci_s, consti_df, induci_df)
        _osmo_consti_s = (
            load_csv("OsmoProdConsti_trait", header=0, index_col=0)
            .squeeze()
            .astype("float32")
        )
        _osmo_induci_s = (
            load_csv("OsmoProdInduci_trait", header=0, index_col=0)
            .squeeze()
            .astype("float32")
        )
        _osmo_consti_df = load_csv("OsmoProdConsti", header=0, index_col=0).astype(
            "float32"
        )
        _osmo_induci_df = load_csv("OsmoProdInduci", header=0, index_col=0).astype(
            "float32"
        )
        microbial_osmolyte_prod_rate = (
            _osmo_consti_s,
            _osmo_induci_s,
            _osmo_consti_df,
            _osmo_induci_df,
        )

        # ...Microbial drought tolerance
        microbial_drought_tol = (
            load_csv("TaxDroughtTol", header=0, index_col=0).squeeze().astype("float32")
        )

        # ...Microbial mortality — must use grid-expanded _fb_grid so outputs
        #   are length n_taxa*gridsize, matching what grid.py expects
        microbial_mortality = Microbes.microbe_mortality(_fb_grid)

    else:
        # ── Original random initialization (unchanged) ────────────────────
        # ...an instance of Microbe class
        Microbes = Microbe(runtime_parameters, parameters)
        # ...Microbial pool initialization
        microbial_pool = Microbes.microbial_pool_initialization()
        # ...Microbial community initialization
        microbial_community = Microbes.microbial_community_initialization(
            microbial_pool[0], microbial_pool[1]
        )
        # ...Microbial minimum ratios
        microbial_min_ratios = Microbes.minimum_cell_quota()
        # ...Microbial enzyme genes
        microbial_enzyme_gene = Microbes.microbe_enzyme_gene()
        # ...Microbial osmolyte genes
        microbial_osmolyte_gene = Microbes.microbe_osmolyte_gene()
        # ...Microbial heat stress genes
        microbial_heat_stress_gene = Microbes.microbe_heat_stress_gene()
        # ...Microbial uptake genes
        microbial_uptake_gene = Microbes.microbe_uptake_gene(
            substrates_req_enzyme, microbial_enzyme_gene, substrates_produced_monomers
        )
        # ...Microbial uptake cost
        microbial_uptake_cost = Microbes.microbe_uptake_cost(microbial_uptake_gene)
        # ...Microbial enzyme production rate
        microbial_enzyme_prod_rate = Microbes.microbe_enzproduction_rate(
            microbial_enzyme_gene, enzymes_attributes
        )
        # ...Microbial heat stress protein production rate
        microbial_heat_stress_protein_init = Microbes.HSP_initialization(
            microbial_community[0], microbial_heat_stress_gene
        )
        # ...Microbial heat stress protein production rate
        microbial_heat_stress_protein_prod_rate = (
            Microbes.microbe_heat_stress_protein_production_rate(
                microbial_heat_stress_gene
            )
        )
        # ...Microbial osmolyte production rate
        microbial_osmolyte_prod_rate = Microbes.microbe_osmoproduction_rate(
            microbial_osmolyte_gene
        )
        # ...Microbial drought tolerance
        microbial_drought_tol = Microbes.microbe_drought_tol(
            microbial_osmolyte_prod_rate[2], microbial_osmolyte_prod_rate[3]
        )
        # ...Microbial mortality
        microbial_mortality = Microbes.microbe_mortality(microbial_pool[1])
        # ...Microbial respiration thermal response parameters
        # microbial_thermal_adaptation = Microbes.microbe_thermal_adaptation_gene()

    # ...Dump all initialized data into a dictionary; NOTE: variables with expand() put on the spatial grid
    gridsize = int(runtime_parameters.loc["gridsize", 1])

    Data_Dictionary = {
        "Substrates": expand(substrates_initial_pool, gridsize),
        "Substrate_index": substrates_index,
        "SubInput": expand(substrates_input_rate, gridsize),
        "ReqEnz": substrates_req_enzyme,
        "MonomersProduced": substrates_produced_monomers,
        "Monomers": expand(monomers_initial_pool, gridsize),
        "Monomer_index": monomers_index,
        "Monomer_ratio": expand(monomer_ratio_initial, gridsize),
        "MonInput": expand(monomers_input_rate, gridsize),
        "MonEnergy_Attrib": monomer_energy_characteristics,
        "Uptake_ReqEnz": expand(monomers_uptake_reqenzyme, gridsize),
        "Enzymes": expand(enzymes_initial_pool, gridsize),
        "Enz_km_sub_spec": expand(enzymes_sub_spec, gridsize),
        "Enz_vmax_sub_spec": expand(enzymes_sub_spec_T, gridsize),
        "Uptake_Vmax_sub_spec": expand(enzymes_uptake_sub_spec, gridsize),
        "Km0": expand(enzymes_Km, gridsize),  # enzyme half-saturation constant
        "Uptake_Km0": expand(
            enzymes_uptake_Km, gridsize
        ),  # transporter half-saturation constant
        "Uptake_Ea": expand(
            enzymes_uptake_Ea, gridsize
        ),  # transporter acitivation energy
        "Uptake_Vmax0": expand(
            enzymes_uptake_Vmax, gridsize
        ),  # transporter reaction rate
        "Ea": expand(enzymes_Ea, gridsize),  # enzyme activation energy
        "Vmax_HT0": expand(enzymes_Vmax_HT0, gridsize),  # enzyme enthalpy
        "Vmax_ST0": expand(enzymes_Vmax_ST0, gridsize),  # enzyme entropy
        "Vmax_Cp": expand(enzymes_Vmax_Cp, gridsize),  # enzyme heat capacity
        "Km_HT0": expand(enzymes_Km_HT0, gridsize),  # enzyme enthalpy
        "Km_ST0": expand(enzymes_Km_ST0, gridsize),  # enzyme entropy
        "Km_Cp": expand(enzymes_Km_Cp, gridsize),  # enzyme heat capacity
        "Uptake_HT0": expand(enzymes_uptake_HT0, gridsize),  #
        "Uptake_ST0": expand(enzymes_uptake_ST0, gridsize),  #
        "Uptake_Cp": expand(enzymes_uptake_Cp, gridsize),  #
        "Vmax0": expand(enzymes_Vmax_T, gridsize),  # enzyme reaction rate
        "EnzAttrib": enzymes_attributes,  # enzyme stoichiometry and energy cost
        "Microbes_pp": microbial_pool[0],  # tuple[0]: microbes preceding placement
        "Microbes": microbial_community[0],  # tuple[1]: initialized spatial microbes
        "fb": microbial_pool[1],  # tuple[0]: fungi index
        "Bac_density": microbial_community[1],  # tuple[1]: bacterial density
        "Fun_density": microbial_community[2],  # tuple[2]: fungi density
        "MinRatios": expand(
            microbial_min_ratios, gridsize
        ),  # microbial cell min. ratios
        "UptakeGenes": expand(
            microbial_uptake_gene, gridsize
        ),  # transporter gene distribution across taxa
        "OsmoGenes": expand(
            microbial_osmolyte_gene, gridsize
        ),  # osmolyte gene distribution across taxa
        "HSPGenes": expand(
            microbial_heat_stress_gene, gridsize
        ),  # heat stress protein gene distribution across taxa
        "EnzGenes": expand(
            microbial_enzyme_gene, gridsize
        ),  # enzyme gene distribution across taxa
        "Initial_HSP_Prop": microbial_heat_stress_protein_init,  # initial heat stress proteins
        "UptakeGenes_trait": expand(
            microbial_uptake_cost[0], gridsize
        ),  # single gene cost of transporter
        "OsmoProdConsti_trait": expand(
            microbial_osmolyte_prod_rate[0], gridsize
        ),  # single gene cost of constitutive osmolyte
        "OsmoProdInduci_trait": expand(
            microbial_osmolyte_prod_rate[1], gridsize
        ),  # single gene cost of inducible osmolyte
        "EnzProdConsti_trait": expand(
            microbial_enzyme_prod_rate[0], gridsize
        ),  # single gene cost of constitutive enzyme
        "EnzProdInduci_trait": expand(
            microbial_enzyme_prod_rate[1], gridsize
        ),  # single gene cost of inducible enzyme
        "HSPProdConsti_trait": expand(
            microbial_heat_stress_protein_prod_rate[0], gridsize
        ),  # single gene cost of constitutive heat stress protein
        "HSPProdInduci_trait": expand(
            microbial_heat_stress_protein_prod_rate[1], gridsize
        ),  # single gene cost of inducible heat stress protein
        "UptakeGenesCost": expand(
            microbial_uptake_cost[1], gridsize
        ),  # distribution of transporter gene cost across taxa
        "UptakeProdCost": expand(
            microbial_uptake_cost[1], gridsize
        ),  # transporter genes by taxon
        "OsmoProdConsti": expand(
            microbial_osmolyte_prod_rate[2], gridsize
        ),  # distribution of consti. osmolyte gene cost across taxa
        "OsmoProdInduci": expand(
            microbial_osmolyte_prod_rate[3], gridsize
        ),  # distribution of induci. osmolyte gene cost across taxa
        "EnzProdConstit": expand(
            microbial_enzyme_prod_rate[2], gridsize
        ),  # distribution of consti. enzyme gene cost across taxa
        "EnzProdInduce": expand(
            microbial_enzyme_prod_rate[3], gridsize
        ),  # distribution of induci. enzyme gene cost across taxa
        "HSPProdConsti": expand(
            microbial_heat_stress_protein_prod_rate[2], gridsize
        ),  # distribution of consti. HSP gene cost across taxa
        "HSPProdInduci": expand(
            microbial_heat_stress_protein_prod_rate[3], gridsize
        ),  # distribution of induci. HSP gene cost across taxa
        # "TaxHSPinit":         expand(microbial_heat_stress_protein_prod_rate[4],gridsize), # initial HSP concentrations
        "TaxDroughtTol": expand(
            microbial_drought_tol, gridsize
        ),  # distribution of taxon-specific drought tol.
        "basal_death_prob": microbial_mortality[0],  # basal death probability
        "death_rate": microbial_mortality[
            1
        ],  # change rate of death prob. agaist mositure
        "Uptake_prod_E_cost": parameters.loc[
            "Uptake_Prod_E_cost", 1
        ],  # constant of transporter maintenance cost
        "ActiveUptakeE_cost": parameters.loc[
            "Active_Uptake_E_cost", 1
        ],  # E cost of active monomer transport
        # 'Cell_Maint_cost': parameters.loc['Cell_maint_C_cost',1], # constant of cell maintenance cost
        "OsmoProd_E_cost": parameters.loc[
            "Osmo_prod_E_cost", 1
        ],  # constant for calculating E cost of osmo prod.
        "OsmoProd_N_cost": parameters.loc[
            "Osmo_prod_N_cost", 1
        ],  # constant for calculating N cost of osmo prod.
        "HSPProd_E_cost": parameters.loc[
            "HSP_prod_E_cost", 1
        ],  # constant for calculating E cost of HSP prod.
        "HSPProd_N_cost": parameters.loc[
            "HSP_prod_N_cost", 1
        ],  # constant for calculating N cost of HSP prod.
        #'Cell_maint_E_cost':          parameters.loc['Cell_maint_E_cost',1],     # maintenance energy cost for microbial biomass
        "Enzyme_Loss_Rate": parameters.loc[
            "Enz_loss_rate", 1
        ],  # enzyme turnover rate (=0.04; Allison 2006)
        "Osmo_Loss_Rate": parameters.loc["Osmo_loss_rate", 1],  # osmolyte turnover rate
        "HSPturnover": parameters.loc["HSP_turnover", 1],  # HSP turnover rate
        # 'Thermal_Traits':    expand(microbial_thermal_adaptation,gridsize),            # Taxon specific parameters for calculating temperature response of respiration costs
        "C_min": parameters.loc["C_min", 1],  # C threshold of cell lysis
        "N_min": parameters.loc["N_min", 1],  # N threshold of cell lysis
        "P_min": parameters.loc["P_min", 1],  # P threshold of cell lysis
        "Erange": parameters.loc["Erange", 1],  # Tolerance range for microbial E:MBC
        "Eratb": parameters.loc["Erat_b", 1],  # Ideal Energy:MBC ratio for bacteria
        "Eratf": parameters.loc["Erat_f", 1],  # Ideal Energy:MBC ratio for fungi
        "Fermoxp": parameters.loc[
            "Ferm_ox_p", 1
        ],  # prop of C converted to CO2 in monomer fermentation
        "FermCond": parameters.loc[
            "Ferm_Cond", 1
        ],  # favorability of fermentation based on environment
        "max_size_b": parameters.loc[
            "max_size_b", 1
        ],  # C quota threshold for bacterial cell division
        "max_size_f": parameters.loc[
            "max_size_f", 1
        ],  # C quota threshold for fungal cell division
        "wp_fc": parameters.loc[
            "wp_fc", 1
        ],  # threshold below which microbes start to respond to drought
        "wp_th": parameters.loc[
            "wp_th", 1
        ],  # threshold below which microbes in full swing to respond to drought
        # optional: moisture sensitivity of the uptake (f_psi = exp(slope * (Psi - wp_fc)) below wp_fc)
        "Uptake_psi_slope": _optional(parameters, "Uptake_psi_slope", 0.5),
        # optional: moisture sensitivity of the enzyme Vmax in the degradation (f_psi = exp(slope * (Psi - wp_fc)) below wp_fc)
        "Enzyme_psi_slope": _optional(parameters, "Enzyme_psi_slope", 0.25),
        # optional: maximum share of biomass (C+N+P) held as HSP C / osmolyte C; production of each compound is scaled by
        # (1 - share / max) and stops at the max (composition feedback, see Grid.metabolism)
        "HSP_max_frac": _optional(parameters, "HSP_max_frac", 0.2),
        "Osmo_max_frac": _optional(parameters, "Osmo_max_frac", 0.2),
        "alpha": parameters.loc[
            "alpha", 1
        ],  # factor delineating curve concavity of microbial response to drought
        "Temp": daily_temp,  # temperature
        "Psi": daily_psi,  # water potential
        "mmrt_Vmax_Km_slope": parameters.loc[
            "mmrt_Vmax_Km_slope", 1
        ],  # slope of Vmax_Km relationship for MMRT model
        "mmrt_Vmax_Km_int": parameters.loc[
            "mmrt_Vmax_Km_int", 1
        ],  # intercept of Vmax_Km relationship for MMRT model
        "uptake_mmrt_Vmax_Km_slope": parameters.loc[
            "Uptake_mmrt_Vmax_Km_slope", 1
        ],  # slope of Uptake Vmax_Km relationship for MMRT model
        "uptake_mmrt_Vmax_Km_int": parameters.loc[
            "Uptake_mmrt_Vmax_Km_int", 1
        ],  # intercept of Uptake Vmax_Km relationship for MMRT model
        #'mmrt_Vmax_Km_error':       parameters.loc['mmrt_Vmax_Km_error', 1],                # error for Vmax_Km relationships for MMRT model
        "mmrt_tref": parameters.loc[
            "MMRT_Tref", 1
        ],  # reference temperature (°K) for kinetic reactions
        "HSP_con_Ea": parameters.loc["HSP_Consti_Ea", 1],
        "HSP_ind_Ea": parameters.loc["HSP_Induci_Ea", 1],
        "HSP_Ea_slope": parameters.loc[
            "Maint_E_HSP_Ea_Slope", 1
        ],  # Slope for calculating maintenance energy rate Ea from heat stress protein abundance
        "HSP_Ea_int": parameters.loc[
            "Maint_E_HSP_Ea_Intercept", 1
        ],  # Slope for calculating maintenance energy rate Ea from heat stress protein abundance
        "HSP_Ea_n": parameters.loc[
            "Maint_E_HSP_Ea_n", 1
        ],  # Slope for calculating maintenance energy rate Ea from heat stress protein abundance
        "maintE_rate_Tref": parameters.loc[
            "Maint_E_Rate_Tref", 1
        ],  # reference temperature (°K) for maintenance energy cost
        "maintE_rate_const": parameters.loc[
            "Maint_E_Rate_Const", 1
        ],  # pre-exponential rate constant for calculating maintenance energy cost
        "Enz_km_temp_sens": parameters.loc[
            "Enz_Km_temp_sens", 1
        ],  # parameter specifying whether Km should vary with temperature
        "Uptake_km_temp_sens": parameters.loc[
            "Uptake_Km_temp_sens", 1
        ],  # parameter specifying whether Uptake Km should vary with temperature
        "km_tref": parameters.loc[
            "Km_Tref", 1
        ],  # Temperature to calculate Km, if Km is temperature independent
        "km_ea": parameters.loc[
            "Km_Ea", 1
        ],  # kj mol-1;activation energy for both enzyme and transporter
        "Arrhenius_tref": parameters.loc[
            "Arrhenius_Tref", 1
        ],  # reference temperature to normalize temperature response (default = 0)
    }

    return Data_Dictionary


# ── Grid-duplication detection/collapse ─────────────────────────────────────
# Every trait/kinetic table wrapped in utility.expand(df, gridsize) inside
# initialize_data() above is pd.concat([df]*gridsize) -- i.e. `gridsize`
# byte-for-byte identical stacked copies of a much smaller per-taxon /
# per-substrate / per-enzyme block. That is the entire source of the
# multi-GB Initialization/ export: none of those copies carry any
# information the first one doesn't already have, since nothing about them
# is spatially stochastic at t=0 (grid.py only ever reads these tables, it
# never mutates them -- see self.ReqEnz = data_init['ReqEnz'] etc. in
# grid.py's __init__).
#
# A handful of Initialization entries (Microbes.csv, and the biomass-derived
# columns of Initial_HSP_Prop.csv) are NOT expand()-wrapped: they come from
# the stochastic per-cell microbial placement in
# microbe.py's microbial_community_initialization(), so they are genuinely
# different from one grid cell to the next and must stay at full resolution.
#
# Rather than hard-coding which dict keys are which (fragile if
# initialize_data() changes), _dedup_grid_blocks() checks the data itself:
# if a table's row count divides evenly by gridsize AND every one of the
# gridsize stacked blocks is byte-identical to the first, only the first
# block is written to disk. Anything that fails either check (including a
# table that isn't grid-expanded at all, e.g. ReqEnz/EnzAttrib) is written
# in full, unchanged from before.
def _dedup_grid_blocks(values: np.ndarray, gridsize):
    """
    Check whether `values` (an array whose first axis is the row axis)
    consists of `gridsize` stacked, identical copies of a smaller block, as
    produced by utility.expand(). Returns the block size (number of rows in
    one copy) if so, otherwise None.
    """
    if gridsize is None or gridsize <= 1:
        return None
    n = values.shape[0]
    if n % gridsize != 0:
        return None
    block = n // gridsize
    if block == 0 or block == n:
        return None
    reshaped = values.reshape((gridsize, block) + values.shape[1:])
    first = reshaped[0]
    equal = reshaped == first
    if not equal.all():
        if np.issubdtype(values.dtype, np.floating):
            both_nan = np.isnan(reshaped) & np.isnan(first)
            equal = equal | both_nan
        if not equal.all():
            return None
    return block


def _write_dedup_manifest(base_path: Path, deduped: list, gridsize) -> None:
    """Write a small manifest + README documenting any grid-deduplicated fields."""
    if not deduped:
        return
    manifest = pd.DataFrame(deduped, columns=["field", "rows_saved", "gridsize"])
    manifest.to_csv(base_path / "grid_deduplicated_fields.csv", index=False)
    with open(base_path / "README_grid_dedup.txt", "w") as f:
        f.write(
            "The fields listed in grid_deduplicated_fields.csv were identical across\n"
            f"every one of the {gridsize} grid cells at the start of this run (they are\n"
            "static trait/kinetic tables broadcast onto the spatial grid by\n"
            "utility.expand()), so only one representative block of rows was written\n"
            f"here instead of {gridsize} duplicate copies. To reconstruct the full\n"
            "grid-expanded table exactly as it existed in memory during the run:\n\n"
            "    import pandas as pd\n"
            "    df = pd.read_csv('<field>.csv', index_col=0)\n"
            f"    df_full = pd.concat([df] * {gridsize})\n\n"
            "Fields NOT listed here (e.g. Microbes.csv, Initial_HSP_Prop.csv) were not\n"
            "uniform across the grid -- they depend on the stochastic microbial\n"
            "placement done in microbe.py's microbial_community_initialization() --\n"
            "and are saved at full spatial resolution as before.\n"
        )


# This is new
def export_initialization_dict_to_csv(
    base_path: Union[Path, str], d: dict, gridsize: int = None
) -> None:
    """Export contents of the initialisation directory to a folder.
    Writes each of the items of a type below to a separate CSV file
      - pandas.DataFrame
      - pandas.Series
      - numpy.ndarray of rank below 2
    All scalar numbers are grouped in a single CSV 'scalars.csv' file.

    When `gridsize` is provided, any table that is a byte-identical
    `gridsize`-fold repeat of a smaller block (as produced by
    utility.expand()) is written as a single representative block instead of
    gridsize duplicate copies -- see _dedup_grid_blocks() above. This is a
    lossless transformation: nothing about these tables varies from one grid
    cell to another at t=0, so the duplicate rows carry no information.
    Tables that are NOT uniform across the grid (e.g. Microbes.csv, which
    reflects stochastic per-cell microbial placement) are left untouched.
    Pass gridsize=None (the default) to fully disable this and reproduce the
    original behaviour.

    Note:
        All other items are ignored following a warning!
        If you need them written you need to add extra entry.
    """

    # Create space for output
    base_path = Path(base_path)
    base_path.mkdir(parents=True, exist_ok=True)

    # Collect all scalar numbers
    scalar_numbers = dict()
    deduped = []

    for name, member in d.items():
        if isinstance(member, (pd.DataFrame, pd.Series)):
            fname = name + ".csv"
            block = _dedup_grid_blocks(member.values, gridsize)
            if block is not None:
                member.iloc[:block].to_csv(base_path / fname)
                deduped.append((name, block, gridsize))
            else:
                member.to_csv(base_path / fname)
        elif isinstance(member, np.ndarray):
            if len(member.shape) <= 2:
                fname = name + ".csv"
                block = _dedup_grid_blocks(member, gridsize)
                if block is not None:
                    np.savetxt(base_path / fname, member[:block], delimiter=",")
                    deduped.append((name, block, gridsize))
                else:
                    np.savetxt(base_path / fname, member, delimiter=",")
            else:
                warnings.warn(
                    f"Member '{name}' of initialisation dictionary could not be saved since "
                    f"it is an array of rank higher than 2 (rank: {len(member.shape)})."
                )
        elif isinstance(member, numbers.Number):
            scalar_numbers[name] = member
        elif isinstance(member, pd.Index):
            pass  # Index objects are already captured as DataFrame row labels; no separate export needed
        else:
            warnings.warn(
                f"Initialisation member '{name}' has unsupported type '{type(member)}'. "
                f"It has not been exported to the output directory '{base_path}'."
            )

    # Print numbers
    pd.Series(scalar_numbers).to_csv(base_path / "scalars.csv")

    # Document any grid-deduplicated fields
    _write_dedup_manifest(base_path, deduped, gridsize)


def export_initialization_dict_to_netcdf(
    base_path: Union[Path, str], d: dict, gridsize: int = None
) -> None:
    """
    Export contents of the initialisation directory to a folder in NetCDF format.
    - pandas.DataFrame, pandas.Series, and numpy.ndarray (rank <= 2) are each
      written to a separate .nc file.
    - All scalar numbers are grouped in a single 'scalars.nc' file.

    When `gridsize` is provided, the same grid-block deduplication described
    in export_initialization_dict_to_csv() is applied before conversion to
    NetCDF -- see _dedup_grid_blocks() above.

    Note:
        All other items are ignored following a warning!
    """
    # Create space for output
    base_path = Path(base_path)
    base_path.mkdir(parents=True, exist_ok=True)

    # Collect all scalar numbers to be processed at the end
    scalar_numbers = dict()
    deduped = []

    for name, member in d.items():
        fname = name + ".nc"  # Use .nc extension for all files

        if isinstance(member, pd.DataFrame):
            try:
                block = _dedup_grid_blocks(member.values, gridsize)
                df_to_save = member.iloc[:block].copy() if block is not None else member.copy()
                if block is not None:
                    deduped.append((name, block, gridsize))
                if isinstance(df_to_save.index, pd.MultiIndex):
                    df_to_save = df_to_save.reset_index()
                ds = xr.Dataset.from_dataframe(df_to_save)
                encoding = {var: {"zlib": True, "complevel": 4} for var in ds.data_vars}
                ds.to_netcdf(base_path / fname, encoding=encoding, format="NETCDF4")
            except Exception as e:
                warnings.warn(
                    f"Could not convert DataFrame '{name}' to NetCDF. Error: {e}"
                )

        elif isinstance(member, pd.Series):
            try:
                block = _dedup_grid_blocks(member.values, gridsize)
                series_to_save = member.iloc[:block] if block is not None else member
                if block is not None:
                    deduped.append((name, block, gridsize))
                da = xr.DataArray.from_series(series_to_save.reset_index(drop=True))
                da.name = name
                da.to_netcdf(
                    base_path / fname,
                    encoding={name: {"zlib": True, "complevel": 4}},
                    format="NETCDF4",
                )
            except Exception as e:
                warnings.warn(
                    f"Could not convert Series '{name}' to NetCDF. Error: {e}"
                )

        elif isinstance(member, pd.Index):
            try:
                # Use the index's name for the dimension and variable, or a default.
                index_name = member.name if member.name is not None else "index_data"
                # Convert the Index to a self-describing DataArray and save.
                da = xr.DataArray(member, dims=[index_name], name=index_name)
                da.to_netcdf(
                    base_path / fname,
                    encoding={index_name: {"zlib": True, "complevel": 4}},
                    format="NETCDF4",
                )
            except Exception as e:
                warnings.warn(f"Could not convert Index '{name}' to NetCDF. Error: {e}")

        elif isinstance(member, np.ndarray):
            if member.ndim <= 2:
                try:
                    block = _dedup_grid_blocks(member, gridsize)
                    array_to_save = member[:block] if block is not None else member
                    if block is not None:
                        deduped.append((name, block, gridsize))
                    # Convert numpy array to an xarray DataArray.
                    # We give it a name and default dimension names.
                    data_array = xr.DataArray(
                        array_to_save, name=name, dims=[f"dim_{i}" for i in range(array_to_save.ndim)]
                    )
                    data_array.to_netcdf(
                        base_path / fname,
                        encoding={name: {"zlib": True, "complevel": 4}},
                        format="NETCDF4",
                    )
                except Exception as e:
                    warnings.warn(
                        f"Could not convert array '{name}' to NetCDF. Error: {e}"
                    )
            else:
                warnings.warn(
                    f"Member '{name}' of initialisation dictionary could not be saved since "
                    f"it is an array of rank higher than 2 (rank: {member.ndim})."
                )

        elif isinstance(member, numbers.Number):
            scalar_numbers[name] = member

        else:
            warnings.warn(
                f"Initialisation member '{name}' has unsupported type '{type(member)}'. "
                f"It has not been exported to the output directory '{base_path}'."
            )

    # Process and save all collected scalars
    if scalar_numbers:
        # Create an xarray Dataset from the dictionary of scalars and save to .nc
        scalars_dataset = xr.Dataset(scalar_numbers)
        encoding = {
            var: {"zlib": True, "complevel": 4} for var in scalars_dataset.data_vars
        }
        scalars_dataset.to_netcdf(
            base_path / "scalars.nc", encoding=encoding, format="NETCDF4"
        )

    # Document any grid-deduplicated fields
    _write_dedup_manifest(base_path, deduped, gridsize)


def export_initialization_dict_to_zarr(store_path: str, d: dict, gridsize: int = None) -> None:
    """Export initialization dictionary to a Zarr store.

    When `gridsize` is provided, the same grid-block deduplication described
    in export_initialization_dict_to_csv() is applied before writing -- see
    _dedup_grid_blocks() above.
    """
    import zarr
    from zarr.codecs import ZstdCodec

    compressor = ZstdCodec(level=3)
    scalar_numbers = dict()
    deduped = []

    for name, member in d.items():
        try:
            if isinstance(member, pd.DataFrame):
                block = _dedup_grid_blocks(member.values, gridsize)
                df = member.iloc[:block].copy() if block is not None else member.copy()
                if block is not None:
                    deduped.append((name, block, gridsize))
                if isinstance(df.index, pd.MultiIndex):
                    df = df.reset_index()
                ds = xr.Dataset.from_dataframe(df)
                encoding = {var: {"compressors": compressor} for var in ds.data_vars}
                ds.to_zarr(store_path, group=name, encoding=encoding, mode="a")

            elif isinstance(member, pd.Series):
                block = _dedup_grid_blocks(member.values, gridsize)
                series = member.iloc[:block] if block is not None else member
                if block is not None:
                    deduped.append((name, block, gridsize))
                da = xr.DataArray(
                    series.reset_index(drop=True).values, dims=["index"], name=name
                )
                ds = da.to_dataset()
                ds.to_zarr(
                    store_path,
                    group=name,
                    encoding={name: {"compressors": compressor}},
                    mode="a",
                )

            elif isinstance(member, pd.Index):
                pass  # captured via DataFrame row labels

            elif isinstance(member, np.ndarray) and member.ndim <= 2:
                block = _dedup_grid_blocks(member, gridsize)
                array_to_save = member[:block] if block is not None else member
                if block is not None:
                    deduped.append((name, block, gridsize))
                da = xr.DataArray(
                    array_to_save, name=name, dims=[f"dim_{i}" for i in range(array_to_save.ndim)]
                )
                ds = da.to_dataset()
                ds.to_zarr(
                    store_path,
                    group=name,
                    encoding={name: {"compressors": compressor}},
                    mode="a",
                )

            elif isinstance(member, numbers.Number):
                scalar_numbers[name] = member

            else:
                warnings.warn(
                    f"Initialisation member '{name}' has unsupported "
                    f"type '{type(member)}' and was not exported."
                )
        except Exception as e:
            warnings.warn(f"Could not export '{name}' to Zarr. Error: {e}")

    if scalar_numbers:
        xr.Dataset(scalar_numbers).to_zarr(store_path, group="scalars", mode="a")

    # Document any grid-deduplicated fields (written next to the store)
    _write_dedup_manifest(Path(str(store_path)).parent, deduped, gridsize)
