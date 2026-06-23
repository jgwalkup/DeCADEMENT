# ------------------
# this module, grid.py, deal with calculations of all microbe-related activites on a spatial grid with a class, Grid().
# Modified to calculate kinetics of substrate degradation, and monomer uptake, and metabolism using MMRT kinetics
# Revised by Jeth Walkup, 2024
# ------------------

import sys
import numpy as np
import pandas as pd
import scipy.constants as const

from microbe import microbe_osmo_psi
from microbe import microbe_mortality_prob as MMP
from utility import expand
from enzyme import Arrhenius
from enzyme import MMRT



# ── Fast groupby-sum replacements ────────────────────────────────────────────
# Replace pd.groupby(level=0).sum() on fixed-structure MultiIndex arrays.
# The MultiIndex is always [taxon_0..taxon_n]*gridsize or [enz_0..enz_n]*gridsize
# so reshape+sum is mathematically identical but ~60x faster.

def _taxa_sum(arr, gridsize, n_taxa, index):
    """Sum a (gridsize*n_taxa, n_cols) DataFrame or (gridsize*n_taxa,) Series
    over the grid dimension, returning a (n_taxa, n_cols) DataFrame or
    (n_taxa,) Series with the taxon index."""
    if hasattr(arr, "values"):
        vals = arr.values
    else:
        vals = np.asarray(arr)
    if vals.ndim == 1:
        result = vals.reshape(gridsize, n_taxa).sum(axis=0)
        return pd.Series(result, index=index, dtype=arr.dtype if hasattr(arr, "dtype") else None)
    else:
        result = vals.reshape(gridsize, n_taxa, vals.shape[1]).sum(axis=0)
        return pd.DataFrame(result, index=index, columns=arr.columns)

def _enz_sum(arr, gridsize, n_enzymes, index):
    """Sum a (gridsize*n_enzymes,) Series over the grid dimension."""
    vals = arr.values if hasattr(arr, "values") else np.asarray(arr)
    result = vals.reshape(gridsize, n_enzymes).sum(axis=0)
    return pd.Series(result, index=index, dtype=arr.dtype if hasattr(arr, "dtype") else None)

def _mon_sum(arr, gridsize, n_monomers, index):
    """Sum a (gridsize*n_monomers, n_taxa) DataFrame over the grid dimension."""
    vals = arr.values if hasattr(arr, "values") else np.asarray(arr)
    if vals.ndim == 1:
        result = vals.reshape(gridsize, n_monomers).sum(axis=0)
        return pd.Series(result, index=index)
    else:
        result = vals.reshape(gridsize, n_monomers, vals.shape[1]).sum(axis=0)
        return pd.DataFrame(result, index=index, columns=arr.columns)

class Grid():
    """
    This class holds all variables related to microbe, substrate, monomer, and enzyme over the spatial grid.

    Accepts returns from the module 'initialization.py' and includes methods as follows:
        1) degrdation():   explicit substrate degradation
        2) uptake():       explicit monomers uptake
        3) metabolism():   cellular processes and emergent CAE, CUE and respiration
        4) mortality():    determine mortality of microbial cells based on mass thresholds
        5) reproduction(): compute cell division and dispersal
        6) repopulation(): resample taxa from the microbial pool and place them on the grid
    Coding philosophy:
        Each method starts with passing some global variables to local ones and creating
        some indices facilitating dataframe index/column processing and ends up with updating
        state variables and passing them back to the global ones. All computation stays in between.
    Reminder:
        Keep a CLOSE EYE on the indexing throughout the matrix/dataframe operations
    """

    def __init__(self, runtime, data_init):
        """
        The constructor of Grid class.

        Parameters:
            runtime:   user-specified parameters
            data_init: dictionary;initialized data from the module 'initialization.py'
        """

        self.cycle = int(runtime.loc['end_time', 1])
        self.gridsize = int(runtime.loc['gridsize', 1])
        self.n_taxa = int(runtime.loc["n_taxa", 1])
        self.n_substrates = int(runtime.loc["n_substrates", 1])
        self.n_enzymes = int(runtime.loc["n_enzymes", 1])
        self.n_monomers = self.n_substrates + 3
        self.kinetics_model = int(
            runtime.loc["Kinetics_Model", 1])  # specifies use of Arrhenius or MMRT to calculate rxn rates
        # self.kinetics_model = runtime.loc["Kinetics_Model", 1]  # Arrhenius or MMRT to calculate rxn rates

        # Degradation
        self.Substrates_init = data_init['Substrates']  # Substrates initialized
        self.Substrates = data_init['Substrates'].copy(deep=True)  # Substrates;df; w/ .copy() avoiding mutation
        self.SubInput = data_init['SubInput']  # Substrate inputs
        self.Enzymes_init = data_init['Enzymes']  # Initial pool of Enzymes
        self.Enzymes = data_init['Enzymes'].copy(deep=True)  # Enzymes
        self.ReqEnz = data_init['ReqEnz']  # Enzymes required by each substrate
        self.Ea = data_init['Ea']  # Enzyme activation energy
        self.HT0 = data_init['Vmax_HT0']  # Enzyme enthalpy
        self.ST0 = data_init['Vmax_ST0']  # Enzyme entropy
        self.Cp = data_init['Vmax_Cp']  # Enzyme heat capacity
        self.Km_HT0 = data_init['Km_HT0']  # Enzyme enthalpy (for Km calculation)
        self.Km_ST0 = data_init['Km_ST0']  # Enzyme entropy (for Km calculation)
        self.Km_Cp = data_init['Km_Cp']  # Enzyme heat capacity (for Km calculation)
        self.Vmax0 = data_init['Vmax0']  # Max. reaction speed # not used in MMRT model (Jeth)
        self.Km0 = data_init['Km0']  # Half-saturation constant # not used in MMRT model (Jeth)
        self.SubstrateRatios = np.float32('nan')  # Substrate stoichiometry
        self.DecayRates = np.float32('nan')  # Substrate decay rate
        self.Vmax_Km_slope = data_init['mmrt_Vmax_Km_slope']  # slope of the enzyme Vmax/Km relationship
        self.Vmax_Km_int = data_init['mmrt_Vmax_Km_int']  # intercept of the Uptake Vmax/Km relationship
        self.Enz_Vmax_sub_spec = data_init[
            'Enz_vmax_sub_spec']  # substrate specificity-efficiency adjustment for degradation Vmax
        self.Enz_Km_sub_spec = data_init[
            'Enz_km_sub_spec']  # substrate specificity-efficiency adjustment for degradation Km
        self.Enz_Km_Temp_Sens = data_init[
            'Enz_km_temp_sens']  # substrate specificity-efficiency adjustment for monomer uptake Vmax
        self.Arrh_Tref = data_init['Arrhenius_tref']

        # Uptake
        self.Microbes_init = data_init['Microbes_pp']  # microbial community before placement
        self.Microbes = data_init['Microbes'].copy(deep=True)  # microbial community after placement
        self.Monomers_init = data_init['Monomers']  # Monomers initialized
        self.Monomers = data_init['Monomers'].copy(deep=True)  # Monomers
        self.MonInput = data_init['MonInput']  # Inputs of monomers
        self.Uptake_Ea = data_init['Uptake_Ea']  # transporter enzyme Ea
        self.Uptake_Vmax0 = data_init['Uptake_Vmax0']  # transporter Vmax
        self.Uptake_Km0 = data_init['Uptake_Km0']  # transporter Km
        self.Uptake_HT0 = data_init['Uptake_HT0']  # transporter HT0
        self.Uptake_ST0 = data_init['Uptake_ST0']  # transporter ST0
        self.Uptake_Cp = data_init['Uptake_Cp']  # transporter Cp
        self.Monomer_ratios = data_init['Monomer_ratio'].copy(deep=True)  # monomer stoichiometry
        self.Monomer_energy = data_init['MonEnergy_Attrib']  # monomer energy characteristics
        self.Mon_energy_resp = expand(
            pd.concat([pd.Series([0, 0], index=['NH4', 'PO4']),
                       self.Monomer_energy['E_resp']], sort=False),
            self.gridsize)
        self.Mon_energy_ferm = expand(
            pd.concat([pd.Series([0, 0], index=['NH4', 'PO4']),
                       self.Monomer_energy['E_ferm']], sort=False),
            self.gridsize)
        self.Uptake_ReqEnz = data_init['Uptake_ReqEnz']  # Enzymes required by monomers
        self.Uptake_Enz_Cost = data_init['UptakeGenesCost']  # Cost of encoding each uptake gene
        self.Uptake_Prod_Cost = data_init['UptakeProdCost']  # C Cost of transporter production
        self.UptakeProd_E_Cost = data_init['Uptake_prod_E_cost']  # Energy cost of transporter production per mg C
        self.Taxon_Uptake_C = np.float32('nan')  # taxon total C uptake
        self.Taxon_Uptake_BC = np.float32('nan')  # taxon uptake of C allocated to biomass
        self.Taxon_Uptake_N = np.float32('nan')  # taxon uptake of N
        self.Taxon_Uptake_P = np.float32('nan')  # taxon uptake of P
        self.Taxon_Uptake_E = np.float32('nan')  # taxon uptake of Energy
        self.Microbe_E_Gain = np.float32('nan')  # taxon net energy gain
        self.Microbe_E_Cost = np.float32('nan')  # taxon net energy costs
        self.Uptake_Vmax_Km_slope = data_init['uptake_mmrt_Vmax_Km_slope']  # slope of the Uptake Vmax/Km relationship
        self.Uptake_Vmax_Km_int = data_init['uptake_mmrt_Vmax_Km_int']  # intercept for the Uptake Vmax/Km relationship
        self.Uptake_sub_spec = data_init[
            'Uptake_Vmax_sub_spec']  # substrate specificity-efficiency adjustment for monomer uptake Vmax
        self.Uptake_Km_Temp_Sens = data_init[
            'Uptake_km_temp_sens']  # substrate specificity-efficiency adjustment for monomer uptake Vmax
        self.km_Tref = data_init['km_tref']  # substrate specificity-efficiency adjustment for monomer uptake Vmax
        self.Pox_Ferm = data_init['Fermoxp']  # proportion of fermented monomer converted to CO2
        # fermentation condition score (0-1): reflects environmental favorability of fermentation vs respiration
        self.ferm_cond = data_init['FermCond']  # currently set to 0.1; can update to vary with environmental conditions
        self.active_uptake_E_cost = data_init['ActiveUptakeE_cost']  # energy cost of transporting monomers

        # Metabolism
        self.Consti_Enzyme_C = data_init["EnzProdConstit"]  # C cost of encoding constitutive enzyme
        self.Induci_Enzyme_C = data_init["EnzProdInduce"]  # C Cost of encoding inducible enzyme
        self.Consti_Osmo_C = data_init['OsmoProdConsti']  # C Cost of encoding constitutive osmolyte
        self.Induci_Osmo_C = data_init['OsmoProdInduci']  # C Cost of encoding inducible osmolyte
        self.Enz_Attrib = data_init['EnzAttrib']  # Enzyme attributes; dataframe
        # self.cell_maint_E_cost = data_init['Cell_maint_E_cost']  # energy cost of biomass maintenance
        self.Osmo_E_cost = data_init['OsmoProd_E_cost']  # E cost per unit of osmo-C production
        self.Osmo_N_cost = data_init['OsmoProd_N_cost']  # N cost per unit of osmo-C production
        self.HSP_E_cost = data_init['HSPProd_E_cost']  # E cost per unit of HSP-C production
        self.HSP_N_cost = data_init['HSPProd_N_cost']  # N cost per unit of HSP-C production
        # self.Osmo_Maint_cost_constant = data_init['Osmo_Maint_cost_constant']       # Respiration cost constant for osmolyte production
        self.Respiration = np.float32('nan')  # Respiration
        self.Taxon_Resp = np.float32('nan')  # Taxon specific Respiration
        self.Taxon_C_Uptake = np.float32('nan')  # Taxon-specific total C uptake
        self.Taxon_E_Uptake = np.float32('nan')  # Taxon-specific E production
        self.Taxon_E_Prod = np.float32('nan')  # Taxon-specific E production
        self.Taxon_BC_Uptake = np.float32('nan')  # Taxon-specific C uptake allocated to biomass
        self.Taxon_Total_Uptake = np.float32('nan')  # Taxon-specific biomass uptake (C+N+P)
        self.Enzyme_Loss_Rate = data_init['Enzyme_Loss_Rate']
        self.OSmo_Loss_Rate = data_init['Osmo_Loss_Rate']
        self.maint_E_HSP_Ea_slope = data_init['HSP_Ea_slope']
        self.maint_E_HSP_Ea_int = data_init['HSP_Ea_int']
        self.maint_E_HSP_Ea_n = data_init['HSP_Ea_n']
        self.maint_E_rate_Tref = data_init['maintE_rate_Tref']
        self.maint_E_rate_const = data_init['maintE_rate_const']
        self.HSP_consti_pre_exp = data_init[
            'HSPProdConsti']  # pre exponential factor for temperature induced constitutive HSP production
        self.HSP_induci_pre_exp = data_init[
            'HSPProdInduci']  # pre exponential factor for temperature induced inducible HSP production
        self.HSP_consti_Ea = data_init[
            'HSP_con_Ea']  # exponential factor for temperature induced constitutive HSP production
        self.HSP_induci_Ea = data_init[
            'HSP_ind_Ea']  # exponential factor for temperature induced inducible HSP production
        # self.Transporters = float('nan')
        # self.Osmolyte_Con = float('nan')
        # self.Osmolyte_Ind = float('nan')
        # self.Enzyme_Con   = float('nan')
        # self.Enzyme_Ind   = float('nan')
        self.HSP_Turnover = data_init['HSPturnover']
        self.Taxon_HSP = data_init['Initial_HSP_Prop']  #

        # Osmolyte pool: same structure as constitutive osmolyte production DataFrame
        # (gridsize*n_taxa rows × osmolyte-gene columns), initialized to zero.
        # Updated in metabolism() each step with production and turnover.
        self.Taxon_Osmolyte = data_init['OsmoProdConsti'].copy().astype('float32')
        self.Taxon_Osmolyte[:] = np.float32(0)

        # Mortality
        self.MinRatios = data_init['MinRatios']  # Minimal cell quotas
        self.C_min = data_init['C_min']  # C threshold value of living cell
        self.N_min = data_init['N_min']  # N threshold value of living cell
        self.P_min = data_init['P_min']  # P threshold value of living cell
        self.E_range = data_init['Erange']  # tolerance on energy:biomass C ratio
        self.E_ratb = data_init['Eratb']  # ideal E:MBC ratio for bacteria
        self.E_ratf = data_init['Eratf']  # ideal E:MBC ratio for fungi
        self.basal_death_prob = data_init['basal_death_prob']  # Basal death probability of microbes
        self.death_rate = data_init['death_rate']  # Change rate of mortality with water potential
        self.tolerance = data_init['TaxDroughtTol']  # taxon drought tolerance
        self.wp_fc = data_init['wp_fc']  # -1.0
        self.wp_th = data_init['wp_th']  # -6.0
        self.alpha = data_init['alpha']  # integer;1
        self.Kill = np.float32('nan')  # total number of cells stochastically killed

        # Reproduction
        self.fb = data_init['fb']  # index of fungal taxa (=1)
        self.max_size_b = data_init['max_size_b']  # threshold of cell division
        self.max_size_f = data_init['max_size_f']  # threshold of cell division
        self.x = int(runtime.loc['x', 1])  # x dimension of grid
        self.y = int(runtime.loc['y', 1])  # y dimension of grid
        self.dist = int(runtime.loc['dist', 1])  # maximum dispersal distance: 1 cell
        self.direct = int(runtime.loc['direct', 1])  # dispersal direction: 0.95

        # Climate data
        self.temp = data_init['Temp']  # series; temperature
        self.psi = data_init['Psi']  # series; water potential

        # Global constants
        self.Km_Ea = np.float32(20)  # kj mol-1;activation energy for both enzyme and transporter
        self.MMRT_Tref = data_init['mmrt_tref']  # reference temperature for calculating rates (20 celcius)
        self.gas_const = np.float32(const.R)  # Universal gas constant
        self.boltz_const = np.float32(const.k)  # Boltzmann's constant
        self.planck_const = np.float32(const.h)  # Planck's constant
        self.Substrates_id = data_init['Substrate_index']
        self.monomer_index = data_init['Monomer_index']
        # self.monomer_index = ["NH4", "PO4", "FermProd", "DeadMic", "DeadEnz"] + ["Mon" + str(i) for i in range(4, self.n_monomers - 2 + 1)]


        # ── Fast groupby helpers ──────────────────────────────────────────────
        # These replace pd.DataFrame.groupby(level=0).sum() calls throughout
        # the model.  The MultiIndex structure is fixed for the entire run:
        #   taxa-indexed arrays:   gridsize repetitions of n_taxa taxon labels
        #   enzyme-indexed arrays: gridsize repetitions of n_enzymes labels
        # Reshaping to (gridsize, n) and summing axis=0 is equivalent to
        # groupby(level=0).sum() but ~60x faster (no index hashing required).

        # Pre-built index arrays for reconstructing output Series/DataFrames
        self._taxa_index = self.Microbes.index[:self.n_taxa]
        self._enz_index  = self.Enzymes.index[:self.n_enzymes]

    def degradation(self, pulse, day):
        """
        Explicit degradation of different substrates.

        Calculation procedure:
        1. Determine substates pool: incl. inputs
        2. Compute Vmax & Km and make them follow the index of Substrates
        3. Follow the Michaelis-Menten equation to compute full degradation rate
        4. Impose the substrate-required enzymes upon the full degradation rate
        5. Adjust cellulose rate with LCI(lignocellulose index)
        """

        # constant of lignocellulose index--LCI
        LCI_slope = np.float32(-0.8)
        # Substrates index by substrate names
        Sub_index = self.Substrates.index

        # total mass of each substrate: C+N+P
        rss = self.Substrates.sum(axis=1)
        # substrate stoichiometry; NOTE:ensure NA(b/c of 0/0 in df) = 0
        SubstrateRatios = self.Substrates.divide(rss, axis=0)
        SubstrateRatios = SubstrateRatios.fillna(0)

        # moisture effects on enzymatic kinetics
        if self.psi[day] >= self.wp_fc:
            f_psi = np.float32(1.0)
        else:
            f_psi = np.exp(0.25 * (self.psi[day] - self.wp_fc)).astype('float32')

        # calculate enzyme Vmax using the user specified kinetics model (MMRT = 1; Arrhenius = 0)
        if self.kinetics_model > 0:
            Vmax_mmrt = MMRT(self.Cp, self.HT0, self.ST0, self.temp[day], self.MMRT_Tref, self.boltz_const,
                             self.planck_const, self.gas_const)
            Vmax_mmrt_sub_spec = Vmax_mmrt * self.Enz_Vmax_sub_spec  # account for substrate specificity-efficiency tradeoff
            Vmax = Vmax_mmrt_sub_spec * f_psi
            # If specified, derive Km_mmrt based on Vmax-Km tradeoff; otherwise, calculate Km assuming
            if self.Enz_Km_Temp_Sens == 1:
                Km_mmrt = MMRT(self.Km_Cp, self.Km_HT0, self.Km_ST0, self.temp[day], self.MMRT_Tref, self.boltz_const,
                               self.planck_const, self.gas_const)
            else:
                Km_mmrt = MMRT(self.Km_Cp, self.Km_HT0, self.Km_ST0, self.km_Tref, self.MMRT_Tref, self.boltz_const,
                               self.planck_const, self.gas_const)
            # account for substrate specifity
            Km_mmrt_sub_spec = Km_mmrt * self.Enz_Km_sub_spec  # account for substrate specificity-efficiency tradeoff
            # # Derive Km based on Vmax-Km tradeoff
            Km = (Km_mmrt_sub_spec * self.Vmax_Km_slope) + self.Vmax_Km_int

        else:
            Vmax = self.Vmax0 * Arrhenius(self.Ea, self.gas_const, self.temp[day],
                                          self.Arrh_Tref) * f_psi  # Vmax: (enz*gridsize) * sub
            Km = self.Km0 * Arrhenius(self.Km_Ea, self.gas_const, self.temp[day],
                                      self.Arrh_Tref)  # Km:   (sub*gridsize) * enz

        # Multiply Vmax by enzyme concentration
        tev_transition = Vmax.mul(self.Enzymes, axis=0)  # (enz*gridsize) * sub
        tev_transition.index = [np.arange(self.gridsize).repeat(self.n_enzymes),
                                tev_transition.index]  # create a MultiIndex
        tev = tev_transition.stack().unstack(1).reset_index(level=0, drop=True)  # (sub*gridsize) * enz
        tev = tev[Km.columns]  # ensure to re-order the columns b/c of python's default alphabetical ordering

        # Michaelis-Menten equation
        Decay = tev.mul(rss, axis=0) / Km.add(rss, axis=0)

        # Pull out each batch of required enzymes and sum across redundant enzymes
        batch1 = (self.ReqEnz.loc['set1'].values * Decay).sum(axis=1)
        # batch2 = (self.ReqEnz.loc['set2'].values * Decay).sum(axis=1)

        # Assess the rate-limiting enzyme and set decay to that rate
        # DecaySums = pd.concat([batch1, batch2],axis=1)
        # DecayRates0 = DecaySums.min(axis=1, skipna=True)

        # Compare to substrate available and take the min, allowing for a tolerance of 1e-9
        DecayRates = pd.concat([batch1, rss], axis=1, sort=False).min(axis=1, skipna=True)

        # Adjust cellulose rate by linking cellulose degradation to lignin concentration (LCI)
        ss7 = self.Substrates.loc[Sub_index == 'Lignin'].sum(axis=1).values
        DecayRates.loc[Sub_index == 'Cellulose'] *= np.float32(1) + (
                ss7 / (ss7 + self.Substrates.loc[Sub_index == 'Cellulose', 'C'])) * LCI_slope

        # Update Substrates Pool by removing decayed C, N, & P. Depending on specific needs, adding inputs of substrates can be done here
        self.Substrates -= SubstrateRatios.mul(DecayRates, axis=0)  # + self.SubInput

        # Pass these two back to the global variables to be used in the next method
        self.SubstrateRatios = SubstrateRatios
        self.DecayRates = DecayRates

    def uptake(self, pulse, day):
        """
        Explicit uptake of different monomers by transporters following the Michaelis-Menten equation.

        Calculaton procedure:
        -> Average monomers across the grid:
        -> Determine pool of monomers: add degradation and input, update stoichimoetry
        -> Maximum uptake:
        -> Uptake by Monomer:
        -> Uptake by Taxon:
        """

        # Every monomer averaged over the grid in each time step
        self.Monomers = expand(
            pd.DataFrame(
                self.Monomers.values.reshape(self.gridsize, self.n_monomers, 3).sum(axis=0) / self.gridsize,
                index=self.monomer_index, columns=self.Monomers.columns),
            self.gridsize)

        # Indices
        is_org = (self.Monomers.index != "NH4") & (self.Monomers.index != "PO4") & (
                self.Monomers.index != "FermProd")  # organic monomers
        # is_mineral = (Monomers.index == "NH4") | (Monomers.index == "PO4")

        # Update monomer ratios in each time step with organic monomers following the substrates
        self.Monomer_ratios[is_org] = self.SubstrateRatios.values

        # Determine monomer pool from decay and input
        # Organic monomers derived from substrate-decomposition
        Decay_Org = self.Monomer_ratios[is_org].mul(self.DecayRates.values, axis=0)
        # inputs of organic and mineral monomers
        # Input_Org = MR_transition[is_org].mul(self.MonInput[is_org].tolist(),axis=0)
        # Input_Mineral = MR_transition[is_mineral].mul((self.MonInput[is_mineral]).tolist(),axis=0)
        # Monomer pool determined
        self.Monomers.loc[is_org] += Decay_Org  # + Input_Org
        # self.Monomers.loc[is_mineral] += Input_Mineral

        # Get the total mass of each monomer: C+N+P
        rsm = self.Monomers.sum(axis=1)
        # Recalculate monomer ratios after updating monomer pool and before uptake calculation
        self.Monomer_ratios.loc[is_org] = self.Monomers.loc[is_org].divide(rsm[is_org], axis=0)
        self.Monomer_ratios = self.Monomer_ratios.fillna(0)

        # Start calculating monomer uptake
        # Moisture impacts on uptake, mimicking the diffusivity implications
        if self.psi[day] >= self.wp_fc:
            f_psi = np.float32(1.0)
        else:
            f_psi = np.exp(0.5 * (self.psi[day] - self.wp_fc)).astype('float32')

        if self.kinetics_model > 0:
            Uptake_Vmax_mmrt = MMRT(self.Uptake_Cp, self.Uptake_HT0, self.Uptake_ST0, self.temp[day], self.MMRT_Tref,
                                    self.boltz_const, self.planck_const, self.gas_const)
            Uptake_Vmax_sub_spec = Uptake_Vmax_mmrt * self.Uptake_sub_spec  # account for substrate specificity-efficiency tradeoff
            Uptake_Vmax = Uptake_Vmax_sub_spec * f_psi  # account for moisture limitation

            # If specified, derive Km_mmrt based on Vmax-Km tradeoff; otherwise, calculate Km assuming
            if self.Uptake_Km_Temp_Sens == 1:
                Uptake_Km = (Uptake_Vmax_sub_spec * self.Uptake_Vmax_Km_slope) + self.Uptake_Vmax_Km_int
            else:
                Uptake_Km_mmrt = MMRT(self.Uptake_Cp, self.Uptake_HT0, self.Uptake_ST0, self.km_Tref, self.MMRT_Tref,
                                      self.boltz_const, self.planck_const, self.gas_const)
                Uptake_Km_sub_spec = Uptake_Km_mmrt * self.Uptake_sub_spec
                Uptake_Km = (Uptake_Km_sub_spec * self.Uptake_Vmax_Km_slope) + self.Uptake_Vmax_Km_int

        else:
            Uptake_Vmax = self.Uptake_Vmax0 * Arrhenius(self.Uptake_Ea, self.gas_const, self.temp[day],
                                                        self.Arrh_Tref) * f_psi
            Uptake_Km = self.Uptake_Km0 * Arrhenius(self.Km_Ea, self.gas_const, self.temp[day], self.Arrh_Tref)

        # Equation for hypothetical potential uptake (per unit of compatible uptake protein)
        Potential_Uptake = (self.Uptake_ReqEnz * Uptake_Vmax).mul(rsm.values, axis=0) / Uptake_Km.add(rsm.values,
                                                                                                      axis=0)

        # Derive the mass of each transporter of each taxon NOTE: transpose the df to Upt*(Taxa*grid)
        MicCXGenes = (self.Uptake_Enz_Cost.mul(self.Microbes[['C', 'N', 'P']].sum(axis=1), axis=0)).T

        # ── Uptake computation: numpy-native to minimise working-set size ──────────
        # All intermediate arrays are kept as (G, M, T) numpy arrays.
        # This avoids materialising 11 separate (G*M, T) Pandas DataFrames
        # simultaneously, cutting the peak working set from ~63 MB (30×30) /
        # ~700 MB (100×100) to ~20 MB / ~220 MB and reducing DRAM pressure.

        G, M, U, T = self.gridsize, self.n_monomers, Potential_Uptake.shape[1], self.n_taxa

        # ── Step 1: Max_Uptake via batched matmul, then scale to available monomer ──
        P_3d  = Potential_Uptake.values.reshape(G, M, U)               # (G, M, U)
        Mc_3d = MicCXGenes.values.reshape(U, G, T).transpose(1, 0, 2)  # (G, U, T)
        # max_up: (G, M, T)  — potential uptake per grid × monomer × taxon
        max_up = (P_3d @ Mc_3d).astype(np.float32)                     # (G, M, T)

        # Scale uptake to available monomer mass (rsm: G*M vector)
        rsm_3d = rsm.values.reshape(G, M, 1).astype(np.float32)        # (G, M, 1)
        csmu   = max_up.sum(axis=2, keepdims=True)                      # (G, M, 1)
        scale  = np.zeros_like(csmu, dtype=np.float32)                  # (G, M, 1)
        np.divide(np.minimum(rsm_3d, csmu), csmu,
                  out=scale, where=(csmu > 0))                           # avoids div/0
        upt    = (max_up * scale).astype(np.float32)                    # (G, M, T)

        # ── Step 2: Update Monomers pool ─────────────────────────────────────────
        # upt.sum(axis=2): (G, M) — total uptake per grid × monomer
        upt_sum_T = upt.sum(axis=2)                                     # (G, M)
        MR_C = self.Monomer_ratios["C"].values.reshape(G, M).astype(np.float32)
        MR_N = self.Monomer_ratios["N"].values.reshape(G, M).astype(np.float32)
        MR_P = self.Monomer_ratios["P"].values.reshape(G, M).astype(np.float32)
        mon_depl = np.stack([MR_C, MR_N, MR_P], axis=2) * upt_sum_T[:, :, None]  # (G,M,3)
        self.Monomers -= mon_depl.reshape(G * M, 3)

        # ── Step 3: C uptake allocation to energy vs biomass ─────────────────────
        self.E_max    = self.E_ratb + self.E_range
        Tax_E_demand  = 1 - ((1 - ((self.E_max - (self.Microbes['E'] /
                                self.Microbes['C'])) / self.E_max)) ** (self.ferm_cond + 1))

        Org_Mon_E_value = self.Monomer_energy['E_value'] ** (
            (2 / np.log10(self.Monomer_energy['E_ferm'])) ** self.ferm_cond)
        Mon_E_value = np.concatenate([[0, 0], Org_Mon_E_value.values]).astype(np.float32)  # (M,)
        # Tax_E_alloc: (M, G*T) → reshape to (G, M, T)
        Tea_3d = np.outer(Mon_E_value,
                          Tax_E_demand.values.astype(np.float32)
                          ).reshape(M, G, T).transpose(1, 0, 2).astype(np.float32)  # (G,M,T)
        Tea_3d = np.nan_to_num(Tea_3d, nan=0.0)

        # Store Tax_E_alloc as DataFrame for any downstream code that reads it
        self.Tax_E_alloc = pd.DataFrame(
            Tea_3d.transpose(1, 0, 2).reshape(M, G * T),
            index=self.monomer_index,
            columns=self.Microbes.index)

        c_upt_3d  = upt * MR_C[:, :, None]                             # (G, M, T) C uptake
        ce_3d     = (c_upt_3d * Tea_3d).astype(np.float32)             # (G, M, T) C→energy
        cb_3d     = (c_upt_3d - ce_3d).astype(np.float32)              # (G, M, T) C→biomass
        n_3d      = upt * MR_N[:, :, None]                             # (G, M, T) N uptake
        p_3d      = upt * MR_P[:, :, None]                             # (G, M, T) P uptake

        # ── Step 4: Energy pathways ───────────────────────────────────────────────
        ME_resp  = self.Mon_energy_resp.values.reshape(G, M).astype(np.float32)  # (G,M)
        ME_ferm  = self.Mon_energy_ferm.values.reshape(G, M).astype(np.float32)
        ce_ferm  = ce_3d * self.ferm_cond                               # (G, M, T)
        ce_resp  = ce_3d * (1 - self.ferm_cond)                        # (G, M, T)
        e_3d     = (ce_resp * ME_resp[:, :, None]
                    + ce_ferm * ME_ferm[:, :, None]).astype(np.float32) # (G, M, T) energy
        co2_3d   = ce_resp + ce_ferm * self.Pox_Ferm                   # (G, M, T) CO2
        ferm_3d  = ce_ferm * (1 - self.Pox_Ferm)                       # (G, M, T) ferm prods

        # ── Step 5: Reduce (G, M, T) → (G, T) then ravel to (G*T,) Series ───────
        mic_idx  = self.Microbes.index

        TUC  = c_upt_3d.sum(axis=1)   # (G, T) total C uptake
        TUCb = cb_3d.sum(axis=1)      # (G, T) biomass C
        TUE  = e_3d.sum(axis=1)       # (G, T) energy
        TUN  = n_3d.sum(axis=1)       # (G, T) N
        TUP  = p_3d.sum(axis=1)       # (G, T) P
        TU   = TUCb + TUN + TUP       # (G, T) total biomass

        self.Taxon_resp_df = pd.DataFrame(co2_3d.sum(axis=1),    # (G, T) CO2
                                          columns=mic_idx[:T])

        # ── Step 6: Fermentation products → Monomers ─────────────────────────────
        FermProd_C = ferm_3d.sum(axis=(1, 2)).astype(np.float32)       # (G,)
        # Add per-cell fermentation C to each grid cell's FermProd row.
        # Monomers has G*M rows; FermProd is the 3rd monomer (index 2).
        self.Monomers.loc['FermProd', 'C'] += FermProd_C

        # ── Reconstruct DataFrames needed by the rest of uptake() ─────────────────
        # TUC_df etc. are (G, T) DataFrames indexed by grid cell with taxon columns.
        # We build them directly from the numpy arrays to avoid extra copies.
        tax_cols = mic_idx[:T]
        TUC_df  = pd.DataFrame(TUC,  columns=tax_cols)
        TUCb_df = pd.DataFrame(TUCb, columns=tax_cols)
        TUE_df  = pd.DataFrame(TUE,  columns=tax_cols)
        TUN_df  = pd.DataFrame(TUN,  columns=tax_cols)
        TUP_df  = pd.DataFrame(TUP,  columns=tax_cols)
        TU_df   = pd.DataFrame(TU,   columns=tax_cols)

        # Update these global variables
        # self.Taxon_Uptake_C = TUC_df.stack().values.astype(np.float32)  # spatial total C uptake
        # self.Taxon_Uptake_BC = TUCb_df.stack().values  # spatial biomass C uptake: array
        # self.Taxon_Uptake_N = TUN_df.stack().values  # spatial N uptake: array
        # self.Taxon_Uptake_P = TUP_df.stack().values  # spatial P uptake: array
        # self.Taxon_Uptake_E = TUE_df.stack().values.astype(np.float32)  # spatial Energy uptake: array
        self.Taxon_Uptake_C = pd.Series(TUC_df.values.ravel(order='C').astype(np.float32), index=self.Microbes.index)
        self.Taxon_Uptake_BC = pd.Series(TUCb_df.values.ravel(order='C').astype(np.float32), index=self.Microbes.index)
        self.Taxon_Uptake_N = pd.Series(TUN_df.values.ravel(order='C').astype(np.float32), index=self.Microbes.index)
        self.Taxon_Uptake_P = pd.Series(TUP_df.values.ravel(order='C').astype(np.float32), index=self.Microbes.index)
        self.Taxon_Uptake_E = pd.Series(TUE_df.values.ravel(order='C').astype(np.float32), index=self.Microbes.index)
        self.Taxon_Uptake_Total = TU_df.values.ravel(order='C').astype(np.float32)  # total biomass (C, N, & P)

        # Taxon-specific C assimilation efficiency (CAE)
        self.CAE_taxon = self.Microbes['C'].copy()  # create a dataframe and set all vals to 0
        self.CAE_taxon[:] = 0
        pos_uptake_index = self.Taxon_Uptake_C > 0
        self.CAE_taxon[pos_uptake_index] = self.Taxon_Uptake_BC[pos_uptake_index] / self.Taxon_Uptake_C[
            pos_uptake_index]

        # System-level CAE
        Taxon_Uptake_C_grid = self.Taxon_Uptake_C.sum()  # Total C Uptake
        Taxon_Uptake_BC_grid = self.Taxon_Uptake_BC.sum()
        if Taxon_Uptake_C_grid == 0:
            self.CAE_system = np.float32(0)
        else:
            self.CAE_system = Taxon_Uptake_BC_grid / Taxon_Uptake_C_grid

        ## Respiration from Constitutive + Inducible(NOTE: missing sum(MicLoss[,"C"]) in the Mortality below)
        # self.Respiration = (
        #        Taxon_Transporter_Maint + Taxon_Growth_Respiration + Taxon_Osmo_Consti_Maint + Taxon_Osmo_Induci_Maint + Taxon_Enzyme_Consti_Maint + Taxon_Enzyme_Induci_Maint).sum(
        #    axis=0)

        # Carbon uptake
        # self.Taxon_C_Uptake = TUC_df
        self.Taxon_C_Uptake = self.CAE_taxon.copy()
        self.Taxon_C_Uptake[:] = self.Taxon_Uptake_C
        self.Taxon_BC_Uptake = self.CAE_taxon.copy()
        self.Taxon_BC_Uptake[:] = self.Taxon_Uptake_BC
        # self.Taxon_BC_Uptake = TUCb_df
        self.Taxon_E_Uptake = pd.Series(
            TUE_df.values.ravel(order='C').astype(np.float32),
            index=self.Microbes.index
        )
        self.Taxon_E_Prod = self.CAE_taxon.copy()
        self.Taxon_E_Prod[:] = self.Taxon_Uptake_E
        self.Taxon_Total_Uptake = pd.Series(
            TU_df.values.ravel(order='C').astype(np.float32),
            index=self.Microbes.index
        )

    def metabolism(self, day):
        """
        Explicitly calculate intra-cellular production of metabolites.

        Handles both constitutive (standing biomass) and inducible (immediate monomers uptake) pathways following:
          1. constitutive enzyme and osmolyte production
          2. inducible enzyme and osmolyte production
          3. emergent CUE & Respiration
          4. update both Enzymes (with production & loss) and Substrates (with dead enzymes)
        """

        # index of dead enzyme in Substrates
        is_deadEnz = self.Substrates.index == "DeadEnz"

        # Scalar of water potential impact: call the function microbe_osmo_psi()
        f_psi = microbe_osmo_psi(self.psi[day], self.alpha, self.wp_fc, self.wp_th)

        # taxon-specific energy cost of active monomer transport (mass of assimilated monomers* active uptake cost)
        self.Tax_active_transport_E_cost = self.Taxon_Total_Uptake * self.active_uptake_E_cost

        # calculate the proportion of microbial biomass that is HSPs
        # HSP concentration in microbial cells
        hsp_cols = [c for c in self.Taxon_HSP.columns if c != 'Total_biomass']
        self.Taxon_HSP[hsp_cols] *= (1 - self.HSP_Turnover)
        hsp_mass = self.Taxon_HSP[hsp_cols].sum(axis=1)  # total HSP mass per microbe
        HSP_prop = hsp_mass.copy()
        HSP_prop[self.Taxon_HSP['Total_biomass'] > 0] = (
                hsp_mass[self.Taxon_HSP['Total_biomass'] > 0] /
                self.Taxon_HSP['Total_biomass'][self.Taxon_HSP['Total_biomass'] > 0]
        )

        # Calculate temperature sensitivity of maintenance energy cost from HSP concentration
        maint_E_cost_Ea = ((HSP_prop / (
                    HSP_prop + self.maint_E_HSP_Ea_n)) * self.maint_E_HSP_Ea_slope) + self.maint_E_HSP_Ea_int
        # calculate temperature dependent maintenance energy cost
        tax_maint_E_cost = self.maint_E_rate_const * Arrhenius(maint_E_cost_Ea, self.gas_const, self.temp[day],
                                                               self.maint_E_rate_Tref)
        self.Tax_cell_maint_E_cost = tax_maint_E_cost.mul(self.Microbes[['C', 'N', 'P']].sum(axis=1), axis=0)

        # Transporter production
        # Taxon Energy cost of transporter production (inducible transporter production)
        Tax_transport_prod_E_cost = (self.Taxon_Total_Uptake * self.Uptake_Prod_Cost.sum(
            axis=1)) * self.UptakeProd_E_Cost

        self.Tax_cell_maint_E_cost.index = self.Microbes.index
        # calculate the pool of energy available for biomass production after accounting for energy produced,
        # and energy spent on active monomer transport, production of new transporters, and cellular maintenance.:
        MicBio_avail_E = (
                self.Microbes.loc[:, 'E']
                + self.Taxon_E_Uptake
                - self.Tax_active_transport_E_cost
                - Tax_transport_prod_E_cost
                - self.Tax_cell_maint_E_cost.squeeze()  # <-- force to Series
        )
        MicBio_avail_E = MicBio_avail_E.squeeze()  # <-- ensure result is Series
        MicBio_avail_E.loc[MicBio_avail_E < 0] = 0

        # ---------------------------------------------------------------------#
        # ......................constitutive processes.........................#
        # ---------------------------------------------------------------------#
        # Acronyms:
        #   OEHCCN  : Total constitutive N cost (HSP + Osmo + Enzyme)
        #   OEHCCE  : Total constitutive E cost (HSP + Osmo + Enzyme)
        #   MNAOEHC : Min of N cost vs N available
        #   MEAOEHC : Min of E cost vs E available
        #   ARRNOEHC: Ratio of available to required N
        #   ARREOEHC: Ratio of available to required E
        #   ARROEHC : Most limiting ratio (N or E)

        # HSP before adjustment
        HSP_consti_prod = self.HSP_consti_pre_exp * Arrhenius(
            self.HSP_consti_Ea, self.gas_const, self.temp[day], self.maint_E_rate_Tref)
        tax_HSP_consti_prod = HSP_consti_prod.mul(self.Microbes[['C', 'N', 'P']].sum(axis=1), axis=0)
        Taxon_HSP_Consti_Cost_N = (tax_HSP_consti_prod * self.HSP_N_cost).sum(axis=1)
        Taxon_HSP_Consti_Cost_E = (tax_HSP_consti_prod * self.HSP_E_cost).sum(axis=1)

        # Osmolyte before adjustment
        Taxon_Osmo_Consti = self.Consti_Osmo_C.mul(self.Microbes['C'], axis=0) * f_psi
        Taxon_Osmo_Consti_Cost_N = (Taxon_Osmo_Consti * self.Osmo_N_cost).sum(axis=1)
        Taxon_Osmo_Consti_Cost_E = (Taxon_Osmo_Consti * self.Osmo_E_cost).sum(axis=1)

        # Enzyme before adjustment
        Taxon_Enzyme_Consti = self.Consti_Enzyme_C.mul(self.Microbes['C'], axis=0)
        Taxon_Enzyme_Consti_Cost_N = (Taxon_Enzyme_Consti.mul(self.Enz_Attrib['N_cost'], axis=1)).sum(axis=1)
        Taxon_Enzyme_Consti_Cost_E = (Taxon_Enzyme_Consti.mul(self.Enz_Attrib['E_cost'], axis=1)).sum(axis=1)

        # N limitation across all three constitutive products
        OEHCCN = Taxon_HSP_Consti_Cost_N + Taxon_Osmo_Consti_Cost_N + Taxon_Enzyme_Consti_Cost_N
        MNAOEHC = pd.concat(
            [OEHCCN[OEHCCN > 0], self.Microbes['N'][OEHCCN > 0]], axis=1
        ).min(axis=1, skipna=True)
        ARRNOEHC = (MNAOEHC / OEHCCN[OEHCCN > 0]).fillna(0)

        # E limitation across all three constitutive products
        OEHCCE = Taxon_HSP_Consti_Cost_E + Taxon_Osmo_Consti_Cost_E + Taxon_Enzyme_Consti_Cost_E
        _slice_e = OEHCCE[OEHCCE > 0]
        _slice_m = MicBio_avail_E[OEHCCE > 0]
        MEAOEHC = pd.concat([_slice_e, _slice_m], axis=1).min(axis=1, skipna=True)
        ARREOEHC = (MEAOEHC / OEHCCE[OEHCCE > 0]).fillna(0)

        # Most limiting factor
        ARROEHC = ARREOEHC.copy()
        ARROEHC.loc[ARREOEHC > ARRNOEHC] = ARRNOEHC.loc[ARREOEHC > ARRNOEHC]

        # HSP adjusted
        tax_HSP_consti_prod[OEHCCN > 0] = tax_HSP_consti_prod[OEHCCN > 0].mul(ARROEHC, axis=0)
        Taxon_HSP_Consti_Cost_N = (tax_HSP_consti_prod * self.HSP_N_cost).sum(axis=1)
        Taxon_HSP_Consti_Cost_C = tax_HSP_consti_prod.sum(axis=1)
        self.Taxon_HSP_Consti_Cost_E = (tax_HSP_consti_prod * self.HSP_E_cost).sum(axis=1)

        # Osmolyte adjusted
        Taxon_Osmo_Consti[OEHCCN > 0] = Taxon_Osmo_Consti[OEHCCN > 0].mul(ARROEHC, axis=0)
        Taxon_Osmo_Consti_Cost_N = (Taxon_Osmo_Consti * self.Osmo_N_cost).sum(axis=1)
        Taxon_Osmo_Consti_Cost_C = Taxon_Osmo_Consti.sum(axis=1)
        Taxon_Osmo_Consti_Cost_E = Taxon_Osmo_Consti_Cost_C * self.Osmo_E_cost

        # Enzyme adjusted
        Taxon_Enzyme_Consti[OEHCCN > 0] = Taxon_Enzyme_Consti[OEHCCN > 0].mul(ARROEHC, axis=0)
        Taxon_Enzyme_Consti_Cost_N = (Taxon_Enzyme_Consti.mul(self.Enz_Attrib['N_cost'], axis=1)).sum(axis=1)
        Taxon_Enzyme_Consti_Cost_P = (Taxon_Enzyme_Consti.mul(self.Enz_Attrib['P_cost'], axis=1)).sum(axis=1)
        Taxon_Enzyme_Consti_Cost_C = Taxon_Enzyme_Consti.sum(axis=1)
        Taxon_Enzyme_Consti_Cost_E = (Taxon_Enzyme_Consti.mul(self.Enz_Attrib['E_cost'], axis=1)).sum(axis=1)

        # Deduct all constitutive E costs from available pool
        MicBio_avail_E -= self.Taxon_HSP_Consti_Cost_E + Taxon_Osmo_Consti_Cost_E + Taxon_Enzyme_Consti_Cost_E
        MicBio_avail_E.loc[MicBio_avail_E < 0] = 0

        # --- Inducible production (HSP, Osmolyte, Enzyme) before resource adjustment ---

        # Inducible HSP: scales with total C uptake, temperature-dependent rate
        tax_HSP_induci_prod_rate = self.HSP_induci_pre_exp * Arrhenius(
            self.HSP_induci_Ea, self.gas_const, self.temp[day], self.maint_E_rate_Tref)
        tax_HSP_induci_prod_rate.index = self.Microbes.index
        tax_HSP_induci_prod = tax_HSP_induci_prod_rate.mul(self.Taxon_Total_Uptake, axis=0)
        Taxon_HSP_Induci_Cost_N = (tax_HSP_induci_prod * self.HSP_N_cost).sum(axis=1)
        Taxon_HSP_Induci_Cost_E = (tax_HSP_induci_prod * self.HSP_E_cost).sum(axis=1)

        # Inducible Osmolyte: scales with biomass C uptake, moisture-dependent
        Taxon_Osmo_Induci = self.Induci_Osmo_C.mul(self.Taxon_Uptake_BC, axis=0) * f_psi
        Taxon_Osmo_Induci_Cost_N = (Taxon_Osmo_Induci * self.Osmo_N_cost).sum(axis=1)
        Taxon_Osmo_Induci_Cost_E = (Taxon_Osmo_Induci * self.Osmo_E_cost).sum(axis=1)

        # Inducible Enzyme: scales with biomass C uptake
        Taxon_Enzyme_Induci = self.Induci_Enzyme_C.mul(self.Taxon_Uptake_BC, axis=0)
        Taxon_Enzyme_Induci_Cost_N = (Taxon_Enzyme_Induci.mul(self.Enz_Attrib['N_cost'], axis=1)).sum(axis=1)
        Taxon_Enzyme_Induci_Cost_E = (Taxon_Enzyme_Induci.mul(self.Enz_Attrib['E_cost'], axis=1)).sum(axis=1)

        # --- N limitation: HSP + Osmolyte + Enzyme inducible ---
        # Acronyms:
        #   OEIHCN : Total inducible N cost (Osmo + Enzyme + HSP)
        #   OEIHAN : N available from uptake
        #   MNAOEIH: Min of N cost vs N available
        #   ARRNOEI: Ratio of available to required N
        OEIHCN = Taxon_HSP_Induci_Cost_N + Taxon_Osmo_Induci_Cost_N + Taxon_Enzyme_Induci_Cost_N
        OEIHAN = self.Taxon_Uptake_N  # already a Series with Microbes.index after your fix
        MNAOEIH = pd.concat(
            [OEIHCN[OEIHCN > 0], OEIHAN[OEIHCN > 0]], axis=1
        ).min(axis=1, skipna=True)
        ARRNOEI = (MNAOEIH / OEIHCN[OEIHCN > 0]).fillna(0)

        # --- E limitation: HSP + Osmolyte + Enzyme inducible ---
        # Acronyms:
        #   OEIHCE : Total inducible E cost (Osmo + Enzyme + HSP)
        #   MEAOEIH: Min of E cost vs E available
        #   ARREOEI: Ratio of available to required E
        OEIHCE = Taxon_HSP_Induci_Cost_E + Taxon_Osmo_Induci_Cost_E + Taxon_Enzyme_Induci_Cost_E
        MEAOEIH = pd.concat(
            [OEIHCE[OEIHCE > 0], MicBio_avail_E[OEIHCE > 0]], axis=1
        ).min(axis=1, skipna=True)
        ARREOEI = (MEAOEIH / OEIHCE[OEIHCE > 0]).fillna(0)

        # Most limiting factor across N and E
        ARROEI = ARREOEI.copy()
        ARROEI.loc[ARREOEI > ARRNOEI] = ARRNOEI.loc[ARREOEI > ARRNOEI]

        # --- Apply adjustment to all three inducible products ---

        # HSP adjusted
        tax_HSP_induci_prod[OEIHCN > 0] = tax_HSP_induci_prod[OEIHCN > 0].mul(ARROEI, axis=0)
        Taxon_HSP_Induci_Cost_N = (tax_HSP_induci_prod * self.HSP_N_cost).sum(axis=1)
        Taxon_HSP_Induci_Cost_C = tax_HSP_induci_prod.sum(axis=1)
        self.Taxon_HSP_Induci_Cost_E = (tax_HSP_induci_prod * self.HSP_E_cost).sum(axis=1)

        # Osmolyte adjusted
        Taxon_Osmo_Induci[OEIHCN > 0] = Taxon_Osmo_Induci.loc[OEIHCN > 0].mul(ARROEI, axis=0)
        Taxon_Osmo_Induci_Cost_N = (Taxon_Osmo_Induci * self.Osmo_N_cost).sum(axis=1)
        Taxon_Osmo_Induci_Cost_C = Taxon_Osmo_Induci.sum(axis=1)
        Taxon_Osmo_Induci_Cost_E = (Taxon_Osmo_Induci * self.Osmo_E_cost).sum(axis=1)

        # Enzyme adjusted
        Taxon_Enzyme_Induci[OEIHCN > 0] = Taxon_Enzyme_Induci.loc[OEIHCN > 0].mul(ARROEI, axis=0)
        Taxon_Enzyme_Induci_Cost_N = (Taxon_Enzyme_Induci.mul(self.Enz_Attrib["N_cost"], axis=1)).sum(axis=1)
        Taxon_Enzyme_Induci_Cost_P = (Taxon_Enzyme_Induci.mul(self.Enz_Attrib["P_cost"], axis=1)).sum(axis=1)
        Taxon_Enzyme_Induci_Cost_C = Taxon_Enzyme_Induci.sum(axis=1)
        Taxon_Enzyme_Induci_Cost_E = (Taxon_Enzyme_Induci.mul(self.Enz_Attrib["E_cost"], axis=1)).sum(axis=1)

        # Derive C, N, & P deposited as biomass from Uptake; ensure no negative values
        # HSP and osmolyte C/N are NOT subtracted here: both are produced within the cell and
        # remain part of microbial biomass (unlike enzymes which are exuded). Only enzyme
        # inducible costs (C/N) are deducted as those molecules leave the cell.
        Microbe_C_Gain = self.Taxon_Uptake_BC - Taxon_Enzyme_Induci_Cost_C
        Microbe_N_Gain = self.Taxon_Uptake_N - Taxon_Enzyme_Induci_Cost_N
        Microbe_P_Gain = self.Taxon_Uptake_P - Taxon_Enzyme_Induci_Cost_P

        # Derive total E costs and gains
        self.Microbe_E_Cost = (self.Tax_active_transport_E_cost + Tax_transport_prod_E_cost
                               + self.Tax_cell_maint_E_cost + self.Taxon_HSP_Consti_Cost_E
                               + self.Taxon_HSP_Induci_Cost_E + Taxon_Enzyme_Induci_Cost_E
                               + Taxon_Osmo_Induci_Cost_E + Taxon_Enzyme_Consti_Cost_E
                               + Taxon_Osmo_Consti_Cost_E)
        self.Microbe_E_Gain = self.Taxon_Uptake_E - self.Microbe_E_Cost

        # ---------------------------------------------------------------------#
        # ...............................Integration...........................#
        # ---------------------------------------------------------------------#
        # Update Microbial pools with GAINS (from uptake) and LOSSES (from constitutive production)
        # HSP and osmolyte C/N costs are excluded: both remain within the cell as part of
        # microbial biomass. Only enzyme constitutive costs (C/N) are deducted (exuded).
        self.Microbes.loc[
            :, 'C'] += Microbe_C_Gain - Taxon_Enzyme_Consti_Cost_C
        self.Microbes.loc[
            :, 'N'] += Microbe_N_Gain - Taxon_Enzyme_Consti_Cost_N
        self.Microbes.loc[:, 'P'] += Microbe_P_Gain - Taxon_Enzyme_Consti_Cost_P
        # Update Microbial E pool
        self.Microbes.loc[:, 'E'] += self.Microbe_E_Gain
        self.Microbes[self.Microbes < 0] = np.float32(0)  # avoid negative values

        # Taxon-specific emergent CUE (accounts for C lost through respiration and enzyme exudation)
        # HSP and osmolyte C are NOT subtracted: both remain in the cell as part of biomass.
        # Only enzyme C (constitutive + inducible, already removed from Microbe_C_Gain) reduces CUE.
        self.CUE_taxon = self.Microbes['C'].copy()  # create a dataframe and set all vals to 0
        self.CUE_taxon[:] = 0
        pos_uptake_index = self.Taxon_Uptake_C > 0
        Mic_BC_Use = Microbe_C_Gain - Taxon_Enzyme_Consti_Cost_C
        self.CUE_taxon[pos_uptake_index] = Mic_BC_Use[pos_uptake_index] / self.Taxon_Uptake_C[pos_uptake_index]

        # System-level emergent CUE
        Taxon_Use_C_grid = self.Taxon_Uptake_C.sum()  # Total C Uptake
        if Taxon_Use_C_grid == 0:
            self.CUE_system = np.float32(0)
        else:
            self.CUE_system = Mic_BC_Use.sum() / Taxon_Use_C_grid

        # Heat stress protein production
        self.Taxon_HSP[tax_HSP_consti_prod.columns] += tax_HSP_consti_prod + tax_HSP_induci_prod
        self.Taxon_HSP['Total_biomass'] = self.Microbes[['C', 'N', 'P']].sum(axis=1).values

        # Osmolyte pool: accumulate this step's production (constitutive + inducible),
        # then apply fractional turnover loss. Osmolytes remain within the cell (like HSP)
        # so turnover represents intracellular degradation, not exudation.
        self.Taxon_Osmolyte += Taxon_Osmo_Consti + Taxon_Osmo_Induci
        self.Taxon_Osmolyte *= np.float32(1) - self.OSmo_Loss_Rate
        # Track total osmolyte C production this step (summed across gene columns) for output
        self.Taxon_Osmo_Cost_C = (Taxon_Osmo_Consti + Taxon_Osmo_Induci).sum(axis=1)
        # Track per-taxon enzyme and HSP C costs (consti + induci) for per-microbe output
        self.Taxon_Enzyme_Cost_C = Taxon_Enzyme_Consti_Cost_C + Taxon_Enzyme_Induci_Cost_C
        self.Taxon_HSP_Cost_C        = Taxon_HSP_Consti_Cost_C + Taxon_HSP_Induci_Cost_C
        self.Taxon_Transport_Cost_C  = Tax_transport_prod_E_cost / self.UptakeProd_E_Cost
        # Track per-pathway E costs (per taxon×gridcell Series) for per-microbe output
        self.Taxon_Transport_E_Cost  = (self.Tax_active_transport_E_cost
                                        + Tax_transport_prod_E_cost)
        self.Taxon_Maint_E_Cost      = self.Tax_cell_maint_E_cost
        self.Taxon_Enzyme_E_Cost     = Taxon_Enzyme_Consti_Cost_E + Taxon_Enzyme_Induci_Cost_E
        self.Taxon_Osmo_E_Cost_total = Taxon_Osmo_Consti_Cost_E  + Taxon_Osmo_Induci_Cost_E
        self.Taxon_HSP_E_Cost        = (self.Taxon_HSP_Consti_Cost_E
                                        + self.Taxon_HSP_Induci_Cost_E)

        # Derive Enzyme production
        self.Taxon_Enzyme_Production = Taxon_Enzyme_Consti + Taxon_Enzyme_Induci  # gene-specific prod of enzyme of each taxon: (taxon*gridsize) * enzyme
        self.Taxon_Enzyme_Production.index = [np.arange(self.gridsize).repeat(self.n_taxa),
                                              self.Taxon_Enzyme_Production.index]  # create a multi-index
        EP_df = pd.DataFrame(
            self.Taxon_Enzyme_Production.values.reshape(self.gridsize, self.n_taxa, self.n_enzymes).sum(axis=1),
            index=np.arange(self.gridsize), columns=self.Taxon_Enzyme_Production.columns
        )  # enzyme-specific production in each grid cell
        # self.Enzyme_Production = EP_df.stack().values  # 1-D array
        self.Enzyme_Production = EP_df.stack()  # 1-D array
        self.Enzyme_Production.index = self.Enzymes.index

        # Derive Enzyme turnover
        self.Enzyme_Loss = self.Enzymes * self.Enzyme_Loss_Rate

        # Update Enzyme pools by adding enzymes produced and substracting the 'dead' enzymes
        self.Enzymes += self.Enzyme_Production - self.Enzyme_Loss
        self.Enzymes_gain = self.Enzyme_Production - self.Enzyme_Loss

        # Update Substrates pools with dead enzymes
        DeadEnz_df = pd.concat(
            [self.Enzyme_Loss, self.Enzyme_Loss.mul(self.Enz_Attrib['N_cost'].tolist() * self.gridsize, axis=0),
             self.Enzyme_Loss.mul(self.Enz_Attrib['P_cost'].tolist() * self.gridsize, axis=0)], axis=1)
        DeadEnz_df.index = [np.arange(self.gridsize).repeat(self.n_enzymes), DeadEnz_df.index]  # create a multi-index
        DeadEnz_gridcell = pd.DataFrame(
            DeadEnz_df.values.reshape(self.gridsize, self.n_enzymes, 3).sum(axis=1),
            index=np.arange(self.gridsize), columns=DeadEnz_df.columns
        )  # total dead enzyme mass per grid cell
        self.Substrates.loc[is_deadEnz] += DeadEnz_gridcell.values.astype(np.float32)


    def mortality(self, day):
        """
        Calculate microbial mortality and update stoichiometry of the alive and microbial pools.

        Kill microbes that are starving deterministically and microbes that are drought intolerant stochastically
        Also update Substrates with input from dead microbes, monomers(with leaching loss), and respiration
        """

        # Constants
        Leaching = np.float32(0.1)  # Abiotic monomer loss rate
        Psi_slope_leach = np.float32(0.5)  # Mositure sensivity of abiotic monomer loss rate
        # Indices
        Mic_index = self.Microbes.index
        is_DeadMic = self.Substrates.index == "DeadMic"
        is_NH4 = self.Monomers.index == "NH4"
        is_PO4 = self.Monomers.index == "PO4"

        # Reset the index to arabic numerals from taxa series
        self.Microbes = self.Microbes.reset_index(drop=True)
        MinRatios = self.MinRatios.reset_index(drop=True)

        # Create a blank dataframe, Death, having the same structure as Microbes
        Death = self.Microbes.copy(deep=True)  # .drop('E',axis=1)
        Death[:] = np.float32(0)
        # Create a series, kill, holding boolean value of False
        kill = pd.Series([False] * self.n_taxa * self.gridsize)

        # Start of calcualtion of mortality first with THRESHOLD
        # Kill microbes deterministically based on threshold values: C_min: 0.086; N_min:0.012; P_min: 0.002
        starve_index = (self.Microbes["C"] > 0) & (
                (self.Microbes["C"] < self.C_min) | (self.Microbes["N"] < self.N_min) | (
                self.Microbes["P"] < self.P_min))
        # Index the dead, put them in Death, and set them to 0 in Microbes
        Death.loc[starve_index] = self.Microbes[starve_index]
        self.Microbes.loc[starve_index] = np.float32(0)
        # Index the locations where microbial cells remain alive
        mic_index = self.Microbes["C"] > 0

        # Kill microbes stochastically based on mortality prob as a function of water potential and drought tolerance
        # call the function, MMP:microbe_mortality_psi()
        r_death = MMP(self.psi[day], self.wp_fc, self.basal_death_prob, self.death_rate, self.tolerance)
        kill.loc[mic_index] = r_death[mic_index] > np.random.uniform(0, 1, sum(mic_index)).astype('float32')
        # Index the dead, put them in Death, and set them to 0 in Microbes
        Death.loc[kill] = self.Microbes[kill]
        self.Microbes.loc[kill] = np.float32(0)
        # Index locations where microbes remain alive
        mic_index = self.Microbes['C'] > 0

        # Calculate the total dead mass (threshold & drought) across taxa in each grid cell
        # Death_gridcell = Death.groupby(Death.index // self.n_taxa).sum(axis=0)
        # Death_gridcell = Death.groupby(Death.index // self.n_taxa).sum()
        Death_gridcell = pd.DataFrame(
            Death.drop('E', axis=1).values.reshape(self.gridsize, self.n_taxa, 3).sum(axis=1))

        # Distinguish between conditions of complete death VS partial death
        # All cells die
        if sum(mic_index) == 0:

            # ...Update Substrates pool by adding dead microbial biomass
            self.Substrates.loc[is_DeadMic] += Death_gridcell.values

        # Partly die and adjust stoichiometry of those remaining alive
        else:

            # Index only those taxa in Microbes that have below-minimum quotas: Mic_subset
            MicrobeRatios = self.Microbes[mic_index].divide(self.Microbes[mic_index].drop('E', axis=1).sum(axis=1),
                                                            axis=0)
            mic_index_sub = (MicrobeRatios["C"] < MinRatios[mic_index]["C"]) | (
                    MicrobeRatios["N"] < MinRatios[mic_index]["N"]) | (
                                    MicrobeRatios["P"] < MinRatios[mic_index]["P"])
            rat_index = self.Microbes.index.map(mic_index_sub).fillna(False)
            # Derive the Microbes wanted
            Mic_subset = self.Microbes[rat_index]
            StartMicrobes = self.Microbes.copy(deep=True)[rat_index]

            # Derive new ratios and Calculate difference between actual and min ratios
            MicrobeRatios = Mic_subset.divide(Mic_subset.sum(axis=1), axis=0)
            MinRat = MinRatios[rat_index]
            Ratio_dif = MicrobeRatios - MinRat
            # Create a df recording the ratio differences < 0
            Ratio_dif_0 = Ratio_dif.copy(deep=True)
            Ratio_dif_0[Ratio_dif > 0] = np.float32(0)
            # Create a df recording the ratio differences > 0
            Excess = Ratio_dif.copy(deep=True)
            Excess[Ratio_dif < 0] = np.float32(0)

            # Determine the limiting nutrient that will be conserved
            Limiting = (-Ratio_dif / MinRat).idxmax(
                axis=1)  # Series of index of the first occurrence of maximum in each row
            # Set all deficient ratios to their minima
            MicrobeRatios[Ratio_dif < 0] = MinRat[Ratio_dif < 0]
            # Reduce the mass fractions for non-deficient elements in proportion to the distance from the minimum
            # ....Partition the total deficit to the excess element(s) in proportion to their distances from their minima
            MicrobeRatios[Ratio_dif > 0] += Excess.mul((Ratio_dif_0.sum(axis=1) / Excess.sum(axis=1)), axis=0)[
                Ratio_dif > 0]

            # Construct hypothetical nutrient quotas for each possible minimum nutrient
            MC = Mic_subset["C"]
            MN = Mic_subset["N"]
            MP = Mic_subset["P"]
            MRC = MicrobeRatios["C"]
            MRN = MicrobeRatios["N"]
            MRP = MicrobeRatios["P"]

            new_C = pd.concat([MC, MN * MRC / MRN, MP * MRC / MRP], axis=1)
            new_C = new_C.fillna(0)
            new_C[np.isinf(new_C)] = np.float32(0)
            new_C.columns = ['C', 'N', 'P']

            new_N = pd.concat([MC * MRN / MRC, MN, MP * MRN / MRP], axis=1)
            new_N = new_N.fillna(0)
            new_N[np.isinf(new_N)] = np.float32(0)
            new_N.columns = ['C', 'N', 'P']

            new_P = pd.concat([MC * MRP / MRC, MN * MRP / MRN, MP], axis=1)
            new_P = new_P.fillna(0)
            new_P[np.isinf(new_P)] = np.float32(0)
            new_P.columns = ['C', 'N', 'P']

            # Insert the appropriate set of nutrient quotas scaled to the minimum nutrient
            C = [new_C.loc[i, Limiting[i]] for i in Limiting.index]  # list
            N = [new_N.loc[i, Limiting[i]] for i in Limiting.index]  # list
            P = [new_P.loc[i, Limiting[i]] for i in Limiting.index]  # list

            # Update Microbes
            self.Microbes.loc[rat_index, 'C'] = C
            self.Microbes.loc[rat_index, 'N'] = N
            self.Microbes.loc[rat_index, 'P'] = P

            # Sum up the element losses from biomass across whole grid and calculate average loss
            MicLoss = StartMicrobes - self.Microbes[rat_index]
            tax_resp_mort = self.Microbes['C'].copy(deep=True)
            tax_resp_mort[:] = self.Taxon_resp_df.values.ravel(order='C')
            # tax_resp_mort.index = self.Microbes.index
            tax_resp_mort.loc[rat_index] += MicLoss['C']
            self.Taxon_Resp = tax_resp_mort
            # Update monomer pools
            self.Monomers.loc[is_NH4, "N"] += sum(MicLoss["N"]) / self.gridsize
            self.Monomers.loc[is_PO4, "P"] += sum(MicLoss["P"]) / self.gridsize

            # Update Substrates pool by adding dead microbial biomass
            self.Substrates.loc[is_DeadMic] += Death_gridcell.values
        # End of if else clause

        # Check microbes energy ratio and set to max (10000) if it is higher
        Tax_E_ratio = self.Microbes['E'] / self.Microbes['C']
        # self.Microbes['E'][Tax_E_ratio > self.E_max] = self.Microbes['C'] * self.E_max
        Tax_E_pos_index = self.Microbes.index.map(Tax_E_ratio > self.E_max)
        self.Microbes.loc[Tax_E_pos_index, 'E'] = self.Microbes['C'] * self.E_max

        # Leaching of monomers
        Leaching_N = self.Monomers.loc[is_NH4, "N"] * Leaching * np.exp(Psi_slope_leach * (self.psi[day] - self.wp_fc))
        Leaching_P = self.Monomers.loc[is_PO4, "P"] * Leaching * np.exp(Psi_slope_leach * (self.psi[day] - self.wp_fc))
        # Update Monomers
        self.Monomers.loc[is_NH4, "N"] -= Leaching_N
        self.Monomers.loc[is_PO4, "P"] -= Leaching_P

        # Restore the index to taxa series
        self.Microbes.index = Mic_index
        self.Taxon_Resp.index = Mic_index
        # system respiration
        self.Respiration = self.Taxon_Resp.sum()

        # Update the death toll of cells
        self.Kill = kill.sum().astype('uint32')

    def reproduction(self, day):
        """
        Calculate reproduction and dispersal, and update microbial composition/distribution on the spatial grid.

        Parameters:
            fb         : index of fungal taxa
            max_size_b : threshold of cell division
            max_size_f : threshold of cell division
            x,y        : x,y dimension of grid
            dist       : maximum dispersal distance: 1 cell
            direct     : dispersal direction: 0.95
        """

        # Microbes' index
        Mic_index = self.Microbes.index
        # Set up the Colonization dataframe: taxon * 3(C,N,&P)
        Colonization = self.Microbes.copy(deep=True)
        Colonization = Colonization.reset_index(drop=True)
        Colonization[:] = np.float32(0)

        # STEP 1: Fungal translocation by calculating average biomass within fungal taxa
        # Count the fungal taxa before cell division
        # Add one or two fungi to the count series based on size
        Fungi_df = pd.Series(data=[0] * self.n_taxa * self.gridsize, index=Mic_index, name='Count', dtype='int8')
        Fungi_df.loc[(self.fb == 1) & (self.Microbes['C'] > 0)] = np.int8(1)
        Fungi_df.loc[(self.fb == 1) & (self.Microbes['C'] > self.max_size_f)] = np.int8(2)
        # dtype is intentionally not inherited from Fungi_df (int8) because
        # the sum over gridsize cells can exceed int8 max (127) on large grids.
        Fungi_count = pd.Series(
            Fungi_df.values.reshape(self.gridsize, self.n_taxa).sum(axis=0),
            index=self._taxa_index)
        # Derive average biomass of fungal taxa
        Microbes_grid = pd.DataFrame(
            self.Microbes.values.reshape(self.gridsize, self.n_taxa, 4).sum(axis=0),
            index=self._taxa_index, columns=self.Microbes.columns)
        Mean_fungi = Microbes_grid.divide(Fungi_count, axis=0)
        Mean_fungi[Fungi_count == 0] = np.float32(0)
        # Expand the fungal average across the grid
        eMF = expand(Mean_fungi, self.gridsize).astype('float32')

        # STEP 2: Cell division & translocate nutrients
        MicrobesBeforeDivision = self.Microbes.copy(deep=True)
        # bacterial cell division
        bac_index = (self.fb == 0) & (self.Microbes['C'] > self.max_size_b)
        self.Microbes[bac_index] = self.Microbes[bac_index] / 2
        # fungal cell division
        fun_index = (self.fb == 1) & (self.Microbes['C'] > self.max_size_f)
        self.Microbes[fun_index] = self.Microbes[fun_index] / 2
        # Put daughter cells into a seperate dataframe, Reprod
        Reprod = MicrobesBeforeDivision - self.Microbes
        # Translocate nutrients within fungal taxa after reproduction
        self.Microbes[(self.fb == 1) & (self.Microbes['C'] > 0)] = eMF[(self.fb == 1) & (self.Microbes['C'] > 0)]
        # Index the daughter cells of fungi vs bacteria
        daughters_b = (Reprod['C'] > 0) & (self.fb == 0)
        daughters_f = (Reprod['C'] > 0) & (self.fb == 1)
        # set all fungi equal to their grid averages for translocation before colonization
        Reprod[daughters_f] = eMF[daughters_f]

        # STEP 3: dispersal calculation
        num_b = sum(daughters_b)
        num_f = sum(daughters_f)
        shift_x = pd.Series(data=[0] * self.gridsize * self.n_taxa, index=Mic_index, dtype='int8')
        shift_y = pd.Series(data=[0] * self.gridsize * self.n_taxa, index=Mic_index, dtype='int8')
        # Bacterial dispersal movements in X & Y direction
        shift_x.loc[daughters_b] = np.random.choice([i for i in range(-self.dist, self.dist + 1)], num_b,
                                                    replace=True).astype('int8')
        shift_y.loc[daughters_b] = np.random.choice([i for i in range(-self.dist, self.dist + 1)], num_b,
                                                    replace=True).astype('int8')
        # Fungi always move positively in x direction, and in y direction constrained to one box away determined by probability "direct"
        shift_x.loc[daughters_f] = np.int8(1)
        shift_y.loc[daughters_f] = np.random.choice([-1, 0, 1], num_f, replace=True,
                                                    p=[0.5 * (1 - self.direct), self.direct,
                                                       0.5 * (1 - self.direct)]).astype('int8')
        # Calculate x,y coordinates of dispersal destinations (% remainder of x/x)
        new_x = (list(np.repeat(range(1, self.x + 1), self.n_taxa)) * self.y + shift_x + self.x) % self.x
        new_y = (list(np.repeat(range(1, self.y + 1), self.n_taxa * self.x)) + shift_y + self.y) % self.y
        new_x[new_x == 0] = self.x  # Substitute coordinates when there is no shift
        new_y[new_y == 0] = self.y  # Substitute coordinates when there is no shift
        # Convert x,y coordinates to a Series of destination locations
        index_series = (self.n_taxa * ((new_y - 1) * self.x + (new_x - 1))) + list(
            range(1, self.n_taxa + 1)) * self.gridsize - 1

        # Step 4: colonization of dispersed microbes
        # Use np.add.at to accumulate biomass, correctly handling multiple daughters landing on the same grid location.
        # The previous assignment (iloc[...] = ...) silently overwrote earlier daughters at colliding destinations.
        col_arr = Colonization.to_numpy(copy=True)
        np.add.at(col_arr, index_series[daughters_b].to_numpy(), Reprod[daughters_b].to_numpy())
        np.add.at(col_arr, index_series[daughters_f].to_numpy(), Reprod[daughters_f].to_numpy())
        Colonization[:] = col_arr
        # Colonization of dispersing microbes
        self.Microbes += Colonization.values

    def repopulation(self, output, day, mic_reinit):
        """
        Reinitialize a microbial community on the spatial grid in a new pulse.

        Meanwhile, start with new substrates, monomers, and enzymes on the grid that are initialized from the very beginning

        Parameters:
            output:     an instance of the Output class, from which the var, MicrobesSeries_repop,
                        referring to taxon-specific total mass over the grid is retrieved
            pulse:      the pulse index
            day:        the day index
            mic_reinit: 0/1; 1 means reinitialization
        Returns:
            update Substrates, Monomers, and Microbes, as well as Enzymes
        """

        # reinitialize substrates, monomers, and enzymes in a new pulse
        self.Substrates = self.Substrates_init.copy(deep=True)
        self.Monomers = self.Monomers_init.copy(deep=True)
        self.Enzymes = self.Enzymes_init.copy(deep=True)

        # reinitialize microbial community in a new pulse if True
        if mic_reinit:

            self.Microbes = self.Microbes_init.copy(deep=True)  # NOTE copy()
            # fb = self.fb[0:self.n_taxa]
            # max_size_b = self.max_size_b
            # max_size_f = self.max_size_f

            # cumulative abundance; note the column index
            # option 1: mass-based.
            # cum_abundance = output.MicrobesSeries_repop.iloc[:,(day+2-self.cycle):(day+2)].sum(axis=1)
            # option 2: abundance-based
            cum_abundance = output.Taxon_count_repop.iloc[:, (day + 2 - self.cycle):(day + 2)].sum(axis=1)

            # account for different cell mass sizes of bacteria and fungi
            # if sum(fb==1) == 0: # no fungal taxa
            #    frequencies = cum_abundance/cum_abundance.sum()
            # else:
            #    cum_abundance.loc[fb==1] = cum_abundance[fb==1]*max_size_b/max_size_f
            #    frequencies = cum_abundance/cum_abundance.sum()

            # Switched to taxon abundance-based, so no more adjustments
            frequencies = cum_abundance / cum_abundance.sum()
            frequencies = frequencies.fillna(0)
            probs = pd.concat([frequencies, 1 - frequencies], axis=1, sort=False)
            # Randomly assign microbes to each grid box based on prior densities
            choose_taxa = np.array([0] * self.gridsize * self.n_taxa).reshape(self.n_taxa, self.gridsize)
            for i in range(self.n_taxa):
                # Alternative 1
                choose_taxa[i, :] = np.random.choice([1, 0], self.gridsize, replace=True, p=probs.iloc[i, :])

            # Note order='F'
            self.Microbes.loc[np.ravel(choose_taxa, order='F') == 0] = np.float32(0)

        # Record the post-repopulation microbial C for mass balance tracking.
        # next_pulse is the 1-based index of the pulse about to begin.
        next_pulse = (day // self.cycle) + 2  # day is last day of current pulse (0-based)
        output.record_pulse_init(self, next_pulse)
