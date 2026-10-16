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

def _limit_ratio(*pairs):
    """Joint resource limitation: for (cost, available) pairs of per-row arrays, return the per-row factor in [0, 1]
    by which production must be scaled so that no cost exceeds what is available (rows without cost: 1)."""
    r = None
    for cost, avail in pairs:
        cost = np.asarray(cost, dtype=np.float32)
        avail = np.maximum(np.asarray(avail, dtype=np.float32), np.float32(0))
        rk = np.ones_like(cost)
        pos = cost > 0
        rk[pos] = np.minimum(avail[pos] / cost[pos], np.float32(1))
        r = rk if r is None else np.minimum(r, rk)
    return r


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
        self.Vmax0 = data_init['Vmax0']  # Max. reaction speed # not used in MMRT model (Jeth)
        self.Km0 = data_init['Km0']  # Half-saturation constant # not used in MMRT model (Jeth)
        self.SubstrateRatios = np.float32('nan')  # Substrate stoichiometry
        self.DecayRates = np.float32('nan')  # Substrate decay rate
        # Vmax-Km relationship (both kinetics models): Km = Vmax x slope x error factor + intercept, at least Km_min
        self.Vmax_Km_slope = np.float32(data_init['Vmax_Km'])
        self.Vmax_Km_int = np.float32(data_init['Vmax_Km_int'])
        self.Km_min = np.float32(data_init['Km_min'])
        self.Km_factor = data_init['Km_factor']  # (G*M, NE) error factor, or None (Km_error = 0)
        self.Enz_Vmax_sub_spec = data_init[
            'Enz_vmax_sub_spec']  # substrate specificity-efficiency adjustment for degradation Vmax
        self.Enz_Km_Temp_Sens = data_init[
            'Enz_km_temp_sens']  # 1: Km follows the daily temperature; 0: Km is evaluated at Km_Tref
        self.Arrh_Tref = data_init['Arrhenius_tref']

        # Uptake
        self.Microbes_init = data_init['Microbes_pp']  # microbial community before placement
        self.Microbes = data_init['Microbes'].copy(deep=True)  # microbial community after placement
        self.Monomers_init = data_init['Monomers']  # Monomers initialized
        self.Monomers = data_init['Monomers'].copy(deep=True)  # Monomers
        self.MonInput = data_init['MonInput']  # Inputs of monomers
        # Ongoing mineral inputs (Input_NH4 as N and Input_PO4 as P per grid cell per day, from parameters.csv); the inputs of the
        # organic monomers and FermProd stay inactive. Zero inputs leave the model unchanged (re-enabled; see uptake()).
        _mi = self.MonInput
        self._nh4_input = np.float32(_mi[_mi.index == "NH4"].iloc[0])
        self._po4_input = np.float32(_mi[_mi.index == "PO4"].iloc[0])
        self._mineral_input_on = bool(self._nh4_input != 0 or self._po4_input != 0)
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
        self.Microbe_E_Gain = np.float32('nan')  # taxon net energy gain
        self.Microbe_E_Cost = np.float32('nan')  # taxon net energy costs
        self.Uptake_Vmax_Km_slope = np.float32(data_init['Uptake_Vmax_Km'])  # slope of the uptake Vmax-Km relationship
        self.Uptake_Vmax_Km_int = np.float32(data_init['Uptake_Vmax_Km_int'])  # intercept of the uptake Vmax-Km relationship
        self.Uptake_Km_min = np.float32(data_init['Uptake_Km_min'])
        self.Uptake_Km_factor = data_init['Uptake_Km_factor']  # (G*M, n_uptake) error factor, or None (Km_error = 0)
        self.Uptake_sub_spec = data_init[
            'Uptake_Vmax_sub_spec']  # substrate specificity-efficiency adjustment for monomer uptake Vmax
        self.Uptake_Km_Temp_Sens = data_init[
            'Uptake_km_temp_sens']  # 1: uptake Km follows the daily temperature; 0: evaluated at Km_Tref
        self.km_Tref = data_init['km_tref']  # temperature (deg C) at which Km is evaluated when not temperature sensitive
        self.Pox_Ferm = data_init['Fermoxp']  # proportion of fermented monomer converted to CO2
        # fermentation condition score (0-1): reflects environmental favorability of fermentation vs respiration
        self.ferm_cond = data_init['FermCond']  # currently set to 0.1; can update to vary with environmental conditions
        self.active_uptake_E_cost = data_init['ActiveUptakeE_cost']  # energy cost of transporting monomers
        # Energy deficit: when on, a cell whose energy costs exceed its energy pays the shortfall by respiring free
        # biomass C (C not bound in HSP/osmolytes) at the respiration energy yield of the necromass monomer (E_resp of
        # Mon_DeadMic; no fermentation, all eroded C leaves as CO2); when off, the shortfall is floored away (E set to 0).
        self.E_deficit_erosion = bool(data_init['E_deficit_erosion'])
        if self.E_deficit_erosion:
            if 'Mon_DeadMic' not in self.Monomer_energy.index:
                raise ValueError("E_deficit_erosion = 1 needs a 'Mon_DeadMic' row in monomer_energy.csv")
            self._erosion_E_yield = np.float32(self.Monomer_energy.loc['Mon_DeadMic', 'E_resp'])  # energy per C eroded

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
        self.Taxon_Resp = np.float32('nan')  # Taxon specific Respiration
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
        self.HSP_Turnover = data_init['HSPturnover']  # intracellular recycling rate of HSP (released C/N stay in the cell)
        self.HSP_max_frac = data_init['HSP_max_frac']  # HSP C / biomass at which HSP production stops
        self.Osmo_max_frac = data_init['Osmo_max_frac']  # osmolyte C / biomass at which osmolyte production stops

        # Intracellular compound pools (C units; N = C x HSP_N_cost / Osmo_N_cost), part of biomass C and N.
        # Rows follow Microbes (gridsize*n_taxa); columns are HSP / osmolyte genes. Start empty; carried through
        # mortality (zeroed for dead cells), reproduction (split, translocated and dispersed with biomass) and
        # repopulation (reset). Exposed as DataFrames through the Taxon_HSP / Taxon_Osmolyte properties.
        self._hsp_cols = [c for c in data_init['Initial_HSP_Prop'].columns if c != 'Total_biomass']
        self._osmo_cols = list(data_init['OsmoProdConsti'].columns)
        n_rows = len(data_init['Microbes'])
        self._hsp_pool = np.zeros((n_rows, len(self._hsp_cols)), dtype=np.float32)
        self._osmo_pool = np.zeros((n_rows, len(self._osmo_cols)), dtype=np.float32)

        # Mortality
        self.MinRatios = data_init['MinRatios']  # Minimal cell quotas
        self.C_min = data_init['C_min']  # C threshold value of living cell
        self.N_min = data_init['N_min']  # N threshold value of living cell
        self.P_min = data_init['P_min']  # P threshold value of living cell
        # Energy:biomass C bounds of each cell (bacteria: E_*_ratio_b, fungi: E_*_ratio_f). E_max: energy above
        # E_max * C is dissipated (mortality step), and the fraction of uptake C routed to energy falls to 0 as E:C
        # approaches it. E_min: at or below it all energy-eligible uptake C is routed to energy.
        _fungal = np.asarray(data_init['fb']) == 1
        self.E_min = np.where(_fungal, data_init['E_min_ratio_f'], data_init['E_min_ratio_b']).astype(np.float32)
        self.E_max = np.where(_fungal, data_init['E_max_ratio_f'], data_init['E_max_ratio_b']).astype(np.float32)
        self.basal_death_prob = data_init['basal_death_prob']  # Basal death probability of microbes
        self.death_rate = data_init['death_rate']  # Change rate of mortality with water potential
        self.Osmo_tol_K = data_init['Osmo_tol_K']  # osmolyte share of biomass giving drought tolerance 0.5
        self.Drought_tol = np.zeros(len(data_init['Microbes']), dtype=np.float32)  # per cell, set in mortality()
        self.wp_fc = data_init['wp_fc']  # -1.0
        self.Leaching = data_init['Leaching']  # daily NH4/PO4 leaching fraction at Psi = wp_fc
        self.Leaching_psi_slope = data_init['Leaching_psi_slope']  # moisture sensitivity of leaching
        self.Leached_N = np.zeros(self.gridsize, dtype=np.float32)  # NH4-N leached per grid cell, set in mortality()
        self.Leached_P = np.zeros(self.gridsize, dtype=np.float32)  # PO4-P leached per grid cell, set in mortality()
        self.uptake_psi_slope = data_init['Uptake_psi_slope']  # moisture sensitivity of uptake
        self.enzyme_psi_slope = data_init['Enzyme_psi_slope']  # moisture sensitivity of enzyme Vmax
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
        self.direct = float(runtime.loc['direct', 1])  # probability that fungal daughters disperse straight (no y shift)

        # Climate data
        self.temp = data_init['Temp']  # series; temperature
        self.psi = data_init['Psi']  # series; water potential

        # Global constants
        # NOTE: self.Km_Ea (formerly a hardcoded 20 kJ/mol shared by both the
        # degradation and uptake Arrhenius-Km calculations) was removed 2026-07-23.
        # It was decoupled from each substrate/enzyme's own Ea, so Km barely
        # shrank relative to Vmax0's raw 1e7-1e8 scale -- silently suppressing
        # both degradation and uptake under Arrhenius kinetics. Both call sites
        # now reuse self.Ea / self.Uptake_Ea instead of a single global constant.
        self.MMRT_Tref = data_init['mmrt_tref']  # reference temperature for calculating rates (20 celcius)
        self.gas_const = np.float32(const.R)  # Universal gas constant
        self.boltz_const = np.float32(const.k)  # Boltzmann's constant
        self.planck_const = np.float32(const.h)  # Planck's constant
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

        # Per-taxon production rates and per-enzyme costs as float32 arrays (fixed for the whole run)
        self._hsp_consti_rate = self.HSP_consti_pre_exp.to_numpy(dtype=np.float32, copy=True)
        self._hsp_induci_rate = self.HSP_induci_pre_exp.to_numpy(dtype=np.float32, copy=True)
        self._osmo_consti_rate = self.Consti_Osmo_C.to_numpy(dtype=np.float32, copy=True)
        self._osmo_induci_rate = self.Induci_Osmo_C.to_numpy(dtype=np.float32, copy=True)
        self._enz_consti_rate = self.Consti_Enzyme_C.to_numpy(dtype=np.float32, copy=True)
        self._enz_induci_rate = self.Induci_Enzyme_C.to_numpy(dtype=np.float32, copy=True)

        def _enz_costs(cols, names):
            return tuple(self.Enz_Attrib[n].reindex(cols).to_numpy(dtype=np.float32, copy=True) for n in names)
        self._enz_cost_consti = _enz_costs(self.Consti_Enzyme_C.columns, ('N_cost', 'P_cost', 'E_cost'))
        self._enz_cost_induci = _enz_costs(self.Induci_Enzyme_C.columns, ('N_cost', 'P_cost', 'E_cost'))
        self._enz_cost_pool = _enz_costs(self._enz_index, ('N_cost', 'P_cost'))

    def degradation(self, pulse, day):
        """
        Explicit degradation of different substrates.

        Calculation procedure:
        1. Determine substates pool: incl. inputs
        2. Compute Vmax & Km and make them follow the index of Substrates
        3. Follow the Michaelis-Menten equation to compute full degradation rate
        4. Impose the substrate-required enzymes upon the full degradation rate
        5. Adjust cellulose rate with LCI(lignocellulose index)

        NOTE: rewritten to numpy-native arithmetic, same strategy used for
        uptake()/metabolism()/mortality()/reproduction() this session. Two
        pandas idioms here were both known-slow AND, in one case, silently
        load-bearing -- both are called out explicitly below because getting
        either wrong changes model output:

        1. `tev_transition.stack().unstack(1).reset_index(level=0, drop=True)`
           reshapes Vmax*Enzymes from (enz*grid, sub) into (sub*grid, enz)
           layout. This is mathematically just `.reshape(G, n_enz, n_sub)
           .transpose(0, 2, 1).reshape(G*n_sub, n_enz)` -- confirmed bit-exact
           against the original via a standalone probe. `.unstack()` also
           resorts the resulting enzyme columns ALPHABETICALLY (e.g. Enz1,
           Enz10, Enz11, Enz12, Enz2, ...), which is why the original
           immediately re-indexes with `tev = tev[Km.columns]` afterward. The
           numpy reshape never touches labels at all, so it lands directly in
           Vmax's own (natural, numeric) enzyme order -- confirmed via a
           second probe (with n_enz=12, deliberately chosen so alphabetical
           and numeric order differ) that this matches Km's column order
           without needing that reindex, since both are built from the same
           'Enz1'..'EnzN' construction order.

        2. `DecayRates = pd.concat([batch1, rss], axis=1).min(axis=1,
           skipna=True)` -- pandas' row-wise min with skipna=True treats NaN
           as "not present" and returns the OTHER (non-NaN) value instead of
           propagating it. `np.minimum(batch1, rss)` does NOT replicate this
           -- it propagates NaN, which (per this session's Consti_Osmo_C bug)
           can silently poison everything downstream forever once introduced.
           `np.fmin(batch1, rss)` is the correct numpy equivalent: it ignores
           NaN in either operand exactly like pandas' skipna=True, and only
           returns NaN if BOTH operands are NaN (matching what skipna=True
           would do in that case too, since there'd be nothing left to pick).
        """

        G   = self.gridsize
        M   = self.n_substrates      # note: local var name mirrors upstream n_substrates convention
        NE  = self.n_enzymes

        def _v(x):
            return x.to_numpy(dtype=np.float32, copy=True) if hasattr(x, "to_numpy") else x

        # constant of lignocellulose index--LCI
        LCI_slope = np.float32(-0.8)
        # Substrates index by substrate names (grid-major, substrate-minor; fixed for the whole run)
        Sub_index_arr = np.asarray(self.Substrates.index)

        Substrates_v = self.Substrates.to_numpy(dtype=np.float32, copy=True)  # (G*M, 3)

        # total mass of each substrate: C+N+P
        rss = Substrates_v.sum(axis=1)  # (G*M,)
        # substrate stoichiometry; NOTE: ensure NaN (b/c of 0/0) = 0, matching the
        # original's `.divide(rss, axis=0).fillna(0)` exactly (not "improving" it)
        SubstrateRatios_v = np.zeros_like(Substrates_v)
        with np.errstate(divide='ignore', invalid='ignore'):
            # np.divide's `where=` guard still evaluates the masked-out branch
            # internally (a known numpy quirk), which would otherwise emit a
            # spurious "invalid value encountered" warning even though those
            # positions are never written to `out` -- harmless, just noisy.
            np.divide(Substrates_v, rss[:, None], out=SubstrateRatios_v, where=(rss[:, None] != 0))

        # moisture effects on enzymatic kinetics
        if self.psi[day] >= self.wp_fc:
            f_psi = np.float32(1.0)
        else:
            f_psi = np.exp(self.enzyme_psi_slope * (self.psi[day] - self.wp_fc)).astype('float32')

        # calculate enzyme Vmax and Km using the user specified kinetics model (MMRT = 1; Arrhenius = 0).
        # Km follows the Vmax-Km relationship in both: Km = Vmax(T_Km) x Vmax_Km x error factor + Vmax_Km_int, where
        # T_Km is the day's temperature (Enz_Km_temp_sens = 1) or Km_Tref (0); Vmax here excludes the moisture factor.
        def _sub_major(x):  # (G*NE, M) -> (G*M, NE)
            return x.reshape(G, NE, M).transpose(0, 2, 1).reshape(G * M, NE)

        if self.kinetics_model > 0:
            Vmax_mmrt = MMRT(_v(self.Cp), _v(self.HT0), _v(self.ST0), self.temp[day], self.MMRT_Tref,
                             self.boltz_const, self.planck_const, self.gas_const)
            Vmax_mmrt_sub_spec = Vmax_mmrt * _v(self.Enz_Vmax_sub_spec)
            Vmax_v = Vmax_mmrt_sub_spec * f_psi  # (G*NE, M)
            if self.Enz_Km_Temp_Sens == 1:
                Km_vmax = _sub_major(np.asarray(Vmax_mmrt_sub_spec, dtype=np.float32))
            else:
                if not hasattr(self, "_km_vmax_tref"):  # constant for the run
                    self._km_vmax_tref = _sub_major(np.asarray(
                        MMRT(_v(self.Cp), _v(self.HT0), _v(self.ST0), self.km_Tref, self.MMRT_Tref,
                             self.boltz_const, self.planck_const, self.gas_const) * _v(self.Enz_Vmax_sub_spec),
                        dtype=np.float32))
                Km_vmax = self._km_vmax_tref
            Km_slope = self.Vmax_Km_slope if self.Km_factor is None else self.Vmax_Km_slope * _v(self.Km_factor)
            Km_v = np.maximum(Km_vmax * Km_slope + self.Vmax_Km_int, self.Km_min)  # (G*M, NE)
        else:
            # NOTE (fix 2026-07-23): Km's Arrhenius shrink factor previously used a
            # hardcoded, model-wide self.Km_Ea = 20 -- completely decoupled from each
            # substrate/enzyme's own Ea (35000-39000). Since Km0 = Vmax0 * Vmax_Km is
            # built from the *raw* Vmax0 (1e7-1e8), it needs the SAME order-of-magnitude
            # shrink that Vmax gets from Arrhenius(Ea~35000+) to land in a sane range;
            # Km_Ea=20 barely shrinks it at all, leaving Km ~1e8 -- orders of magnitude
            # above any realistic substrate mass -- which silently suppressed decay to
            # ~0 under Arrhenius kinetics (confirmed numerically: ~1.8e-7 %/day).
            # Fixed by reusing each substrate/enzyme's own Ea instead of a single global
            # constant. self.Ea is (G*NE, M)-shaped (enzyme-major) while self.Km0 is
            # (G*M, NE)-shaped (substrate-major); cache the reshape/transpose (same
            # operation already used below for tev_transition_v -> tev_v) since Ea
            # itself never changes across days.
            if not hasattr(self, "_ea_v"):
                self._ea_v = _v(self.Ea)  # (G*NE, M) -- fixed for the whole run
                self._ea_km_v = self._ea_v.reshape(G, NE, M).transpose(0, 2, 1).reshape(G * M, NE)  # (G*M, NE)
            Vmax_v = _v(self.Vmax0) * Arrhenius(self._ea_v, self.gas_const, self.temp[day],
                                                self.Arrh_Tref) * f_psi          # (G*NE, M)
            T_km = self.temp[day] if self.Enz_Km_Temp_Sens == 1 else self.km_Tref
            Km_v = np.maximum(_v(self.Km0) * Arrhenius(self._ea_km_v, self.gas_const, T_km, self.Arrh_Tref),
                              self.Km_min)                                        # (G*M, NE)

        Vmax_v = np.asarray(Vmax_v, dtype=np.float32)
        Km_v = np.asarray(Km_v, dtype=np.float32)

        # Multiply Vmax by enzyme concentration, then reshape (enz*grid, sub) -> (sub*grid, enz).
        # See docstring point 1: plain reshape/transpose, bit-exact vs. the original's
        # stack/unstack + explicit column reindex, and far cheaper.
        Enzymes_v = self.Enzymes.to_numpy(dtype=np.float32, copy=True)  # (G*NE,)
        tev_transition_v = Vmax_v * Enzymes_v[:, None]                  # (G*NE, M)
        tev_v = tev_transition_v.reshape(G, NE, M).transpose(0, 2, 1).reshape(G * M, NE)  # (G*M, NE)

        # Michaelis-Menten equation
        Decay_v = (tev_v * rss[:, None]) / (Km_v + rss[:, None])  # (G*M, NE)

        if getattr(self, "debug_degradation_days", 0) and day < self.debug_degradation_days:
            print(f"  -- degradation() [numpy] shape probe (day {day}) --")
            print(f"    Vmax_v.shape={Vmax_v.shape}  Km_v.shape={Km_v.shape}  "
                  f"tev_v.shape={tev_v.shape}  Decay_v.shape={Decay_v.shape}")

        # Pull out each batch of required enzymes and sum across redundant enzymes.
        # self.ReqEnz never changes across days (fixed biochemical requirement table),
        # so the 'set1' slice + array extraction is cached on first use instead of
        # repeating a MultiIndex `.loc[]` selection every single day.
        if not hasattr(self, "_reqenz_set1_v"):
            self._reqenz_set1_v = self.ReqEnz.loc['set1'].to_numpy(dtype=np.float32, copy=True)  # (G*M, NE)
        batch1 = (self._reqenz_set1_v * Decay_v).sum(axis=1)  # (G*M,)

        # Compare to substrate available and take the min -- np.fmin (not np.minimum!)
        # to replicate pandas' `.min(axis=1, skipna=True)` NaN-as-missing semantics;
        # see docstring point 2.
        DecayRates = np.fmin(batch1, rss)  # (G*M,)

        # Adjust cellulose rate by linking cellulose degradation to lignin concentration (LCI)
        is_lignin = Sub_index_arr == 'Lignin'
        is_cellulose = Sub_index_arr == 'Cellulose'
        if is_lignin.any() and is_cellulose.any():
            ss7 = Substrates_v[is_lignin].sum(axis=1)  # (G,) -- one Lignin row per grid cell
            cellulose_C = Substrates_v[is_cellulose, 0]  # (G,) -- C column of the Cellulose row per cell
            with np.errstate(divide='ignore', invalid='ignore'):
                # Matches the original exactly: no zero-denominator guard here in
                # either version, so a genuine 0/0 (ss7 + cellulose_C == 0) still
                # produces NaN in both -- only the warning is suppressed, not the
                # underlying (unchanged, original) behavior.
                DecayRates[is_cellulose] *= np.float32(1) + (ss7 / (ss7 + cellulose_C)) * LCI_slope

        
        # Update Substrates Pool by removing decayed C, N, & P. Depending on specific needs, adding inputs of substrates can be done here
        Substrates_v -= SubstrateRatios_v * DecayRates[:, None]
        self.Substrates[:] = Substrates_v

        # Pass these two back to the global variables to be used in the next method
        self.SubstrateRatios = pd.DataFrame(SubstrateRatios_v, index=self.Substrates.index,
                                            columns=self.Substrates.columns)
        self.DecayRates = pd.Series(DecayRates, index=self.Substrates.index)

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
        # Ongoing NH4 and PO4 inputs (N and P per grid cell per day), added before uptake; none if both are 0
        if self._mineral_input_on:
            self.Monomers.loc[self.Monomers.index == "NH4", "N"] += self._nh4_input
            self.Monomers.loc[self.Monomers.index == "PO4", "P"] += self._po4_input

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
            f_psi = np.exp(self.uptake_psi_slope * (self.psi[day] - self.wp_fc)).astype('float32')

        if self.kinetics_model > 0:
            Uptake_Vmax_mmrt = MMRT(self.Uptake_Cp, self.Uptake_HT0, self.Uptake_ST0, self.temp[day], self.MMRT_Tref,
                                    self.boltz_const, self.planck_const, self.gas_const)
            Uptake_Vmax_sub_spec = Uptake_Vmax_mmrt * self.Uptake_sub_spec  # account for substrate specificity-efficiency tradeoff
            Uptake_Vmax = Uptake_Vmax_sub_spec * f_psi  # account for moisture limitation

            # Km from the Vmax-Km relationship, at the day's temperature or at Km_Tref (see degradation())
            if self.Uptake_Km_Temp_Sens == 1:
                Uptake_Km_vmax = Uptake_Vmax_sub_spec
            else:
                if not hasattr(self, "_uptake_km_vmax_tref"):  # constant for the run
                    self._uptake_km_vmax_tref = MMRT(self.Uptake_Cp, self.Uptake_HT0, self.Uptake_ST0, self.km_Tref,
                                                     self.MMRT_Tref, self.boltz_const, self.planck_const,
                                                     self.gas_const) * self.Uptake_sub_spec
                Uptake_Km_vmax = self._uptake_km_vmax_tref
            Uptake_Km_slope = (self.Uptake_Vmax_Km_slope if self.Uptake_Km_factor is None
                               else self.Uptake_Km_factor.to_numpy(dtype=np.float32) * self.Uptake_Vmax_Km_slope)
            Uptake_Km = ((Uptake_Km_vmax * Uptake_Km_slope) + self.Uptake_Vmax_Km_int).clip(lower=self.Uptake_Km_min)

        else:
            # NOTE (fix 2026-07-23): same fix as degradation()'s Arrhenius branch --
            # Uptake_Km used to be shrunk by the hardcoded self.Km_Ea = 20 instead of
            # Uptake_Ea, leaving Uptake_Km ~1e8 (barely shrunk from Uptake_Vmax0's raw
            # 1e7-1e8 scale) vs. realistic monomer masses -- silently suppressing
            # uptake to a negligible-but-technically-nonzero level. self.Uptake_Km0
            # and self.Uptake_Ea are both (monomer x transporter)-shaped already, so
            # (unlike degradation's Ea/Km0, which are transposed relative to each
            # other) no reshape is needed here.
            Uptake_Vmax = self.Uptake_Vmax0 * Arrhenius(self.Uptake_Ea, self.gas_const, self.temp[day],
                                                        self.Arrh_Tref) * f_psi
            T_km = self.temp[day] if self.Uptake_Km_Temp_Sens == 1 else self.km_Tref
            Uptake_Km = (self.Uptake_Km0 * Arrhenius(self.Uptake_Ea, self.gas_const, T_km, self.Arrh_Tref)
                         ).clip(lower=self.Uptake_Km_min)

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
        Mic_C_v = self.Microbes['C'].to_numpy(dtype=np.float32)
        Mic_E_v = self.Microbes['E'].to_numpy(dtype=np.float32)
        with np.errstate(divide='ignore', invalid='ignore'):
            # distance of E:C below E_max, as a fraction of the E_min..E_max range (1 at or below E_min, 0 at E_max)
            E_room = np.clip((self.E_max - Mic_E_v / Mic_C_v) / (self.E_max - self.E_min), 0, 1)
            Tax_E_demand = 1 - ((1 - E_room) ** (self.ferm_cond + 1))
        self.Tax_E_demand = Tax_E_demand  # fraction of uptake C routed to energy (per cell; NaN in empty cells)

        Org_Mon_E_value = self.Monomer_energy['E_value'] ** (
            (2 / np.log10(self.Monomer_energy['E_ferm'])) ** self.ferm_cond)
        Mon_E_value = np.concatenate([[0, 0], Org_Mon_E_value.values]).astype(np.float32)  # (M,)
        # fraction of each monomer's C routed to energy, per cell: (M, G*T) -> (G, M, T)
        Tea_3d = np.outer(Mon_E_value,
                          Tax_E_demand.astype(np.float32)
                          ).reshape(M, G, T).transpose(1, 0, 2).astype(np.float32)  # (G,M,T)
        Tea_3d = np.nan_to_num(Tea_3d, nan=0.0)

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

        # Taxon-cell uptake (G*T,) in Microbes' row order (grid outer, taxon inner)
        def _cell(x):
            return pd.Series(x.ravel(order='C').astype(np.float32), index=self.Microbes.index)
        self.Taxon_Uptake_C = _cell(TUC)      # total C taken up
        self.Taxon_Uptake_BC = _cell(TUCb)    # C routed to biomass
        self.Taxon_Uptake_N = _cell(TUN)
        self.Taxon_Uptake_P = _cell(TUP)
        self.Taxon_Uptake_E = _cell(TUE)      # energy produced from the C routed to energy
        self.Taxon_Total_Uptake = _cell(TU)   # biomass C + N + P taken up

        # Taxon-specific C assimilation efficiency (CAE)
        #
        # NOTE: this used to be `self.CAE_taxon[pos_uptake_index] = A[pos_uptake_index] /
        # B[pos_uptake_index]` -- a boolean-Series-indexed assignment on Microbes.index,
        # which repeats each taxon label gridsize times (e.g. 10,000x on a 100x100 grid).
        # Under pandas 3.0, that pattern can fall back to an expensive non-unique-index
        # JOIN (pandas.core.reshape.merge.get_join_indexers_non_unique) instead of a cheap
        # positional mask, and on a heavily-duplicated index that join's memory use exploded
        # to 50+ GB in testing on a 100x100/100-taxon grid -- not a leak, but a pathological
        # alignment cost triggered by this specific indexing pattern. Doing the same math on
        # raw numpy arrays (positional boolean masking, no index/labels involved at all)
        # avoids the join entirely. output.py's _taxa_sum() already accepts plain arrays via
        # np.asarray(), so CAE_taxon doesn't need to stay a pandas Series.
        Taxon_Uptake_C_arr = self.Taxon_Uptake_C.to_numpy(dtype=np.float32)
        Taxon_Uptake_BC_arr = self.Taxon_Uptake_BC.to_numpy(dtype=np.float32)
        CAE_taxon_v = np.zeros(len(self.Microbes), dtype=np.float32)
        pos_uptake_index = Taxon_Uptake_C_arr > 0
        CAE_taxon_v[pos_uptake_index] = Taxon_Uptake_BC_arr[pos_uptake_index] / Taxon_Uptake_C_arr[pos_uptake_index]
        self.CAE_taxon = CAE_taxon_v

        # System-level CAE
        Taxon_Uptake_C_grid = self.Taxon_Uptake_C.sum()  # Total C Uptake
        Taxon_Uptake_BC_grid = self.Taxon_Uptake_BC.sum()
        if Taxon_Uptake_C_grid == 0:
            self.CAE_system = np.float32(0)
        else:
            self.CAE_system = Taxon_Uptake_BC_grid / Taxon_Uptake_C_grid


    def metabolism(self, day):
        """
        Intracellular production of enzymes, heat shock proteins (HSP) and osmolytes; CUE (finished in mortality()).

        Compound classes
          Enzymes are exported: their C, N and P leave the cell (to the Enzymes pool).
          HSP and osmolytes are intracellular and part of biomass: their C and N stay in Microbes['C'/'N'] but are
          bound in the compound pools (self._hsp_pool, self._osmo_pool; C units, N = C x *_N_cost) and so are not
          available to other syntheses until turnover returns them to the free pool.
        Turnover (HSP_turnover, Osmo_loss_rate) is intracellular recycling: each day a fraction of each pool is
          released back to the free C and N of the same cell (no loss from the cell).
        Composition feedback: production of each intracellular compound X is scaled by
          g_X = clip(1 - f_X / X_max_frac, 0, 1),   f_X = X pool C / biomass (C + N + P),
          so a cell cannot be made up of more than X_max_frac of X (HSP_max_frac, Osmo_max_frac).
        Resource limitation (joint, all three compound classes, same factor for each class):
          constitutive: free biomass C and N (total minus the C and N bound in HSP and osmolytes) and the energy left
                        after transport and maintenance;
          inducible:    the C routed to biomass and the N taken up this step, and the energy left after constitutive
                        production.
        Constitutive osmolyte production is not gated by water potential; inducible osmolyte production is
          (microbe_osmo_psi).
        Energy costs of production are paid from the energy pool; maintenance and transport costs are charged in full.
        """
        G, T = self.gridsize, self.n_taxa
        is_deadEnz = self.Substrates.index == "DeadEnz"

        # Moisture modifier of inducible osmolyte production
        f_psi = microbe_osmo_psi(self.psi[day], self.alpha, self.wp_fc, self.wp_th)

        Microbes_C = self.Microbes['C'].to_numpy(dtype=np.float32, copy=True)
        Microbes_N = self.Microbes['N'].to_numpy(dtype=np.float32, copy=True)
        Microbes_P = self.Microbes['P'].to_numpy(dtype=np.float32, copy=True)
        Microbes_E = self.Microbes['E'].to_numpy(dtype=np.float32, copy=True)
        Microbes_CNP = Microbes_C + Microbes_N + Microbes_P

        Taxon_Total_Uptake_v = np.asarray(self.Taxon_Total_Uptake, dtype=np.float32)
        Taxon_Uptake_BC_v = np.asarray(self.Taxon_Uptake_BC, dtype=np.float32)
        Taxon_Uptake_N_v = np.asarray(self.Taxon_Uptake_N, dtype=np.float32)
        Taxon_Uptake_P_v = np.asarray(self.Taxon_Uptake_P, dtype=np.float32)
        Taxon_Uptake_E_v = np.asarray(self.Taxon_Uptake_E, dtype=np.float32)
        Taxon_Uptake_C_v = np.asarray(self.Taxon_Uptake_C, dtype=np.float32)

        # ── Intracellular pools: turnover (recycling to the free pool), composition, free C and N ──────────
        self._hsp_pool *= np.float32(1) - np.float32(self.HSP_Turnover)
        self._osmo_pool *= np.float32(1) - np.float32(self.OSmo_Loss_Rate)
        hsp_C = self._hsp_pool.sum(axis=1)
        osmo_C = self._osmo_pool.sum(axis=1)
        alive_cnp = Microbes_CNP > 0
        HSP_prop = np.zeros_like(Microbes_CNP)
        Osmo_prop = np.zeros_like(Microbes_CNP)
        HSP_prop[alive_cnp] = hsp_C[alive_cnp] / Microbes_CNP[alive_cnp]
        Osmo_prop[alive_cnp] = osmo_C[alive_cnp] / Microbes_CNP[alive_cnp]
        g_hsp = np.clip(np.float32(1) - HSP_prop / np.float32(self.HSP_max_frac), 0, 1).astype(np.float32)
        g_osmo = np.clip(np.float32(1) - Osmo_prop / np.float32(self.Osmo_max_frac), 0, 1).astype(np.float32)
        free_C = np.maximum(Microbes_C - hsp_C - osmo_C, 0)
        free_N = np.maximum(Microbes_N - np.float32(self.HSP_N_cost) * hsp_C - np.float32(self.Osmo_N_cost) * osmo_C, 0)

        # ── Transport and maintenance energy (charged in full) ───────────────────────────────────────────────
        Tax_active_transport_E_cost = Taxon_Total_Uptake_v * np.float32(self.active_uptake_E_cost)
        maint_E_cost_Ea = ((HSP_prop / (HSP_prop + self.maint_E_HSP_Ea_n)) * self.maint_E_HSP_Ea_slope
                           ) + self.maint_E_HSP_Ea_int
        tax_maint_E_cost = self.maint_E_rate_const * Arrhenius(
            maint_E_cost_Ea, self.gas_const, self.temp[day], self.maint_E_rate_Tref)
        Tax_cell_maint_E_cost = tax_maint_E_cost * Microbes_CNP
        Uptake_Prod_Cost_sum = self.Uptake_Prod_Cost.to_numpy(dtype=np.float32, copy=True).sum(axis=1)
        Tax_transport_prod_E_cost = (Taxon_Total_Uptake_v * Uptake_Prod_Cost_sum) * np.float32(self.UptakeProd_E_Cost)

        # energy available for biosynthesis: the pool plus this step's energy, minus transport and maintenance
        MicBio_avail_E = (Microbes_E + Taxon_Uptake_E_v - Tax_active_transport_E_cost
                          - Tax_transport_prod_E_cost - Tax_cell_maint_E_cost)
        MicBio_avail_E = np.maximum(MicBio_avail_E, np.float32(0))

        Enz_N_cost_v, Enz_P_cost_v, Enz_E_cost_v = self._enz_cost_consti
        Enz_N_cost_i, Enz_P_cost_i, Enz_E_cost_i = self._enz_cost_induci
        HSP_N, HSP_E = np.float32(self.HSP_N_cost), np.float32(self.HSP_E_cost)
        Osmo_N, Osmo_E = np.float32(self.Osmo_N_cost), np.float32(self.Osmo_E_cost)
        Arr_hsp_consti = Arrhenius(self.HSP_consti_Ea, self.gas_const, self.temp[day], self.maint_E_rate_Tref)
        Arr_hsp_induci = Arrhenius(self.HSP_induci_Ea, self.gas_const, self.temp[day], self.maint_E_rate_Tref)

        # ── Constitutive production: potential, then joint C / N / E limitation on free resources ──────────────
        tax_HSP_consti_prod = (self._hsp_consti_rate * Arr_hsp_consti
                               * (Microbes_CNP * g_hsp)[:, None])                       # (G*T, n_HSP), C units
        Taxon_Osmo_Consti = self._osmo_consti_rate * (Microbes_C * g_osmo)[:, None]     # (G*T, n_osmo), no Psi gate
        Taxon_Enzyme_Consti = self._enz_consti_rate * Microbes_C[:, None]               # (G*T, n_enz)

        cC = tax_HSP_consti_prod.sum(axis=1) + Taxon_Osmo_Consti.sum(axis=1) + Taxon_Enzyme_Consti.sum(axis=1)
        cN = (tax_HSP_consti_prod.sum(axis=1) * HSP_N + Taxon_Osmo_Consti.sum(axis=1) * Osmo_N
              + (Taxon_Enzyme_Consti * Enz_N_cost_v).sum(axis=1))
        cE = (tax_HSP_consti_prod.sum(axis=1) * HSP_E + Taxon_Osmo_Consti.sum(axis=1) * Osmo_E
              + (Taxon_Enzyme_Consti * Enz_E_cost_v).sum(axis=1))
        r_consti = _limit_ratio((cC, free_C), (cN, free_N), (cE, MicBio_avail_E))

        tax_HSP_consti_prod *= r_consti[:, None]
        Taxon_Osmo_Consti *= r_consti[:, None]
        Taxon_Enzyme_Consti *= r_consti[:, None]

        Taxon_HSP_Consti_Cost_C = tax_HSP_consti_prod.sum(axis=1)
        Taxon_HSP_Consti_Cost_E = Taxon_HSP_Consti_Cost_C * HSP_E
        Taxon_Osmo_Consti_Cost_C = Taxon_Osmo_Consti.sum(axis=1)
        Taxon_Osmo_Consti_Cost_E = Taxon_Osmo_Consti_Cost_C * Osmo_E
        Taxon_Enzyme_Consti_Cost_C = Taxon_Enzyme_Consti.sum(axis=1)
        Taxon_Enzyme_Consti_Cost_N = (Taxon_Enzyme_Consti * Enz_N_cost_v).sum(axis=1)
        Taxon_Enzyme_Consti_Cost_P = (Taxon_Enzyme_Consti * Enz_P_cost_v).sum(axis=1)
        Taxon_Enzyme_Consti_Cost_E = (Taxon_Enzyme_Consti * Enz_E_cost_v).sum(axis=1)

        MicBio_avail_E = np.maximum(
            MicBio_avail_E - (Taxon_HSP_Consti_Cost_E + Taxon_Osmo_Consti_Cost_E + Taxon_Enzyme_Consti_Cost_E),
            np.float32(0))

        # ── Inducible production: potential, then joint C / N / E limitation on this step's uptake ──────────────
        tax_HSP_induci_prod = (self._hsp_induci_rate * Arr_hsp_induci
                               * (Taxon_Total_Uptake_v * g_hsp)[:, None])
        Taxon_Osmo_Induci = self._osmo_induci_rate * (Taxon_Uptake_BC_v * g_osmo * f_psi)[:, None]
        Taxon_Enzyme_Induci = self._enz_induci_rate * Taxon_Uptake_BC_v[:, None]

        iC = tax_HSP_induci_prod.sum(axis=1) + Taxon_Osmo_Induci.sum(axis=1) + Taxon_Enzyme_Induci.sum(axis=1)
        iN = (tax_HSP_induci_prod.sum(axis=1) * HSP_N + Taxon_Osmo_Induci.sum(axis=1) * Osmo_N
              + (Taxon_Enzyme_Induci * Enz_N_cost_i).sum(axis=1))
        iE = (tax_HSP_induci_prod.sum(axis=1) * HSP_E + Taxon_Osmo_Induci.sum(axis=1) * Osmo_E
              + (Taxon_Enzyme_Induci * Enz_E_cost_i).sum(axis=1))
        r_induci = _limit_ratio((iC, Taxon_Uptake_BC_v), (iN, Taxon_Uptake_N_v), (iE, MicBio_avail_E))

        tax_HSP_induci_prod *= r_induci[:, None]
        Taxon_Osmo_Induci *= r_induci[:, None]
        Taxon_Enzyme_Induci *= r_induci[:, None]

        Taxon_HSP_Induci_Cost_C = tax_HSP_induci_prod.sum(axis=1)
        Taxon_HSP_Induci_Cost_E = Taxon_HSP_Induci_Cost_C * HSP_E
        Taxon_Osmo_Induci_Cost_C = Taxon_Osmo_Induci.sum(axis=1)
        Taxon_Osmo_Induci_Cost_E = Taxon_Osmo_Induci_Cost_C * Osmo_E
        Taxon_Enzyme_Induci_Cost_C = Taxon_Enzyme_Induci.sum(axis=1)
        Taxon_Enzyme_Induci_Cost_N = (Taxon_Enzyme_Induci * Enz_N_cost_i).sum(axis=1)
        Taxon_Enzyme_Induci_Cost_P = (Taxon_Enzyme_Induci * Enz_P_cost_i).sum(axis=1)
        Taxon_Enzyme_Induci_Cost_E = (Taxon_Enzyme_Induci * Enz_E_cost_i).sum(axis=1)

        # ── Integration ─────────────────────────────────────────────────────────────────────────────────────────
        # Only exported enzymes leave the cell; HSP and osmolyte C and N stay in biomass (moved to the bound pools).
        Microbe_C_Gain = Taxon_Uptake_BC_v - Taxon_Enzyme_Induci_Cost_C
        Microbe_N_Gain = Taxon_Uptake_N_v - Taxon_Enzyme_Induci_Cost_N
        Microbe_P_Gain = Taxon_Uptake_P_v - Taxon_Enzyme_Induci_Cost_P

        Microbe_E_Cost = (Tax_active_transport_E_cost + Tax_transport_prod_E_cost + Tax_cell_maint_E_cost
                          + Taxon_HSP_Consti_Cost_E + Taxon_HSP_Induci_Cost_E
                          + Taxon_Osmo_Consti_Cost_E + Taxon_Osmo_Induci_Cost_E
                          + Taxon_Enzyme_Consti_Cost_E + Taxon_Enzyme_Induci_Cost_E)
        Microbe_E_Gain = Taxon_Uptake_E_v - Microbe_E_Cost

        Microbes_C = np.maximum(Microbes_C + Microbe_C_Gain - Taxon_Enzyme_Consti_Cost_C, np.float32(0))
        Microbes_N = np.maximum(Microbes_N + Microbe_N_Gain - Taxon_Enzyme_Consti_Cost_N, np.float32(0))
        Microbes_P = np.maximum(Microbes_P + Microbe_P_Gain - Taxon_Enzyme_Consti_Cost_P, np.float32(0))
        E_raw = Microbes_E + Microbe_E_Gain
        Microbes_E = np.maximum(E_raw, np.float32(0))
        self.Microbes['C'] = Microbes_C
        self.Microbes['N'] = Microbes_N
        self.Microbes['P'] = Microbes_P
        self.Microbes['E'] = Microbes_E

        self._hsp_pool += tax_HSP_consti_prod + tax_HSP_induci_prod
        self._osmo_pool += Taxon_Osmo_Consti + Taxon_Osmo_Induci
        self._clip_pools(Microbes_C, Microbes_N)

        # ── Energy deficit: energy costs (transport, maintenance) above the energy available ──────────────────────
        # Production is limited to the energy available, so a shortfall comes from transport and maintenance only.
        E_deficit = np.maximum(-E_raw, np.float32(0))
        C_eroded = np.zeros_like(E_deficit)
        if self.E_deficit_erosion and E_deficit.any():
            free_C_now = np.maximum(Microbes_C - self._hsp_pool.sum(axis=1) - self._osmo_pool.sum(axis=1), 0)
            C_eroded = np.minimum(E_deficit / self._erosion_E_yield, free_C_now).astype(np.float32)
            Microbes_C = Microbes_C - C_eroded           # E stays at 0: the shortfall is paid, not stored
            self.Microbes['C'] = Microbes_C
            self.Taxon_resp_df += C_eroded.reshape(G, T)  # all eroded C is respired (with this day's uptake CO2)
        # energy that was not paid for (with erosion off, all of the shortfall; on, what free C could not cover)
        E_unpaid = np.maximum(E_deficit - C_eroded * (self._erosion_E_yield if self.E_deficit_erosion else 0), 0)

        # ── CUE, part 1: C retained so far (HSP and osmolyte C stays in biomass). Finished in mortality(), after the
        # C respired by stoichiometric rebalancing (see _update_cue) ─────────────────────────────────────────────────
        Mic_BC_Use = Microbe_C_Gain - Taxon_Enzyme_Consti_Cost_C
        if self.E_deficit_erosion:
            Mic_BC_Use = Mic_BC_Use - C_eroded           # biomass respired to pay energy deficits
        self._cue_C_retained = np.asarray(Mic_BC_Use, dtype=np.float32)
        self._cue_C_uptake = np.asarray(Taxon_Uptake_C_v, dtype=np.float32)
        self._update_cue()                                # provisional until mortality() subtracts rebalancing losses

        # ── Publish taxon-level results (plain float32 arrays) ───────────────────────────────────────────────────
        self.Tax_active_transport_E_cost = Tax_active_transport_E_cost
        self.Tax_cell_maint_E_cost = Tax_cell_maint_E_cost
        self.Taxon_HSP_Consti_Cost_E = Taxon_HSP_Consti_Cost_E
        self.Taxon_HSP_Induci_Cost_E = Taxon_HSP_Induci_Cost_E
        self.Microbe_E_Cost = Microbe_E_Cost
        self.Microbe_E_Gain = Microbe_E_Gain
        self.Taxon_Osmo_Cost_C = Taxon_Osmo_Consti_Cost_C + Taxon_Osmo_Induci_Cost_C     # C moved into the pool
        self.Taxon_Enzyme_Cost_C = Taxon_Enzyme_Consti_Cost_C + Taxon_Enzyme_Induci_Cost_C
        self.Taxon_HSP_Cost_C = Taxon_HSP_Consti_Cost_C + Taxon_HSP_Induci_Cost_C          # C moved into the pool
        self.Taxon_Transport_Cost_C = Tax_transport_prod_E_cost / np.float32(self.UptakeProd_E_Cost)
        self.Taxon_Transport_E_Cost = Tax_active_transport_E_cost + Tax_transport_prod_E_cost
        self.Taxon_Maint_E_Cost = Tax_cell_maint_E_cost
        self.Taxon_Enzyme_E_Cost = Taxon_Enzyme_Consti_Cost_E + Taxon_Enzyme_Induci_Cost_E
        self.Taxon_Osmo_E_Cost_total = Taxon_Osmo_Consti_Cost_E + Taxon_Osmo_Induci_Cost_E
        self.Taxon_HSP_E_Cost = Taxon_HSP_Consti_Cost_E + Taxon_HSP_Induci_Cost_E
        self.HSP_prop = HSP_prop                  # HSP C / biomass (C+N+P), start of step after turnover
        self.Osmo_prop = Osmo_prop                # osmolyte C / biomass (C+N+P), start of step after turnover
        self.HSP_regulation = g_hsp               # composition feedback on HSP production (1 = none)
        self.Osmo_regulation = g_osmo             # composition feedback on osmolyte production (1 = none)
        self.Consti_limitation = r_consti         # joint C/N/E limitation of constitutive production (1 = none)
        self.Induci_limitation = r_induci         # joint C/N/E limitation of inducible production (1 = none)
        self.maint_E_cost_Ea = maint_E_cost_Ea
        self.Taxon_E_deficit = E_deficit          # energy shortfall before any payment from biomass
        self.Taxon_Erosion_C = C_eroded           # biomass C respired to pay the shortfall (0 when the switch is off)
        self.Taxon_E_unpaid = E_unpaid            # shortfall left unpaid (floored away: energy not accounted for)
        self.tax_maint_E_cost = tax_maint_E_cost

        # ── Enzyme pools: production (summed over taxa per cell), turnover, dead enzymes to DeadEnz ──────────────
        Taxon_Enzyme_Production_v = Taxon_Enzyme_Consti + Taxon_Enzyme_Induci                # (G*T, n_enz)
        Enzyme_Production_v = Taxon_Enzyme_Production_v.reshape(G, T, self.n_enzymes).sum(axis=1).reshape(-1)
        Enzymes_v = self.Enzymes.to_numpy(dtype=np.float32, copy=True)
        Enzyme_Loss_v = Enzymes_v * np.float32(self.Enzyme_Loss_Rate)
        self.Enzymes[:] = Enzymes_v + Enzyme_Production_v - Enzyme_Loss_v
        self.Enzymes_gain = Enzyme_Production_v - Enzyme_Loss_v
        self.Enzyme_Production = Enzyme_Production_v
        self.Enzyme_Loss = Enzyme_Loss_v

        EnzLoss_2d = Enzyme_Loss_v.reshape(G, self.n_enzymes)
        N_cost_enz, P_cost_enz = self._enz_cost_pool
        DeadEnz_C = EnzLoss_2d.sum(axis=1)
        DeadEnz_N = (EnzLoss_2d * N_cost_enz[None, :]).sum(axis=1)
        DeadEnz_P = (EnzLoss_2d * P_cost_enz[None, :]).sum(axis=1)
        dead_cols = list(self.Substrates.columns[:3])
        self.Substrates.loc[is_deadEnz, dead_cols] = (
            self.Substrates.loc[is_deadEnz, dead_cols].to_numpy(dtype=np.float32, copy=True)
            + np.stack([DeadEnz_C, DeadEnz_N, DeadEnz_P], axis=1)
        )

    def _update_cue(self, rebalance_C_loss=None):
        """CUE = (C taken up - C respired for energy - enzyme C - C eroded for energy deficits - C respired by
        stoichiometric rebalancing) / C taken up, per cell (CUE_taxon) and for the grid (CUE_system); values below 0
        are set to 0. Cells without C uptake get 0. Dead-cell C (necromass) is not a loss here."""
        retained = self._cue_C_retained if rebalance_C_loss is None else self._cue_C_retained - rebalance_C_loss
        up = self._cue_C_uptake
        cue = np.zeros_like(up)
        pos = up > 0
        cue[pos] = retained[pos] / up[pos]
        self.CUE_taxon = np.maximum(cue, np.float32(0))
        tot = float(up.sum(dtype=np.float64))
        self.CUE_system = np.float32(0) if tot == 0 else np.float32(max(float(retained.sum(dtype=np.float64)) / tot, 0.0))

    def _clip_pools(self, C, N):
        """Keep the bound pools inside the cell: bound C <= biomass C and bound N <= biomass N (pools of empty cells -> 0).

        Production never violates this by construction; it is needed after mortality, where stoichiometric
        rebalancing can remove biomass C or N from a living cell."""
        hsp_C = self._hsp_pool.sum(axis=1)
        osmo_C = self._osmo_pool.sum(axis=1)
        bC = hsp_C + osmo_C
        bN = np.float32(self.HSP_N_cost) * hsp_C + np.float32(self.Osmo_N_cost) * osmo_C
        f = np.ones_like(bC)
        over_C = bC > C
        over_N = bN > N
        f[over_C] = C[over_C] / bC[over_C]
        f[over_N] = np.minimum(f[over_N], N[over_N] / bN[over_N])
        f[C <= 0] = 0
        if (f < 1).any():
            self._hsp_pool *= f[:, None]
            self._osmo_pool *= f[:, None]

    @property
    def Taxon_HSP(self):
        """HSP pools (C units) per cell and gene, with the cell's biomass (C+N+P) -- same layout output.py reads."""
        df = pd.DataFrame(self._hsp_pool, index=self.Microbes.index, columns=self._hsp_cols)
        df['Total_biomass'] = self.Microbes[['C', 'N', 'P']].to_numpy(dtype=np.float32).sum(axis=1)
        return df

    @property
    def Taxon_Osmolyte(self):
        """Osmolyte pools (C units) per cell and gene."""
        return pd.DataFrame(self._osmo_pool, index=self.Microbes.index, columns=self._osmo_cols)

    def mortality(self, day):
        """
        Calculate microbial mortality and update stoichiometry of the alive and microbial pools.

        Kill microbes that are starving deterministically and microbes that are drought intolerant stochastically
        Also update Substrates with input from dead microbes, monomers(with leaching loss), and respiration
        """

        # Leaching (abiotic NH4 and PO4 loss) parameters
        Leaching = np.float32(self.Leaching)
        Psi_slope_leach = np.float32(self.Leaching_psi_slope)
        # Indices
        Mic_index = self.Microbes.index
        is_DeadMic = self.Substrates.index == "DeadMic"
        is_NH4 = self.Monomers.index == "NH4"
        is_PO4 = self.Monomers.index == "PO4"

        n_taxa, gridsize = self.n_taxa, self.gridsize

        # NOTE: this rewrite avoids pandas boolean-masked assignment on the
        # (gridsize*n_taxa)-length Microbes.index (which repeats each taxon
        # label gridsize times) by pulling the large per-taxon-per-cell
        # structures into plain numpy arrays up front and doing everything
        # with positional numpy indexing -- the same pattern used to fix the
        # catastrophic non-unique-index join in uptake()'s CAE_taxon line,
        # and to speed up metabolism(). self.Microbes/self.Taxon_Resp are
        # reconstructed as pandas objects with the ORIGINAL Mic_index only at
        # the very end, so nothing downstream needs to change.
        Mic = self.Microbes.to_numpy(dtype=np.float32, copy=True)              # (N,4): C,N,P,E
        MinRat_cnp = self.MinRatios[['C', 'N', 'P']].to_numpy(dtype=np.float32, copy=True)   # (N,3)
        Death = np.zeros_like(Mic)
        kill = np.zeros(n_taxa * gridsize, dtype=bool)

        # Start of calcualtion of mortality first with THRESHOLD
        # Kill microbes deterministically based on threshold values: C_min; N_min; P_min
        starve_index = (Mic[:, 0] > 0) & (
                (Mic[:, 0] < self.C_min) | (Mic[:, 1] < self.N_min) | (Mic[:, 2] < self.P_min))
        # Index the dead, put them in Death, and set them to 0 in Microbes
        Death[starve_index] = Mic[starve_index]
        Mic[starve_index] = np.float32(0)
        self._hsp_pool[starve_index] = np.float32(0)   # pool C/N is part of the dead biomass (-> DeadMic)
        self._osmo_pool[starve_index] = np.float32(0)
        # Index the locations where microbial cells remain alive
        mic_index = Mic[:, 0] > 0

        # Kill microbes stochastically based on mortality prob as a function of water potential and drought tolerance
        # call the function, MMP:microbe_mortality_psi()
        # Drought tolerance of each cell from its current osmolyte composition (saturating in the osmolyte share
        # of biomass, same form as the HSP effect on maintenance Ea): tol = s / (s + Osmo_tol_K), s = osmolyte C / (C+N+P)
        cnp = Mic[:, :3].sum(axis=1)
        osmo_share = np.zeros_like(cnp)
        has_mass = cnp > 0
        osmo_share[has_mass] = self._osmo_pool.sum(axis=1)[has_mass] / cnp[has_mass]
        self.Drought_tol = osmo_share / (osmo_share + np.float32(self.Osmo_tol_K))
        r_death = np.asarray(MMP(self.psi[day], self.wp_fc, self.basal_death_prob, self.death_rate, self.Drought_tol))
        kill[mic_index] = r_death[mic_index] > np.random.uniform(0, 1, int(mic_index.sum())).astype('float32')
        # Index the dead, put them in Death, and set them to 0 in Microbes
        Death[kill] = Mic[kill]
        Mic[kill] = np.float32(0)
        self._hsp_pool[kill] = np.float32(0)
        self._osmo_pool[kill] = np.float32(0)
        # Index locations where microbes remain alive
        mic_index = Mic[:, 0] > 0

        # Calculate the total dead mass (threshold & drought) across taxa in each grid cell
        Death_gridcell = Death[:, :3].reshape(gridsize, n_taxa, 3).sum(axis=1)  # (gridsize,3): C,N,P

        # Distinguish between conditions of complete death VS partial death
        # All cells die
        if not mic_index.any():

            # ...Update Substrates pool by adding dead microbial biomass
            self.Substrates.loc[is_DeadMic] += Death_gridcell
            # Respiration of the day (uptake CO2 and any deficit-erosion CO2; no rebalancing losses). Previously left
            # untouched here, so the last value was repeated every day once all cells were dead.
            self.Taxon_Resp = pd.Series(self.Taxon_resp_df.to_numpy(copy=True).ravel(order='C').astype(np.float32),
                                        index=Mic_index)

        # Partly die and adjust stoichiometry of those remaining alive
        else:

            # Index only those taxa in Microbes that have below-minimum quotas
            idx = np.nonzero(mic_index)[0]
            Mic_alive = Mic[idx]                                    # (k,4)
            row_sum_cnp = Mic_alive[:, :3].sum(axis=1)              # C+N+P (E excluded), matches .drop('E',axis=1).sum(axis=1)
            with np.errstate(invalid='ignore', divide='ignore'):
                ratios_cnp = Mic_alive[:, :3] / row_sum_cnp[:, None]
            MinRat_alive = MinRat_cnp[idx]                           # (k,3)
            below_min = (ratios_cnp[:, 0] < MinRat_alive[:, 0]) | \
                        (ratios_cnp[:, 1] < MinRat_alive[:, 1]) | \
                        (ratios_cnp[:, 2] < MinRat_alive[:, 2])

            # Rebalance the cells whose C, N or P share of biomass (C + N + P) is below its minimum. Energy is not part of
            # cellular stoichiometry and is left unchanged. As in DEMENTpy: deficient elements are set to their minimum
            # share, the non-deficient elements give up that deficit in proportion to how far they are above their
            # minimum, and the cell keeps all of its most limiting element (largest relative deficit); the other two
            # are reduced to the new ratios relative to it (losses: C respired, N to NH4, P to PO4).
            idx2 = idx[below_min]                                    # positions needing rebalancing
            Mic_subset = Mic[idx2, :3]                               # (m,3): C,N,P
            StartMicrobes = Mic[idx2].copy()                         # (m,4), for the losses below
            MinRat3 = MinRat_cnp[idx2]                               # (m,3)
            with np.errstate(invalid='ignore', divide='ignore'):
                MicrobeRatios3 = Mic_subset / Mic_subset.sum(axis=1)[:, None]
            Ratio_dif3 = MicrobeRatios3 - MinRat3
            Ratio_dif_0_3 = np.where(Ratio_dif3 > 0, np.float32(0), Ratio_dif3)   # deficits (<= 0)
            Excess3 = np.where(Ratio_dif3 < 0, np.float32(0), Ratio_dif3)         # excesses (>= 0)

            # Limiting element: largest relative deficit (first occurrence of the maximum, as pandas idxmax)
            with np.errstate(invalid='ignore', divide='ignore'):
                metric3 = -Ratio_dif3 / MinRat3
            limiting_pos = np.argmax(np.where(np.isnan(metric3), -np.inf, metric3), axis=1)   # 0=C, 1=N, 2=P

            # Deficient shares to their minima; excess shares reduced by the total deficit, in proportion to the excess
            mask_deficient3 = Ratio_dif3 < 0
            mask_excess3 = Ratio_dif3 > 0
            MicrobeRatios3 = MicrobeRatios3.copy()
            MicrobeRatios3[mask_deficient3] = MinRat3[mask_deficient3]
            with np.errstate(invalid='ignore', divide='ignore'):
                scale = Ratio_dif_0_3.sum(axis=1) / Excess3.sum(axis=1)
                MicrobeRatios3[mask_excess3] += (Excess3 * scale[:, None])[mask_excess3]

            # Nutrient quotas scaled to each possible limiting element, then pick the limiting one
            MC, MN, MP = Mic_subset[:, 0], Mic_subset[:, 1], Mic_subset[:, 2]
            MRC, MRN, MRP = MicrobeRatios3[:, 0], MicrobeRatios3[:, 1], MicrobeRatios3[:, 2]
            with np.errstate(invalid='ignore', divide='ignore'):
                new_C = np.stack([MC, MN * MRC / MRN, MP * MRC / MRP], axis=1)
                new_N = np.stack([MC * MRN / MRC, MN, MP * MRN / MRP], axis=1)
                new_P = np.stack([MC * MRP / MRC, MN * MRP / MRN, MP], axis=1)
            for arr in (new_C, new_N, new_P):
                arr[~np.isfinite(arr)] = np.float32(0)
            rows = np.arange(len(idx2))
            C_final = new_C[rows, limiting_pos]
            N_final = new_N[rows, limiting_pos]
            P_final = new_P[rows, limiting_pos]

            # Update Microbes (E unchanged)
            Mic[idx2, 0] = C_final
            Mic[idx2, 1] = N_final
            Mic[idx2, 2] = P_final

            # Sum up the element losses from biomass across whole grid and calculate average loss
            MicLoss = StartMicrobes - Mic[idx2]
            tax_resp_mort = self.Taxon_resp_df.to_numpy(copy=True).ravel(order='C').astype(np.float32)
            tax_resp_mort[idx2] += MicLoss[:, 0]
            rebalance_C_loss = np.zeros(n_taxa * gridsize, dtype=np.float32)
            rebalance_C_loss[idx2] = MicLoss[:, 0]
            self._update_cue(rebalance_C_loss)            # CUE after the C respired by rebalancing
            self.Taxon_Resp = pd.Series(tax_resp_mort, index=Mic_index)
            # Update monomer pools
            self.Monomers.loc[is_NH4, "N"] += np.float32(MicLoss[:, 1].sum() / gridsize)
            self.Monomers.loc[is_PO4, "P"] += np.float32(MicLoss[:, 2].sum() / gridsize)

            # Update Substrates pool by adding dead microbial biomass
            self.Substrates.loc[is_DeadMic] += Death_gridcell
        # End of if else clause

        # Cap the energy:C ratio at E_max_ratio_b / _f (energy above the cap is dissipated)
        with np.errstate(invalid='ignore', divide='ignore'):
            Tax_E_ratio = Mic[:, 3] / Mic[:, 0]
        e_cap_index = Tax_E_ratio > self.E_max
        Mic[e_cap_index, 3] = Mic[e_cap_index, 0] * self.E_max[e_cap_index]

        # Leaching of monomers
        Leaching_N = self.Monomers.loc[is_NH4, "N"] * Leaching * np.exp(Psi_slope_leach * (self.psi[day] - self.wp_fc))
        Leaching_P = self.Monomers.loc[is_PO4, "P"] * Leaching * np.exp(Psi_slope_leach * (self.psi[day] - self.wp_fc))
        # Update Monomers
        self.Monomers.loc[is_NH4, "N"] -= Leaching_N
        self.Monomers.loc[is_PO4, "P"] -= Leaching_P
        self.Leached_N = np.asarray(Leaching_N, dtype=np.float32)  # per grid cell, this day
        self.Leached_P = np.asarray(Leaching_P, dtype=np.float32)

        # Rebalancing may have removed C or N from living cells: keep the bound pools inside the cell
        self._clip_pools(Mic[:, 0], Mic[:, 1])

        # Restore Microbes to a DataFrame with the original taxon-repeated index
        self.Microbes = pd.DataFrame(Mic, index=Mic_index, columns=['C', 'N', 'P', 'E'])

        # Update the death toll of cells
        self.Kill = np.uint32(kill.sum())

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
        n_taxa, gridsize = self.n_taxa, self.gridsize
        N = n_taxa * gridsize

        # NOTE: rewritten to operate on plain numpy arrays throughout instead
        # of pandas Series/DataFrames indexed by the (gridsize*n_taxa)-length
        # Mic_index (which repeats each taxon label gridsize times). Boolean
        # masked assignment on that duplicated index is both slow (general
        # pandas per-op overhead) and, under pandas 3.0, can fall into an
        # expensive non-unique-index join in some cases (see uptake()'s
        # CAE_taxon fix) -- positional numpy indexing sidesteps both issues.
        # Intracellular pools travel with biomass: appended as extra columns so that fungal translocation
        # (averaging), cell division (halving) and dispersal/colonization (adding) treat them like C, N, P and E.
        n_hsp = self._hsp_pool.shape[1]
        Mic = np.hstack([self.Microbes.to_numpy(dtype=np.float32, copy=True),
                         self._hsp_pool, self._osmo_pool])          # (N, 4 + n_HSP + n_osmo)
        K = Mic.shape[1]
        fb = np.asarray(self.fb)                                    # (N,) 0/1

        # STEP 1: Fungal translocation by calculating average biomass within fungal taxa
        # Count the fungal taxa before cell division
        # Add one or two fungi to the count series based on size
        C0 = Mic[:, 0]
        is_fungal_alive = (fb == 1) & (C0 > 0)
        is_fungal_big = (fb == 1) & (C0 > self.max_size_f)
        fungi_flags = np.zeros(N, dtype=np.int64)
        fungi_flags[is_fungal_alive] = 1
        fungi_flags[is_fungal_big] = 2
        # dtype is intentionally not int8 here because the sum over gridsize
        # cells can exceed int8 max (127) on large grids.
        Fungi_count = fungi_flags.reshape(gridsize, n_taxa).sum(axis=0)   # (n_taxa,)

        # Derive average biomass of fungal taxa
        Microbes_grid = Mic.reshape(gridsize, n_taxa, K).sum(axis=0)      # (n_taxa,K)
        with np.errstate(invalid='ignore', divide='ignore'):
            Mean_fungi = Microbes_grid / Fungi_count[:, None]
        Mean_fungi[Fungi_count == 0] = np.float32(0)
        # Expand the fungal average across the grid
        eMF = np.tile(Mean_fungi, (gridsize, 1)).astype(np.float32)      # (N,K)

        # STEP 2: Cell division & translocate nutrients
        MicrobesBeforeDivision = Mic.copy()
        # bacterial cell division
        bac_index = (fb == 0) & (Mic[:, 0] > self.max_size_b)
        Mic[bac_index] = Mic[bac_index] / 2
        # fungal cell division
        fun_index = (fb == 1) & (Mic[:, 0] > self.max_size_f)
        Mic[fun_index] = Mic[fun_index] / 2
        # Put daughter cells into a seperate array, Reprod
        Reprod = MicrobesBeforeDivision - Mic
        # Translocate nutrients within fungal taxa after reproduction
        fungal_alive_now = (fb == 1) & (Mic[:, 0] > 0)
        Mic[fungal_alive_now] = eMF[fungal_alive_now]
        # Index the daughter cells of fungi vs bacteria
        daughters_b = (Reprod[:, 0] > 0) & (fb == 0)
        daughters_f = (Reprod[:, 0] > 0) & (fb == 1)
        # set all fungi equal to their grid averages for translocation before colonization
        Reprod[daughters_f] = eMF[daughters_f]

        # STEP 3: dispersal calculation
        num_b = int(daughters_b.sum())
        num_f = int(daughters_f.sum())
        shift_x = np.zeros(N, dtype=np.int8)
        shift_y = np.zeros(N, dtype=np.int8)
        # Bacterial dispersal movements in X & Y direction
        shift_x[daughters_b] = np.random.choice(
            np.arange(-self.dist, self.dist + 1), num_b, replace=True).astype('int8')
        shift_y[daughters_b] = np.random.choice(
            np.arange(-self.dist, self.dist + 1), num_b, replace=True).astype('int8')
        # Fungi always move positively in x direction, and in y direction constrained to one box away determined by probability "direct"
        shift_x[daughters_f] = np.int8(1)
        shift_y[daughters_f] = np.random.choice(
            [-1, 0, 1], num_f, replace=True,
            p=[0.5 * (1 - self.direct), self.direct, 0.5 * (1 - self.direct)]).astype('int8')
        # Calculate x,y coordinates of dispersal destinations (% remainder of x/x)
        x_coords = np.tile(np.repeat(np.arange(1, self.x + 1), n_taxa), self.y)
        y_coords = np.repeat(np.arange(1, self.y + 1), n_taxa * self.x)
        new_x = (x_coords.astype(np.int64) + shift_x + self.x) % self.x
        new_y = (y_coords.astype(np.int64) + shift_y + self.y) % self.y
        new_x[new_x == 0] = self.x  # Substitute coordinates when there is no shift
        new_y[new_y == 0] = self.y  # Substitute coordinates when there is no shift
        # Convert x,y coordinates to destination locations
        taxa_offset = np.tile(np.arange(1, n_taxa + 1), gridsize)
        index_series = (n_taxa * ((new_y - 1) * self.x + (new_x - 1))) + taxa_offset - 1
        index_series = index_series.astype(np.int64)

        # Step 4: colonization of dispersed microbes
        # Use np.add.at to accumulate biomass, correctly handling multiple daughters landing on the same grid location.
        Colonization = np.zeros((N, K), dtype=np.float32)
        np.add.at(Colonization, index_series[daughters_b], Reprod[daughters_b])
        np.add.at(Colonization, index_series[daughters_f], Reprod[daughters_f])
        # Colonization of dispersing microbes
        Mic += Colonization

        self.Microbes = pd.DataFrame(Mic[:, :4], index=Mic_index, columns=['C', 'N', 'P', 'E'])
        self._hsp_pool = np.ascontiguousarray(Mic[:, 4:4 + n_hsp])
        self._osmo_pool = np.ascontiguousarray(Mic[:, 4 + n_hsp:])

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
        self._hsp_pool[:] = np.float32(0)
        self._osmo_pool[:] = np.float32(0)

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
