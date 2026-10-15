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
        return pd.Series(
            result, index=index, dtype=arr.dtype if hasattr(arr, "dtype") else None
        )
    else:
        result = vals.reshape(gridsize, n_taxa, vals.shape[1]).sum(axis=0)
        return pd.DataFrame(result, index=index, columns=arr.columns)


def _enz_sum(arr, gridsize, n_enzymes, index):
    """Sum a (gridsize*n_enzymes,) Series over the grid dimension."""
    vals = arr.values if hasattr(arr, "values") else np.asarray(arr)
    result = vals.reshape(gridsize, n_enzymes).sum(axis=0)
    return pd.Series(
        result, index=index, dtype=arr.dtype if hasattr(arr, "dtype") else None
    )


def _mon_sum(arr, gridsize, n_monomers, index):
    """Sum a (gridsize*n_monomers, n_taxa) DataFrame over the grid dimension."""
    vals = arr.values if hasattr(arr, "values") else np.asarray(arr)
    if vals.ndim == 1:
        result = vals.reshape(gridsize, n_monomers).sum(axis=0)
        return pd.Series(result, index=index)
    else:
        result = vals.reshape(gridsize, n_monomers, vals.shape[1]).sum(axis=0)
        return pd.DataFrame(result, index=index, columns=arr.columns)


class Grid:
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

        self.cycle = int(runtime.loc["end_time", 1])
        self.gridsize = int(runtime.loc["gridsize", 1])
        self.n_taxa = int(runtime.loc["n_taxa", 1])
        self.n_substrates = int(runtime.loc["n_substrates", 1])
        self.n_enzymes = int(runtime.loc["n_enzymes", 1])
        self.n_monomers = self.n_substrates + 3
        self.kinetics_model = int(
            runtime.loc["Kinetics_Model", 1]
        )  # specifies use of Arrhenius or MMRT to calculate rxn rates
        # self.kinetics_model = runtime.loc["Kinetics_Model", 1]  # Arrhenius or MMRT to calculate rxn rates

        # Degradation
        self.Substrates_init = data_init["Substrates"]  # Substrates initialized
        self.Substrates = data_init["Substrates"].copy(
            deep=True
        )  # Substrates;df; w/ .copy() avoiding mutation
        self.SubInput = data_init["SubInput"]  # Substrate inputs
        self.Enzymes_init = data_init["Enzymes"]  # Initial pool of Enzymes
        self.Enzymes = data_init["Enzymes"].copy(deep=True)  # Enzymes
        self.ReqEnz = data_init["ReqEnz"]  # Enzymes required by each substrate
        self.Ea = data_init["Ea"]  # Enzyme activation energy
        self.HT0 = data_init["Vmax_HT0"]  # Enzyme enthalpy
        self.ST0 = data_init["Vmax_ST0"]  # Enzyme entropy
        self.Cp = data_init["Vmax_Cp"]  # Enzyme heat capacity
        self.Km_HT0 = data_init["Km_HT0"]  # Enzyme enthalpy (for Km calculation)
        self.Km_ST0 = data_init["Km_ST0"]  # Enzyme entropy (for Km calculation)
        self.Km_Cp = data_init["Km_Cp"]  # Enzyme heat capacity (for Km calculation)
        self.Vmax0 = data_init[
            "Vmax0"
        ]  # Max. reaction speed # not used in MMRT model (Jeth)
        self.Km0 = data_init[
            "Km0"
        ]  # Half-saturation constant # not used in MMRT model (Jeth)
        self.SubstrateRatios = np.float32("nan")  # Substrate stoichiometry
        self.DecayRates = np.float32("nan")  # Substrate decay rate
        self.Vmax_Km_slope = data_init[
            "mmrt_Vmax_Km_slope"
        ]  # slope of the enzyme Vmax/Km relationship
        self.Vmax_Km_int = data_init[
            "mmrt_Vmax_Km_int"
        ]  # intercept of the Uptake Vmax/Km relationship
        self.Enz_Vmax_sub_spec = data_init[
            "Enz_vmax_sub_spec"
        ]  # substrate specificity-efficiency adjustment for degradation Vmax
        self.Enz_Km_sub_spec = data_init[
            "Enz_km_sub_spec"
        ]  # substrate specificity-efficiency adjustment for degradation Km
        self.Enz_Km_Temp_Sens = data_init[
            "Enz_km_temp_sens"
        ]  # substrate specificity-efficiency adjustment for monomer uptake Vmax
        self.Arrh_Tref = data_init["Arrhenius_tref"]

        # Uptake
        self.Microbes_init = data_init[
            "Microbes_pp"
        ]  # microbial community before placement
        self.Microbes = data_init["Microbes"].copy(
            deep=True
        )  # microbial community after placement
        self.Monomers_init = data_init["Monomers"]  # Monomers initialized
        self.Monomers = data_init["Monomers"].copy(deep=True)  # Monomers
        self.MonInput = data_init["MonInput"]  # Inputs of monomers
        self.Uptake_Ea = data_init["Uptake_Ea"]  # transporter enzyme Ea
        self.Uptake_Vmax0 = data_init["Uptake_Vmax0"]  # transporter Vmax
        self.Uptake_Km0 = data_init["Uptake_Km0"]  # transporter Km
        self.Uptake_HT0 = data_init["Uptake_HT0"]  # transporter HT0
        self.Uptake_ST0 = data_init["Uptake_ST0"]  # transporter ST0
        self.Uptake_Cp = data_init["Uptake_Cp"]  # transporter Cp
        self.Monomer_ratios = data_init["Monomer_ratio"].copy(
            deep=True
        )  # monomer stoichiometry
        self.Monomer_energy = data_init[
            "MonEnergy_Attrib"
        ]  # monomer energy characteristics
        self.Mon_energy_resp = expand(
            pd.concat(
                [
                    pd.Series([0, 0], index=["NH4", "PO4"]),
                    self.Monomer_energy["E_resp"],
                ],
                sort=False,
            ),
            self.gridsize,
        )
        self.Mon_energy_ferm = expand(
            pd.concat(
                [
                    pd.Series([0, 0], index=["NH4", "PO4"]),
                    self.Monomer_energy["E_ferm"],
                ],
                sort=False,
            ),
            self.gridsize,
        )
        self.Uptake_ReqEnz = data_init["Uptake_ReqEnz"]  # Enzymes required by monomers
        self.Uptake_Enz_Cost = data_init[
            "UptakeGenesCost"
        ]  # Cost of encoding each uptake gene
        self.Uptake_Prod_Cost = data_init[
            "UptakeProdCost"
        ]  # C Cost of transporter production
        self.UptakeProd_E_Cost = data_init[
            "Uptake_prod_E_cost"
        ]  # Energy cost of transporter production per mg C
        self.Taxon_Uptake_C = np.float32("nan")  # taxon total C uptake
        self.Taxon_Uptake_BC = np.float32(
            "nan"
        )  # taxon uptake of C allocated to biomass
        self.Taxon_Uptake_N = np.float32("nan")  # taxon uptake of N
        self.Taxon_Uptake_P = np.float32("nan")  # taxon uptake of P
        self.Taxon_Uptake_E = np.float32("nan")  # taxon uptake of Energy
        self.Microbe_E_Gain = np.float32("nan")  # taxon net energy gain
        self.Microbe_E_Cost = np.float32("nan")  # taxon net energy costs
        self.Uptake_Vmax_Km_slope = data_init[
            "uptake_mmrt_Vmax_Km_slope"
        ]  # slope of the Uptake Vmax/Km relationship
        self.Uptake_Vmax_Km_int = data_init[
            "uptake_mmrt_Vmax_Km_int"
        ]  # intercept for the Uptake Vmax/Km relationship
        self.Uptake_sub_spec = data_init[
            "Uptake_Vmax_sub_spec"
        ]  # substrate specificity-efficiency adjustment for monomer uptake Vmax
        self.Uptake_Km_Temp_Sens = data_init[
            "Uptake_km_temp_sens"
        ]  # substrate specificity-efficiency adjustment for monomer uptake Vmax
        self.km_Tref = data_init[
            "km_tref"
        ]  # substrate specificity-efficiency adjustment for monomer uptake Vmax
        self.Pox_Ferm = data_init[
            "Fermoxp"
        ]  # proportion of fermented monomer converted to CO2
        # fermentation condition score (0-1): reflects environmental favorability of fermentation vs respiration
        self.ferm_cond = data_init[
            "FermCond"
        ]  # currently set to 0.1; can update to vary with environmental conditions
        self.active_uptake_E_cost = data_init[
            "ActiveUptakeE_cost"
        ]  # energy cost of transporting monomers

        # Metabolism
        self.Consti_Enzyme_C = data_init[
            "EnzProdConstit"
        ]  # C cost of encoding constitutive enzyme
        self.Induci_Enzyme_C = data_init[
            "EnzProdInduce"
        ]  # C Cost of encoding inducible enzyme
        self.Consti_Osmo_C = data_init[
            "OsmoProdConsti"
        ]  # C Cost of encoding constitutive osmolyte
        self.Induci_Osmo_C = data_init[
            "OsmoProdInduci"
        ]  # C Cost of encoding inducible osmolyte
        self.Enz_Attrib = data_init["EnzAttrib"]  # Enzyme attributes; dataframe
        # self.cell_maint_E_cost = data_init['Cell_maint_E_cost']  # energy cost of biomass maintenance
        self.Osmo_E_cost = data_init[
            "OsmoProd_E_cost"
        ]  # E cost per unit of osmo-C production
        self.Osmo_N_cost = data_init[
            "OsmoProd_N_cost"
        ]  # N cost per unit of osmo-C production
        self.HSP_E_cost = data_init[
            "HSPProd_E_cost"
        ]  # E cost per unit of HSP-C production
        self.HSP_N_cost = data_init[
            "HSPProd_N_cost"
        ]  # N cost per unit of HSP-C production
        # self.Osmo_Maint_cost_constant = data_init['Osmo_Maint_cost_constant']       # Respiration cost constant for osmolyte production
        self.Respiration = np.float32("nan")  # Respiration
        self.Taxon_Resp = np.float32("nan")  # Taxon specific Respiration
        self.Taxon_C_Uptake = np.float32("nan")  # Taxon-specific total C uptake
        self.Taxon_E_Uptake = np.float32("nan")  # Taxon-specific E production
        self.Taxon_E_Prod = np.float32("nan")  # Taxon-specific E production
        self.Taxon_BC_Uptake = np.float32(
            "nan"
        )  # Taxon-specific C uptake allocated to biomass
        self.Taxon_Total_Uptake = np.float32(
            "nan"
        )  # Taxon-specific biomass uptake (C+N+P)
        self.Enzyme_Loss_Rate = data_init["Enzyme_Loss_Rate"]
        self.OSmo_Loss_Rate = data_init["Osmo_Loss_Rate"]
        self.maint_E_HSP_Ea_slope = data_init["HSP_Ea_slope"]
        self.maint_E_HSP_Ea_int = data_init["HSP_Ea_int"]
        self.maint_E_HSP_Ea_n = data_init["HSP_Ea_n"]
        self.maint_E_rate_Tref = data_init["maintE_rate_Tref"]
        self.maint_E_rate_const = data_init["maintE_rate_const"]
        self.HSP_consti_pre_exp = data_init[
            "HSPProdConsti"
        ]  # pre exponential factor for temperature induced constitutive HSP production
        self.HSP_induci_pre_exp = data_init[
            "HSPProdInduci"
        ]  # pre exponential factor for temperature induced inducible HSP production
        self.HSP_consti_Ea = data_init[
            "HSP_con_Ea"
        ]  # exponential factor for temperature induced constitutive HSP production
        self.HSP_induci_Ea = data_init[
            "HSP_ind_Ea"
        ]  # exponential factor for temperature induced inducible HSP production
        # self.Transporters = float('nan')
        # self.Osmolyte_Con = float('nan')
        # self.Osmolyte_Ind = float('nan')
        # self.Enzyme_Con   = float('nan')
        # self.Enzyme_Ind   = float('nan')
        self.HSP_Turnover = data_init["HSPturnover"]
        self.Taxon_HSP = data_init["Initial_HSP_Prop"]  #

        # Osmolyte pool: same structure as constitutive osmolyte production DataFrame
        # (gridsize*n_taxa rows × osmolyte-gene columns), initialized to zero.
        # Updated in metabolism() each step with production and turnover.
        self.Taxon_Osmolyte = data_init["OsmoProdConsti"].copy().astype("float32")
        self.Taxon_Osmolyte[:] = np.float32(0)

        # Mortality
        self.MinRatios = data_init["MinRatios"]  # Minimal cell quotas
        self.C_min = data_init["C_min"]  # C threshold value of living cell
        self.N_min = data_init["N_min"]  # N threshold value of living cell
        self.P_min = data_init["P_min"]  # P threshold value of living cell
        self.E_range = data_init["Erange"]  # tolerance on energy:biomass C ratio
        self.E_ratb = data_init["Eratb"]  # ideal E:MBC ratio for bacteria
        self.E_ratf = data_init["Eratf"]  # ideal E:MBC ratio for fungi
        self.basal_death_prob = data_init[
            "basal_death_prob"
        ]  # Basal death probability of microbes
        self.death_rate = data_init[
            "death_rate"
        ]  # Change rate of mortality with water potential
        self.tolerance = data_init["TaxDroughtTol"]  # taxon drought tolerance
        self.wp_fc = data_init["wp_fc"]  # -1.0
        self.wp_th = data_init["wp_th"]  # -6.0
        self.alpha = data_init["alpha"]  # integer;1
        self.Kill = np.float32("nan")  # total number of cells stochastically killed

        # Reproduction
        self.fb = data_init["fb"]  # index of fungal taxa (=1)
        self.max_size_b = data_init["max_size_b"]  # threshold of cell division
        self.max_size_f = data_init["max_size_f"]  # threshold of cell division
        self.x = int(runtime.loc["x", 1])  # x dimension of grid
        self.y = int(runtime.loc["y", 1])  # y dimension of grid
        self.dist = int(runtime.loc["dist", 1])  # maximum dispersal distance: 1 cell
        self.direct = int(runtime.loc["direct", 1])  # dispersal direction: 0.95

        # Climate data
        self.temp = data_init["Temp"]  # series; temperature
        self.psi = data_init["Psi"]  # series; water potential

        # Global constants
        self.Km_Ea = data_init[
            "km_ea"
        ]  # kj mol-1;activation energy for both enzyme and transporter
        self.MMRT_Tref = data_init[
            "mmrt_tref"
        ]  # reference temperature for calculating rates (20 celcius)
        self.gas_const = np.float32(const.R)  # Universal gas constant
        self.boltz_const = np.float32(const.k)  # Boltzmann's constant
        self.planck_const = np.float32(const.h)  # Planck's constant
        self.Substrates_id = data_init["Substrate_index"]
        self.monomer_index = data_init["Monomer_index"]
        # self.monomer_index = ["NH4", "PO4", "FermProd", "DeadMic", "DeadEnz"] + ["Mon" + str(i) for i in range(4, self.n_monomers - 2 + 1)]

        # ── Fast groupby helpers ──────────────────────────────────────────────
        # These replace pd.DataFrame.groupby(level=0).sum() calls throughout
        # the model.  The MultiIndex structure is fixed for the entire run:
        #   taxa-indexed arrays:   gridsize repetitions of n_taxa taxon labels
        #   enzyme-indexed arrays: gridsize repetitions of n_enzymes labels
        # Reshaping to (gridsize, n) and summing axis=0 is equivalent to
        # groupby(level=0).sum() but ~60x faster (no index hashing required).

        # Pre-built index arrays for reconstructing output Series/DataFrames
        self._taxa_index = self.Microbes.index[: self.n_taxa]
        self._enz_index = self.Enzymes.index[: self.n_enzymes]

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

        G = self.gridsize
        M = (
            self.n_substrates
        )  # note: local var name mirrors upstream n_substrates convention
        NE = self.n_enzymes

        def _v(x):
            return (
                x.to_numpy(dtype=np.float32, copy=True) if hasattr(x, "to_numpy") else x
            )

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
        with np.errstate(divide="ignore", invalid="ignore"):
            # np.divide's `where=` guard still evaluates the masked-out branch
            # internally (a known numpy quirk), which would otherwise emit a
            # spurious "invalid value encountered" warning even though those
            # positions are never written to `out` -- harmless, just noisy.
            np.divide(
                Substrates_v,
                rss[:, None],
                out=SubstrateRatios_v,
                where=(rss[:, None] != 0),
            )

        # moisture effects on enzymatic kinetics
        if self.psi[day] >= self.wp_fc:
            f_psi = np.float32(1.0)
        else:
            f_psi = np.exp(0.25 * (self.psi[day] - self.wp_fc)).astype("float32")

        # calculate enzyme Vmax using the user specified kinetics model (MMRT = 1; Arrhenius = 0)
        if self.kinetics_model > 0:
            Vmax_mmrt = MMRT(
                _v(self.Cp),
                _v(self.HT0),
                _v(self.ST0),
                self.temp[day],
                self.MMRT_Tref,
                self.boltz_const,
                self.planck_const,
                self.gas_const,
            )
            Vmax_mmrt_sub_spec = Vmax_mmrt * _v(self.Enz_Vmax_sub_spec)
            Vmax_v = Vmax_mmrt_sub_spec * f_psi  # (G*NE, M)
            if self.Enz_Km_Temp_Sens == 1:
                Km_mmrt = MMRT(
                    _v(self.Km_Cp),
                    _v(self.Km_HT0),
                    _v(self.Km_ST0),
                    self.temp[day],
                    self.MMRT_Tref,
                    self.boltz_const,
                    self.planck_const,
                    self.gas_const,
                )
            else:
                Km_mmrt = MMRT(
                    _v(self.Km_Cp),
                    _v(self.Km_HT0),
                    _v(self.Km_ST0),
                    self.km_Tref,
                    self.MMRT_Tref,
                    self.boltz_const,
                    self.planck_const,
                    self.gas_const,
                )
            Km_mmrt_sub_spec = Km_mmrt * _v(self.Enz_Km_sub_spec)
            Km_v = (
                Km_mmrt_sub_spec * self.Vmax_Km_slope
            ) + self.Vmax_Km_int  # (G*M, NE)
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
                self._ea_km_v = (
                    self._ea_v.reshape(G, NE, M).transpose(0, 2, 1).reshape(G * M, NE)
                )  # (G*M, NE)
            Vmax_v = (
                _v(self.Vmax0)
                * Arrhenius(self._ea_v, self.gas_const, self.temp[day], self.Arrh_Tref)
                * f_psi
            )  # (G*NE, M)
            Km_v = _v(self.Km0) * Arrhenius(
                self.Km_Ea, self.gas_const, self.temp[day], self.Arrh_Tref
            )  # (G*M, NE)

        Vmax_v = np.asarray(Vmax_v, dtype=np.float32)
        Km_v = np.asarray(Km_v, dtype=np.float32)

        # Multiply Vmax by enzyme concentration, then reshape (enz*grid, sub) -> (sub*grid, enz).
        # See docstring point 1: plain reshape/transpose, bit-exact vs. the original's
        # stack/unstack + explicit column reindex, and far cheaper.
        Enzymes_v = self.Enzymes.to_numpy(dtype=np.float32, copy=True)  # (G*NE,)
        tev_transition_v = Vmax_v * Enzymes_v[:, None]  # (G*NE, M)
        tev_v = (
            tev_transition_v.reshape(G, NE, M).transpose(0, 2, 1).reshape(G * M, NE)
        )  # (G*M, NE)

        # Michaelis-Menten equation
        Decay_v = (tev_v * rss[:, None]) / (Km_v + rss[:, None])  # (G*M, NE)

        if (
            getattr(self, "debug_degradation_days", 0)
            and day < self.debug_degradation_days
        ):
            print(f"  -- degradation() [numpy] shape probe (day {day}) --")
            print(
                f"    Vmax_v.shape={Vmax_v.shape}  Km_v.shape={Km_v.shape}  "
                f"tev_v.shape={tev_v.shape}  Decay_v.shape={Decay_v.shape}"
            )

        # Pull out each batch of required enzymes and sum across redundant enzymes.
        # self.ReqEnz never changes across days (fixed biochemical requirement table),
        # so the 'set1' slice + array extraction is cached on first use instead of
        # repeating a MultiIndex `.loc[]` selection every single day.
        if not hasattr(self, "_reqenz_set1_v"):
            self._reqenz_set1_v = self.ReqEnz.loc["set1"].to_numpy(
                dtype=np.float32, copy=True
            )  # (G*M, NE)
        batch1 = (self._reqenz_set1_v * Decay_v).sum(axis=1)  # (G*M,)

        # Compare to substrate available and take the min -- np.fmin (not np.minimum!)
        # to replicate pandas' `.min(axis=1, skipna=True)` NaN-as-missing semantics;
        # see docstring point 2.
        DecayRates = np.fmin(batch1, rss)  # (G*M,)

        # Adjust cellulose rate by linking cellulose degradation to lignin concentration (LCI)
        is_lignin = Sub_index_arr == "Lignin"
        is_cellulose = Sub_index_arr == "Cellulose"
        ss7 = Substrates_v[is_lignin].sum(
            axis=1
        )  # (G,) -- one Lignin row per grid cell
        cellulose_C = Substrates_v[
            is_cellulose, 0
        ]  # (G,) -- C column of the Cellulose row per cell
        with np.errstate(divide="ignore", invalid="ignore"):
            # Matches the original exactly: no zero-denominator guard here in
            # either version, so a genuine 0/0 (ss7 + cellulose_C == 0) still
            # produces NaN in both -- only the warning is suppressed, not the
            # underlying (unchanged, original) behavior.
            DecayRates[is_cellulose] *= (
                np.float32(1) + (ss7 / (ss7 + cellulose_C)) * LCI_slope
            )

        # Update Substrates Pool by removing decayed C, N, & P. Depending on specific needs, adding inputs of substrates can be done here
        Substrates_v -= SubstrateRatios_v * DecayRates[:, None]
        self.Substrates[:] = Substrates_v

        # Pass these two back to the global variables to be used in the next method
        self.SubstrateRatios = pd.DataFrame(
            SubstrateRatios_v,
            index=self.Substrates.index,
            columns=self.Substrates.columns,
        )
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
                self.Monomers.values.reshape(self.gridsize, self.n_monomers, 3).sum(
                    axis=0
                )
                / self.gridsize,
                index=self.monomer_index,
                columns=self.Monomers.columns,
            ),
            self.gridsize,
        )

        # Indices
        is_org = (
            (self.Monomers.index != "NH4")
            & (self.Monomers.index != "PO4")
            & (self.Monomers.index != "FermProd")
        )  # organic monomers
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
        self.Monomer_ratios.loc[is_org] = self.Monomers.loc[is_org].divide(
            rsm[is_org], axis=0
        )
        self.Monomer_ratios = self.Monomer_ratios.fillna(0)

        # Start calculating monomer uptake
        # Moisture impacts on uptake, mimicking the diffusivity implications
        if self.psi[day] >= self.wp_fc:
            f_psi = np.float32(1.0)
        else:
            f_psi = np.exp(0.5 * (self.psi[day] - self.wp_fc)).astype("float32")

        if self.kinetics_model > 0:
            Uptake_Vmax_mmrt = MMRT(
                self.Uptake_Cp,
                self.Uptake_HT0,
                self.Uptake_ST0,
                self.temp[day],
                self.MMRT_Tref,
                self.boltz_const,
                self.planck_const,
                self.gas_const,
            )
            Uptake_Vmax_sub_spec = (
                Uptake_Vmax_mmrt * self.Uptake_sub_spec
            )  # account for substrate specificity-efficiency tradeoff
            Uptake_Vmax = (
                Uptake_Vmax_sub_spec * f_psi
            )  # account for moisture limitation

            # If specified, derive Km_mmrt based on Vmax-Km tradeoff; otherwise, calculate Km assuming
            if self.Uptake_Km_Temp_Sens == 1:
                Uptake_Km = (
                    Uptake_Vmax_sub_spec * self.Uptake_Vmax_Km_slope
                ) + self.Uptake_Vmax_Km_int
            else:
                Uptake_Km_mmrt = MMRT(
                    self.Uptake_Cp,
                    self.Uptake_HT0,
                    self.Uptake_ST0,
                    self.km_Tref,
                    self.MMRT_Tref,
                    self.boltz_const,
                    self.planck_const,
                    self.gas_const,
                )
                Uptake_Km_sub_spec = Uptake_Km_mmrt * self.Uptake_sub_spec
                Uptake_Km = (
                    Uptake_Km_sub_spec * self.Uptake_Vmax_Km_slope
                ) + self.Uptake_Vmax_Km_int

        else:
            # NOTE (fix 2026-07-23): same fix as degradation()'s Arrhenius branch --
            # Uptake_Km used to be shrunk by the hardcoded self.Km_Ea = 20 instead of
            # Uptake_Ea, leaving Uptake_Km ~1e8 (barely shrunk from Uptake_Vmax0's raw
            # 1e7-1e8 scale) vs. realistic monomer masses -- silently suppressing
            # uptake to a negligible-but-technically-nonzero level. self.Uptake_Km0
            # and self.Uptake_Ea are both (monomer x transporter)-shaped already, so
            # (unlike degradation's Ea/Km0, which are transposed relative to each
            # other) no reshape is needed here.
            Uptake_Vmax = (
                self.Uptake_Vmax0
                * Arrhenius(
                    self.Uptake_Ea, self.gas_const, self.temp[day], self.Arrh_Tref
                )
                * f_psi
            )
            Uptake_Km = self.Uptake_Km0 * Arrhenius(
                self.Km_Ea, self.gas_const, self.temp[day], self.Arrh_Tref
            )

        # Equation for hypothetical potential uptake (per unit of compatible uptake protein)
        Potential_Uptake = (self.Uptake_ReqEnz * Uptake_Vmax).mul(
            rsm.values, axis=0
        ) / Uptake_Km.add(rsm.values, axis=0)

        # Derive the mass of each transporter of each taxon NOTE: transpose the df to Upt*(Taxa*grid)
        MicCXGenes = (
            self.Uptake_Enz_Cost.mul(self.Microbes[["C", "N", "P"]].sum(axis=1), axis=0)
        ).T

        # ── Uptake computation: numpy-native to minimise working-set size ──────────
        # All intermediate arrays are kept as (G, M, T) numpy arrays.
        # This avoids materialising 11 separate (G*M, T) Pandas DataFrames
        # simultaneously, cutting the peak working set from ~63 MB (30×30) /
        # ~700 MB (100×100) to ~20 MB / ~220 MB and reducing DRAM pressure.

        G, M, U, T = (
            self.gridsize,
            self.n_monomers,
            Potential_Uptake.shape[1],
            self.n_taxa,
        )

        # ── Step 1: Max_Uptake via batched matmul, then scale to available monomer ──
        P_3d = Potential_Uptake.values.reshape(G, M, U)  # (G, M, U)
        Mc_3d = MicCXGenes.values.reshape(U, G, T).transpose(1, 0, 2)  # (G, U, T)
        # max_up: (G, M, T)  — potential uptake per grid × monomer × taxon
        max_up = (P_3d @ Mc_3d).astype(np.float32)  # (G, M, T)

        # Scale uptake to available monomer mass (rsm: G*M vector)
        rsm_3d = rsm.values.reshape(G, M, 1).astype(np.float32)  # (G, M, 1)
        csmu = max_up.sum(axis=2, keepdims=True)  # (G, M, 1)
        scale = np.zeros_like(csmu, dtype=np.float32)  # (G, M, 1)
        np.divide(
            np.minimum(rsm_3d, csmu), csmu, out=scale, where=(csmu > 0)
        )  # avoids div/0
        upt = (max_up * scale).astype(np.float32)  # (G, M, T)

        # ── Step 2: Update Monomers pool ─────────────────────────────────────────
        # upt.sum(axis=2): (G, M) — total uptake per grid × monomer
        upt_sum_T = upt.sum(axis=2)  # (G, M)
        MR_C = self.Monomer_ratios["C"].values.reshape(G, M).astype(np.float32)
        MR_N = self.Monomer_ratios["N"].values.reshape(G, M).astype(np.float32)
        MR_P = self.Monomer_ratios["P"].values.reshape(G, M).astype(np.float32)
        mon_depl = (
            np.stack([MR_C, MR_N, MR_P], axis=2) * upt_sum_T[:, :, None]
        )  # (G,M,3)
        self.Monomers -= mon_depl.reshape(G * M, 3)

        # ── Step 3: C uptake allocation to energy vs biomass ─────────────────────
        self.E_max = self.E_ratb + self.E_range
        Tax_E_demand = 1 - (
            (
                1
                - (
                    (self.E_max - (self.Microbes["E"] / self.Microbes["C"]))
                    / self.E_max
                )
            )
            ** (self.ferm_cond + 1)
        )

        Org_Mon_E_value = self.Monomer_energy["E_value"] ** (
            (2 / np.log10(self.Monomer_energy["E_ferm"])) ** self.ferm_cond
        )
        Mon_E_value = np.concatenate([[0, 0], Org_Mon_E_value.values]).astype(
            np.float32
        )  # (M,)
        # Tax_E_alloc: (M, G*T) → reshape to (G, M, T)
        Tea_3d = (
            np.outer(Mon_E_value, Tax_E_demand.values.astype(np.float32))
            .reshape(M, G, T)
            .transpose(1, 0, 2)
            .astype(np.float32)
        )  # (G,M,T)
        Tea_3d = np.nan_to_num(Tea_3d, nan=0.0)

        # Store Tax_E_alloc as DataFrame for any downstream code that reads it
        self.Tax_E_alloc = pd.DataFrame(
            Tea_3d.transpose(1, 0, 2).reshape(M, G * T),
            index=self.monomer_index,
            columns=self.Microbes.index,
        )

        c_upt_3d = upt * MR_C[:, :, None]  # (G, M, T) C uptake
        ce_3d = (c_upt_3d * Tea_3d).astype(np.float32)  # (G, M, T) C→energy
        cb_3d = (c_upt_3d - ce_3d).astype(np.float32)  # (G, M, T) C→biomass
        n_3d = upt * MR_N[:, :, None]  # (G, M, T) N uptake
        p_3d = upt * MR_P[:, :, None]  # (G, M, T) P uptake

        # ── Step 4: Energy pathways ───────────────────────────────────────────────
        ME_resp = self.Mon_energy_resp.values.reshape(G, M).astype(np.float32)  # (G,M)
        ME_ferm = self.Mon_energy_ferm.values.reshape(G, M).astype(np.float32)
        ce_ferm = ce_3d * self.ferm_cond  # (G, M, T)
        ce_resp = ce_3d * (1 - self.ferm_cond)  # (G, M, T)
        e_3d = (ce_resp * ME_resp[:, :, None] + ce_ferm * ME_ferm[:, :, None]).astype(
            np.float32
        )  # (G, M, T) energy
        co2_3d = ce_resp + ce_ferm * self.Pox_Ferm  # (G, M, T) CO2
        ferm_3d = ce_ferm * (1 - self.Pox_Ferm)  # (G, M, T) ferm prods

        # ── Step 5: Reduce (G, M, T) → (G, T) then ravel to (G*T,) Series ───────
        mic_idx = self.Microbes.index

        TUC = c_upt_3d.sum(axis=1)  # (G, T) total C uptake
        TUCb = cb_3d.sum(axis=1)  # (G, T) biomass C
        TUE = e_3d.sum(axis=1)  # (G, T) energy
        TUN = n_3d.sum(axis=1)  # (G, T) N
        TUP = p_3d.sum(axis=1)  # (G, T) P
        TU = TUCb + TUN + TUP  # (G, T) total biomass

        self.Taxon_resp_df = pd.DataFrame(
            co2_3d.sum(axis=1),  # (G, T) CO2
            columns=mic_idx[:T],
        )

        # ── Step 6: Fermentation products → Monomers ─────────────────────────────
        FermProd_C = ferm_3d.sum(axis=(1, 2)).astype(np.float32)  # (G,)
        # Add per-cell fermentation C to each grid cell's FermProd row.
        # Monomers has G*M rows; FermProd is the 3rd monomer (index 2).
        self.Monomers.loc["FermProd", "C"] += FermProd_C

        # ── Reconstruct DataFrames needed by the rest of uptake() ─────────────────
        # TUC_df etc. are (G, T) DataFrames indexed by grid cell with taxon columns.
        # We build them directly from the numpy arrays to avoid extra copies.
        tax_cols = mic_idx[:T]
        TUC_df = pd.DataFrame(TUC, columns=tax_cols)
        TUCb_df = pd.DataFrame(TUCb, columns=tax_cols)
        TUE_df = pd.DataFrame(TUE, columns=tax_cols)
        TUN_df = pd.DataFrame(TUN, columns=tax_cols)
        TUP_df = pd.DataFrame(TUP, columns=tax_cols)
        TU_df = pd.DataFrame(TU, columns=tax_cols)

        # Update these global variables
        # self.Taxon_Uptake_C = TUC_df.stack().values.astype(np.float32)  # spatial total C uptake
        # self.Taxon_Uptake_BC = TUCb_df.stack().values  # spatial biomass C uptake: array
        # self.Taxon_Uptake_N = TUN_df.stack().values  # spatial N uptake: array
        # self.Taxon_Uptake_P = TUP_df.stack().values  # spatial P uptake: array
        # self.Taxon_Uptake_E = TUE_df.stack().values.astype(np.float32)  # spatial Energy uptake: array
        self.Taxon_Uptake_C = pd.Series(
            TUC_df.values.ravel(order="C").astype(np.float32), index=self.Microbes.index
        )
        self.Taxon_Uptake_BC = pd.Series(
            TUCb_df.values.ravel(order="C").astype(np.float32),
            index=self.Microbes.index,
        )
        self.Taxon_Uptake_N = pd.Series(
            TUN_df.values.ravel(order="C").astype(np.float32), index=self.Microbes.index
        )
        self.Taxon_Uptake_P = pd.Series(
            TUP_df.values.ravel(order="C").astype(np.float32), index=self.Microbes.index
        )
        self.Taxon_Uptake_E = pd.Series(
            TUE_df.values.ravel(order="C").astype(np.float32), index=self.Microbes.index
        )
        self.Taxon_Uptake_Total = TU_df.values.ravel(order="C").astype(
            np.float32
        )  # total biomass (C, N, & P)

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
        CAE_taxon_v[pos_uptake_index] = (
            Taxon_Uptake_BC_arr[pos_uptake_index] / Taxon_Uptake_C_arr[pos_uptake_index]
        )
        self.CAE_taxon = CAE_taxon_v

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
            TUE_df.values.ravel(order="C").astype(np.float32), index=self.Microbes.index
        )
        self.Taxon_E_Prod = self.CAE_taxon.copy()
        self.Taxon_E_Prod[:] = self.Taxon_Uptake_E
        self.Taxon_Total_Uptake = pd.Series(
            TU_df.values.ravel(order="C").astype(np.float32), index=self.Microbes.index
        )

    def metabolism(self, day):
        """
        Explicitly calculate intra-cellular production of metabolites.

        Handles both constitutive (standing biomass) and inducible (immediate monomers uptake) pathways following:
          1. constitutive enzyme and osmolyte production
          2. inducible enzyme and osmolyte production
          3. emergent CUE & Respiration
          4. update both Enzymes (with production & loss) and Substrates (with dead enzymes)

        NOTE: rewritten to numpy-native (G*T,)/(G*n_enz,)-flat array arithmetic, mirroring the
        optimization already applied to uptake(). Unlike uptake(), nothing here needs to reduce
        *across* taxa or grid cells within a row, so no 3-D (G,T,...) reshape is needed -- every
        operation here is row-wise over the same (gridsize*n_taxa)-length axis Microbes.index
        already uses. The original pandas version paid index-alignment/boolean-mask overhead on
        every one of ~40 chained Series/DataFrame ops per call; this version does the identical
        math on raw .values arrays and only re-wraps results in pandas objects where something
        downstream actually needs one (self.Microbes, self.Enzymes, self.Substrates, self.Taxon_HSP,
        self.Taxon_Osmolyte -- everything else metabolism() sets is read downstream only via
        output.py's _taxa_sum()/_enz_sum(), which already fall back to np.asarray() for plain
        arrays, so those are left as plain float32 numpy arrays).

        BUGFIX (2026-07-22): the N/E "most limiting resource" ratio (ARROEHC / ARROEI)
        used to be left undefined (NaN) for any row where its own E-cost was exactly 0
        while the N-cost was positive -- mirroring the original pandas code, which built
        this ratio only over the boolean subset where the E-cost was positive, then
        aligned it against the (larger) N-cost-positive subset, silently introducing NaN
        wherever a row was in the second subset but not the first. This is a genuine bug
        in the original model, not an intentional design choice: neither pandas' nor
        numpy's `x < 0` clipping catches NaN (comparisons against NaN are always False),
        so a NaN in Microbes['E'] introduced this way propagates forward permanently once
        introduced. Confirmed against a real run: energy production was legitimately
        nonzero on day 0, then exactly 0/NaN every day after, because some cell first hit
        this maskN-without-maskE condition during day 0's metabolism() call. Fixed by
        falling back to the N-only ratio (ARRNOEHC / ARRNOEI, always defined wherever
        N-cost is positive) instead of NaN when there's no E-cost to weigh against --
        see the inline comment where ARROEHC/ARROEI are computed below. This does not
        change behavior at all in the (typical) case where both costs are positive.
        """

        G, T = self.gridsize, self.n_taxa
        idx = (
            self.Microbes.index
        )  # unchanged; row order for every (G*T,)-length quantity below

        # index of dead enzyme in Substrates
        is_deadEnz = self.Substrates.index == "DeadEnz"

        # Scalar of water potential impact: call the function microbe_osmo_psi()
        f_psi = microbe_osmo_psi(self.psi[day], self.alpha, self.wp_fc, self.wp_th)

        # ── Cache raw .values arrays once (avoid repeated pandas attribute access/overhead) ──
        Microbes_C = self.Microbes["C"].to_numpy(dtype=np.float32, copy=True)
        Microbes_N = self.Microbes["N"].to_numpy(dtype=np.float32, copy=True)
        Microbes_P = self.Microbes["P"].to_numpy(dtype=np.float32, copy=True)
        Microbes_E = self.Microbes["E"].to_numpy(dtype=np.float32, copy=True)
        Microbes_CNP = (
            Microbes_C + Microbes_N + Microbes_P
        )  # == Microbes[['C','N','P']].sum(axis=1)

        Taxon_Total_Uptake_v = np.asarray(self.Taxon_Total_Uptake, dtype=np.float32)
        Taxon_Uptake_BC_v = np.asarray(self.Taxon_Uptake_BC, dtype=np.float32)
        Taxon_Uptake_N_v = np.asarray(self.Taxon_Uptake_N, dtype=np.float32)
        Taxon_Uptake_P_v = np.asarray(self.Taxon_Uptake_P, dtype=np.float32)
        Taxon_Uptake_E_v = np.asarray(self.Taxon_Uptake_E, dtype=np.float32)
        Taxon_Uptake_C_v = np.asarray(self.Taxon_Uptake_C, dtype=np.float32)

        # taxon-specific energy cost of active monomer transport (mass of assimilated monomers * active uptake cost)
        Tax_active_transport_E_cost = Taxon_Total_Uptake_v * np.float32(
            self.active_uptake_E_cost
        )

        # calculate the proportion of microbial biomass that is HSPs
        hsp_cols = [c for c in self.Taxon_HSP.columns if c != "Total_biomass"]
        HSP_v = self.Taxon_HSP[hsp_cols].to_numpy(
            dtype=np.float32, copy=True
        )  # (G*T, n_HSP)
        HSP_v *= np.float32(1) - np.float32(self.HSP_Turnover)
        self.Taxon_HSP[hsp_cols] = HSP_v

        Total_biomass_v = self.Taxon_HSP["Total_biomass"].to_numpy(
            dtype=np.float32, copy=True
        )
        hsp_mass = HSP_v.sum(axis=1)  # total HSP mass per microbe
        HSP_prop = hsp_mass.copy()
        pos_tb = Total_biomass_v > 0
        HSP_prop[pos_tb] = hsp_mass[pos_tb] / Total_biomass_v[pos_tb]

        # Calculate temperature sensitivity of maintenance energy cost from HSP concentration
        maint_E_cost_Ea = (
            (HSP_prop / (HSP_prop + self.maint_E_HSP_Ea_n)) * self.maint_E_HSP_Ea_slope
        ) + self.maint_E_HSP_Ea_int
        tax_maint_E_cost = self.maint_E_rate_const * Arrhenius(
            maint_E_cost_Ea, self.gas_const, self.temp[day], self.maint_E_rate_Tref
        )
        Tax_cell_maint_E_cost = tax_maint_E_cost * Microbes_CNP

        # Transporter production
        Uptake_Prod_Cost_sum = self.Uptake_Prod_Cost.to_numpy(
            dtype=np.float32, copy=True
        ).sum(axis=1)
        Tax_transport_prod_E_cost = (
            Taxon_Total_Uptake_v * Uptake_Prod_Cost_sum
        ) * np.float32(self.UptakeProd_E_Cost)

        # Pool of energy available for biomass production after uptake/transport/maintenance costs
        MicBio_avail_E = (
            Microbes_E
            + Taxon_Uptake_E_v
            - Tax_active_transport_E_cost
            - Tax_transport_prod_E_cost
            - Tax_cell_maint_E_cost
        )
        MicBio_avail_E = np.where(MicBio_avail_E < 0, np.float32(0), MicBio_avail_E)

        if (
            getattr(self, "debug_metabolism_days", 0)
            and day < self.debug_metabolism_days
        ):

            def _rep(name, arr):
                a = np.asarray(arr, dtype=np.float64)
                n_nan = np.isnan(a).sum()
                n_inf = np.isinf(a).sum()
                print(
                    f"    [metab day {day}] {name}: n_nan={n_nan} n_inf={n_inf} "
                    f"shape={a.shape} sum={np.nansum(a):.4g}"
                )

            def _rep_df(name, df):
                a = np.asarray(df, dtype=np.float64)
                n_nan = np.isnan(a).sum()
                cols = getattr(df, "columns", None)
                print(f"    [metab day {day}] {name}: n_nan={n_nan} shape={a.shape}")
                if n_nan and cols is not None:
                    col_nan = np.isnan(a).any(axis=0)
                    bad_cols = [str(c) for c, bad in zip(cols, col_nan) if bad]
                    print(f"        -> NaN-containing columns: {bad_cols}")
                    row_nan = np.isnan(a).any(axis=1)
                    print(
                        f"        -> n_rows_with_any_nan={row_nan.sum()} "
                        f"(of {a.shape[0]})"
                    )

            print(f"  -- metabolism() checkpoint A (day {day}) --")
            _rep_df("self.Consti_Enzyme_C (raw input)", self.Consti_Enzyme_C)
            _rep_df("self.Induci_Enzyme_C (raw input)", self.Induci_Enzyme_C)
            _rep_df("self.Consti_Osmo_C (raw input)", self.Consti_Osmo_C)
            _rep_df("self.Induci_Osmo_C (raw input)", self.Induci_Osmo_C)
            _rep_df("self.HSP_consti_pre_exp (raw input)", self.HSP_consti_pre_exp)
            _rep_df("self.HSP_induci_pre_exp (raw input)", self.HSP_induci_pre_exp)
            _rep_df("self.Uptake_Prod_Cost (raw input)", self.Uptake_Prod_Cost)
            _rep("HSP_prop", HSP_prop)
            _rep("maint_E_cost_Ea", maint_E_cost_Ea)
            _rep("tax_maint_E_cost", tax_maint_E_cost)
            _rep("Tax_cell_maint_E_cost", Tax_cell_maint_E_cost)
            _rep("Uptake_Prod_Cost_sum", Uptake_Prod_Cost_sum)
            _rep("Tax_transport_prod_E_cost", Tax_transport_prod_E_cost)
            _rep("MicBio_avail_E (initial)", MicBio_avail_E)

        # ---------------------------------------------------------------------#
        # ......................constitutive processes.........................#
        # ---------------------------------------------------------------------#

        # Reindex the small (n_enzymes,)-length per-enzyme cost columns to match column order,
        # once, defensively (cheap -- avoids relying on positional ordering happening to match).
        consti_enz_cols = self.Consti_Enzyme_C.columns
        Enz_N_cost_v = (
            self.Enz_Attrib["N_cost"]
            .reindex(consti_enz_cols)
            .to_numpy(dtype=np.float32, copy=True)
        )
        Enz_P_cost_v = (
            self.Enz_Attrib["P_cost"]
            .reindex(consti_enz_cols)
            .to_numpy(dtype=np.float32, copy=True)
        )
        Enz_E_cost_v = (
            self.Enz_Attrib["E_cost"]
            .reindex(consti_enz_cols)
            .to_numpy(dtype=np.float32, copy=True)
        )

        # HSP before adjustment
        HSP_consti_prod = self.HSP_consti_pre_exp.to_numpy(
            dtype=np.float32, copy=True
        ) * Arrhenius(
            self.HSP_consti_Ea, self.gas_const, self.temp[day], self.maint_E_rate_Tref
        )
        tax_HSP_consti_prod = HSP_consti_prod * Microbes_CNP[:, None]
        Taxon_HSP_Consti_Cost_N = (tax_HSP_consti_prod * self.HSP_N_cost).sum(axis=1)
        Taxon_HSP_Consti_Cost_E = (tax_HSP_consti_prod * self.HSP_E_cost).sum(axis=1)

        # Osmolyte before adjustment
        Taxon_Osmo_Consti = (
            self.Consti_Osmo_C.to_numpy(dtype=np.float32, copy=True)
            * Microbes_C[:, None]
            * f_psi
        )
        Taxon_Osmo_Consti_Cost_N = (Taxon_Osmo_Consti * self.Osmo_N_cost).sum(axis=1)
        Taxon_Osmo_Consti_Cost_E = (Taxon_Osmo_Consti * self.Osmo_E_cost).sum(axis=1)

        # Enzyme before adjustment
        Taxon_Enzyme_Consti = (
            self.Consti_Enzyme_C.to_numpy(dtype=np.float32, copy=True)
            * Microbes_C[:, None]
        )
        Taxon_Enzyme_Consti_Cost_N = (Taxon_Enzyme_Consti * Enz_N_cost_v).sum(axis=1)
        Taxon_Enzyme_Consti_Cost_E = (Taxon_Enzyme_Consti * Enz_E_cost_v).sum(axis=1)

        # N limitation across all three constitutive products
        OEHCCN = (
            Taxon_HSP_Consti_Cost_N
            + Taxon_Osmo_Consti_Cost_N
            + Taxon_Enzyme_Consti_Cost_N
        )
        maskN = OEHCCN > 0
        MNAOEHC = np.minimum(OEHCCN, Microbes_N)
        ARRNOEHC = np.where(maskN, MNAOEHC / np.where(maskN, OEHCCN, 1), np.nan)

        # E limitation across all three constitutive products
        OEHCCE = (
            Taxon_HSP_Consti_Cost_E
            + Taxon_Osmo_Consti_Cost_E
            + Taxon_Enzyme_Consti_Cost_E
        )
        maskE = OEHCCE > 0
        MEAOEHC = np.minimum(OEHCCE, MicBio_avail_E)
        ARREOEHC = np.where(maskE, MEAOEHC / np.where(maskE, OEHCCE, 1), np.nan)

        # Most limiting factor: use whichever of the N-ratio / E-ratio is more
        # restrictive when BOTH costs are positive; fall back to the N-ratio
        # alone when there's no E-cost to weigh against (maskE False).
        #
        # BUGFIX (was NaN-producing in both the original pandas code and an
        # earlier version of this rewrite): the original computed ARROEHC as
        # `ARREOEHC.copy()` where ARREOEHC is a Series indexed ONLY over the
        # maskE==True subset (built via `OEHCCE[OEHCCE > 0]`-style boolean
        # selection, which drops rows rather than keeping them as NaN). It
        # then did `tax_HSP_consti_prod[OEHCCN > 0].mul(ARROEHC, axis=0)` --
        # a pandas alignment between the maskN subset and the maskE subset.
        # For any row where maskN is True but maskE is False (an N-cost with
        # no corresponding E-cost that step), ARROEHC has no entry for that
        # row, so the aligned multiply silently produces NaN, which then
        # gets written into tax_HSP_consti_prod for that row. Since neither
        # pandas' nor numpy's `x < 0` clipping catches NaN (NaN compares
        # False against everything, so `Microbes['E'][Microbes['E']<0] = 0`
        # never clears it), that NaN persists in self.Microbes['E'] from
        # that day forward -- confirmed against a real run (E production
        # legitimately nonzero on day 0, then NaN forever after day 0's
        # metabolism() call once some cell first hit this maskN-without-
        # maskE condition). Falling back to ARRNOEHC (the N-only ratio,
        # always defined wherever maskN is True) instead of NaN fixes this
        # without changing behavior in the maskE==True case at all.
        ARROEHC = np.where(maskE, ARREOEHC, ARRNOEHC)
        swap = maskE & maskN & (ARREOEHC > ARRNOEHC)
        ARROEHC = np.where(swap, ARRNOEHC, ARROEHC)

        if (
            getattr(self, "debug_metabolism_days", 0)
            and day < self.debug_metabolism_days
        ):
            _rep("Enz_N_cost_v (consti, reindexed)", Enz_N_cost_v)
            _rep("Enz_E_cost_v (consti, reindexed)", Enz_E_cost_v)
            _rep("OEHCCN", OEHCCN)
            _rep("OEHCCE", OEHCCE)
            _rep("ARRNOEHC", ARRNOEHC)
            _rep("ARREOEHC", ARREOEHC)
            _rep("ARROEHC (final)", ARROEHC)
            print(
                f"    [metab day {day}] maskN sum={maskN.sum()} maskE sum={maskE.sum()} "
                f"maskN&~maskE={int((maskN & ~maskE).sum())} ~maskN&maskE={int((~maskN & maskE).sum())}"
            )

        # HSP adjusted
        tax_HSP_consti_prod = np.where(
            maskN[:, None], tax_HSP_consti_prod * ARROEHC[:, None], tax_HSP_consti_prod
        )
        Taxon_HSP_Consti_Cost_N = (tax_HSP_consti_prod * self.HSP_N_cost).sum(axis=1)
        Taxon_HSP_Consti_Cost_C = tax_HSP_consti_prod.sum(axis=1)
        Taxon_HSP_Consti_Cost_E = (tax_HSP_consti_prod * self.HSP_E_cost).sum(axis=1)

        # Osmolyte adjusted
        Taxon_Osmo_Consti = np.where(
            maskN[:, None], Taxon_Osmo_Consti * ARROEHC[:, None], Taxon_Osmo_Consti
        )
        Taxon_Osmo_Consti_Cost_N = (Taxon_Osmo_Consti * self.Osmo_N_cost).sum(axis=1)
        Taxon_Osmo_Consti_Cost_C = Taxon_Osmo_Consti.sum(axis=1)
        Taxon_Osmo_Consti_Cost_E = Taxon_Osmo_Consti_Cost_C * self.Osmo_E_cost

        # Enzyme adjusted
        Taxon_Enzyme_Consti = np.where(
            maskN[:, None], Taxon_Enzyme_Consti * ARROEHC[:, None], Taxon_Enzyme_Consti
        )
        Taxon_Enzyme_Consti_Cost_N = (Taxon_Enzyme_Consti * Enz_N_cost_v).sum(axis=1)
        Taxon_Enzyme_Consti_Cost_P = (Taxon_Enzyme_Consti * Enz_P_cost_v).sum(axis=1)
        Taxon_Enzyme_Consti_Cost_C = Taxon_Enzyme_Consti.sum(axis=1)
        Taxon_Enzyme_Consti_Cost_E = (Taxon_Enzyme_Consti * Enz_E_cost_v).sum(axis=1)

        # Deduct all constitutive E costs from available pool
        MicBio_avail_E = MicBio_avail_E - (
            Taxon_HSP_Consti_Cost_E
            + Taxon_Osmo_Consti_Cost_E
            + Taxon_Enzyme_Consti_Cost_E
        )
        MicBio_avail_E = np.where(MicBio_avail_E < 0, np.float32(0), MicBio_avail_E)

        if (
            getattr(self, "debug_metabolism_days", 0)
            and day < self.debug_metabolism_days
        ):
            _rep("Taxon_HSP_Consti_Cost_E (adj)", Taxon_HSP_Consti_Cost_E)
            _rep("Taxon_Osmo_Consti_Cost_E (adj)", Taxon_Osmo_Consti_Cost_E)
            _rep("Taxon_Enzyme_Consti_Cost_E (adj)", Taxon_Enzyme_Consti_Cost_E)
            _rep("MicBio_avail_E (post-consti)", MicBio_avail_E)

        # --- Inducible production (HSP, Osmolyte, Enzyme) before resource adjustment ---

        induci_enz_cols = self.Induci_Enzyme_C.columns
        Enz_N_cost_i = (
            self.Enz_Attrib["N_cost"]
            .reindex(induci_enz_cols)
            .to_numpy(dtype=np.float32, copy=True)
        )
        Enz_P_cost_i = (
            self.Enz_Attrib["P_cost"]
            .reindex(induci_enz_cols)
            .to_numpy(dtype=np.float32, copy=True)
        )
        Enz_E_cost_i = (
            self.Enz_Attrib["E_cost"]
            .reindex(induci_enz_cols)
            .to_numpy(dtype=np.float32, copy=True)
        )

        tax_HSP_induci_prod_rate = self.HSP_induci_pre_exp.to_numpy(
            dtype=np.float32, copy=True
        ) * Arrhenius(
            self.HSP_induci_Ea, self.gas_const, self.temp[day], self.maint_E_rate_Tref
        )
        tax_HSP_induci_prod = tax_HSP_induci_prod_rate * Taxon_Total_Uptake_v[:, None]
        Taxon_HSP_Induci_Cost_N = (tax_HSP_induci_prod * self.HSP_N_cost).sum(axis=1)
        Taxon_HSP_Induci_Cost_E = (tax_HSP_induci_prod * self.HSP_E_cost).sum(axis=1)

        Taxon_Osmo_Induci = (
            self.Induci_Osmo_C.to_numpy(dtype=np.float32, copy=True)
            * Taxon_Uptake_BC_v[:, None]
            * f_psi
        )
        Taxon_Osmo_Induci_Cost_N = (Taxon_Osmo_Induci * self.Osmo_N_cost).sum(axis=1)
        Taxon_Osmo_Induci_Cost_E = (Taxon_Osmo_Induci * self.Osmo_E_cost).sum(axis=1)

        Taxon_Enzyme_Induci = (
            self.Induci_Enzyme_C.to_numpy(dtype=np.float32, copy=True)
            * Taxon_Uptake_BC_v[:, None]
        )
        Taxon_Enzyme_Induci_Cost_N = (Taxon_Enzyme_Induci * Enz_N_cost_i).sum(axis=1)
        Taxon_Enzyme_Induci_Cost_E = (Taxon_Enzyme_Induci * Enz_E_cost_i).sum(axis=1)

        # --- N limitation: HSP + Osmolyte + Enzyme inducible ---
        OEIHCN = (
            Taxon_HSP_Induci_Cost_N
            + Taxon_Osmo_Induci_Cost_N
            + Taxon_Enzyme_Induci_Cost_N
        )
        OEIHAN = Taxon_Uptake_N_v
        maskNi = OEIHCN > 0
        MNAOEIH = np.minimum(OEIHCN, OEIHAN)
        ARRNOEI = np.where(maskNi, MNAOEIH / np.where(maskNi, OEIHCN, 1), np.nan)

        # --- E limitation: HSP + Osmolyte + Enzyme inducible ---
        OEIHCE = (
            Taxon_HSP_Induci_Cost_E
            + Taxon_Osmo_Induci_Cost_E
            + Taxon_Enzyme_Induci_Cost_E
        )
        maskEi = OEIHCE > 0
        MEAOEIH = np.minimum(OEIHCE, MicBio_avail_E)
        ARREOEI = np.where(maskEi, MEAOEIH / np.where(maskEi, OEIHCE, 1), np.nan)

        # Most limiting factor across N and E -- same bugfix as the constitutive
        # block above: fall back to the N-only ratio (ARRNOEI) when there's no
        # E-cost to weigh against, instead of the original's NaN.
        ARROEI = np.where(maskEi, ARREOEI, ARRNOEI)
        swap_i = maskEi & maskNi & (ARREOEI > ARRNOEI)
        ARROEI = np.where(swap_i, ARRNOEI, ARROEI)

        # --- Apply adjustment to all three inducible products ---

        tax_HSP_induci_prod = np.where(
            maskNi[:, None], tax_HSP_induci_prod * ARROEI[:, None], tax_HSP_induci_prod
        )
        Taxon_HSP_Induci_Cost_N = (tax_HSP_induci_prod * self.HSP_N_cost).sum(axis=1)
        Taxon_HSP_Induci_Cost_C = tax_HSP_induci_prod.sum(axis=1)
        Taxon_HSP_Induci_Cost_E = (tax_HSP_induci_prod * self.HSP_E_cost).sum(axis=1)

        Taxon_Osmo_Induci = np.where(
            maskNi[:, None], Taxon_Osmo_Induci * ARROEI[:, None], Taxon_Osmo_Induci
        )
        Taxon_Osmo_Induci_Cost_N = (Taxon_Osmo_Induci * self.Osmo_N_cost).sum(axis=1)
        Taxon_Osmo_Induci_Cost_C = Taxon_Osmo_Induci.sum(axis=1)
        Taxon_Osmo_Induci_Cost_E = (Taxon_Osmo_Induci * self.Osmo_E_cost).sum(axis=1)

        Taxon_Enzyme_Induci = np.where(
            maskNi[:, None], Taxon_Enzyme_Induci * ARROEI[:, None], Taxon_Enzyme_Induci
        )
        Taxon_Enzyme_Induci_Cost_N = (Taxon_Enzyme_Induci * Enz_N_cost_i).sum(axis=1)
        Taxon_Enzyme_Induci_Cost_P = (Taxon_Enzyme_Induci * Enz_P_cost_i).sum(axis=1)
        Taxon_Enzyme_Induci_Cost_C = Taxon_Enzyme_Induci.sum(axis=1)
        Taxon_Enzyme_Induci_Cost_E = (Taxon_Enzyme_Induci * Enz_E_cost_i).sum(axis=1)

        # Derive C, N, & P deposited as biomass from Uptake; ensure no negative values
        Microbe_C_Gain = Taxon_Uptake_BC_v - Taxon_Enzyme_Induci_Cost_C
        Microbe_N_Gain = Taxon_Uptake_N_v - Taxon_Enzyme_Induci_Cost_N
        Microbe_P_Gain = Taxon_Uptake_P_v - Taxon_Enzyme_Induci_Cost_P

        # Derive total E costs and gains
        Microbe_E_Cost = (
            Tax_active_transport_E_cost
            + Tax_transport_prod_E_cost
            + Tax_cell_maint_E_cost
            + Taxon_HSP_Consti_Cost_E
            + Taxon_HSP_Induci_Cost_E
            + Taxon_Enzyme_Induci_Cost_E
            + Taxon_Osmo_Induci_Cost_E
            + Taxon_Enzyme_Consti_Cost_E
            + Taxon_Osmo_Consti_Cost_E
        )
        Microbe_E_Gain = Taxon_Uptake_E_v - Microbe_E_Cost

        if (
            getattr(self, "debug_metabolism_days", 0)
            and day < self.debug_metabolism_days
        ):
            _rep("Enz_N_cost_i (induci, reindexed)", Enz_N_cost_i)
            _rep("Enz_E_cost_i (induci, reindexed)", Enz_E_cost_i)
            _rep("OEIHCN", OEIHCN)
            _rep("OEIHCE", OEIHCE)
            _rep("ARROEI (final)", ARROEI)
            print(
                f"    [metab day {day}] maskNi sum={maskNi.sum()} maskEi sum={maskEi.sum()} "
                f"maskNi&~maskEi={int((maskNi & ~maskEi).sum())} ~maskNi&maskEi={int((~maskNi & maskEi).sum())}"
            )
            _rep("Microbe_E_Cost", Microbe_E_Cost)
            _rep("Microbe_E_Gain", Microbe_E_Gain)
            _rep("Microbes_E (before integration)", Microbes_E)

        # ---------------------------------------------------------------------#
        # ...............................Integration...........................#
        # ---------------------------------------------------------------------#
        Microbes_C = Microbes_C + Microbe_C_Gain - Taxon_Enzyme_Consti_Cost_C
        Microbes_N = Microbes_N + Microbe_N_Gain - Taxon_Enzyme_Consti_Cost_N
        Microbes_P = Microbes_P + Microbe_P_Gain - Taxon_Enzyme_Consti_Cost_P
        Microbes_E = Microbes_E + Microbe_E_Gain
        Microbes_C = np.where(Microbes_C < 0, np.float32(0), Microbes_C)
        Microbes_N = np.where(Microbes_N < 0, np.float32(0), Microbes_N)
        Microbes_P = np.where(Microbes_P < 0, np.float32(0), Microbes_P)
        Microbes_E = np.where(Microbes_E < 0, np.float32(0), Microbes_E)
        self.Microbes["C"] = Microbes_C
        self.Microbes["N"] = Microbes_N
        self.Microbes["P"] = Microbes_P
        self.Microbes["E"] = Microbes_E

        if (
            getattr(self, "debug_metabolism_days", 0)
            and day < self.debug_metabolism_days
        ):
            _rep("Microbes_E (final, post-clip)", Microbes_E)
            n_alive = int(((Microbes_C + Microbes_N + Microbes_P) > 0).sum())
            print(f"    [metab day {day}] n_alive (CNP>0) = {n_alive}")

        # Taxon-specific emergent CUE
        CUE_taxon = np.zeros(G * T, dtype=np.float32)
        pos_uptake_index = Taxon_Uptake_C_v > 0
        Mic_BC_Use = Microbe_C_Gain - Taxon_Enzyme_Consti_Cost_C
        CUE_taxon[pos_uptake_index] = (
            Mic_BC_Use[pos_uptake_index] / Taxon_Uptake_C_v[pos_uptake_index]
        )

        # System-level emergent CUE
        Taxon_Use_C_grid = Taxon_Uptake_C_v.sum()
        if Taxon_Use_C_grid == 0:
            CUE_system = np.float32(0)
        else:
            CUE_system = Mic_BC_Use.sum() / Taxon_Use_C_grid

        # Heat stress protein production
        HSP_v = HSP_v + tax_HSP_consti_prod + tax_HSP_induci_prod
        self.Taxon_HSP[hsp_cols] = HSP_v
        self.Taxon_HSP["Total_biomass"] = Microbes_C + Microbes_N + Microbes_P

        # Osmolyte pool
        Osmo_v = self.Taxon_Osmolyte.to_numpy(dtype=np.float32, copy=True)
        Osmo_v = Osmo_v + Taxon_Osmo_Consti + Taxon_Osmo_Induci
        Osmo_v = Osmo_v * (np.float32(1) - np.float32(self.OSmo_Loss_Rate))
        self.Taxon_Osmolyte[:] = Osmo_v

        Taxon_Osmo_Cost_C = (Taxon_Osmo_Consti + Taxon_Osmo_Induci).sum(axis=1)
        Taxon_Enzyme_Cost_C = Taxon_Enzyme_Consti_Cost_C + Taxon_Enzyme_Induci_Cost_C
        Taxon_HSP_Cost_C = Taxon_HSP_Consti_Cost_C + Taxon_HSP_Induci_Cost_C
        Taxon_Transport_Cost_C = Tax_transport_prod_E_cost / np.float32(
            self.UptakeProd_E_Cost
        )
        Taxon_Transport_E_Cost = Tax_active_transport_E_cost + Tax_transport_prod_E_cost
        Taxon_Maint_E_Cost = Tax_cell_maint_E_cost
        Taxon_Enzyme_E_Cost = Taxon_Enzyme_Consti_Cost_E + Taxon_Enzyme_Induci_Cost_E
        Taxon_Osmo_E_Cost_total = Taxon_Osmo_Consti_Cost_E + Taxon_Osmo_Induci_Cost_E
        Taxon_HSP_E_Cost = Taxon_HSP_Consti_Cost_E + Taxon_HSP_Induci_Cost_E

        # ── Publish taxon-level results (plain float32 arrays; output.py's _taxa_sum()/
        #    _enz_sum() already accept either pandas or plain-array input) ──────────────
        self.Tax_active_transport_E_cost = Tax_active_transport_E_cost
        self.Tax_cell_maint_E_cost = Tax_cell_maint_E_cost
        self.Taxon_HSP_Consti_Cost_E = Taxon_HSP_Consti_Cost_E
        self.Taxon_HSP_Induci_Cost_E = Taxon_HSP_Induci_Cost_E
        self.Microbe_E_Cost = Microbe_E_Cost
        self.Microbe_E_Gain = Microbe_E_Gain
        self.CUE_taxon = CUE_taxon
        self.CUE_system = CUE_system
        self.Taxon_Osmo_Cost_C = Taxon_Osmo_Cost_C
        self.Taxon_Enzyme_Cost_C = Taxon_Enzyme_Cost_C
        self.Taxon_HSP_Cost_C = Taxon_HSP_Cost_C
        self.Taxon_Transport_Cost_C = Taxon_Transport_Cost_C
        self.Taxon_Transport_E_Cost = Taxon_Transport_E_Cost
        self.Taxon_Maint_E_Cost = Taxon_Maint_E_Cost
        self.Taxon_Enzyme_E_Cost = Taxon_Enzyme_E_Cost
        self.Taxon_Osmo_E_Cost_total = Taxon_Osmo_E_Cost_total
        self.Taxon_HSP_E_Cost = Taxon_HSP_E_Cost

        # Derive Enzyme production: gene-specific production of each taxon (G*T, n_enzymes),
        # reduced across taxa within each grid cell -> (G, n_enzymes) -> flattened to match
        # self.Enzymes' (G*n_enzymes,) index order (grid outer, enzyme inner -- same convention
        # _enz_sum()/_taxa_sum() rely on elsewhere). No MultiIndex is built (nothing downstream
        # reads self.Taxon_Enzyme_Production's index -- only output.py's _taxa_sum(), which
        # works off .values regardless of index).
        Taxon_Enzyme_Production_v = (
            Taxon_Enzyme_Consti + Taxon_Enzyme_Induci
        )  # (G*T, n_enzymes)
        self.Taxon_Enzyme_Production = pd.DataFrame(
            Taxon_Enzyme_Production_v, index=idx, columns=consti_enz_cols
        )
        Enzyme_Production_v = (
            Taxon_Enzyme_Production_v.reshape(G, T, self.n_enzymes)
            .sum(axis=1)
            .reshape(-1)
        )

        # Derive Enzyme turnover
        Enzymes_v = self.Enzymes.to_numpy(dtype=np.float32, copy=True)
        Enzyme_Loss_v = Enzymes_v * np.float32(self.Enzyme_Loss_Rate)

        # Update Enzyme pools by adding enzymes produced and subtracting the 'dead' enzymes
        self.Enzymes[:] = Enzymes_v + Enzyme_Production_v - Enzyme_Loss_v
        self.Enzymes_gain = Enzyme_Production_v - Enzyme_Loss_v
        self.Enzyme_Production = Enzyme_Production_v
        self.Enzyme_Loss = Enzyme_Loss_v

        # Update Substrates pools with dead enzymes: reshape (G*n_enzymes,) -> (G, n_enzymes),
        # sum across enzymes per grid cell to get total dead-enzyme C/N/P mass added to each
        # cell's 'DeadEnz' substrate row.
        n_enz = self.n_enzymes
        EnzLoss_2d = Enzyme_Loss_v.reshape(G, n_enz)
        # Reindex defensively to self._enz_index (the per-cell enzyme order Enzymes/Enzyme_Loss's
        # flat (G*n_enz,) layout uses), rather than assuming Enz_Attrib's row order already matches.
        N_cost_enz = (
            self.Enz_Attrib["N_cost"]
            .reindex(self._enz_index)
            .to_numpy(dtype=np.float32, copy=True)
        )
        P_cost_enz = (
            self.Enz_Attrib["P_cost"]
            .reindex(self._enz_index)
            .to_numpy(dtype=np.float32, copy=True)
        )
        DeadEnz_C = EnzLoss_2d.sum(axis=1)
        DeadEnz_N = (EnzLoss_2d * N_cost_enz[None, :]).sum(axis=1)
        DeadEnz_P = (EnzLoss_2d * P_cost_enz[None, :]).sum(axis=1)
        dead_cols = list(
            self.Substrates.columns[:3]
        )  # original relied on positional [C,N,P] order too
        self.Substrates.loc[is_deadEnz, dead_cols] = self.Substrates.loc[
            is_deadEnz, dead_cols
        ].to_numpy(dtype=np.float32, copy=True) + np.stack(
            [DeadEnz_C, DeadEnz_N, DeadEnz_P], axis=1
        )

    def mortality(self, day):
        """
        Calculate microbial mortality and update stoichiometry of the alive and microbial pools.

        Kill microbes that are starving deterministically and microbes that are drought intolerant stochastically
        Also update Substrates with input from dead microbes, monomers(with leaching loss), and respiration
        """

        # Constants
        Leaching = np.float32(0.1)  # Abiotic monomer loss rate
        Psi_slope_leach = np.float32(
            0.5
        )  # Mositure sensivity of abiotic monomer loss rate
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
        Mic = self.Microbes.to_numpy(dtype=np.float32, copy=True)  # (N,4): C,N,P,E
        MinRat_cnp = self.MinRatios[["C", "N", "P"]].to_numpy(
            dtype=np.float32, copy=True
        )  # (N,3)
        MinRat_full4 = self.MinRatios[["C", "N", "P", "E"]].to_numpy(
            dtype=np.float32, copy=True
        )  # (N,4)

        Death = np.zeros_like(Mic)
        kill = np.zeros(n_taxa * gridsize, dtype=bool)

        # Start of calcualtion of mortality first with THRESHOLD
        # Kill microbes deterministically based on threshold values: C_min; N_min; P_min
        starve_index = (Mic[:, 0] > 0) & (
            (Mic[:, 0] < self.C_min)
            | (Mic[:, 1] < self.N_min)
            | (Mic[:, 2] < self.P_min)
        )
        # Index the dead, put them in Death, and set them to 0 in Microbes
        Death[starve_index] = Mic[starve_index]
        Mic[starve_index] = np.float32(0)
        # Index the locations where microbial cells remain alive
        mic_index = Mic[:, 0] > 0

        # Kill microbes stochastically based on mortality prob as a function of water potential and drought tolerance
        # call the function, MMP:microbe_mortality_psi()
        r_death = np.asarray(
            MMP(
                self.psi[day],
                self.wp_fc,
                self.basal_death_prob,
                self.death_rate,
                self.tolerance,
            )
        )
        kill[mic_index] = r_death[mic_index] > np.random.uniform(
            0, 1, int(mic_index.sum())
        ).astype("float32")
        # Index the dead, put them in Death, and set them to 0 in Microbes
        Death[kill] = Mic[kill]
        Mic[kill] = np.float32(0)
        # Index locations where microbes remain alive
        mic_index = Mic[:, 0] > 0

        # Calculate the total dead mass (threshold & drought) across taxa in each grid cell
        Death_gridcell = (
            Death[:, :3].reshape(gridsize, n_taxa, 3).sum(axis=1)
        )  # (gridsize,3): C,N,P

        # Distinguish between conditions of complete death VS partial death
        # All cells die
        if not mic_index.any():
            # ...Update Substrates pool by adding dead microbial biomass
            self.Substrates.loc[is_DeadMic] += Death_gridcell
            # NOTE: self.Taxon_Resp is deliberately left untouched here, matching
            # the original -- it retains whatever value metabolism() last set.

        # Partly die and adjust stoichiometry of those remaining alive
        else:
            # Index only those taxa in Microbes that have below-minimum quotas
            idx = np.nonzero(mic_index)[0]
            Mic_alive = Mic[idx]  # (k,4)
            row_sum_cnp = Mic_alive[:, :3].sum(
                axis=1
            )  # C+N+P (E excluded), matches .drop('E',axis=1).sum(axis=1)
            with np.errstate(invalid="ignore", divide="ignore"):
                ratios_cnp = Mic_alive[:, :3] / row_sum_cnp[:, None]
            MinRat_alive = MinRat_cnp[idx]  # (k,3)
            below_min = (
                (ratios_cnp[:, 0] < MinRat_alive[:, 0])
                | (ratios_cnp[:, 1] < MinRat_alive[:, 1])
                | (ratios_cnp[:, 2] < MinRat_alive[:, 2])
            )

            # Derive the Microbes wanted
            idx2 = idx[below_min]  # positions needing rebalancing
            Mic_subset = Mic[idx2]  # (m,4)
            StartMicrobes = Mic_subset.copy()
            MinRat4 = MinRat_full4[idx2]  # (m,4): C,N,P,E mins

            # Derive new ratios and Calculate difference between actual and min ratios.
            # NOTE: this operates on all 4 columns (C,N,P,E), matching the original's
            # `Mic_subset.divide(Mic_subset.sum(axis=1), axis=0)` -- Mic_subset there is
            # the full 4-column DataFrame, so the denominator is C+N+P+E (not C+N+P as
            # in the first ratio check above), and MinRatios itself has an 'E' column
            # (the taxon's E:biomass-C ratio, ~5000 -- see microbe.py's
            # minimum_cell_quota()). That E term dominates Ratio_dif_0.sum(axis=1)/
            # Excess.sum(axis=1) below (the "excess redistribution" scale factor), so
            # dropping E here (an earlier version of this rewrite did exactly that)
            # changes the rebalancing outcome substantially, not just cosmetically --
            # confirmed via a side-by-side day-by-day biomass comparison against the
            # original on a real run, where the E-dropped version caused runaway
            # population collapse within ~10 days that doesn't occur with E included.
            row_sum_all4 = Mic_subset.sum(axis=1)
            with np.errstate(invalid="ignore", divide="ignore"):
                MicrobeRatios4 = Mic_subset / row_sum_all4[:, None]  # (m,4)
            Ratio_dif4 = MicrobeRatios4 - MinRat4  # (m,4)
            # Create an array recording the ratio differences < 0
            Ratio_dif_0_4 = np.where(Ratio_dif4 > 0, np.float32(0), Ratio_dif4)
            # Create an array recording the ratio differences > 0
            Excess4 = np.where(Ratio_dif4 < 0, np.float32(0), Ratio_dif4)

            # Determine the limiting nutrient that will be conserved
            # (first occurrence of max in each row, skipping NaNs -- matches
            # pandas idxmax(axis=1) with skipna=True). Considers all 4 columns,
            # matching the original's `(-Ratio_dif / MinRat).idxmax(axis=1)` where
            # Ratio_dif/MinRat are the full 4-column (C,N,P,E) frames -- empirically
            # (confirmed against the original) 'E' is never actually the argmax in
            # practice, since MinRat['E'] (~5000) dwarfs any plausible C/N/P deficit
            # ratio, so index 3 is guarded against below purely as a safety net.
            with np.errstate(invalid="ignore", divide="ignore"):
                metric4 = -Ratio_dif4 / MinRat4
            metric4_safe = np.where(np.isnan(metric4), -np.inf, metric4)
            limiting_pos = np.argmax(
                metric4_safe, axis=1
            )  # 0=C, 1=N, 2=P, (3=E, guarded below)

            # Set all deficient ratios to their minima
            mask_deficient4 = Ratio_dif4 < 0
            mask_excess4 = Ratio_dif4 > 0
            MicrobeRatios4 = MicrobeRatios4.copy()
            MicrobeRatios4[mask_deficient4] = MinRat4[mask_deficient4]
            # Reduce the mass fractions for non-deficient elements in proportion to the distance from the minimum
            sum_deficit = Ratio_dif_0_4.sum(axis=1)
            sum_excess = Excess4.sum(axis=1)
            with np.errstate(invalid="ignore", divide="ignore"):
                scale = sum_deficit / sum_excess
                MicrobeRatios4[mask_excess4] += (Excess4 * scale[:, None])[mask_excess4]

            # Construct hypothetical nutrient quotas for each possible minimum nutrient
            MC, MN, MP = Mic_subset[:, 0], Mic_subset[:, 1], Mic_subset[:, 2]
            MRC, MRN, MRP = (
                MicrobeRatios4[:, 0],
                MicrobeRatios4[:, 1],
                MicrobeRatios4[:, 2],
            )

            with np.errstate(invalid="ignore", divide="ignore"):
                new_C = np.stack([MC, MN * MRC / MRN, MP * MRC / MRP], axis=1)
                new_N = np.stack([MC * MRN / MRC, MN, MP * MRN / MRP], axis=1)
                new_P = np.stack([MC * MRP / MRC, MN * MRP / MRN, MP], axis=1)
            for arr in (new_C, new_N, new_P):
                arr[~np.isfinite(arr)] = np.float32(0)

            # Insert the appropriate set of nutrient quotas scaled to the minimum nutrient.
            # Safety net: if 'E' (index 3) were ever the limiting column -- which never
            # happens in practice (see note above) -- new_C/new_N/new_P have no 4th
            # column to index, so fall back to the argmax restricted to C/N/P only for
            # those (should-never-occur) rows rather than raising.
            if np.any(limiting_pos > 2):
                bad = limiting_pos > 2
                limiting_pos = limiting_pos.copy()
                limiting_pos[bad] = np.argmax(metric4_safe[bad][:, :3], axis=1)
            rows = np.arange(len(idx2))
            C_final = new_C[rows, limiting_pos]
            N_final = new_N[rows, limiting_pos]
            P_final = new_P[rows, limiting_pos]

            # Update Microbes
            Mic[idx2, 0] = C_final
            Mic[idx2, 1] = N_final
            Mic[idx2, 2] = P_final

            # Sum up the element losses from biomass across whole grid and calculate average loss
            MicLoss = StartMicrobes - Mic[idx2]
            tax_resp_mort = (
                self.Taxon_resp_df.to_numpy(copy=True)
                .ravel(order="C")
                .astype(np.float32)
            )
            tax_resp_mort[idx2] += MicLoss[:, 0]
            self.Taxon_Resp = pd.Series(tax_resp_mort, index=Mic_index)
            # Update monomer pools
            self.Monomers.loc[is_NH4, "N"] += np.float32(MicLoss[:, 1].sum() / gridsize)
            self.Monomers.loc[is_PO4, "P"] += np.float32(MicLoss[:, 2].sum() / gridsize)

            # Update Substrates pool by adding dead microbial biomass
            self.Substrates.loc[is_DeadMic] += Death_gridcell
        # End of if else clause

        # Check microbes energy ratio and set to max (10000) if it is higher
        with np.errstate(invalid="ignore", divide="ignore"):
            Tax_E_ratio = Mic[:, 3] / Mic[:, 0]
        e_cap_index = Tax_E_ratio > self.E_max
        Mic[e_cap_index, 3] = Mic[e_cap_index, 0] * self.E_max

        # Leaching of monomers
        Leaching_N = (
            self.Monomers.loc[is_NH4, "N"]
            * Leaching
            * np.exp(Psi_slope_leach * (self.psi[day] - self.wp_fc))
        )
        Leaching_P = (
            self.Monomers.loc[is_PO4, "P"]
            * Leaching
            * np.exp(Psi_slope_leach * (self.psi[day] - self.wp_fc))
        )
        # Update Monomers
        self.Monomers.loc[is_NH4, "N"] -= Leaching_N
        self.Monomers.loc[is_PO4, "P"] -= Leaching_P

        # Restore Microbes to a DataFrame with the original taxon-repeated index
        self.Microbes = pd.DataFrame(Mic, index=Mic_index, columns=["C", "N", "P", "E"])
        # system respiration
        self.Respiration = self.Taxon_Resp.sum()

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
        Mic = self.Microbes.to_numpy(dtype=np.float32, copy=True)  # (N,4): C,N,P,E
        fb = np.asarray(self.fb)  # (N,) 0/1

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
        Fungi_count = fungi_flags.reshape(gridsize, n_taxa).sum(axis=0)  # (n_taxa,)

        # Derive average biomass of fungal taxa
        Microbes_grid = Mic.reshape(gridsize, n_taxa, 4).sum(axis=0)  # (n_taxa,4)
        with np.errstate(invalid="ignore", divide="ignore"):
            Mean_fungi = Microbes_grid / Fungi_count[:, None]
        Mean_fungi[Fungi_count == 0] = np.float32(0)
        # Expand the fungal average across the grid
        eMF = np.tile(Mean_fungi, (gridsize, 1)).astype(np.float32)  # (N,4)

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
            np.arange(-self.dist, self.dist + 1), num_b, replace=True
        ).astype("int8")
        shift_y[daughters_b] = np.random.choice(
            np.arange(-self.dist, self.dist + 1), num_b, replace=True
        ).astype("int8")
        # Fungi always move positively in x direction, and in y direction constrained to one box away determined by probability "direct"
        shift_x[daughters_f] = np.int8(1)
        shift_y[daughters_f] = np.random.choice(
            [-1, 0, 1],
            num_f,
            replace=True,
            p=[0.5 * (1 - self.direct), self.direct, 0.5 * (1 - self.direct)],
        ).astype("int8")
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
        Colonization = np.zeros((N, 4), dtype=np.float32)
        np.add.at(Colonization, index_series[daughters_b], Reprod[daughters_b])
        np.add.at(Colonization, index_series[daughters_f], Reprod[daughters_f])
        # Colonization of dispersing microbes
        Mic += Colonization

        self.Microbes = pd.DataFrame(Mic, index=Mic_index, columns=["C", "N", "P", "E"])

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
            cum_abundance = output.Taxon_count_repop.iloc[
                :, (day + 2 - self.cycle) : (day + 2)
            ].sum(axis=1)

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
            choose_taxa = np.array([0] * self.gridsize * self.n_taxa).reshape(
                self.n_taxa, self.gridsize
            )
            for i in range(self.n_taxa):
                # Alternative 1
                choose_taxa[i, :] = np.random.choice(
                    [1, 0], self.gridsize, replace=True, p=probs.iloc[i, :]
                )

            # Note order='F'
            self.Microbes.loc[np.ravel(choose_taxa, order="F") == 0] = np.float32(0)

        # Record the post-repopulation microbial C for mass balance tracking.
        # next_pulse is the 1-based index of the pulse about to begin.
        next_pulse = (
            day // self.cycle
        ) + 2  # day is last day of current pulse (0-based)
        output.record_pulse_init(self, next_pulse)
