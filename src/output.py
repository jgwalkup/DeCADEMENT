# output.py module dealing with outputs of DEMENTpy.
# Bin Wang in January, 2020
# modified by Jeth Walkup September, 2024

import numbers
import warnings
# new
from pathlib import Path

# to here
import numpy as np
import pandas as pd
import xarray as xr

from initialization import export_initialization_dict_to_csv
from initialization import export_initialization_dict_to_netcdf
from initialization import export_initialization_dict_to_zarr

# ── Fast groupby-sum helpers (mirrors grid.py) ────────────────────────────────
def _taxa_sum(arr, gridsize, n_taxa, index):
    """Sum (gridsize*n_taxa, ...) over the grid axis → (n_taxa, ...)."""
    vals = arr.values if hasattr(arr, "values") else np.asarray(arr)
    if vals.ndim == 1:
        result = vals.reshape(gridsize, n_taxa).sum(axis=0)
        return pd.Series(result, index=index)
    else:
        result = vals.reshape(gridsize, n_taxa, vals.shape[1]).sum(axis=0)
        return pd.DataFrame(result, index=index, columns=arr.columns)

def _enz_sum(arr, gridsize, n_enzymes, index):
    """Sum (gridsize*n_enzymes,) over the grid axis → (n_enzymes,)."""
    vals = arr.values if hasattr(arr, "values") else np.asarray(arr)
    return pd.Series(vals.reshape(gridsize, n_enzymes).sum(axis=0), index=index)


class Output():
    """
    This class deals with outputs.

    Accepts data derived from the initialization.py and grid.py modules, and have two methods:
        output():             stores all time series
        microbes_abundance(): a special method
    """

    def __init__(self,runtime,data_init):
        """
        The constructor of Output class.

        Parameters:
            runtime:   user-specified parameters when running the model
            data_init: data dictionary;all initialized data from the 'initialization.py' module
        Returns:
            Initialization:   all data initialized preceding the execution of grid.py: dictionary
            Microbial_traits: microbial traits only pulled out from Initialization: dataframe
            SubstratesSeries: substrate-specific total mass over the grid: Substrate * C
            Substrates_Sum:   total substrates over the grid: day * C
            MonomersSeries:   monomoer-specific total mass over the grid: Monomer * C
            Monomers_Sum:     total monomers over the grid: day * C
            MicrobesSeries:   taxon-specific total biomass over the grid: Taxon * (C,N,P)
            Microbes_Sum:     total biomass over the grid: day * (C,N,P)
            TransporterSeries: taxon-specific total transporter production over the grid
            EnzymesSeries:    enzyme-specific total mass over the grid: Enzyme * C
            Enzymes_Sum:      total enzyme summed up over the grid: list
            OsmolyteSeries:   taxon-specific total osmolyte production over the grid
            RespSeries:       total respiration over the grid
            CUE_system:       Carbon Use Efficiency (CUE) for the system
            CAE_system:       Carbon Assimilation Efficiency (CAE) for the system
            CAE_TaxonSeries:  taxon specific Carbon Assimilation Efficiency (CAE)
            CUE_TaxonSeries:  taxon specific Carbon Use Efficiency (CUE)
            Growth_Yield:     growth yield
        Notes:
            Variables that are not tracked but instead by DEMENT(R version):
                NH4Series:
                PO4Series:
        """

        # A couple of vars used in processing outputs
        n_taxa       = int(runtime.loc['n_taxa',    1])            # number of taxa
        n_enzymes    = int(runtime.loc['n_enzymes', 1])            # number of enzymes
        n_substrates = int(runtime.loc['n_substrates', 1])         # number of substrates
        n_monomers   = n_substrates + 3                            # number of monomers
        gridsize     = int(runtime.loc['gridsize',  1])            # grid cells
        # Use the true taxon identities preserved from the (pool/community)
        # inputs rather than positional labels, so every per-taxon output carries
        # the global super-community ID.  Mirrors the index[:n] pattern used for
        # substrates/monomers below; falls back to positional if unavailable.
        try:
            Mic_index = list(data_init['Microbes_pp'].index[:n_taxa])
            if len(Mic_index) != n_taxa:
                raise ValueError("Microbes_pp index shorter than n_taxa")
        except Exception:
            Mic_index = ["Tax" + str(i) for i in range(1, n_taxa + 1)]   # fallback
        Enz_index = ["Enz" + str(i) for i in range(1, n_enzymes + 1)]  # enzyme index

        # Store dimensions and pre-built indices for fast reshape-based groupby
        self.gridsize     = gridsize
        self.n_taxa       = n_taxa
        self.n_enzymes    = n_enzymes
        self.n_substrates = n_substrates
        self.n_monomers   = n_monomers
        self.n_pulses     = int(runtime.loc['pulse', 1])
        self.cycle        = int(runtime.loc['end_time', 1])
        self._taxa_index  = Mic_index
        self._enz_index   = Enz_index
        self._sub_index   = data_init['Substrates'].index[:n_substrates]
        self._mon_index   = data_init['Monomers'].index[:n_monomers]

        # Pass all runtime parameters to Runtime
        self.Runtime = runtime
        # Pass all initialized data (a dictionary) to 'Initialization'
        self.Initialization = data_init

        # Pull out microbial traits data and put them in a dataframe:Microbial_traits
        data = np.concatenate((data_init['fb'][0:n_taxa],
                               data_init['Microbes_pp']['C'][0:n_taxa],
                               data_init['Microbes_pp']['N'][0:n_taxa],
                               data_init['Microbes_pp']['P'][0:n_taxa],
                               data_init['UptakeGenes'].sum(axis=1)[0:n_taxa],
                               data_init['EnzGenes'].sum(axis=1)[0:n_taxa],
                               data_init['OsmoGenes'].sum(axis=1)[0:n_taxa],
                               data_init['UptakeGenes_trait'][0:n_taxa],
                               data_init['EnzProdConsti_trait'][0:n_taxa],
                               data_init['EnzProdInduci_trait'][0:n_taxa],
                               data_init['OsmoProdConsti_trait'][0:n_taxa],
                               data_init['OsmoProdInduci_trait'][0:n_taxa],
                               data_init['TaxDroughtTol'][0:n_taxa])).reshape(13,n_taxa).T
        columns = ['ForB','C','N','P','Uptake_Gene','Enz_Gene','Osmo_Gene','Uptake_Cost',
                   'Enz_Consti_Cost','Enz_Induci_Cost','Osmo_Consti_Cost','Osmo_Induci_Cost',
                   'Drought_tolerance']
        self.Microbial_traits = pd.DataFrame(data=data, index=Mic_index, columns=columns)


        # Account for inputs in mass balance
        # self.Cum_Substrate      = (data_init['Substrates'].groupby(level=0,sort=False).sum()).sum(axis=0)
        # self.Cum_Monomer        = (data_init['Monomers'].groupby(level=0,sort=False).sum()).sum(axis=0)
        # self.Cum_Monomer_ratios = ecosystem.Cum_Monomer_ratios

        # Degradation rates
        #DecayRates_list       = [0] * n_substrates * gridsize
        #DecayRates_array      = np.array(DecayRates_list).reshape(n_substrates * gridsize,1)
        #DecayRates_df         = pd.DataFrame(data = DecayRates_array, index= data_init['Substrates'].index, columns= ['C'])
        #self.DecayRatesSeries = DecayRates_df.groupby(level=0,sort=False).sum()

        # Leaching
        #self.Cum_Leaching_N   = pd.Series([0],index=[0])
        #self.Cum_Leaching_P   = pd.Series([0],index=[0])

        # Substrates
        Substrates_grid = pd.DataFrame(
            data_init['Substrates'].values.reshape(gridsize, n_substrates, 3).sum(axis=0),
            index=data_init['Substrates'].index[:n_substrates],
            columns=data_init['Substrates'].columns)
        # (name set via .rename(0) below)
        self.SubstratesSeries     = Substrates_grid['C'].rename(0)
        #self.Substrates_Sum      = pd.Series([Substrates_grid['C'].sum()],index=[0])

        # Monomers
        Monomers_grid = pd.DataFrame(
            data_init['Monomers'].values.reshape(gridsize, n_monomers, 3).sum(axis=0),
            index=self._mon_index, columns=data_init['Monomers'].columns)
        # (name set via .rename(0) below)
        self.MonomersSeries     = Monomers_grid["C"].rename(0)
        #self.Monomers_Sum      = pd.Series([sum(Monomers_grid["C"])],index=[0])
        #self.NH4Series         = pd.Series([Monomers_grid.loc["NH4","N"]],index=[0])
        #self.PO4Series         = pd.Series([Monomers_grid.loc["PO4","P"]],index=[0])

        # Microbes
        Microbes_grid = pd.DataFrame(
            data_init['Microbes'].values.reshape(gridsize, n_taxa, 4).sum(axis=0),
            index=self._taxa_index, columns=data_init['Microbes'].columns)
        # (name set via .rename(0) below)                                               # Rename the series name to 0
        self.MicrobesSeries = Microbes_grid['C'].rename(0)  # Taxon-specific total biomass summed over the grid
        self.EnergySeries = Microbes_grid['E'].rename(0)  # Taxon-specific total biomass summed over the grid

        # Taxon HSP: initial state summed over grid, stacked to (taxon, HSP_type) × day
        hsp0 = pd.DataFrame(
            data_init['Initial_HSP_Prop'].drop(columns='Total_biomass').values
            .reshape(gridsize, n_taxa, -1).sum(axis=0),
            index=self._taxa_index,
            columns=data_init['Initial_HSP_Prop'].drop(columns='Total_biomass').columns)
        hsp0_stacked = hsp0.stack()
        hsp0_stacked.name = 0
        self.Taxon_HSP_Series = hsp0_stacked

        # Number of individuals of Taxon count
        #self.Taxon_index            = data_init['Microbes']['C'] > 0
        #self.Taxon_index.name       = 0
        #self.Taxon_count       = self.Taxon_index.groupby(level=0,sort=False).sum().astype('uint32')
        #self.Taxon_count_repop = self.Taxon_index.groupby(level=0,sort=False).sum().astype('uint32')
        # Taxon _count for taxon specific outputs
        #self.taxon_count = pd.Series(data=[0] * n_taxa, index=Mic_index)

        # Taxon count
        self.Taxon_count = pd.Series(data=[0] * n_taxa, index=Mic_index)
        self.Taxon_count_repop = pd.Series(data=[0] * n_taxa, index=Mic_index)

        # Transporters: taxon-specific production summed over the grid by taxon
        #self.TransporterSeries = pd.Series(data=[0]*n_taxa,index=Mic_index)

        # Enyzmes: taxon-specific production summed over the grid by taxon
        #self.EnzymeConSeries    = pd.Series(data=[0]*n_taxa,index=Mic_index)
        #self.EnzymeIndSeries    = pd.Series(data=[0]*n_taxa,index=Mic_index)
        #self.Enzyme_TaxonSeries = pd.Series(data=[0]*n_taxa,index=Mic_index)

        # Osmolytes: taxon-specific production summed over the grid by taxon
        #self.OsmolyteConSeries  = pd.Series(data=[0]*n_taxa,index=Mic_index)
        #self.OsmolyteIndSeries  = pd.Series(data=[0]*n_taxa,index=Mic_index)
        self.Osmolyte_TaxonSeries = pd.Series(data=[0]*n_taxa, index=Mic_index)

        # Growth yield
        self.Growth_Yield = pd.Series(data=[0]*n_taxa,index=Mic_index)

        # Growth respiration
        self.Growth_resp = pd.Series(data=[0]*n_taxa,index=Mic_index)

        # Taxon-specific C assimilation efficiency (CAE)
        self.CAE_TaxonSeries = pd.Series(data=[0] * n_taxa, index=Mic_index)

        # Taxon-specific Carbon Use Efficiency (CUE)
        self.CUE_TaxonSeries = pd.Series(data=[0] * n_taxa, index=Mic_index)

        # Taxon-specific Respiration
        self.Respiration_Taxon = pd.Series(data=[0] * n_taxa, index=Mic_index)

        # Taxon-specific total Carbon Uptake
        self.C_Uptake_Taxon = pd.Series(data=[0] * n_taxa, index=Mic_index)

        # Taxon-specific Carbon Uptake allocated to biomass
        self.BC_Uptake_Taxon = pd.Series(data=[0] * n_taxa, index=Mic_index)

        #Taxon-specific Energy production
        self.E_Prod_Taxon = pd.Series(data=[0] * n_taxa, index=Mic_index)

        # taxon net energy gain
        self.Microbe_E_Gain = pd.Series(data=[0] * n_taxa, index=Mic_index)

        # taxon net energy costs
        self.Microbe_E_Cost = pd.Series(data=[0] * n_taxa, index=Mic_index)

        # Taxon-specific transporter maintenance cost
        #self.Taxon_Transporter_Maint = pd.Series(data=[0] * n_taxa, index=Mic_index)

        # Taxon-specific C assimilation cost
        # self.Taxon_AE = pd.Series(data=[0] * n_taxa, index=Mic_index)

        # Taxon-specific assimilation respiration
        # self.Taxon_AE_resp = pd.Series(data=[0] * n_taxa, index=Mic_index)

        # Taxon-specific cell maintenance cost
        self.Taxon_Cell_Maint = pd.Series(data=[0] * n_taxa, index=Mic_index)

        # Per-taxon C allocation tracking (per alive microbe, time series)
        self.Enzyme_C_Taxon     = pd.Series(data=[0] * n_taxa, index=Mic_index)
        self.HSP_C_Taxon        = pd.Series(data=[0] * n_taxa, index=Mic_index)
        self.Osmo_C_Taxon       = pd.Series(data=[0] * n_taxa, index=Mic_index)
        self.Transport_C_Taxon  = pd.Series(data=[0] * n_taxa, index=Mic_index)

        # Per-taxon E allocation tracking by pathway (per alive microbe, time series)
        self.Transport_E_Taxon = pd.Series(data=[0] * n_taxa, index=Mic_index)
        self.Maint_E_Taxon     = pd.Series(data=[0] * n_taxa, index=Mic_index)
        self.Enzyme_E_Taxon    = pd.Series(data=[0] * n_taxa, index=Mic_index)
        self.Osmo_E_Taxon      = pd.Series(data=[0] * n_taxa, index=Mic_index)
        self.HSP_E_Taxon       = pd.Series(data=[0] * n_taxa, index=Mic_index)


        # Enzymes
        Enzymes_grid = pd.Series(
            data_init['Enzymes'].values.reshape(gridsize, n_enzymes).sum(axis=0),
            index=data_init['Enzymes'].index[:n_enzymes])  # total of each enzyme summed over the spatial grid
        Enzymes_grid.name  = 0
        self.EnzymesSeries = Enzymes_grid
        #self.Enzymes_Sum  = pd.Series([Enzymes_grid.sum()],index=[0])

        # all enzyme production
        #Enzymes_prod = data_init['Enzymes'].groupby(level=0,sort=False).sum()
        #Enzymes_prod.name = 0
        self.EnzymesLoss = pd.Series(data=[0]*n_enzymes, index = Enz_index)
        self.EnzymesGain = pd.Series(data=[0]*n_enzymes, index = Enz_index)
        self.EnzymesProduction = pd.Series(data=[0]*n_enzymes, index = Enz_index)
        # self.EnzymesProduction = pd.Series([0], dtype='float32')


        # Emergent properties over the grid
        self.RespSeries = pd.Series([0],index=[0], dtype='float32')  # respiration
        self.CAE_system = pd.Series([0],index=[0], dtype='float32')  # emergent CAE (Carbon Assimilation Efficiency)
        self.CUE_system = pd.Series([0],index=[0], dtype='float32')  # emergent CUE (Carbon Use Efficiency)
        self.Kill       = pd.Series([0],index=[0], dtype='uint32')   # stochastic death toll

        # Per-pulse cumulative respiration (for mass balance checks)
        # CumRespByPulse: Series indexed by pulse number (1-based), value = total C respired in that pulse.
        # Note: accumulation happens inside output(), which is called every day when interval=1 (the
        # standard setting).  If interval > 1, only days that fall on an output step are accumulated,
        # so CumRespByPulse will under-count.  Use interval=1 for mass balance runs.
        self.n_pulses   = int(runtime.loc['pulse', 1])
        self.cycle      = int(runtime.loc['end_time', 1])   # days per pulse
        self.CumRespByPulse      = pd.Series(data=0.0, index=range(1, self.n_pulses + 1), dtype='float64')
        self._pulse_resp_accum   = 0.0   # running accumulator; reset at start of each pulse

        # Per-pulse initial microbial C (for mass balance checks).
        # Substrates, Monomers, and Enzymes are reset to their t=0 values each pulse, so
        # the t=0 snapshot serves as the initial state for all pulses.  Microbes are
        # different: repopulation() resets to Microbes_pp (pre-placement) and then
        # re-places stochastically, producing a different total each pulse.  We record
        # the actual post-repopulation microbial C here for each pulse so the mass
        # balance checker can use the true initial state.  Pulse 1 uses data_init['Microbes']
        # (the t=0 post-placement community); subsequent pulses are filled by
        # record_pulse_init(), called from grid.repopulation().
        mic_C_init = float(
            data_init['Microbes'].values.reshape(gridsize, n_taxa, 4).sum(axis=0)[:, 0].sum())
        self.MicC_PulseInit = pd.Series(data=0.0, index=range(1, self.n_pulses + 1), dtype='float64')
        self.MicC_PulseInit[1] = mic_C_init


        # ── Per-day accumulation buffers (list-of-columns → single concat at export) ──
        # Replaces per-day pd.concat (O(n²)) with O(1) list.append.
        # _assemble_series() converts buffers to DataFrames once at export.
        self._buf = {}
        for _s in [
            "SubstratesSeries", "MonomersSeries", "Taxon_count", "MicrobesSeries",
            "EnergySeries", "Taxon_HSP_Series", "C_Uptake_Taxon", "BC_Uptake_Taxon",
            "Enzyme_C_Taxon", "HSP_C_Taxon", "Osmo_C_Taxon", "Transport_C_Taxon",
            "Transport_E_Taxon", "Maint_E_Taxon", "Enzyme_E_Taxon", "Osmo_E_Taxon",
            "HSP_E_Taxon", "E_Prod_Taxon", "Microbe_E_Gain", "Microbe_E_Cost",
            "CAE_TaxonSeries", "CAE_system", "CUE_TaxonSeries", "CUE_system",
            "Respiration_Taxon", "RespSeries", "EnzymesSeries", "EnzymesProduction",
            "EnzymesGain", "EnzymesLoss", "Kill", "Taxon_count_repop",
        ]:
            self._buf[_s] = []

    def output(self,ecosystem,day):
        """
        Records outputs in various variables of each iteration.

        Parameters:
            ecosystem: object from the grid.py module
            day:       the day to record to outputs
        Returns:
        """

        # Account for inputs in mass balance
        # self.Cum_Substrate = self.Cum_Substrate + ecosystem.SubInput.sum(axis = 0)
        # self.Cum_Monomer = self.Cum_Monomer + (ecosystem.MonInput.mul(ecosystem.Cum_Monomer_ratios,axis=0)).sum(axis=0)

        # DecayRates
        # DecayRates_grid = ecosystem.DecayRates.groupby(level=0,sort=False).sum()
        # self.DecayRatesSeries = pd.concat([self.DecayRatesSeries,DecayRates_grid],axis=1,sort=False)

        # Substrates
        Substrates_grid = pd.DataFrame(
            ecosystem.Substrates.values.reshape(self.gridsize, self.n_substrates, 3).sum(axis=0),
            index=self._sub_index, columns=ecosystem.Substrates.columns)
        # (name set via .rename(day+1) below)  # index the output by day
        self._buf["SubstratesSeries"].append(Substrates_grid['C'].rename(day + 1))
        #self.Substrates_Sum   = pd.concat([self.Substrates_Sum,pd.Series([Substrates_grid['C'].sum()],index=[day+1])],axis=0,sort=False)

        # Monomers
        Monomers_grid = pd.DataFrame(
            ecosystem.Monomers.values.reshape(self.gridsize, self.n_monomers, 3).sum(axis=0),
            index=self._mon_index, columns=ecosystem.Monomers.columns)
        # (name set via .rename(day+1) below)
        self._buf["MonomersSeries"].append(Monomers_grid['C'].rename(day + 1))
        #self.Monomers_Sum   = pd.concat([self.Monomers_Sum, pd.Series([sum(Monomers_grid["C"])],  index=[day+1])], axis=0, sort=False)
        #self.NH4Series      = pd.concat([self.NH4Series, pd.Series([Monomers_grid.loc["NH4","N"]],index=[day+1])], axis=0, sort=False)
        #self.PO4Series      = pd.concat([self.PO4Series, pd.Series([Monomers_grid.loc["PO4","P"]],index=[day+1])], axis=0, sort=False)

        # Interim Microbes
        #Microbes_interim_grid = ecosystem.Microbes_interim.groupby(level=0,sort=False).sum()
        #self.Microbes_Interim = pd.concat([self.MicrobesSeries,Microbes_interim_grid['C']],axis=1,sort=False)
        # Count taxon for averaging taxon CAE
        #Taxon_index = (ecosystem.Microbes_interim)['C'] > 0
        #Taxon_index.name = day + 1
        #taxon_count = Taxon_index.groupby(level=0,sort=False).sum()
        #self.Taxon_count = pd.concat([self.Taxon_count,taxon_count],axis=1,sort=False)

        # Microbe
        ## Taxon abundance
        self.Taxon_index      = ecosystem.Microbes['C'] > 0
        self.Taxon_index.name = day + 1
        self._buf["Taxon_count"].append(
            pd.Series(self.Taxon_index.values.reshape(self.gridsize, self.n_taxa).sum(axis=0).astype('uint32'),
                      index=self._taxa_index, name=day + 1))
        self.taxon_count = pd.Series(
            self.Taxon_index.values.reshape(self.gridsize, self.n_taxa).sum(axis=0),
            index=self._taxa_index)

        ## Taxon biomass
        Microbes_grid = pd.DataFrame(
            ecosystem.Microbes.values.reshape(self.gridsize, self.n_taxa, 4).sum(axis=0),
            index=self._taxa_index, columns=ecosystem.Microbes.columns)
        # (name set via .rename(day+1) below)
        self._buf["MicrobesSeries"].append(Microbes_grid['C'].rename(day + 1))
        #self.Microbes_Sum   = pd.concat([self.Microbes_Sum,pd.Series([Microbes_grid['C'].sum()],index=[day+1])],axis=0,sort=False)
        self._buf["EnergySeries"].append(Microbes_grid['E'].rename(day + 1))

        # Taxon HSP: sum over grid, stack (taxon, HSP_type), append as new day column
        hsp_grid = pd.DataFrame(
            ecosystem.Taxon_HSP.drop(columns='Total_biomass').values
            .reshape(self.gridsize, self.n_taxa, -1).sum(axis=0),
            index=self._taxa_index,
            columns=ecosystem.Taxon_HSP.drop(columns='Total_biomass').columns)
        hsp_stacked = hsp_grid.stack()
        hsp_stacked.name = day + 1
        self._buf["Taxon_HSP_Series"].append(hsp_stacked)

        # Microbe-Taxon-specific C Uptake
        #C_Uptake_Taxon_grid = _taxa_sum(ecosystem.Taxon_C_Uptake, self.gridsize, self.n_taxa, self._taxa_index)
        C_Uptake_Taxon_grid = _taxa_sum(ecosystem.Taxon_C_Uptake, self.gridsize, self.n_taxa, self._taxa_index)
        C_Uptake_Taxon_grid = C_Uptake_Taxon_grid/self.taxon_count
        C_Uptake_Taxon_grid[np.isinf(C_Uptake_Taxon_grid)] = np.nan
        C_Uptake_Taxon_grid.name = day + 1
        self._buf["C_Uptake_Taxon"].append(C_Uptake_Taxon_grid)

        #C_Uptake_Taxon_grid = ecosystem.Taxon_C_Uptake.groupby(level=0,sort=False).sum()
        BC_Uptake_Taxon_grid = _taxa_sum(ecosystem.Taxon_BC_Uptake, self.gridsize, self.n_taxa, self._taxa_index)
        BC_Uptake_Taxon_grid = BC_Uptake_Taxon_grid/self.taxon_count
        BC_Uptake_Taxon_grid[np.isinf(BC_Uptake_Taxon_grid)] = np.nan
        BC_Uptake_Taxon_grid.name = day + 1
        self._buf["BC_Uptake_Taxon"].append(BC_Uptake_Taxon_grid)

        # Per-taxon C allocation: enzyme, HSP, osmolyte, transporter production (per alive microbe)
        for attr, eco_attr in [("Enzyme_C_Taxon",    "Taxon_Enzyme_Cost_C"),
                                ("HSP_C_Taxon",       "Taxon_HSP_Cost_C"),
                                ("Osmo_C_Taxon",      "Taxon_Osmo_Cost_C"),
                                ("Transport_C_Taxon", "Taxon_Transport_Cost_C")]:
            grid_val = _taxa_sum(getattr(ecosystem, eco_attr), self.gridsize, self.n_taxa, self._taxa_index)
            grid_val = grid_val / self.taxon_count
            grid_val[np.isinf(grid_val)] = np.nan
            grid_val.name = day + 1
            self._buf[attr].append(grid_val)

        # Per-taxon E allocation by pathway: transport, maintenance, enzyme, osmolyte, HSP
        for attr, eco_attr in [("Transport_E_Taxon", "Taxon_Transport_E_Cost"),
                                ("Maint_E_Taxon",     "Taxon_Maint_E_Cost"),
                                ("Enzyme_E_Taxon",    "Taxon_Enzyme_E_Cost"),
                                ("Osmo_E_Taxon",      "Taxon_Osmo_E_Cost_total"),
                                ("HSP_E_Taxon",       "Taxon_HSP_E_Cost")]:
            grid_val = _taxa_sum(getattr(ecosystem, eco_attr), self.gridsize, self.n_taxa, self._taxa_index)
            grid_val = grid_val / self.taxon_count
            grid_val[np.isinf(grid_val)] = np.nan
            grid_val.name = day + 1
            self._buf[attr].append(grid_val)

        # Microbe-Taxon-specific C Uptake
        E_Prod_Taxon_grid = _taxa_sum(ecosystem.Taxon_E_Prod, self.gridsize, self.n_taxa, self._taxa_index)
        E_Prod_Taxon_grid = E_Prod_Taxon_grid/self.taxon_count
        E_Prod_Taxon_grid[np.isinf(E_Prod_Taxon_grid)] = np.nan
        E_Prod_Taxon_grid.name = day + 1
        self._buf["E_Prod_Taxon"].append(E_Prod_Taxon_grid)

        # Microbe-Taxon-specific net energy gain
        Microbe_E_Gain_grid = _taxa_sum(ecosystem.Microbe_E_Gain, self.gridsize, self.n_taxa, self._taxa_index)
        Microbe_E_Gain_grid = Microbe_E_Gain_grid/self.taxon_count
        Microbe_E_Gain_grid[np.isinf(Microbe_E_Gain_grid)] = np.nan
        Microbe_E_Gain_grid.name = day + 1
        self._buf["Microbe_E_Gain"].append(Microbe_E_Gain_grid)

        # Microbe-Taxon-specific net energy cost
        Microbe_E_Cost_grid = _taxa_sum(ecosystem.Microbe_E_Cost, self.gridsize, self.n_taxa, self._taxa_index)
        Microbe_E_Cost_grid = Microbe_E_Cost_grid/self.taxon_count
        Microbe_E_Cost_grid[np.isinf(Microbe_E_Cost_grid)] = np.nan
        Microbe_E_Cost_grid.name = day + 1
        self._buf["Microbe_E_Cost"].append(Microbe_E_Cost_grid)

        # Microbe-Taxon-specific Monomer transporter maintenance
        # Taxon_Transporter_Maint_grid = ecosystem.Taxon_Transporter_Maint.groupby(level=0,sort=False).sum()
        # Taxon_Transporter_Maint_grid = Taxon_Transporter_Maint_grid
        # Taxon_Transporter_Maint_grid.name = day + 1
        # self.Taxon_Transporter_Maint = pd.concat([self.Taxon_Transporter_Maint, Taxon_Transporter_Maint_grid],axis=1,sort=False)

        # Microbe-Taxon-specific C assimilation efficiency
        # Taxon_AE_grid = ecosystem.Taxon_AE.groupby(level=0,sort=False).sum()
        # Taxon_AE_grid = Taxon_AE_grid
        # Taxon_AE_grid.name = day + 1
        # self.Taxon_AE = pd.concat([self.Taxon_AE, Taxon_AE_grid],axis=1,sort=False)

        # Microbe-Taxon-specific C assimilation efficiency respiration (?)
        # Taxon_AE_resp_grid = ecosystem.Taxon_AE_Resp.groupby(level=0,sort=False).sum()
        # Taxon_AE_resp_grid = Taxon_AE_resp_grid
        # Taxon_AE_resp_grid.name = day + 1
        # self.Taxon_AE_resp = pd.concat([self.Taxon_AE_resp, Taxon_AE_resp_grid],axis=1,sort=False)


        # Microbe-Taxon-specific C Uptake
        # Taxon_Cell_Maint_grid = ecosystem.Taxon_Cell_Maint.groupby(level=0,sort=False).sum()
        # Taxon_Cell_Maint_grid = Taxon_Cell_Maint_grid
        # Taxon_Cell_Maint_grid.name = day + 1
        # self.Taxon_Cell_Maint = pd.concat([self.Taxon_Cell_Maint, Taxon_Cell_Maint_grid],axis=1,sort=False)


        # Microbe-Taxon-specific C Uptake
        # Respiration_Taxon_grid = ecosystem.Taxon_Respiration.groupby(level=0,sort=False).sum()
        # Respiration_Taxon_grid = Respiration_Taxon_grid
        # Respiration_Taxon_grid.name = day + 1
        # self.Respiration_Taxon = pd.concat([self.Respiration_Taxon,Respiration_Taxon_grid],axis=1,sort=False)

        # Microbe-Transporters
        #Transporter_grid = ecosystem.Transporters.groupby(level=0,sort=False).sum() #taxon-specific transpotors summed over the grid by taxon
        #Transporter_grid.name = day + 1
        #self.TransporterSeries = pd.concat([self.TransporterSeries,Transporter_grid],axis=1,sort=False)

        # Microbe-Enzymes
        # Constitutive
        #Enzyme_Con_grid = ecosystem.Enzyme_Con.groupby(level=0,sort=False).sum() # Taxon-specific osmolytes summed over the grid by taxon
        #Enzyme_Con_grid.name = day + 1
        #self.EnzymeConSeries = pd.concat([self.EnzymeConSeries,Enzyme_Con_grid],axis=1,sort=False)
        # Inducible
        #Enzyme_Ind_grid = ecosystem.Enzyme_Ind.groupby(level=0,sort=False).sum()
        #Enzyme_Ind_grid.name = day + 1
        #self.EnzymeIndSeries = pd.concat([self.EnzymeIndSeries,Enzyme_Ind_grid],axis=1,sort=False)
        # Total
        #self.Enzyme_TaxonSeries = pd.concat([self.Enzyme_TaxonSeries,Enzyme_Con_grid+Enzyme_Ind_grid],axis=1,sort=False)

        # Microbe-Osmolytes
        # Constitutive
        #Osmolyte_Con_grid = ecosystem.Osmolyte_Con.groupby(level=0,sort=False).sum()
        #Osmolyte_Con_grid.name = day + 1
        #self.OsmolyteConSeries = pd.concat([self.OsmolyteConSeries,Osmolyte_Con_grid],axis=1,sort=False)
        # Inducible
        #Osmolyte_Ind_grid = ecosystem.Osmolyte_Ind.groupby(level=0,sort=False).sum()
        #Osmolyte_Ind_grid.name = day + 1
        #self.OsmolyteIndSeries = pd.concat([self.OsmolyteIndSeries,Osmolyte_Ind_grid],axis=1,sort=False)
        # Total
        #self.Osmolyte_TaxonSeries = pd.concat([self.Osmolyte_TaxonSeries,Osmolyte_Con_grid+Osmolyte_Ind_grid],axis=1,sort=False)

        # Growth yield by Taxon
        # GY_grid = ecosystem.Growth_Yield.groupby(level=0,sort=False).sum()
        # GY_grid = GY_grid/self.taxon_count
        # GY_grid.name = day + 1
        # self.Growth_Yield = pd.concat([self.Growth_Yield,GY_grid],axis=1,sort=False)

        # Growth respiration by Taxon
        # Growth_resp_grid = ecosystem.Taxon_Growth_Respiration.groupby(level=0,sort=False).sum()
        # Growth_resp_grid = Growth_resp_grid/self.taxon_count
        # Growth_resp_grid.name = day + 1
        # self.Growth_resp = pd.concat([self.Growth_resp,Growth_resp_grid],axis=1,sort=False)

        # Microbe-Taxon-specific CAE
        CAE_Taxon_grid = _taxa_sum(ecosystem.CAE_taxon, self.gridsize, self.n_taxa, self._taxa_index)
        CAE_Taxon_grid = CAE_Taxon_grid/self.taxon_count
        CAE_Taxon_grid.name = day + 1
        self._buf["CAE_TaxonSeries"].append(CAE_Taxon_grid)

        # Carbon assimilation efficiency
        CAE_grid = pd.Series(data=[ecosystem.CAE_system])
        CAE_grid.name = day + 1
        self._buf["CAE_system"].append(CAE_grid)

        # Microbe-Taxon-specific CUE
        CUE_Taxon_grid = _taxa_sum(ecosystem.CUE_taxon, self.gridsize, self.n_taxa, self._taxa_index)
        CUE_Taxon_grid = CUE_Taxon_grid/self.taxon_count
        CUE_Taxon_grid.name = day + 1
        self._buf["CUE_TaxonSeries"].append(CUE_Taxon_grid)

        # Carbon use efficiency
        CUE_grid = pd.Series(data = [ecosystem.CUE_system])
        CUE_grid.name = day + 1
        self._buf["CUE_system"].append(CUE_grid)

        # Respiration by Taxon
        Respiration_grid = _taxa_sum(ecosystem.Taxon_Resp, self.gridsize, self.n_taxa, self._taxa_index)
        Respiration_grid = Respiration_grid/self.taxon_count
        Respiration_grid.name = day + 1
        self._buf["Respiration_Taxon"].append(Respiration_grid)

        # Respiration
        Resp_grid = pd.Series(data = [ecosystem.Taxon_Resp.sum()])
        Resp_grid.name = day + 1
        self._buf["RespSeries"].append(Resp_grid)

        # Per-pulse cumulative respiration accumulation
        # 'day' is the global 0-based index across all pulses.
        # The last day of pulse p (1-based) has global index: p * self.cycle - 1.
        # We detect end-of-pulse when (day + 1) is a multiple of self.cycle.
        self._pulse_resp_accum += float(ecosystem.Taxon_Resp.sum())
        _day_in_pulse = (day + 1) % self.cycle   # 0 means last day of a pulse
        if _day_in_pulse == 0:
            _pulse_num = (day + 1) // self.cycle  # which pulse just ended (1-based)
            self.CumRespByPulse[_pulse_num] = self._pulse_resp_accum
            self._pulse_resp_accum = 0.0          # reset for next pulse

        # Enzymes over the grid
        # Derive the enzyme-specific enzyme production summed over the grid by enzyme
        Enzymes_grid = _enz_sum(ecosystem.Enzymes, self.gridsize, self.n_enzymes, self._enz_index)
        Enzymes_grid.name  = day + 1
        self._buf["EnzymesSeries"].append(Enzymes_grid)
        #self.Enzymes_Sum  = pd.concat([self.Enzymes_Sum,pd.Series([Enzymes_grid.sum()],index=[day+1])],axis=0,sort=False)

        Enzymes_prod = _enz_sum(ecosystem.Enzyme_Production, self.gridsize, self.n_enzymes, self._enz_index)
        Enzymes_prod.name = day + 1
        self._buf["EnzymesProduction"].append(Enzymes_prod)

        Enzymes_gain = _enz_sum(ecosystem.Enzymes_gain, self.gridsize, self.n_enzymes, self._enz_index)
        Enzymes_gain.name = day + 1
        self._buf["EnzymesGain"].append(Enzymes_gain)

        Enzymes_loss = _enz_sum(ecosystem.Enzyme_Loss, self.gridsize, self.n_enzymes, self._enz_index)
        Enzymes_loss.name = day + 1
        self._buf["EnzymesLoss"].append(Enzymes_loss)

        # Death toll of stochasticity origin
        #Kill_data = ecosystem.Kill
        #Kill_data.name = day + 1
        #self.Kill       = pd.concat([self.Kill, Kill_data], axis=0, sort=False)
        self._buf["Kill"].append(pd.Series([ecosystem.Kill], index=[day + 1], dtype='uint32'))

    def record_pulse_init(self, ecosystem, next_pulse):
        """Record the post-repopulation microbial C for mass balance tracking.

        Called from grid.repopulation() immediately after the microbial community
        has been reset and re-placed for the next pulse.  Stores the actual total
        microbial C into MicC_PulseInit[next_pulse] so that mass_balance.py can
        use the true (stochastic) initial microbial state for each pulse rather
        than the t=0 snapshot, which differs because repopulation() resets to
        Microbes_pp (pre-placement) before re-placing, not to the t=0 post-
        placement community.

        Parameters:
            ecosystem:   the Grid object, after repopulation has completed
            next_pulse:  1-based index of the pulse that is about to begin
        """
        if next_pulse > self.n_pulses:
            return
        mic_C = float(ecosystem.Microbes.values.reshape(
            ecosystem.gridsize, ecosystem.n_taxa, 4).sum(axis=0)[:, 0].sum())
        self.MicC_PulseInit[next_pulse] = mic_C

    def microbes_abundance(self,ecosystem,day):
        """
        Seperately output Microbes from each time step and put them in a dataframe.

        Aims to eventually deal with reinitializing microbial community on the
        grid in a new pulse via calculating the cumulative frequency of different
        taxa over each cycle/pulse. The reason to have this separate method instead of using the method
        above is b/c of the need to track outputs in every iteration, whereas the
        method above tracks outputs with a certain time interval (unless the interval set to 1)

        Parameters:
            ecosystem: an instance of the Grid object, in which only the Microbes is used.
            day:       the iteration index
        Returns:
            MicrobesSeries_repop: dataframe; tracking total biomass of differing taxa
            Taxon_count_repop:    dataframe; tracking abundances of differing taxa
        """
        # Track biomass of every taxon
        #Microbes_grid             = ecosystem.Microbes.groupby(level=0,sort=False).sum()
        #Microbes_grid['C'].name   = day + 1
        #self.MicrobesSeries_repop = pd.concat([self.MicrobesSeries_repop,Microbes_grid['C']],axis=1,sort=False)

        # Track abundance of every taxon
        Taxon_index      = ecosystem.Microbes['C'] > 0
        Taxon_index.name = day + 1
        taxon_count      = pd.Series(
            Taxon_index.values.reshape(self.gridsize, self.n_taxa).sum(axis=0).astype('uint32'),
            index=self._taxa_index)
        self._buf["Taxon_count_repop"].append(taxon_count)

    def microbes_tradeoff(self, ecosystem, year, day):
        """
        Track allocation of C uptake to different components.

        These components include: osmolytes, enzymes, and carbon gain
        """

        # Microbe-Transporters
        # Transporter_grid = ecosystem.Transporters.groupby(level=0,sort=False).sum() #taxon-specific transpotors summed over the grid by taxon
        # Transporter_grid.name = day + 1
        # self.TransporterSeries = pd.concat([self.TransporterSeries,Transporter_grid],axis=1,sort=False)

        # Microbe-Enzymes
        # Constitutive
        # Enzyme_Con_grid = ecosystem.Enzyme_Con.groupby(level=0,sort=False).sum() # Taxon-specific osmolytes summed over the grid by taxon
        # Enzyme_Con_grid.name = day + 1
        # self.EnzymeConSeries = pd.concat([self.EnzymeConSeries,Enzyme_Con_grid],axis=1,sort=False)
        # Inducible
        Enzyme_Ind_grid = _taxa_sum(ecosystem.Taxon_Enzyme_Cost_C, self.gridsize, self.n_taxa, self._taxa_index)
        Enzyme_Ind_grid.name = self.cycle * year + (day + 1)
        # self.EnzymeIndSeries = pd.concat([self.EnzymeIndSeries,Enzyme_Ind_grid],axis=1,sort=False)
        # Total
        self.Enzyme_TaxonSeries = pd.concat([self.Enzyme_TaxonSeries, Enzyme_Ind_grid], axis=1, sort=False)

        # Microbe-Osmolytes
        # Constitutive
        # Osmolyte_Con_grid = ecosystem.Osmolyte_Con.groupby(level=0,sort=False).sum()
        # Osmolyte_Con_grid.name = day + 1
        # self.OsmolyteConSeries = pd.concat([self.OsmolyteConSeries,Osmolyte_Con_grid],axis=1,sort=False)
        # Inducible
        Osmolyte_Ind_grid = _taxa_sum(ecosystem.Taxon_Osmo_Cost_C, self.gridsize, self.n_taxa, self._taxa_index)
        Osmolyte_Ind_grid.name = self.cycle * year + (day + 1)
        # self.OsmolyteIndSeries = pd.concat([self.OsmolyteIndSeries,Osmolyte_Ind_grid],axis=1,sort=False)
        # Total
        self.Osmolyte_TaxonSeries = pd.concat([self.Osmolyte_TaxonSeries, Osmolyte_Ind_grid], axis=1, sort=False)

        # Growth yield by Taxon
        GY_grid = _taxa_sum(ecosystem.Microbe_C_Gain, self.gridsize, self.n_taxa, self._taxa_index)
        GY_grid.name = self.cycle * year + (day + 1)
        self.Growth_yield = pd.concat([self.Growth_yield, GY_grid], axis=1, sort=False)





    def flush_pulse_buffers(self):
        """Assemble Taxon_count_repop from its list buffer into a real DataFrame.

        Call at the end of each pulse, after the last output() / microbes_abundance()
        call but before repopulation() reads Taxon_count_repop via .iloc[].
        Keeps the series readable mid-run while eliminating quadratic concat growth.
        All other series remain buffered until export.
        """
        buf = self._buf.get("Taxon_count_repop", [])
        if not buf:
            return
        new_cols = pd.concat(buf, axis=1, sort=False)
        existing = self.Taxon_count_repop
        if existing is not None and len(getattr(existing, "columns", [])) > 0:
            self.Taxon_count_repop = pd.concat([existing, new_cols],
                                               axis=1, sort=False)
        else:
            self.Taxon_count_repop = new_cols
        self._buf["Taxon_count_repop"] = []

    def _assemble_series(self):
        """Convert per-day list buffers to DataFrames (called once per export).

        Replaces the incremental pd.concat pattern (O(n²)) with a single
        pd.concat per series (O(n)).  The initial-state column stored in
        self.NAME during __init__ is prepended automatically.
        """
        for name, buf in self._buf.items():
            if not buf:
                continue
            axis = 0 if name == "Kill" else 1
            new_data = pd.concat(buf, axis=axis, sort=False)
            existing = getattr(self, name, None)
            if existing is not None and not (
                    hasattr(existing, "__len__") and len(existing) == 0):
                try:
                    setattr(self, name, pd.concat([existing, new_data],
                                                  axis=axis, sort=False))
                except Exception:
                    setattr(self, name, new_data)
            else:
                setattr(self, name, new_data)
            self._buf[name] = []

    def export_to_csv(self, base_path: Path | str) -> None:
        """Export contents of the output file to a directory.
        Exports each class member of type pandas.DataFrame to a separate CSV file.
        All pandas.Series members are combined in a DataFrame and printed dto 'series.csv' file.
        Similarly all scalar numerical members are grouped in 'scalars.csv'.
        Parameters:
          base_path : Path
            A path that names the root directory where contents will be exported.
            If the directory does not exist it will be created.
        """
        self._assemble_series()
        # Create space for output
        base_path = Path(base_path)
        base_path.mkdir(parents=True, exist_ok=True)

        # Collect all series and scalar data
        # We will dump them at the end
        series_data = dict()
        scalar_numbers = dict()

        for name, member in vars(self).items():
            if name.startswith('_'):
                continue  # skip private helper attributes (indices, caches)
            if isinstance(member, pd.DataFrame):
                fname = name + ".csv"
                member.to_csv(base_path / fname)
            elif isinstance(member, pd.Series):
                series_data[name] = member
            elif isinstance(member, numbers.Number):
                scalar_numbers[name] = member
            elif name == "Initialization":
                # Special case - Initialization dictionary
                # Serialise it to a subfolder
                path = base_path / name
                export_initialization_dict_to_csv(path, member, gridsize=self.gridsize)
            else:
                warnings.warn(
                    f"Output member '{name}' has unsupported type '{type(member)}'. "
                    f"It has not been exported to the output directory '{base_path}'."
                )

        for name, series in series_data.items():
            try:
                series.to_csv(base_path / f"{name}.csv", header=[name])
            except Exception as e:
                warnings.warn(f"Could not export series '{name}': {e}")

        # Print numbers
        pd.Series(scalar_numbers).to_csv(base_path / "scalars.csv")

    def export_to_netcdf(self, base_path: Path | str) -> None:
        """Export contents of the output file to a directory in NetCDF format.
        - Each pandas.DataFrame member is saved to a separate .nc file.
        - All pandas.Series members are combined and saved to a single 'series.nc' file.
        - All scalar numerical members are grouped and saved to 'scalars.nc'.

        Parameters:
          base_path : Path
            A path that names the root directory where contents will be exported.
            If the directory does not exist it will be created.
        """
        self._assemble_series()
        # Create space for output
        base_path = Path(base_path)
        base_path.mkdir(parents=True, exist_ok=True)

        # Collect all series and scalar data
        series_data = dict()
        scalar_numbers = dict()

        for name, member in vars(self).items():
            if name.startswith('_'):
                continue  # skip private helper attributes (indices, caches)
            # convert each DataFrame to an xarray Dataset and save to .nc
            if isinstance(member, pd.DataFrame):
                # Ensure column names are strings
                member.columns = member.columns.astype(str)
                if member.index.name is not None:
                    member.index.name = str(member.index.name)
                fname = name + ".nc" # use the .nc extension
                try:
                    xarray_member = xr.Dataset.from_dataframe(member)
                    encoding = {var: {"zlib": True, "complevel": 4} for var in xarray_member.data_vars}
                    xarray_member.to_netcdf(base_path / fname, encoding=encoding, format="NETCDF4")
                except Exception as e:
                    warnings.warn(
                        f"Could not export DataFrame '{name}' to NetCDF. Error: {e}"
                    )

            elif isinstance(member, pd.Series):
                series_data[name] = member

            elif isinstance(member, numbers.Number):
                scalar_numbers[name] = member

            elif name == "Initialization":
                # Special case - Initialization dictionary
                # Serialise it to a subfolder
                path = base_path / name
                export_initialization_dict_to_netcdf(path, member, gridsize=self.gridsize)

            else:
                warnings.warn(
                    f"Output member '{name}' has unsupported type '{type(member)}'. "
                    f"It has not been exported to the output directory '{base_path}'."
                )

        # process and save Series — each written to its own .nc file
        for name, series in series_data.items():
            try:
                if isinstance(series.index, pd.MultiIndex):
                    da = xr.DataArray.from_series(series)
                else:
                    da = xr.DataArray(series.values, dims=["index"], name=name)
                    da = da.assign_coords(index=series.index.values)
                da.name = name
                da.to_netcdf(base_path / f"{name}.nc",
                             encoding={name: {"zlib": True, "complevel": 4}},
                             format="NETCDF4")
            except Exception as e:
                warnings.warn(f"Could not export Series '{name}' to NetCDF. Error: {e}")


        if scalar_numbers:
            # Create an xarray Dataset directly from the dictionary of scalars.
            # Each key will become a variable in the NetCDF file.
            scalars_dataset = xr.Dataset(scalar_numbers)
            # Save the scalars Dataset to a NetCDF file.
            encoding = {var: {"zlib": True, "complevel": 4} for var in scalars_dataset.data_vars}
            scalars_dataset.to_netcdf(base_path / "scalars.nc",
                                      encoding=encoding, format="NETCDF4")

    def export_to_zarr(self, base_path: Path | str) -> None:
        """Export contents of the output object to a single Zarr store.
        - Each DataFrame is saved as a separate Zarr group.
        - Each Series is saved as a separate Zarr group.
        - Scalar numbers are saved together in a 'scalars' group.
        - Initialization dict is saved in an 'Initialization' sub-store.

        Parameters:
          base_path : Path
            Path for the output Zarr store directory (e.g. 'run_001.zarr').
            If it exists it will be overwritten.
        """
        self._assemble_series()
        import zarr
        from zarr.codecs import ZstdCodec
        compressor = ZstdCodec(level=3)

        base_path = Path(base_path)
        if base_path.suffix != ".zarr":
            base_path = base_path.with_suffix(".zarr")

        series_data = dict()
        scalar_numbers = dict()

        for name, member in vars(self).items():
            if name.startswith('_'):
                continue  # skip private helper attributes (indices, caches)
            if isinstance(member, pd.DataFrame):
                member.columns = member.columns.astype(str)
                if member.index.name is not None:
                    member.index.name = str(member.index.name)
                try:
                    if isinstance(member.index, pd.MultiIndex):
                        member = member.copy().reset_index()
                    ds = xr.Dataset.from_dataframe(member)
                    encoding = {var: {"compressors": compressor} for var in ds.data_vars}
                    ds.to_zarr(str(base_path), group=name, encoding=encoding, mode="a")
                except Exception as e:
                    warnings.warn(f"Could not export DataFrame '{name}' to Zarr. Error: {e}")

            elif isinstance(member, pd.Series):
                series_data[name] = member

            elif isinstance(member, numbers.Number):
                scalar_numbers[name] = member

            elif name == "Initialization":
                init_store = str(base_path) + ".init.zarr"
                export_initialization_dict_to_zarr(init_store, member, gridsize=self.gridsize)

            else:
                warnings.warn(
                    f"Output member '{name}' has unsupported type '{type(member)}'. "
                    f"It has not been exported to '{base_path}'."
                )

        # Save each Series as its own group
        for name, series in series_data.items():
            try:
                if isinstance(series.index, pd.MultiIndex):
                    da = xr.DataArray.from_series(series)
                else:
                    da = xr.DataArray(series.values, dims=["index"], name=name)
                    da = da.assign_coords(index=series.index.values)
                da.name = name
                ds = da.to_dataset()
                encoding = {name: {"compressors": compressor}}
                ds.to_zarr(str(base_path), group=name, encoding=encoding, mode="a")
            except Exception as e:
                warnings.warn(f"Could not export Series '{name}' to Zarr. Error: {e}")

        # Save scalars as a single group
        if scalar_numbers:
            scalars_dataset = xr.Dataset(scalar_numbers)
            scalars_dataset.to_zarr(str(base_path), group="scalars", mode="a")
