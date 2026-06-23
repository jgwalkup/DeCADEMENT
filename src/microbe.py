"""
This microbe.py module has one class and two functions.
    
    Microbe():                class
    microbe_osmo_psi():       function; moisture modifier of inducible osmolyte production efficiency
    microbe_mortality_prob(): function; cell mortality probability
"""

import numpy as np
import pandas as pd
from utility import LHS
# import scipy.constants as const

class Microbe():
    """
    This class holds all variables related to microbes.

    Methods involving composition,stoichiometry, enzyme and gene production, as well as responses to environmental factors.
    These methods include:
        1) microbial_community_initialization(): initialize microbial community on the spatial grid
        2) minimum_cell_quota():          get minimum ratios
        3) microbe_enzyme_gene():         derive taxon-specific genes for enzyme
        4) microbe_osmolyte_gene():       derive genes encoding for osmolytes
        5) microbe_uptake_gene():         derive transporter genes
        6) microbe_uptake_cost():         metabolic cost of producing transporter
        7) microbe_enzproduction_rate():  cost of producing enzymes
        8) microbe_osmoproduction_rate(): cost of producing osmolytes
        9) microbe_drought_tol():         microbial drougth tolerance
       10) microbe_mortality():           paramters pertaining to microbial mortality
    """
    
    def __init__(self,runtime,parameters):
        """
        The constructor of Microbe class.

        Parameters:
            runtime:     user-specified runtime parameters
            parameters: input dataframe of all parameters
        """
        
        # system setup parameters
        self.gridsize    = int(runtime.loc['gridsize',1])
        self.n_taxa      = int(runtime.loc['n_taxa',1])
        self.n_enzymes   = int(runtime.loc['n_enzymes',1])
        self.n_substrates= int(runtime.loc['n_substrates',1])
        self.n_monomers  = int(runtime.loc['n_substrates',1])+3  # +3 b/c of two inorganic monomers & ferm. exudates
        self.n_uptake    = int(runtime.loc['n_uptake',1])        # number of transporters for each taxon
        self.n_osmolyte  = int(runtime.loc['n_osmolytes',1])     # system-allowed number of osmotic compound
        self.n_HSP       = int(runtime.loc['n_HSP',1])          # number of heat stress proteins genes
        self.taxa_per_box = runtime.loc['taxa_per_box', 1]      # Probability of each taxon entering a grid cell
        # index of fungal taxa in a microbial pool (1);1-d array
        self.fb = np.random.choice([1,0], self.n_taxa, replace=True, p=[runtime.loc['fb',1],(1-runtime.loc['fb',1])]).astype('int32')
        self.index  =  ["Tax" + str(i) for i in range(1,self.n_taxa+1)]    # Taxon index

        # microbial cell size
        self.Cfrac_b    = parameters.loc['Cfrac_b',1]     # Bacterial C fraction: 0.825 mg mg-1
        self.Nfrac_b    = parameters.loc['Nfrac_b',1]     # Bacterial N fraction: 0.16  mg mg-1
        self.Pfrac_b    = parameters.loc['Pfrac_b',1]     # Bacterial P fraction: 0.015 mg mg-1
        self.Erat_b     = parameters.loc['Erat_b', 1]         # Bacterial energy:biomass C ratio:   5000 ATP mg-1
        self.Cfrac_f    = parameters.loc['Cfrac_f',1]     # Fungal C fraction: 0.9	mg mg-1
        self.Nfrac_f    = parameters.loc['Nfrac_f',1]     # Fungal N fraction: 0.09	mg mg-1
        self.Pfrac_f    = parameters.loc['Pfrac_f',1]     # Fungal P fraction: 0.01	mg mg-1
        self.Erat_f     = parameters.loc['Erat_f',1]      # Fungal energy:biomass C ratio:   5000 ATP mg-1
        self.max_size_b = parameters.loc['max_size_b',1]  # C quota threshold for bacterial cell division: 2 mg cm-3
        self.max_size_f = parameters.loc['max_size_f',1]  # C quota threshold for fungal cell division: 50 mg cm-3
        self.Crange     = parameters.loc['Crange',1]      # Tolerance on C fraction: 0.09 mg mg-1
        self.Nrange     = parameters.loc['Nrange',1]      # Tolerance on N fraction: 0.04 mg mg-1
        self.Prange     = parameters.loc['Prange',1]      # Tolerance on P fraction: 0.005 mg mg-1
        self.Erange     = parameters.loc['Erange',1]      # tolerance on energy:biomass C ratio

        # transporter
        self.Uptake_C_cost_min = parameters.loc['Uptake_C_cost_min',1]      # Minimum C cost of one transporter gene being encoded 
        self.Uptake_C_cost_max = parameters.loc['Uptake_C_cost_max',1]      # Maximum C cost of one transporter gene being encoded
        self.NormalizeUptake   = parameters.loc['NormalizeUptake',1]        # Normalize uptake investment for the number of uptake genes;default:0
        # enzyme
        self.Enz_per_taxon_min = int(parameters.loc['Enz_per_taxon_min',1]) # =0;Minimum number of enzyme genes a taxon can hold
        self.Enz_per_taxon_max = int(parameters.loc['Enz_per_taxon_max',1]) # =40;Maximum number of enzyme genes a taxon can hold
        self.Constit_Prod_min  = parameters.loc['Constit_Prod_min',1]       # =0.00001;Minimum per enzyme production cost as a fraction of biomass
        self.Constit_Prod_max  = parameters.loc['Constit_Prod_max',1]       # =0.0001; Maximum per enzyme production cost as a fraction of biomass
        self.Enz_Prod_min      = parameters.loc['Enz_Prod_min',1]           # =0.00001;Minimum per enzyme production cost as a fraction of C uptake rate
        self.Enz_Prod_max      = parameters.loc['Enz_Prod_max',1]           # =0.0001; Maximum ...
        self.NormalizeProd     = parameters.loc['NormalizeProd',1]          # Normalize enzyme production for the number of enzyme genes;default:0
        # osmolyte
        self.Osmo_per_taxon_min = int(parameters.loc['Osmo_per_taxon_min',1]) # Minimum number of osmotic gene 
        self.Osmo_per_taxon_max = int(parameters.loc['Osmo_per_taxon_max',1]) # Max. of osmotic gene
        self.Osmo_Consti_Prod_min = parameters.loc['Osmo_Consti_Prod_min',1]  # constitutive cost min
        self.Osmo_Consti_Prod_max = parameters.loc['Osmo_Consti_Prod_max',1]  # constitutive cost max
        self.Osmo_Induci_Prod_min = parameters.loc['Osmo_Induci_Prod_min',1]  # inducible cost min
        self.Osmo_Induci_Prod_max = parameters.loc['Osmo_Induci_Prod_max',1]  # inducible cost max
        # heat stress proteins
        self.HSP_per_taxon_min = int(parameters.loc['HSP_per_taxon_min',1]) # Minimum number of heat stress protein gene
        self.HSP_per_taxon_max = int(parameters.loc['HSP_per_taxon_max',1]) # Max. of HSP genes
        self.HSP_Consti_Prod_min = parameters.loc['HSP_Consti_Prod_min',1]  # constitutive cost min
        self.HSP_Consti_Prod_max = parameters.loc['HSP_Consti_Prod_max',1]  # constitutive cost max
        self.HSP_Induci_Prod_min = parameters.loc['HSP_Induci_Prod_min',1]  # inducible cost min
        self.HSP_Induci_Prod_max = parameters.loc['HSP_Induci_Prod_max',1]  # inducible cost max

        # mortality
        self.death_rate_bac = parameters.loc['death_rate_bac',1]              # Bacterial basal mortality prob.
        self.death_rate_fun = parameters.loc['death_rate_fun',1]              # Fungal basal mortality probability
        self.beta_bac       = parameters.loc['beta_bac',1]                    # Bacterial mortality coefficient
        self.beta_fun       = parameters.loc['beta_fun',1]                    # Fungal mortality coefficient
        # thermal response
        #self.TR_mag_Cp_Slope = parameters.loc['TR_mag_Cp_slope',1]         # Slope for calculating Cp of maintenance respiration to adjust temperature sensitivity by thermal adaptation magnitude

    def microbial_pool_initialization(self):
        """
        Initialize a microbial pool of bacteria and/or fungi.

        Calculation procedure:
            Create a dataframe with only bacteria
            Substitute part of bacteria with fungi
        Parameters:
            Cfrac_b: Bacterial C fraction:    0.825 mg mg-1
            Nfrac_b: Bacterial N fraction:    0.160 mg mg-1
            Pfrac_b: Bacterial P fraction:    0.015 mg mg-1
            Erat_b:  Bacterial E:MBC ratio:   5000 ATP mg-1
            Cfrac_f: Fungal C fraction:       0.90  mg mg-1
            Nfrac_f: Fungal N fraction:       0.09  mg mg-1
            Pfrac_f: Fungal P fraction:       0.01  mg mg-1
            Erat_f:  Fungal E:MBC ratio:      5000 ATP mg-1
            max_size_b: C quota threshold for bacterial cell division
            max_size_f: C quota threshold for fungal cell division
        Returns:
            microbes_pp: dataframe; (taxa*gridsize)*3; pool of microbes
            fb_grid:     1D array;  index of fungi taxa
        """

        BacC = 0.5 * self.max_size_b
        BacN = BacC * self.Nfrac_b / self.Cfrac_b
        BacP = BacC * self.Pfrac_b / self.Cfrac_b
        BacE = BacC * self.Erat_b  # energy pool relative to biomass C
        microbes_array = np.tile([BacC, BacN, BacP, BacE], (self.n_taxa * self.gridsize, 1))
        microbes_pp = pd.DataFrame(data=microbes_array, index=self.index * self.gridsize, columns=["C", "N", "P", "E"],
                                   dtype='float32')

        # Derive the Fungi index by expanding the fb to the spatial grid
        fb_grid = np.tile(self.fb, self.gridsize)

        # Fungal pool sizes for all elements
        FunC = 0.5 * self.max_size_f
        FunN = FunC * self.Nfrac_f / self.Cfrac_f
        FunP = FunC * self.Pfrac_f / self.Cfrac_f
        FunE = FunC * self.Erat_f

        # Substitute with fungal pools of elements: fb_grid == 1
        microbes_pp.loc[fb_grid == 1, "C"] = FunC
        microbes_pp.loc[fb_grid == 1, "N"] = FunN
        microbes_pp.loc[fb_grid == 1, "P"] = FunP
        microbes_pp.loc[fb_grid == 1, "E"] = FunE

        # Derive the number of fungi and bacteria taxa
        Bac_index = microbes_pp['C'] == BacC
        Fun_index = microbes_pp['C'] == FunC
        Bac_taxa = int(sum(Bac_index) / self.gridsize)
        Fun_taxa = int(sum(Fun_index) / self.gridsize)
        print('Before placement--', 'Bac_taxa=', Bac_taxa, 'Fun_taxa=', Fun_taxa)

        return microbes_pp, fb_grid

    def microbial_community_initialization(self, microbes_pp, fb_grid):
        """
        Initialize a spatially explicit microbial community based on microbial pool.

        Randomly place the microbail pool on the sptatial grid and perform stats on the initialized community
            and output them for record.

        Parameters:
            microbes_pp: microbial pool
            fb_grid:     index of fungal taxa
        Returns:
            microbes_df: df; (n_taxa*gridsize)*4; initialized spatially explicit community
            Bac_density: scalar; bacterial density over the grid
            Fun_density: scalar; fungal density over the grid
        """

        # Export the microbial community preceding placement on the grid
        microbes_df = microbes_pp.copy(deep=True)

        # Randomly place the microbial community created above on the spatial grid to initialize a spatially explicit microbial community
        pb = self.taxa_per_box
        choose_taxa = np.random.choice([1, 0], self.n_taxa * self.gridsize, replace=True, p=[pb, (1 - pb)])
        pf = pb * self.max_size_b / self.max_size_f
        #choose_taxa[fb_grid == 1] = np.random.choice([1, 0], sum(fb_grid), replace=True, p=[pf, (1 - pf)])
        choose_taxa[fb_grid == 1] = np.random.choice([1, 0], len(choose_taxa[fb_grid == 1]), replace=True, p=[pf, (1 - pf)])
        microbes_df.loc[choose_taxa == 0] = np.float32(0)

        BacC = 0.5 * self.max_size_b
        FunC = 0.5 * self.max_size_f
        # Derive the number of fungi and bacteria taxa
        Bac_index = microbes_df['C'] == BacC
        Fun_index = microbes_df['C'] == FunC
        Bac_taxa = microbes_df[Bac_index].groupby(level=0).sum().shape[0]
        Fun_taxa = microbes_df[Fun_index].groupby(level=0).sum().shape[0]
        print('After placement--', 'Bac_taxa=', Bac_taxa, 'Fun_taxa=', Fun_taxa)

        Bac_density = sum(Bac_index) * BacC / self.gridsize
        Fun_density = sum(Fun_index) * FunC / self.gridsize
        print('After placement--', 'Bac_density=', Bac_density, 'Fun_density=', Fun_density)

        return microbes_df, Bac_density, Fun_density

    def minimum_cell_quota(self):
        """
        This will be used in the mortality calculation.
        
        Parameters:
            Cfrac_b: Bacterial C fraction:    0.825 mg mg-1
            Nfrac_b: Bacterial N fraction:    0.160 mg mg-1
            Pfrac_b: Bacterial P fraction:    0.015 mg mg-1
            Erat_b:  Bacterial E:MBC ratio:   5000 ATP mg-1
            Cfrac_f: Fungal C fraction:       0.90  mg mg-1
            Nfrac_f: Fungal N fraction:       0.09  mg mg-1
            Pfrac_f: Fungal P fraction:       0.01  mg mg-1
            Erat_f:  Fungal E:MBC ratio:      5000 ATP mg-1
            Crange:  Tolerance on C fraction: 0.090 mg mg-1
            Nrange:  Tolerance on N fraction: 0.040 mg mg-1
            Prange:  Tolerance on P fraction: 0.005 mg mg-1
            Erange:  Tolerance on E:MBC ratio
        Return:
            MinRatios: dataframe   
        """

        # First derive the optimal stoichiometry of bacterial taxa
        OptimalRatios_array = np.tile([self.Cfrac_b,self.Nfrac_b,self.Pfrac_b,self.Erat_b],(self.n_taxa,1))
        # index = ["Tax" + str(i) for i in range(1,self.n_taxa + 1)]
        OptimalRatios_df = pd.DataFrame(data=OptimalRatios_array, index=self.index, columns=["C","N","P","E"],dtype='float32')
        # Then substitute with fungal stoichiometry
        OptimalRatios_df[self.fb==1] = np.tile([self.Cfrac_f,self.Nfrac_f,self.Pfrac_f,self.Erat_f],(sum(self.fb),1))
        # Calcualte ratio range
        RangeRatios_array = np.tile([self.Crange,self.Nrange,self.Prange,self.Erange],(self.n_taxa,1))
        RangeRatios_df    = pd.DataFrame(data=RangeRatios_array, index=self.index, columns=["C","N","P","E"], dtype='float32')
        # Calculate minimum cell quotas regarding C, N, & P
        MinRatios = OptimalRatios_df - RangeRatios_df
        
        return MinRatios
    

    def microbe_enzyme_gene(self):
        """ 
        Derive taxon-specific enzyme genes.
        
        This method draws the number of genes per taxon from a uniform distribution; Not all taxa need to have an enzyme gene.

        Parameters:
            Enz_per_taxon_min: minimum # genes a taxon can have (0)
            Enz_per_taxon_max: maximum # genes a taxon can have (40)   
        Return:
            EnzGenes: Rows are taxa; cols are genes;values: 0/1
        """

        # gene pool producing the number of enzymes that system requires 
        n_genes = self.n_enzymes

        # taxon-specific number of genes
        genes_per_taxon = np.random.choice(range(self.Enz_per_taxon_min,self.Enz_per_taxon_max+1),self.n_taxa,replace=True)
        
        # randomly assign genes from the gene pool to each individual taxon
        def assign_gene(i_taxon):
            """Randomly assign a specific number of different genes to a taxon."""

            probability_list = [0]*n_genes
            probability_list[0:genes_per_taxon[i_taxon]] = [1] * genes_per_taxon[i_taxon]
            taxon = np.random.choice(probability_list,n_genes,replace=False)
            return taxon
        
        EnzGenes_list = [assign_gene(i) for i in range(self.n_taxa)]    # list of 1D array
        columns       = ['Enz' + str(i) for i in range(1,n_genes+1)]
        EnzGenes      = pd.DataFrame(data=np.vstack(EnzGenes_list), index=self.index, columns=columns, dtype='int8')
        
        return EnzGenes

    def microbe_osmolyte_gene(self):
        """
        Derive the taxon-specific number of genes encoding osmolytes.

        This method draws the number of genes per taxon from a uniform distribution; every taxon must have an osmotic gene.
        
        Parameters:
            n_osmolyte:         total number of different genes allowed by a system
            Osmo_per_taxon_min: min of osmotic genes in a taxon
            Osmo_per_taxon_max: max of osmotic genes in a taxon
        Return:
            OsmoGenes: dataframe; row: taxon; col: genes
        """

        # gene pools producing the system-required number of osmolytes
        n_genes = self.n_osmolyte

        # derive the number of osmolyte genes each taxon can have randomly
        genes_per_taxon = np.random.choice(range(self.Osmo_per_taxon_min,self.Osmo_per_taxon_max+1),self.n_taxa,replace=True)
        
        # randomly determine different genes within a taxon
        def assign_gene(i_taxon):
            """Randomly assign a specific number of different genes to a taxon"""

            probability_list = [0]*n_genes
            probability_list[0:genes_per_taxon[i_taxon]] = [1] * genes_per_taxon[i_taxon]
            taxon = np.random.choice(probability_list, n_genes, replace=False)
            return taxon
        
        OsmoGenes_list = [assign_gene(i) for i in range(self.n_taxa)]
        # index          = ["Tax" + str(i) for i in range(1,self.n_taxa+1)]
        columns        = ['Osmo'+ str(i) for i in range(1,n_genes+1)]
        OsmoGenes      = pd.DataFrame(data=np.vstack(OsmoGenes_list), index=self.index, columns=columns, dtype='int8')
        
        return OsmoGenes


    def microbe_heat_stress_gene(self):
        """
        Derive the taxon-specific number of genes encoding heat stress proteins (HSP).
        This method draws the number of genes per taxon from a uniform distribution; every taxon must have an HSPtic gene.

        Parameters:
            n_HSP:         total number of different genes allowed by a system
            HSP_per_taxon_min: min of HSP genes in a taxon
            HSP_per_taxon_max: max of HSP genes in a taxon
        Return:
            HSPGenes: dataframe; row: taxon; col: genes
        """

        # gene pools producing the system-required number of HSPs
        n_genes = self.n_HSP

        # derive the number of HSP genes each taxon can have randomly
        genes_per_taxon = np.random.choice(range(self.HSP_per_taxon_min, self.HSP_per_taxon_max + 1), self.n_taxa,
                                           replace=True)

        # randomly determine different genes within a taxon
        def assign_gene(i_taxon):
            """Randomly assign a specific number of different genes to a taxon"""

            probability_list = [0] * n_genes
            probability_list[0:genes_per_taxon[i_taxon]] = [1] * genes_per_taxon[i_taxon]
            taxon = np.random.choice(probability_list, n_genes, replace=False)
            return taxon

        HSPGenes_list = [assign_gene(i) for i in range(self.n_taxa)]
        # index          = ["Tax" + str(i) for i in range(1,self.n_taxa+1)]
        columns = ['HSP' + str(i) for i in range(1, n_genes + 1)]
        HSPGenes = pd.DataFrame(data=np.vstack(HSPGenes_list), index=self.index, columns=columns, dtype='int8')

        return HSPGenes

     
    def microbe_uptake_gene(self,ReqEnz,EnzGenes,MonomersProduced):
        """
        Derive the genes encoding uptake enzymes(transporters).
        
        All taxa must have at least one uptake gene.
        
        Parameters:
            ReqEnz:           3D df (2 * n_substrates * n_enzymes); substrate-required enzymes; from the Substrate class
            EnzGenes:         2D df(taxon*gene); genes encoding enzymes degrading substrates; from the above 'microbe_enzyme_gene' method
            MonomersProduced: dataframe; sub * monomers; from the Substrate class
        Return:
            UptakeGenes: dataframe:Rows-taxa; Cols-genes; Values: 0 (no gene) or 1 (yes)
        """
        
        # substrate-required enzymes
        Sub_enz = ReqEnz.loc['set1'].iloc[0:self.n_substrates,:] + ReqEnz.loc['set2'].iloc[0:self.n_substrates,:]
        # Matrix multiplication to relate taxa to the monomers they can generate with their enzymes
        UptakeGenes = EnzGenes @ Sub_enz.T @ MonomersProduced
        UptakeGenes.iloc[:,0:3] = np.int8(1)
        UptakeGenes[lambda df: df>0] = np.int8(1)  # == UptakeGenes[UptakeGenes > 0] = 1

        # Ensure every taxon is likely to have an uptake enzyme for at least 1 organic monomer
        # Not guaranteed unless uptake_prob = 1 (refer to the R version)
        probability_list    = [0]*(self.n_monomers-3)
        probability_list[0] = 1
        for i in range(self.n_taxa):
            if sum(UptakeGenes.iloc[i,3:self.n_monomers]) == 0:
                UptakeGenes.iloc[i,3:self.n_monomers] = np.random.choice(probability_list,self.n_monomers-3,replace=False).astype('int8')

        # Give each taxon a random number of additional uptake genes between the number they have and n_upgenes
        for i in range(self.n_taxa):
            n_zero = sum(UptakeGenes.iloc[i,:][UptakeGenes.iloc[i,:]==0])
            if n_zero == 0: # has all genes
                continue
            probability_list = [0]*n_zero
            locator = np.random.choice(range(1, n_zero+1), 1, replace=True) # derive n_zero >= locater >= 1 (NOTE)
            probability_list[0:int(locator)] = [1]*int(locator)
            UptakeGenes.iloc[i,:][UptakeGenes.iloc[i,:]==0] = np.random.choice(probability_list,n_zero,replace=False).astype('int8')
           
        return UptakeGenes
    
    
    def microbe_uptake_cost(self,UptakeGenes):                         
        """
        Derive the taxon-specific cost of every single gene of uptake transporter.
        
        Note this cost (in terms of Fraction of Biomass C) is same among different genes from the same taxon
        
        Parameters:
            UptakeGenes:       dataframe; taxon-specific transporter genes; from the above method, microbe_uptake_gene()
            NormalizeUptake:   Normalize uptake investment for the number of uptake genes (0: No, 1: Yes)
            Uptake_C_cost_min: 0.001	transporter mg-1 biomass C
            Uptake_C_cost_max: 0.01	transporter mg-1 biomass C
        Returns:
            UptakeProd_series: series; taxon-specific transporter cost
            UptakeGenes_Cost:  dataframe; taxon- and gene-specific transporter cost
        """                         
        
        # LHS sampling of transporter production efficiency for each taxon
        #index = ["Tax" + str(i) for i in range(1,self.n_taxa+1)]
        UptakeProd_array  = LHS(self.n_taxa, self.Uptake_C_cost_min, self.Uptake_C_cost_max, 'uniform')
        UptakeProd_series = pd.Series(data=UptakeProd_array, index=self.index, dtype='float32')

        # If true (== 1), then the uptake potential is normalized to the number of uptake genes (n_uptake)
        if self.NormalizeUptake == 1:
            Normalize   = UptakeGenes.sum(axis=1)/self.n_uptake
            UptakeGenes = UptakeGenes.divide(Normalize,axis=0)
        
        # Derive gene-specific production efficiency of each taxon
        UptakeGenes_Cost = UptakeGenes.mul(UptakeProd_series,axis=0)
        
        return UptakeProd_series, UptakeGenes_Cost
    
    
    def microbe_enzproduction_rate(self,EnzGenes,EnzAttrib):
        """
        Derive the taxon-specific fraction of 'available C' as enzymes: note that.

        This fraction only varies with taxon, which is independent of gene within a taxon
        
        Parameters:
            EnzGenes:         taxon-specific available genes for enzyme;dataframe (taxon * genes); from the above microbe_enzyme_gene() method
            EnzAttrib:        dataframe [enzyme * 4 (C,N,P,maintenence cost)]; basic enzyme stoichiometry
            Constit_Prod_min: min of constitutive enzyme production efficiency
            Constit_Prod_max: max of constitutive enzyme production efficiency
            Enz_Prod_min:     min of inducible enzyme production efficiency
            Enz_Prod_max:     max of inducible enzyme production efficiency
            NormalizeProd:    scalar; normalize enzyme production for the number of enzyme genes (0: No, 1: Yes)  
        Returns:
            Tax_EnzProdConsti:   series; taxon-specific fraction of available C as enzyme for each taxon
            Tax_EnzProdInduce:   series; taxon-specific fraction of available C as enzyme for each taxon
            Tax_Induce_Enzyme_C: dataframe; taxon-specific fraction of available C as enzyme for each enzyme
            Tax_Consti_Enzyme_C: dataframe; taxon-specific fraction of available C as enzyme for each enzyme
        """
        
        # LHS sampling of constitutive and inducible enzyme production efficiency for each taxon
        # index = ["Tax" + str(i) for i in range(1,self.n_taxa+1)]
        Tax_EnzProdConsti_array = LHS(self.n_taxa, self.Constit_Prod_min, self.Constit_Prod_max, 'uniform')
        Tax_EnzProdInduce_array = LHS(self.n_taxa, self.Enz_Prod_min,     self.Enz_Prod_max,     'uniform')
        Tax_EnzProdConsti       = pd.Series(data=Tax_EnzProdConsti_array, index=self.index, dtype='float32')
        Tax_EnzProdInduce       = pd.Series(data=Tax_EnzProdInduce_array, index=self.index, dtype='float32')
        
        # Derive the production rate of every single gene of each taxon
        EnzProdConsti = EnzGenes.mul(Tax_EnzProdConsti, axis=0)
        EnzProdInduce = EnzGenes.mul(Tax_EnzProdInduce, axis=0)

        # Normalize constituitive and inducible enzyme production
        if self.NormalizeProd == 1:
            Normalize     = EnzGenes.sum(axis=1)/self.Enz_per_taxon_max
            EnzProdConsti = EnzProdConsti.divide(Normalize, axis=0)
            EnzProdInduce = EnzProdInduce.divide(Normalize, axis=0)
            EnzProdConsti[Normalize==0] = 0
            EnzProdInduce[Normalize==0] = 0

        # Relative enzyme carbon cost for each enzyme gene of each taxon
        Tax_Consti_Enzyme_C = EnzProdConsti.mul(EnzAttrib["C_cost"], axis=1) #C_cost = 1; so it doesn't matter
        Tax_Induce_Enzyme_C = EnzProdInduce.mul(EnzAttrib["C_cost"], axis=1)  
    
        return Tax_EnzProdConsti, Tax_EnzProdInduce, Tax_Consti_Enzyme_C, Tax_Induce_Enzyme_C 
    
    
    def microbe_osmoproduction_rate(self,OsmoGenes):
        """
        Distribution of osmolyte production rate (i.e.,proportion of available C as osmolytes).

        Parameters:
            OsmoGenes:             dataframe(taxon*genes); generated by the above microbe_osmolyte_gene() method
            Osmo_Constit_Prod_min: min of constitutive osmolyte production efficiency
            Osmo_Constit_Prod_max: max of constitutive osmolyte production efficiency
            Osmo_Induci_Prod_min:  min of inducible osmolyte production efficiency
            Osmo_Induci_Prod_max:  max of inducible osmolyte production efficiency
            NormalizeProd:         0   
        Returns:
            Tax_OsmoProd_Consti_series:
            Tax_OsmoProd_Induci_series:
            Tax_Consti_Osmo_C: taxon-specific fraction of available C as osmolytes for each gene
            Tax_Induce_Osmo_C: taxon-specific fraction of available C as osmolytes for each gene
        """
        # LHS sampling of osmolyte production efficiency for every taxon
        Tax_OsmoProd_Consti = LHS(self.n_taxa, self.Osmo_Consti_Prod_min, self.Osmo_Consti_Prod_max, 'uniform')
        Tax_OsmoProd_Induci = LHS(self.n_taxa, self.Osmo_Induci_Prod_min, self.Osmo_Induci_Prod_max, 'uniform')
        
        #index = ["Tax" + str(i) for i in range(1,self.n_taxa+1)]
        Tax_OsmoProd_Consti_series = pd.Series(data=Tax_OsmoProd_Consti, index=self.index, dtype='float32')
        Tax_OsmoProd_Induci_series = pd.Series(data=Tax_OsmoProd_Induci, index=self.index, dtype='float32')
                                              
        # Derive the production rate of every single gene of each taxon
        Tax_Consti_Osmo_C = OsmoGenes.mul(Tax_OsmoProd_Consti_series,axis=0)
        Tax_Induci_Osmo_C = OsmoGenes.mul(Tax_OsmoProd_Induci_series,axis=0)
        
        return Tax_OsmoProd_Consti_series, Tax_OsmoProd_Induci_series, Tax_Consti_Osmo_C, Tax_Induci_Osmo_C

    def microbe_heat_stress_protein_production_rate(self, HSPGenes):
        """
        Distribution of heat stress protein production rate (i.e.,proportion of available C as HSP).

        Parameters:
            HSPGenes:             dataframe(taxon*genes); generated by the above microbe_heat_stress_gene() method
            HSP_Constit_Prod_min: min of constitutive HSP production efficiency
            HSP_Constit_Prod_max: max of constitutive HSP production efficiency
            HSP_Induci_Prod_min:  min of inducible HSP production efficiency
            HSP_Induci_Prod_max:  max of inducible HSP production efficiency
            NormalizeProd:         0
        Returns:
            Tax_HSP_Prod_Consti_series:
            Tax_HSP_Prod_Induci_series:
            Tax_Consti_HSP_C: taxon-specific fraction of available C as HSP for each gene
            Tax_Induce_HSP_C: taxon-specific fraction of available C as HSP for each gene

        """
        # LHS sampling of HSP production efficiency for every taxon
        Tax_HSP_Prod_Consti = LHS(self.n_taxa, self.HSP_Consti_Prod_min, self.HSP_Consti_Prod_max, 'uniform')
        Tax_HSP_Prod_Induci = LHS(self.n_taxa, self.HSP_Induci_Prod_min, self.HSP_Induci_Prod_max, 'uniform')

        # index = ["Tax" + str(i) for i in range(1,self.n_taxa+1)]
        Tax_HSP_Prod_Consti_series = pd.Series(data=Tax_HSP_Prod_Consti, index=self.index, dtype='float32')
        Tax_HSP_Prod_Induci_series = pd.Series(data=Tax_HSP_Prod_Induci, index=self.index, dtype='float32')

        # Derive the production rate of every single gene of each taxon
        Tax_Consti_HSP_C = HSPGenes.mul(Tax_HSP_Prod_Consti_series, axis=0)
        Tax_Induci_HSP_C = HSPGenes.mul(Tax_HSP_Prod_Induci_series, axis=0)

        return Tax_HSP_Prod_Consti_series, Tax_HSP_Prod_Induci_series, Tax_Consti_HSP_C, Tax_Induci_HSP_C

    def HSP_initialization(self, microbes_df, HSPGenes):
        """
        Initialize the per-gene HSP mass for each taxon on the spatial grid.

        Parameters:
            microbes_df: DataFrame; the initialized microbial community (gridsize*n_taxa rows,
                         columns C/N/P/E), used to set initial Total_biomass.
            HSPGenes:    DataFrame; taxon * HSP genes (0/1), NOT yet expanded to the grid.
                         Columns are ['HSP1', 'HSP2', ...].
        Returns:
            HSP_init:    DataFrame; rows = microbes_df.index (gridsize*n_taxa),
                         columns = HSP gene columns + ['Total_biomass'].
                         All HSP gene columns initialised to 0.
        """
        # Expand gene structure to full grid so rows match Microbes
        from utility import expand
        gridsize = len(microbes_df) // len(HSPGenes)
        HSP_init = expand(HSPGenes.astype('float32'), gridsize).copy()
        HSP_init.index = microbes_df.index

        # Zero out all gene columns — no HSP mass at t=0
        HSP_init[:] = np.float32(0)

        # Add a Total_biomass column for tracking HSP proportion
        HSP_init['Total_biomass'] = microbes_df[['C', 'N', 'P']].sum(axis=1).values

        return HSP_init

    def microbe_drought_tol(self,Tax_Consti_Osmo_C,Tax_Induci_Osmo_C):
        """
        Derive taxon-specific drought tolerance value.
        
        Drought tolerance is postively correlated with taxon-specific inducible osmotic allocation efficiency.
        
        Parameter:
            Tax_Induci_Osmo_C: dataframe
        Return:
            Tax_tolerance: series;float32
        """

        Tax_Osmo_Alloc = Tax_Induci_Osmo_C.sum(axis=1) + Tax_Consti_Osmo_C.sum(axis=1)
        #Tax_tolerance = Tax_Osmo_Alloc.rank(axis=0,method='min')/self.n_taxa
        Tax_tolerance = (Tax_Osmo_Alloc - Tax_Osmo_Alloc.min())/(Tax_Osmo_Alloc.max()-Tax_Osmo_Alloc.min())
        Tax_tolerance = Tax_tolerance.fillna(0)
        
        return Tax_tolerance
    
     
    def microbe_mortality(self,fb_grid):
        """
        Derive taxon-specific microbial mortality parameters.
        
        Parameters:
            fb_grid: 1D array;index of bacteria (0) vs fungi (1) over the grid    
        Returns:
            basal_death_prop: spatial taxon-specific basal death probability; 1D array
            death_rate:       spatial death_rate; 1D array   
        """

        basal_death_prob = ((1-fb_grid) * self.death_rate_bac + fb_grid * self.death_rate_fun).astype('float32')
        death_rate       = ((1-fb_grid) * self.beta_bac + fb_grid * self.beta_fun).astype('float32')

        return basal_death_prob, death_rate 


def microbe_osmo_psi(wp,alfa,wp_fc,wp_th):
    """
    Derive water potential modifier of inducible osmolyte production.

    Inducible production of osmolytes triggered when PSI declines to a **threshold** value,wp_fc,
    below which the production increases and reaches maxima at water potential of wp_th
    ---------------------------------------------------------------------------
    Parameters:
        wp:    scalar;water potential at a daiy step 
        alfa:  scalar;shape factor quantifying curve concavity; could be distinguished btw bacteria and fungi
        wp_fc: scalar;water potential at field capacity
        wp_th: scalar;water potential threshold   
    Returns:
        f_osmo:scalar; modifier of inducible production of osmoylte   
    References:
        Based on Manzoni et al. 2012 Ecology, a synthesis study. 
    """
    
    if wp >= wp_fc:
        f_osmo = 0.0
    #elif wp <= wp_th:
    #    f_osmo = 1.0
    else:
        x = np.log(wp/wp_fc)
        y = np.log(wp_th/wp_fc)
        
        f_osmo = (x/y)**alfa

    return np.float32(f_osmo)


def microbe_mortality_prob(wp,wp_fc,basal_death_prob,death_rate,Tax_tolerance):
    """
    Microbial mortality probability as a function of water potential and drought tolerance.

    Paramters:
        wp:               scalar; water potential
        wp_fc:            scalar; field capacity
        basal_death_prob: array; basal mortality prob. distinguished btw fungi and bacteria
        death_rate:       array; mortality change rate with moisture
        Tax_tolerance:    dataframe; taxon-specific drought tolerance   
    Returns:
        mortality_rate:   taxon-specific mortality probability
    References:
        Allison and Goulden,2017,Soil Biology and Biochemistry
    """
    
    if wp >= wp_fc:
        mortality_rate = basal_death_prob
    else:
        tolerance = Tax_tolerance.to_numpy()
        # option 1
        mortality_rate = basal_death_prob * (1 - (1-tolerance)*(wp-wp_fc)*death_rate)
        # option 2
        #mortality_rate = death_rate * (1/np.exp(tolerance)) * (1 - beta*(wp-wp_fc))
    
    return mortality_rate.astype('float32')


