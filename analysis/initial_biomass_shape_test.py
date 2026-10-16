"""Diagnostic: minimal model (ambient climate, first year) with the starting biomass multiplied by MULT (C, N, P and E
of every cell, after initialisation; cells above max_size divide on the first step). Not a model option.

    python analysis/initial_biomass_shape_test.py MULT VMAX_KM SEED   (run from the repository root)
"""
import sys, types, numpy as np, pandas as pd
sys.modules.setdefault('xarray', types.ModuleType('xarray'))
sys.path.insert(0, 'analysis'); import sweep_minimal as sm
mult, km, seed = float(sys.argv[1]), float(sys.argv[2]), int(sys.argv[3])
Base = sm.Grid
class G2(Base):
    def __init__(self, *a, **k):
        super().__init__(*a, **k)
        self.Microbes[['C', 'N', 'P', 'E']] *= np.float32(mult)   # diagnostic only: larger starting biomass
sm.Grid = G2
s, ser = sm.run_one('input_minimal_e1500_crich_km25_in', {'R_endo': 1, 'Vmax_Km': km}, 20, 365, seed,
                    climate='input/climate_data/climate_ambient_temp.csv')
d = pd.DataFrame(ser).assign(mult=mult, Vmax_Km=km, seed=seed)
d.to_csv(f'initbio_{mult:g}_{km:g}_{seed}.csv', index=False)
