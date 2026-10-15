# input_minimal_e1500

Minimal scenario (one taxon, one substrate IdealSub, one enzyme, one transporter slot; traits pegged) with the
parameter set settled in the sensitivity analyses (see `analysis_output/SENSITIVITY_SUMMARY.md`).
Copy of `input_minimal/` with the changes below. Grid is 20x20 (400 cells), single pulse, `end_time` 365,
flat 15 C climate (replace `climate.csv` for other scenarios), `Kinetics_Model` 1 (MMRT).

| File | Parameter | input_minimal | input_minimal_e1500 |
|---|---|---|---|
| monomer_energy.csv | `Mon_IdealSub` E_resp | 730 | **1500** (E_ferm 230, E_value 1 unchanged) |
| parameters.csv | `Maint_E_Rate_Const` | 2.0 | **2.5** |
| parameters.csv | `Maint_E_HSP_Ea_Intercept` | 30000 | **45000** |
| parameters.csv | `Maint_E_HSP_Ea_Slope` | -30000 | **-90000** |
| parameters.csv | `HSP_prod_E_cost` | 600 | **300** |
| parameters.csv | `HSP_Induci_Prod_min/max` | 0.003 | **0.05** |
| parameters.csv | `HSP_Consti_Prod_min/max` | 0.0015 | 0.0015 (already the settled value) |
| runtime.txt | `x`, `y`, `gridsize` | 100, 100, 10000 | **20, 20, 400** (grid size used in all analyses) |

Unchanged but part of the settled set: `Maint_E_HSP_Ea_n` 0.2, `Maint_E_Rate_Tref` 273, `HSP_Consti_Ea` / `HSP_Induci_Ea` 30000,
`Erat_b` 700, `Erange` 700, `Enz_Prod` 0.002, `Constit_Prod` 0.0002, `Uptake_C_cost` 0.0375, `Enz_Prod_E_cost` / `Uptake_Prod_E_cost` 600,
`Active_Uptake_E_cost` 10, `HSP_prod_N_cost` 0.3, `HSP_turnover` 0.1, `Ferm_Cond` 0.1, `death_rate_bac` 0.001.
Initial substrate IdealSub C:N:P = 100:5:0.5 (200:10:1). Alternative tested set: `E_resp` 1700 (everything else the same).

Notes
- The deterministic identity monomer-transporter mapping (`src/monomer.py`, applies when `n_uptake == n_monomers` and `Uptake_per_monomer == 1`) is required to reproduce the results.
- With a single taxon, drought tolerance was 0 by construction (normalised across taxa; with the current `src/microbe.py` the relative scale starts at `Osmo_tol_lower`, default 0, so it is now 1 unless `Osmo_tol_lower` >= the allocation), so a drying Psi kills the population; the analyses with real temperature used constant Psi -0.5.
- `parameters_minimal.csv` (a backup copy of the old parameters in `input_minimal/`) was deliberately not copied here.
