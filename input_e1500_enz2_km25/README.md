# input_e1500_enz2_km25

Copy of `input_e1500_enz2/` with all four Vmax-Km slopes set to 2.5 (built by `analysis/build_full_input_km.py`):

| parameter | kinetics | role | before | now |
|---|---|---|---|---|
| `mmrt_Vmax_Km_slope` | MMRT | enzyme (degradation) Km = MMRT term x slope | 7 | 2.5 |
| `Vmax_Km` | Arrhenius | enzyme Km0 = Vmax0 x slope | 7 | 2.5 |
| `Uptake_mmrt_Vmax_Km_slope` | MMRT | transporter Km = Vmax x slope | 5 | 2.5 |
| `Uptake_Vmax_Km` | Arrhenius | transporter Km0 = Vmax0 x slope | 5 | 2.5 |

Everything else (energy yields, maintenance/HSP block, enzyme and transporter production, initial enzyme pool, `taxa_per_box` 0.05, `Kinetics_Model` 1) is unchanged from `input_e1500_enz2`; the climate file is the real climate. Drought tolerance rows added to `parameters.csv`: `Osmo_tol_lower` 0 and `Osmo_tol_ref` 0. With `Osmo_tol_ref` 0, tolerance is relative: (allocation - `Osmo_tol_lower`) / (maximum allocation - `Osmo_tol_lower`), clipped to 0-1, so the scale runs from 0 to the maximum osmolyte allocation of the community (`Osmo_tol_lower` < 0 would use the community minimum, the original method). All taxa have the same allocation here (0.005 constitutive + 0.05 inducible = 0.055), so every taxon has tolerance 1 (no drought mortality increase); it only differs between taxa if `Osmo_Consti_Prod` / `Osmo_Induci_Prod` are given a min-max range.

Why 2.5: the enzyme and uptake slope scans (`analysis_output/km_sweep_1to3/`, `analysis_output/uptake_sweep/`) showed degradation was Km-limited at the old values. The combination enzyme 2.5 / uptake 2.5 itself was not part of those scans (the crossed scan used enzyme slopes 2, 3, 4); at uptake 2.5, MMRT 24 C, the enzyme-slope-3 and enzyme-slope-2 results bracket it (72 % and 81 % degraded, energy:C 175 and 71, CUE 0.27 and 0.24). Measured for the combination itself (`analysis_output/full_gradient_resupply_km25/`, flat wet climate, 5 communities, HSP, MMRT, pulse 1, day 365): substrate degraded 42 / 59 / 72 / 78 / 81 / 81 / 81 / 77 / 72 / 62 % at 4 / 8 / 12 / 16 / 20 / 24 / 28 / 32 / 36 / 40 C; biomass C 55-82 (28.6 at the old slopes in the flat reference); energy:C 149 at 20 C, 82 at 24 C, 55 at 28 C, 53 at 32 C; CUE 0.30 at 16 C, 0.24 at 24 C, 0.21-0.23 at 28-40 C. Energy:C in the no-HSP control is 0-2 from 24 C upward. Details in `analysis_output/SENSITIVITY_SUMMARY.md`.

Not recalibrated at the new slopes: the enzyme-production multiplier (2), initial enzyme pool (0.004-0.005) and energy yields were set at the old slopes.
Run with the `deca_dement` conda environment.


Mineral inputs (set after this folder was first used): `Input_NH4` 0.02 (N per grid cell per day) and `Input_PO4` 0.002 (P per grid cell per day), applied in `Grid.uptake`. Earlier results obtained with this folder before this change used inputs of 0. Everything else is unchanged.

Moisture sensitivity (set after this folder was first used): `Uptake_psi_slope` 0.25 (was 0.5, the former hardcoded value) and `Enzyme_psi_slope` 0.25 (equal to the former hardcoded value): below `wp_fc` the uptake and enzyme Vmax are multiplied by exp(slope x (Psi - wp_fc)). Earlier results from this folder used an uptake slope of 0.5.
