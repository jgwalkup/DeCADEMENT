# input_minimal_e1500_km25

Copy of `input_minimal_e1500/` with the four Vmax-Km slopes set to 2.5 (as in `input_e1500_enz2_km25`): `mmrt_Vmax_Km_slope` and `Vmax_Km` 7 -> 2.5,
`Uptake_mmrt_Vmax_Km_slope` and `Uptake_Vmax_Km` 5 -> 2.5 (built by `analysis/build_full_input_km.py`). Substrate unchanged: IdealSub C 100, N 5, P 0.5.


Mineral inputs (set after this folder was first used): `Input_NH4` 0.02 (N per grid cell per day) and `Input_PO4` 0.002 (P per grid cell per day), applied in `Grid.uptake`. Earlier results obtained with this folder before this change used inputs of 0. Everything else is unchanged.

Moisture sensitivity (set after this folder was first used): `Uptake_psi_slope` 0.25 (was 0.5, the former hardcoded value) and `Enzyme_psi_slope` 0.25 (equal to the former hardcoded value): below `wp_fc` the uptake and enzyme Vmax are multiplied by exp(slope x (Psi - wp_fc)). Earlier results from this folder used an uptake slope of 0.5.
