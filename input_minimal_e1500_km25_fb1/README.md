# input_minimal_e1500_km25_fb1

Copy of `input_minimal_e1500_km25/` with `fb` set to 1 in `runtime.txt` (was 0): the single taxon is fungal (initial cell C 25, division above 50, C:N:P 0.9:0.09:0.01,
translocation of biomass between its cells, +x dispersal). With `taxa_per_box` 1 it is placed in a cell with probability 1 x max_size_b / max_size_f = 0.04.
Everything else is unchanged.


Mineral inputs (set after this folder was first used): `Input_NH4` 0.02 (N per grid cell per day) and `Input_PO4` 0.002 (P per grid cell per day), applied in `Grid.uptake`. Earlier results obtained with this folder before this change used inputs of 0. Everything else is unchanged.

Moisture sensitivity (set after this folder was first used): `Uptake_psi_slope` 0.25 (was 0.5, the former hardcoded value) and `Enzyme_psi_slope` 0.25 (equal to the former hardcoded value): below `wp_fc` the uptake and enzyme Vmax are multiplied by exp(slope x (Psi - wp_fc)). Earlier results from this folder used an uptake slope of 0.5.
