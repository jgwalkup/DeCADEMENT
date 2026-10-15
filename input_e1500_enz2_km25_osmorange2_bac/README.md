# input_e1500_enz2_km25_osmorange2_bac

Copy of `input_e1500_enz2_km25_osmorange2/` with `fb` set to 0 in `runtime.txt` (was 0.2): `fb` is the probability that a taxon is fungal, so every taxon is a bacterium
(initial cell C = 1 instead of 25 for fungi; bacterial death-rate, size and stoichiometry parameters apply to all taxa).
Everything else is unchanged: slopes 2.5/2.5, osmolyte ranges (constitutive 0.002-0.005, inducible 0.04-0.08), `Osmo_tol_lower` 0, `Osmo_tol_ref` 0.
Because `fb` changes the random draws made at initialisation, the taxon traits of a given seed differ from those of the folder with fungi.


Mineral inputs (set after this folder was first used): `Input_NH4` 0.02 (N per grid cell per day) and `Input_PO4` 0.002 (P per grid cell per day), applied in `Grid.uptake`. Earlier results obtained with this folder before this change used inputs of 0. Everything else is unchanged.

Moisture sensitivity (set after this folder was first used): `Uptake_psi_slope` 0.25 (was 0.5, the former hardcoded value) and `Enzyme_psi_slope` 0.25 (equal to the former hardcoded value): below `wp_fc` the uptake and enzyme Vmax are multiplied by exp(slope x (Psi - wp_fc)). Earlier results from this folder used an uptake slope of 0.5.
