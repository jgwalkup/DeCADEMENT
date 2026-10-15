# input_e1500_enz2_km25_osmorange

Copy of `input_e1500_enz2_km25/` (enzyme and uptake Vmax-Km slopes 2.5; `Osmo_tol_lower` 0, `Osmo_tol_ref` 0) with a range of osmolyte production across taxa
(built by `analysis/build_full_input_osmo_range.py`):

| parameter | before | now |
|---|---|---|
| `Osmo_Induci_Prod_min` / `_max` | 0.05 / 0.05 | 0.01 / 0.09 |
| `Osmo_Consti_Prod_min` / `_max` | 0.005 / 0.005 | unchanged (0.005 / 0.005) |

Each taxon has one osmolyte gene (`Osmo_per_taxon` 1) and draws its own rates, so total allocation (constitutive + inducible) runs 0.015-0.095 with mean 0.055 (the earlier settled value).
With `Osmo_tol_ref` 0 and `Osmo_tol_lower` 0, relative drought tolerance = allocation / maximum allocation, i.e. 0.16-1 across taxa (mean 0.58).
Everything else is unchanged. Run with the `deca_dement` conda environment.


Mineral inputs (set after this folder was first used): `Input_NH4` 0.02 (N per grid cell per day) and `Input_PO4` 0.002 (P per grid cell per day), applied in `Grid.uptake`. Earlier results obtained with this folder before this change used inputs of 0. Everything else is unchanged.

Moisture sensitivity (set after this folder was first used): `Uptake_psi_slope` 0.25 (was 0.5, the former hardcoded value) and `Enzyme_psi_slope` 0.25 (equal to the former hardcoded value): below `wp_fc` the uptake and enzyme Vmax are multiplied by exp(slope x (Psi - wp_fc)). Earlier results from this folder used an uptake slope of 0.5.
