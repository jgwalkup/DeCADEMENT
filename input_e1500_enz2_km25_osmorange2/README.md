# input_e1500_enz2_km25_osmorange2

Copy of `input_e1500_enz2_km25/` (enzyme and uptake Vmax-Km slopes 2.5; `Osmo_tol_lower` 0, `Osmo_tol_ref` 0) with ranges for both constitutive and inducible osmolyte production
(built by `analysis/build_full_input_osmo_range.py --induci 0.04 0.08 --consti 0.002 0.005`):

| parameter | `input_e1500_enz2_km25` | now |
|---|---|---|
| `Osmo_Consti_Prod_min` / `_max` | 0.005 / 0.005 | 0.002 / 0.005 |
| `Osmo_Induci_Prod_min` / `_max` | 0.05 / 0.05 | 0.04 / 0.08 |

Each taxon has one osmolyte gene and draws its own constitutive and inducible rates independently, so total allocation (constitutive + inducible) runs 0.042-0.085 (mean about 0.0635).
With `Osmo_tol_ref` 0 and `Osmo_tol_lower` 0, relative drought tolerance = allocation / maximum allocation (about 0.5-1 across taxa).
Everything else is unchanged. Run with the `deca_dement` conda environment.


Mineral inputs (set after this folder was first used): `Input_NH4` 0.02 (N per grid cell per day) and `Input_PO4` 0.002 (P per grid cell per day), applied in `Grid.uptake`. Earlier results obtained with this folder before this change used inputs of 0. Everything else is unchanged.

Moisture sensitivity (set after this folder was first used): `Uptake_psi_slope` 0.25 (was 0.5, the former hardcoded value) and `Enzyme_psi_slope` 0.25 (equal to the former hardcoded value): below `wp_fc` the uptake and enzyme Vmax are multiplied by exp(slope x (Psi - wp_fc)). Earlier results from this folder used an uptake slope of 0.5.
