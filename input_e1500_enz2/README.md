# input_e1500_enz2

Full multi-taxon, multi-substrate model inputs adopting the parameters settled in the minimal-model analyses
(see `analysis_output/SENSITIVITY_SUMMARY.md`, sections 3.10-3.11). Based on `input/`; only `parameters.csv`,
`monomer_energy.csv` and `runtime.txt` differ. `input/` and `input_e1500/` are unchanged.

1. Energy yields: every `E_resp` and `E_ferm` x2, rounded to the nearest ten (FermProd `E_ferm` 2 -> 4); `E_value` unchanged.
2. Maintenance / HSP block: `Maint_E_Rate_Const` 2.5, `Maint_E_HSP_Ea_Intercept` 45000, `Maint_E_HSP_Ea_Slope` -90000,
   `HSP_prod_E_cost` 300; `HSP_Consti_Prod` 0.0015, `HSP_Induci_Prod` 0.05 (pegged).
3. Enzyme / transporter production: per-cell totals matched to the minimal model (`input_e1500`), then the enzyme rate ranges
   (`Enz_Prod`, `Constit_Prod`) doubled (calibration multiplier 2): `Enz_Prod` 9.24e-5 to 3.08e-4, `Constit_Prod` 1e-5 to 3e-5;
   `Uptake_C_cost` 0.01221 to 0.02441. Community-mean per-cell totals about 0.00396 (inducible enzymes), 0.00039 (constitutive), 0.225 (transporters).
4. `runtime.txt`: `taxa_per_box` 0.05 (was 0.20) to reduce the initial biomass-to-substrate mismatch.
5. Initial enzyme pool: `Enz_min` 0.004 and `Enz_max` 0.005 (was 0.04 and 0.05), which removes the early degradation burst and
   early energy exhaustion that the full community showed with the original pool.

Rebuild: `python analysis/build_full_input_e1500.py --dst ../input_e1500_enz2 --enz-mult 2 --enz-init 0.004 0.005 --taxa-per-box 0.05`
(verified identical to this folder's `parameters.csv`, `monomer_energy.csv` and `runtime.txt`).
Run with the `deca_dement` conda environment.

Check against the minimal model (20x20, flat wet climate, 365 days, mean of 3 communities; minimal model in brackets), degraded at days 100 / 200 / 365
and energy:C at day 365:
  15 C: 5.6 / 11.9 / 26.8% and 393 (16 C: 6.2 / 11.6 / 29.0% and 662);  28 C: 6.2 / 12.9 / 28.2% and 249 (7.3 / 14.5 / 39.4% and 568);
  36 C: 5.5 / 11.1 / 22.2% and 110 (5.7 / 10.5 / 25.6% and 466).

Known issues / open choices
- Drought tolerance was 0 for all taxa when this setup was run (identical osmolyte allocation) (with the current `src/microbe.py`, relative tolerance has a lower bound of `Osmo_tol_lower`, default 0, so identical allocation now gives tolerance 1; add `Osmo_tol_lower,<allocation>` to reproduce the earlier tolerance-0 behaviour); `Osmo_tol_ref` (opt-in, see `input_minimal_e1500_osmo`) not applied here.
- The full community's degradation is less temperature-sensitive than the minimal model's (8 C: 20% vs 11%; 28 C: 28% vs 39% at day 365).
- At 36 C energy:C falls below 100 around day 285.
