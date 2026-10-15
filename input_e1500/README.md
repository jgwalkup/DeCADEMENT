# input_e1500

Full multi-taxon, multi-substrate model inputs with the parameters settled in the minimal-model analyses
(`analysis_output/SENSITIVITY_SUMMARY.md`). Copy of `input/` built by `analysis/build_full_input_e1500.py`; only
`parameters.csv` and `monomer_energy.csv` differ. `input/` is unchanged.

1. Energy yields (`monomer_energy.csv`): every `E_resp` and `E_ferm` multiplied by 2, rounded to the nearest ten
   (a non-zero yield is never rounded to 0, so FermProd `E_ferm` 2 -> 4). `E_value` unchanged.
2. Settled maintenance / HSP block: `Maint_E_Rate_Const` 2.5, `Maint_E_HSP_Ea_Intercept` 45000, `Maint_E_HSP_Ea_Slope` -90000,
   `HSP_prod_E_cost` 300; `HSP_Consti_Prod` 0.0015 and `HSP_Induci_Prod` 0.05 pegged (min = max).
3. Enzyme / transporter production matched to the minimal model's per-cell totals (rate x number of genes), keeping the
   relative width of each range: `Enz_Prod` 4.62e-5 to 1.54e-4 (inducible total 0.002 at 20 genes), `Constit_Prod` 5e-6 to 1.5e-5
   (constitutive total 0.0002), `Uptake_C_cost` 0.01221 to 0.02441 (transporter total 0.225 at 12.29 genes).
   Measured community means over 6 random communities: 0.00193 / 0.00019 / 0.22489 (targets 0.002 / 0.0002 / 0.225).

Not changed: osmolyte parameters (drought tolerance was 0 for all taxa when built, because all have the same allocation (with the current `src/microbe.py`, relative tolerance has a lower bound of `Osmo_tol_lower`, default 0, so identical allocation now gives tolerance 1; add `Osmo_tol_lower,<allocation>` to reproduce the earlier tolerance-0 behaviour); see `Osmo_tol_ref` in
`input_minimal_e1500_osmo/README.md`), `n` (0.2), `Maint_E_Rate_Tref`, kinetics ranges, runtime.
Run the model with the `deca_dement` conda environment (the base environment's scipy 1.15.2 no longer imports on macOS 27).
