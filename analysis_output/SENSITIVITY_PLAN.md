# Sensitivity-analysis plan: minimal and full model

Written 2026-10-19 after reviewing the analyses run since the last behaviour-changing code update (2026-10-16: repopulation
at the initial density, psi-slope defaults 0.1, starvation and dead-C tracking; later changes were diagnostics and the
`Decay_MM_Form` switch, off by default). Results are in `SENSITIVITY_SUMMARY.md`; this file is the plan and the
decisions that frame it. Tick items off here as they are done.

## Decisions (2026-10-19)

| Item | Decision |
|---|---|
| Drought-mortality slope | `beta_bac` = `beta_fun` = 1 (all input folders) |
| Endogenous respiration | On by default: `R_endo,1` in every input folder; the code default for a missing row is now 1 |
| HSP | HSP should carry a modest cost at cold temperatures (resources diverted from enzymes and growth), not an extreme loss of efficiency |
| Sensitivity design | Morris elementary effects (catches interactions and parameters not considered yet); a few one-at-a-time curves only where a shape is needed |
| Fungi | Skipped: bacterial minimal model, and the full model run without fungi (`fb` 0). Fungal distribution is a future project |
| CUE | Emergent CUE 0.26-0.35 is accepted; higher CUE is a secondary benchmark when choosing between otherwise equivalent parameter sets |
| Mineral uptake (2026-10-19) | NH4 and PO4 enter by passive diffusion, no transport energy costs (`Mineral_Uptake_Passive` 1) |
| Rebalancing energy (2026-10-19) | C respired in stoichiometric rebalancing yields energy (`R_stoich` 1, independent of `R_endo`) |
| Nutrient inputs (2026-10-19) | `Input_NH4` 0.05, `Input_PO4` 0.005 (minimal km25 folders): no N or P limitation |
| HSP cost (2026-10-19) | `HSP_prod_E_cost` 600: HSP pays off above about 10 C, -7 to -14 % biomass at 5-7.5 C |
| Energy-limitation metric (2026-10-19) | Endogenous-respiration C while > 5 % substrate is left (or its share of respiration); a count of days with any cell in deficit grows with grid size |

Calibration targets carried forward (minimal and full model): about 80 % or more of the substrate degraded with a steady
decline (not front-loaded); ambient close to, or barely, energy-limited; +5 C tipping into energy limitation; organisms
able to grow somewhat in the dry season.

## Phase 1: close out the minimal model (bacterial, `input_minimal_e1500_crich_km25_in`)

- [x] **1.1 Re-baseline the final set.** 5 seeds, three climates, 5 pulses, R_endo on and off; degradation curve, CUE,
  energy-limited days, minimum E:C, kills, biomass. Grid-size check on the baseline: 20x20 vs 50x50, plus one 100x100.
  Decide the replication (seeds) and grid size used for the Morris runs from this.
  - 2026-10-19: done (section 34). Pulses reproducible, seed SD small, outputs grid-size invariant: Morris at 20x20,
    one pulse, 2 seeds.
- [x] **1.2 HSP cold cost.** HSP vs control (HSP production off) on a constant-temperature gradient (0-40 C, wet Psi),
  MMRT. Target: a modest penalty at low temperature (earlier window: 5-10 % lower biomass) with a benefit at warm
  temperatures. If the cold penalty is missing, tune HSP production or cost (`HSP_Consti_Prod`, `HSP_prod_E_cost`) or
  the Ea slope; keep `Maint_E_Rate_Tref` and the maintenance rate fixed, as they set the energy-limitation calibration.
  - 2026-10-19: lever screen done (`SENSITIVITY_SUMMARY.md` section 28). The current crossover is about 5 C (no cold
    cost in any climate); HSP_prod_E_cost 600 gives about 9 C, Tref 283 K with a rescaled rate constant about 12.5 C.
    Choose the lever, then re-baseline (1.1).
  - 2026-10-19: done (sections 29-33): N limitation removed, HSP_prod_E_cost 600 chosen.
- [ ] **1.3 MMRT vs Arrhenius at the final set.** Gradient and real climate; re-match `Vmax0` / `Uptake_Vmax0` in the
  15-20 C range if needed (the current match was made at maintenance 2.5).
- [ ] **1.4 Morris screening.** Three climates analysed separately, plus the warming response (+5 C minus ambient, and
  +10 C minus ambient, from the same parameter points). Outputs: % degraded at day 365, day 50 % and 95 % degraded,
  annual CUE, energy-limited days and minimum E:C while > 5 % substrate is left, drought kills, maximum and mean biomass.
  Report mu* (overall influence) and sigma (interactions and non-linearity) per output, with bootstrap confidence intervals.
  - Design: p = 4 levels, r = 10 trajectories for a first pass (about 25 factors -> 260 points x 3 climates x 2 seeds
    ~ 1,560 runs, about 2 h here); r = 20 for the reported version if rankings are not yet stable. Factors are scaled on
    a log range (x/÷ 1.5, or a stated range for slopes and costs) around the final set.
  - Candidate factors:
    - Energy and maintenance: `Maint_E_Rate_Const`, `Maint_E_HSP_Ea_Intercept`, `Maint_E_HSP_Ea_Slope`,
      `Mon_IdealSub` E_resp (needs `set_params` support for `monomer_energy.csv`), `Ferm_Cond`.
    - HSP: `HSP_Consti_Prod`, `HSP_Induci_Prod`, `HSP_prod_E_cost`, `HSP_turnover`.
    - Enzymes: `Constit_Prod`, `Enz_Prod`, `Enz_Prod_E_cost`, `Enz_loss_rate`, `Vmax_Km`.
    - Uptake: `Uptake_Vmax_Km`, `Uptake_C_cost`, `Active_Uptake_E_cost`, `Uptake_Prod_E_cost`.
    - Water and drought: `beta_bac`, `Osmo_tol_K`, `Osmo_Consti_Prod`, `Osmo_prod_E_cost`, `Osmo_loss_rate`,
      `Enzyme_psi_slope`, `Uptake_psi_slope`.
    - Mortality: `death_rate_bac`.
  - Follow-up: one-at-a-time curves (with seeds) for the top 4-6 factors, to show response shapes and thresholds.
- [ ] **1.5 Write-up.** Consolidate `SENSITIVITY_SUMMARY.md`: final parameter table, Morris table, and earlier sections
  marked as superseded where the final set changed their numbers.

## Phase 2: port to the full model (no fungi)

- [ ] **2.1 Translate the parameters.** Enzyme production is per gene (optionally normalized by gene count,
  `NormalizeProd`), so match each taxon's total production, not the per-gene value. Carry over maintenance 8,
  `Vmax_Km` 6, `Uptake_Vmax_Km` 0.25, osmolyte cost 150, beta 1, psi slopes 0.1, R_endo 1, and the Arrhenius
  `Vmax0` / `Uptake_Vmax0` (matched on the main substrate only; check across the 12 substrates, which have their own Ea).
  Build the input folder with a script, as for the earlier full-model folders.
- [ ] **2.2 Single-taxon check.** A minimal-equivalent taxon on the full substrate set: does it reproduce the minimal
  behaviour?
- [ ] **2.3 Community calibration** (100 taxa, small grid, `fb` 0). Targets as above. Main levers: centre and width of
  the production ranges, `taxa_per_box`.
- [ ] **2.4 Full-model Morris.** Community-level factors (trait-range widths, share of HSP taxa, `n_taxa`,
  `taxa_per_box`) plus the top 4-6 factors from 1.4; gradient and climates; MMRT and Arrhenius.
- [ ] **2.5 Production-size check on HPC3** (100x100, 5 pulses) before the community experiments.
