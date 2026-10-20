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
| HSP cost (2026-10-19) | `HSP_prod_E_cost` 600: HSP pays off above about 10 C, -7 to -14 % biomass at 5-7.5 C (at maintenance 8). Kept at 600 after the maintenance change: crossover about 12.5 C, -15 to -17 % biomass at 2.5-7.5 C (section 39) |
| Maintenance and enzyme turnover (2026-10-19) | `Maint_E_Rate_Const` 4 (was 8), `Enz_loss_rate` 0.05 (was 0.04) in the minimal km25 folders: steady ambient decomposition through the year, CUE 0.36-0.37, warm runs finish the substrate and then become energy-limited (SENSITIVITY_SUMMARY 38-38b). Confirmed over 5 pulses x 5 seeds (section 39); minimal baseline final |
| Negligible parameters (2026-10-19) | Kept at their baseline values as defaults (Morris, sections 36-37): `Enz_Prod_E_cost` 600, `Active_Uptake_E_cost` 10, `Uptake_Prod_E_cost` 600, `Uptake_C_cost` 0.0375, `Uptake_Vmax_Km` 0.25, `Uptake_psi_slope` 0.1; HSP and osmolyte production rates and osmolyte parameters unchanged |
| Energy-limitation metric (2026-10-19) | Whole-year endogenous-respiration C and its share of respiration, substrate left when energy limitation starts, substrate left at day 365 (section 35). Not a count of days with any cell in deficit (grows with grid size), and not only the > 5 %-substrate window (misleads between kinetics) |
| Full-model port (2026-10-19) | `input_full_port_bac/` (`analysis/build_full_input_port.py`): litter E_resp x1.245 (e1500 table), necromass / FermProd rows at the minimal values, enzyme production 4.5x the minimal totals per taxon (NormalizeProd 1; Constit 5.0625e-5, Enz_Prod 9e-4 per gene), Input_NH4 0.1. Single taxon (ambient, MMRT): 91 % degraded, day 50 % 130, day 95 % 313, CUE 0.37, endogenous respiration 34 %. Accepted as is; no change to cellulose / lignin yields or to labile-substrate decay (sections 41-41f) |
| Full-model community settings (2026-10-20) | `input_full_port_bac/` rebuilt: `Input_NH4` 0.15, `Input_PO4` 0.01, `Enzymes_per_sub` 2, `Enz_loss_rate` 0.02, Constit_Prod 1e-5-2e-5 and Enz_Prod 2e-4-4e-4 per gene (NormalizeProd 1), other trait ranges x/1.25, 30-40 enzyme genes. Ten seeds (ambient, 20x20): 94 % degraded (85-99), day 50 % 143, day 95 % 306, net CUE 0.17, 94 taxa alive on day 365; constant over five pulses (sections 43-48k). Net CUE (C retained / C taken up over the year and all taxa, unfloored) is the CUE to report. Not adopted: daily substrate inputs, labile-substrate Km, repopulation window (deferred) |

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
- [x] **1.3 MMRT vs Arrhenius at the final set.** Gradient and real climate; re-match `Vmax0` / `Uptake_Vmax0` in the
  15-20 C range if needed (the current match was made at maintenance 2.5).
  - 2026-10-19: done (section 35). The 12.5-20 C match holds; Vmax0 / Uptake_Vmax0 unchanged.
- [x] **1.4 Morris screening** (second pass done 2026-10-19, section 40; set up: `analysis/morris_factors_minimal_v2.csv`, 14 factors around the final baseline (the 6 negligible and 6 borderline factors held at baseline), narrowed ranges for Maint_E_Rate_Const 2.5-6.4, Enz_loss_rate 0.035-0.07, HSP_prod_E_cost 300-900, beta 0.5-2; r = 20 from a pool of 80, seed 2, 300 points, 3,600 runs; design in `analysis_output/morris_minimal_v2/`. First pass done 2026-10-19, sections 36-37: 3,240 runs; open: 20 trajectories if the middle of the ranking matters; Maint_E_Rate_Const x Enz_loss_rate grid (3-8 x 0.02-0.06, `analysis/grid2d_sweep.py`) to choose both; `analysis/morris_minimal.py`, factor table `analysis/morris_factors_minimal.csv`, HPC3 array `analysis/submit_morris_minimal.sh`; both kinetics, MMRT and Arrhenius; 26 factors (Ferm_Cond excluded), `Maint_E_Rate_Const` 2.5-10). Three climates analysed separately, plus the warming response (+5 C minus ambient, and
  +10 C minus ambient, from the same parameter points). Outputs: % degraded at day 365, day 50 % and 95 % degraded,
  annual CUE, whole-year endogenous respiration and its share of respiration, substrate left at the onset of energy
  limitation and at day 365, drought kills, maximum and mean biomass.
  Report mu* (overall influence) and sigma (interactions and non-linearity) per output, with bootstrap confidence intervals.
  - Design: p = 4 levels, r = 10 trajectories for a first pass (about 25 factors -> 260 points x 3 climates x 2 seeds
    ~ 1,560 runs, about 2 h here); r = 20 for the reported version if rankings are not yet stable. Factors are scaled on
    a log range (x/÷ 1.5, or a stated range for slopes and costs) around the final set.
  - Candidate factors:
    - Energy and maintenance: `Maint_E_Rate_Const`, `Maint_E_HSP_Ea_Intercept`, `Maint_E_HSP_Ea_Slope`,
      `Mon_IdealSub` E_resp.
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

- [x] **2.1 Translate the parameters** (2026-10-19, section 41: `analysis/build_full_input_port.py` -> `input_full_port_bac/`; litter E_resp x1.245, necromass rows at the minimal values, NormalizeProd 1 with the minimal production totals). Original plan: Enzyme production is per gene (optionally normalized by gene count,
  `NormalizeProd`), so match each taxon's total production, not the per-gene value. Carry over maintenance 8,
  `Vmax_Km` 6, `Uptake_Vmax_Km` 0.25, osmolyte cost 150, beta 1, psi slopes 0.1, R_endo 1, and the Arrhenius
  `Vmax0` / `Uptake_Vmax0` (matched on the main substrate only; check across the 12 substrates, which have their own Ea).
  Build the input folder with a script, as for the earlier full-model folders.
- [x] **2.2 Single-taxon check** (done 2026-10-19, sections 41-41f; settings accepted, see the decisions table. First pass, section 41: at the minimal production totals the taxon degrades 9 % of the litter; enzyme specificity cuts decay capacity per unit production about 7x; production x3.5-7 restores 85-98 % with a fast labile phase and a cellulose/lignin tail; open: production factor, realised energy yield (x1.45), early CUE drop). A minimal-equivalent taxon on the full substrate set: does it reproduce the minimal
  behaviour?
- [x] **2.3 Community calibration** (done 2026-10-20, sections 42-48k; settings in the decisions table. First runs 2026-10-19, section 42: narrow trait ranges x/1.25, 30-40 genes; 97-98 % degraded, CUE 0.33-0.35, N limitation back at NH4 0.1) (100 taxa, small grid, `fb` 0). Targets as above. Main levers: centre and width of
  the production ranges, `taxa_per_box`.
- [ ] **2.4 Full-model Morris.** Community-level factors (trait-range widths, share of HSP taxa, `n_taxa`,
  `taxa_per_box`) plus the top 4-6 factors from 1.4; gradient and climates; MMRT and Arrhenius.
- [ ] **2.5 Production-size check on HPC3** (100x100, 5 pulses) before the community experiments.
