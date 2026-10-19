# Open model issues

Items from the 2026-10-15 code review that are decided but not yet implemented, or deliberately deferred. Update this
file as items are resolved (the git log has the details of each change).

## Decided, to implement next

- **N limitation and unregulated mineral uptake** (found 2026-10-19; SENSITIVITY_SUMMARY section 29). With the C-rich
  substrate (C:N 40) and `Input_NH4` 0.02, the warm runs hit the minimum N quota during peak growth (rebalancing respires
  C; CUE drops). Raising N supply removes that, but NH4 uptake is not regulated by N demand: cells hoard N up to the C-share
  minimum, rebalancing returns the excess to NH4, and every unit taken up is charged the active-transport and
  transporter-production energy costs, so the cycle drains energy (collapse at `Input_NH4` >= 0.1). The intent is that N
  never limits. Options: (1) cap NH4 (and PO4) uptake at each cell's need (deficit to its optimal ratio plus the day's
  growth), then raise `Input_NH4` until N never limits; (2) exempt mineral uptake from the transport energy costs;
  (3) parameters only (lower substrate C:N or more NH4), which leaves the hoarding cycle in place. Decision pending.

## Deferred (revisit later)

- **Organic monomer and substrate inputs.** `Input_NH4`/`Input_PO4` act in both the minimal and the full model (checked
  2026-10-15). Bin Wang's version also routed organic monomer inputs (C:N:P from `elem_ratios.csv`) and added
  substrate inputs (`SubInput`) daily; the `Sub`/`Mon` columns of `sub_mon_inputs.csv` are inert. Deferred.

- **Fungal distribution in the minimal model.** Fungal placement probability is scaled by `max_size_b/max_size_f`
  (about 4% of cells occupied at t0 with `taxa_per_box` 1), and fungal translocation averages every living fungal
  cell each day. Use the bacterial minimal variant for sensitivity analyses until this is revisited.
  Repopulation now applies the same fungal scaling as the initial placement (fixed 2026-10-16; before, pulses 2-6
  started with every cell occupied).
- **Repopulation density (resolved 2026-10-16).** DEMENTpy's repopulation placed on average one
  (bacterial-equivalent) taxon per grid cell after each pulse, while the initial placement places
  n_taxa x taxa_per_box (5 in the full model), so later pulses started about 5 times thinner (20x20 mixed test: 1.0 vs
  4.9 C per grid cell). Repopulation now restores the t0 density: shares (fungi in bacterial-equivalent units) are
  multiplied by n_taxa x taxa_per_box, a taxon is capped at one cell per grid cell, and the excess of capped taxa goes
  to the others in proportion to their shares (`Grid._repop_probabilities`). One-taxon minimal runs are unchanged
  (bit-identical); full-model multi-pulse results from before this change (e.g. the start-biomass test and the
  pulse-to-pulse energy decline) used the old density.
- **Transporters as tracked cellular compounds.** The two transport energy costs work as designed
  (transporter production per transporter produced; active transport per monomer taken up). Treating transporters
  like HSP and osmolytes (tracked pools) would need their production computed before `uptake()`.
- **`Ferm_Cond` plays three roles**: fraction of energy-C fermented, exponent of the energy-demand curve, exponent of
  the monomer energy-value transform.
- **Leaching shape and scope** (left as is by decision, 2026-10-15). The rate
  `Leaching * exp(Leaching_psi_slope * (Psi - wp_fc))` also applies in dry soil; only NH4 and PO4 leach. About a third
  of the NH4 input leaches over a year. Daily and per-pulse leached N and P are in the output.
- **Inducible osmolyte production is nearly inert under drought** (tied to uptake and leftover energy, both of which
  collapse); drought tolerance currently behaves as a constitutive trait.
- **Endogenous respiration (R_endo) default.** The switch `R_endo` is 0 in every input folder (shortfalls floored away; the output
  `E_unpaid_Taxon` shows how much). The maintenance sensitivity runs use both settings.
- **Drought response sets decomposition under the real climate.** Drought mortality cuts first-year degradation from
  about 43 % to 8 % in the minimal model (SENSITIVITY_SUMMARY 6); tolerance (`Osmo_tol_K`, osmolyte production) and
  `death_rate`/`beta` therefore need calibrating before the climate runs are interpreted.
- **Calibration**: `HSP_max_frac`, `Osmo_max_frac` (defaults 0.2) and `Osmo_tol_K` (default 0.05) are placeholders;
  the linear composition feedback acts as a constant tax well below the cap.

## Housekeeping

- One shared run harness for the analysis scripts (consistent parameter overrides and defaults; e.g. `n_tracker.py`
  overwrites osmolyte parameters unless `--keep-osmo`).
- Obsolete since the osmolyte-content tolerance: `analysis/full_climate_osmo_range.py`,
  `plot_full_climate_osmo_range.py`, `plot_full_climate.py`, `plot_full_climate_km25.py`,
  `build_full_input_osmo_range.py` (they set or label `Osmo_tol_ref`/`Osmo_tol_lower`, which no longer exist).
- Parameter rows that are never read: `Uptake_H_S_error`, `HSP_prod_C_cost`, `HSP_prod_P_cost`, `Km_Ea`.
- Older module copies in `src/` (`grid_vJW.py`, `initialization_big_files.py`) still use the retired Km rows and
  tables; they are not imported by the model.
- Performance: stop expanding grid-uniform tables to every cell; cache temperature-dependent kinetics.
- **Decay rate-law switch (2026-10-16).** `Decay_MM_Form` (runtime.txt) and `Enz_Sat_K` (parameters.csv) are provisional
  names (the runtime row `Kinetics_Model` already selects MMRT vs Arrhenius). To decide: final names; whether enzyme
  saturation should also apply to uptake (currently degradation only); calibration of Enz_Sat_K and enzyme production
  for form 1 (SENSITIVITY_SUMMARY 13).
