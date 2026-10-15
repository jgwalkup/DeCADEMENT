# Open model issues

Items from the 2026-10-15 code review that are decided but not yet implemented, or deliberately deferred. Update this
file as items are resolved (the git log has the details of each change).

## Decided, to implement next

- **Energy deficit (optional mechanism).** When energy costs exceed the energy available, pay the deficit by respiring
  free biomass C, with the energy yield of the necromass monomer (`Mon_DeadMic` in `monomer_energy.csv`). Off by
  default; when off, the model behaves as now (the energy pool is floored at zero, so the shortfall is created from
  nothing and an empty energy pool is never lethal).
- **One set of Km parameters for both kinetics models.** Merge `Vmax_Km`/`mmrt_Vmax_Km_slope`,
  `Vmax_Km_int`/`mmrt_Vmax_Km_int`, `Uptake_Vmax_Km`/`Uptake_mmrt_Vmax_Km_slope`,
  `Uptake_Vmax_Km_int`/`Uptake_mmrt_Vmax_Km_int`, `Km_error`/`mmrt_Km_error`. Reuse the Vmax thermodynamic tables
  for Km instead of a second random draw (`Km_HT0`, `Km_ST0`, `Km_Cp`).
- **Split `Erange`** into `E_min_ratio` (current value 0) and `E_max_ratio` (current value 1400).
- **NH4 and PO4 inputs in both the minimal and the full model.** `Input_NH4`/`Input_PO4` already act in `Grid.uptake`
  for any input folder. Bin Wang's version also routed organic monomer inputs (with C:N:P from `elem_ratios.csv`)
  and added substrate inputs (`SubInput`) each day; decide whether to restore those too. The `Sub`/`Mon` columns of
  `sub_mon_inputs.csv` are currently inert.

## Deferred (revisit later)

- **Fungal distribution in the minimal model.** Fungal placement probability is scaled by `max_size_b/max_size_f`
  (about 4% of cells occupied at t0 with `taxa_per_box` 1), and fungal translocation averages every living fungal
  cell each day. Use the bacterial minimal variant for sensitivity analyses until this is revisited.
- **Transporters as tracked cellular compounds.** The two transport energy costs work as designed
  (transporter production per transporter produced; active transport per monomer taken up). Treating transporters
  like HSP and osmolytes (tracked pools) would need their production computed before `uptake()`.
- **`Ferm_Cond` plays three roles**: fraction of energy-C fermented, exponent of the energy-demand curve, exponent of
  the monomer energy-value transform.
- **Leaching shape and scope.** The rate `Leaching * exp(Leaching_psi_slope * (Psi - wp_fc))` (parameters since
  2026-10-15, defaults 0.1 and 0.5) still applies in dry soil (1.4%/day at -5 MPa), where there is no drainage; only
  NH4 and PO4 leach (no dissolved organic C/N). In a 60-day minimal run about 74% of the NH4 input leached, so
  leaching largely sets how much of `Input_NH4` microbes receive. Daily and per-pulse leached N and P are now in the
  output (`LeachingSeries`, `CumLeachByPulse`).
- **Inducible osmolyte production is nearly inert under drought** (tied to uptake and leftover energy, both of which
  collapse); drought tolerance currently behaves as a constitutive trait.
- **Calibration**: `HSP_max_frac`, `Osmo_max_frac` (defaults 0.2) and `Osmo_tol_K` (default 0.05) are placeholders;
  the linear composition feedback acts as a constant tax well below the cap.

## Housekeeping

- One shared run harness for the analysis scripts (consistent parameter overrides and defaults; e.g. `n_tracker.py`
  overwrites osmolyte parameters unless `--keep-osmo`).
- Obsolete since the osmolyte-content tolerance: `analysis/full_climate_osmo_range.py`,
  `plot_full_climate_osmo_range.py`, `plot_full_climate.py`, `plot_full_climate_km25.py`,
  `build_full_input_osmo_range.py` (they set or label `Osmo_tol_ref`/`Osmo_tol_lower`, which no longer exist).
- Parameter rows that are never read: `mmrt_Km_error`, `Uptake_H_S_error`, `HSP_prod_C_cost`, `HSP_prod_P_cost`,
  `Km_Ea`.
- Performance: stop expanding grid-uniform tables to every cell; cache temperature-dependent kinetics.
