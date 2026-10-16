# Open model issues

Items from the 2026-10-15 code review that are decided but not yet implemented, or deliberately deferred. Update this
file as items are resolved (the git log has the details of each change).

## Decided, to implement next

(none at the moment)

## Deferred (revisit later)

- **Organic monomer and substrate inputs.** `Input_NH4`/`Input_PO4` act in both the minimal and the full model (checked
  2026-10-15). Bin Wang's version also routed organic monomer inputs (C:N:P from `elem_ratios.csv`) and added
  substrate inputs (`SubInput`) daily; the `Sub`/`Mon` columns of `sub_mon_inputs.csv` are inert. Deferred.

- **Fungal distribution in the minimal model.** Fungal placement probability is scaled by `max_size_b/max_size_f`
  (about 4% of cells occupied at t0 with `taxa_per_box` 1), and fungal translocation averages every living fungal
  cell each day. Use the bacterial minimal variant for sensitivity analyses until this is revisited.
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
- **Energy deficit default.** `E_deficit_erosion` is 0 in every input folder (shortfalls floored away; the output
  `E_unpaid_Taxon` shows how much). The maintenance sensitivity runs use both settings.
- **`E_min_ratio_b/_f` biosynthesis reserve.** Besides anchoring the energy-demand curve (its intended role), energy
  below `E_min * C` is withheld from enzyme/HSP/osmolyte production. Inert at the current value 0; confirm or remove.
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
