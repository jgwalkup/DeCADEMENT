# input_minimal_e1500_osmo

Copy of `input_minimal_e1500/` for runs with drought stress (real Psi). Adds one parameter row, `Osmo_tol_ref` (0.055):
if present and > 0, the taxon's drought tolerance = min(1, total osmolyte allocation / Osmo_tol_ref), where allocation =
`Osmo_Consti_Prod` + `Osmo_Induci_Prod` (x genes per taxon). Without it the model normalises allocation between the community's
minimum and maximum, which gives tolerance 0 for a single taxon (and for any community whose taxa all have the same allocation).
Requires the opt-in change in `src/microbe.py` (backup: `src/microbe.py.bak_before_osmo_tol_ref`).
Mortality under drought: basal x (1 + (1 - tolerance) x (wp_fc - Psi) x beta), beta = 10, basal = 0.001/day.
