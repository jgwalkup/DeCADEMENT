# input_minimal_e1500_crich_km25_in

Copy of `input_minimal_e1500_crich_km25/` with the ongoing mineral inputs switched on: `Input_NH4` 0.02 (N per grid cell per day) and `Input_PO4` 0.002 (P per grid cell per day).
The inputs are applied in `Grid.uptake` (see `src/grid.py`); everything else is unchanged (Km slopes 2.5/2.5).
