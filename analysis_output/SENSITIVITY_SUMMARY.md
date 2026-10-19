# Minimal-scenario sensitivity analyses: consolidated summary

Scope: single taxon, single substrate (IdealSub, C:N:P 200:10:1), one enzyme, one transporter slot; traits pegged
(`input_minimal/`); 20x20 grid of 400 identical clone cells; 365-day single pulse unless stated; `Kinetics_Model` 1 (MMRT).
Metrics: substrate degraded, biomass C, CUE (mean of last 30 days), energy:C of living biomass, cells alive.
"Control" = same run with `HSP_Consti_Prod = HSP_Induci_Prod = 0`.

## 1. Settled parameter set (E_resp 1500)

Applied as overrides in the analysis scripts; `input_minimal/parameters.csv` still holds the original baseline (column 2).

| Parameter | Original baseline | Settled | Notes |
|---|---|---|---|
| `Mon_IdealSub` E_resp (monomer_energy.csv) | 730 | **1500** (alt. 1700) | E_ferm 230, E_value 1 unchanged |
| `Maint_E_Rate_Const` | 2.0 | **2.5** | |
| `Maint_E_HSP_Ea_Intercept` | 30,000 | **45,000** | |
| `Maint_E_HSP_Ea_Slope` | -30,000 | **-90,000** | near practical limit (Ea -> 0 near -150k at HSP fraction ~0.2) |
| `Maint_E_HSP_Ea_n` | 0.2 | 0.2 | not re-varied around the final set |
| `HSP_prod_E_cost` | 600 | **300** | |
| `HSP_Induci_Prod` | 0.003 | **0.05** (was 0.3 until the fine-grid search in 3.8) | content saturates by ~0.3; above ~0.15 it throttles enzymes |
| `HSP_Consti_Prod` | 0.0015 | **0.0015** (was 0.0003) | this sets HSP content |
| `HSP_Consti_Ea`, `HSP_Induci_Ea`, `Maint_E_Rate_Tref` | 30,000 / 30,000 / 273 | unchanged | Tref deliberately not varied |
| Others (`Erat_b` 700, `Erange` 700, `Enz_Prod` 0.002, `Constit_Prod` 0.0002, `Uptake_C_cost` 0.0375, `Enz_Prod_E_cost`/`Uptake_Prod_E_cost` 600, `Active_Uptake_E_cost` 10, `Ferm_Cond` 0.1, `death_rate_bac` 0.001) | | unchanged | |

Performance of the settled set (5 seeds, flat-temperature gradient, day 365): biomass change vs control -6.5% / -7.6% at 4 / 8 C,
-3% to -8% at 12-20 C, +5.5% at 24 C, +50% (SD 5.5) at 28 C, +195% at 32 C, +305% at 36 C; mean CUE (4-32 C) 0.395, 0.40 at 4 C.
E_resp 1700 alternative: cold -8.4 / -6.8%, crossover 28 C (+12% at 28 C), CUE 0.44.
Under the real climate with constant wet Psi (-0.5 MPa): HSP -2% biomass, +0.02 CUE; with +2/+4/+6 C warming: -8.0 / -2.8 / +9.0% biomass,
CUE gain +0.05 / +0.06 / +0.07, energy:C gain +63 / +86 / +100.

The settled set is saved as files in `input_minimal_e1500/` (copy of `input_minimal/` with these changes and a README; 20x20 grid); `input_minimal/` is unchanged.

## 2. Analysis stages and the baseline each used

| Stage | Baseline | Horizon / replication |
|---|---|---|
| A. One-parameter sweeps (maintenance, enzyme/HSP/transporter production and energy costs) | original, flat 15 C | 365 d, seed 0 |
| B1. Thermal: maintenance reference; slope, n; HSP_Consti_Ea; HSP_Induci_Ea | original (constitutive-dominated HSP: 0.0015 / 0.003), no control | 365 d, seed 0, 4-40 C |
| B2. Thermal: HSP cost, slope, n, intercept, Maint_E_Rate_Const | inducible-dominated HSP (0.3 / 0.0003), E_resp 730, no control | 365 d, seed 0 |
| C. HSP-vs-control window design (cost, slope, intercept, Const, E_max, E_resp) | variable; controls included | 365 d, 5-8 seeds for final screens |
| D. 730-day gradient; real-climate and warming runs | settled set | 730 d; 6 pulses x 5 seeds |

## 3. Results

### 3.1 One-parameter sweeps at flat 15 C (stage A, original baseline)
- Collapse boundary (energy:C -> 0, expressed as multiples of baseline): Maint_E_Rate_Const ~2.2x, HSP_Consti_Prod ~3.7x, HSP_prod_E_cost ~3.3x,
  Uptake_Prod_E_cost ~27x, Enz_Prod_E_cost ~47x, Enz_Prod ~19x. Uptake_C_cost is unimodal and fragile at both ends.
- HSP production rate and HSP energy cost behave almost identically (energy effect depends mainly on rate x cost); both monotonically reduce
  degradation (24.4 -> 12.0% over the rate range) at 15 C with no interior optimum.
- Collapse mechanism: energy exhaustion (Microbes E -> 0 shuts off enzyme/HSP production); not N limitation (an earlier N-limitation reading was wrong).

**Stage A re-run at the settled set (`stageA_rerun/`, `input_minimal_e1500`, 4 / 15 / 28 / 36 C, HSP vs control, 3 seeds, 15 multipliers of baseline per parameter):**
Collapse = first multiple of baseline with >= 100 of 400 cells lost (HSP arm, 15 C; other temperatures in `stageA_summary.csv`).

| Parameter (baseline) | Old collapse (15 C) | New collapse, 15 C | Other temperatures (HSP arm) |
|---|---|---|---|
| Maint_E_Rate_Const (2.5) | 2.2x | 4.5x (control 3.5x) | 4 C 2.5x, 28 C 4.0x, 36 C 3.0x; control collapses at 1.0x at 36 C |
| Enz_Prod (0.002) | 19x (E:C < 1) | none to 40x; E:C < 1 at 18x (control 14x) = substrate exhausted | - |
| HSP_Consti_Prod (0.0015) | 3.7x | none to 16x | 4 C 12x, 36 C 12x |
| HSP_prod_E_cost (300) | 3.3x | 20x | 4 C 10x, 28 C 16x, 36 C 8x |
| Uptake_C_cost (0.0375) | unimodal | none to 10x; low end collapses at 0.1x (0.25x at 4 and 36 C) | - |
| Enz_Prod_E_cost (600) | 47x | none to 128x | 4 C 48x |
| Uptake_Prod_E_cost (600) | 27x | 48x | 4 C 24x |

- Margins to collapse widened because the settled energy yield is ~2x; the fragile ends are now maintenance (4.5x) and the low end of Uptake_C_cost (<= 0.25x).
- At 15 C degradation is insensitive (within ~5%) to the HSP and energy-cost parameters over 0-5x baseline; it is most sensitive to Enz_Prod (2x: +12%, 4x: +53%) and Uptake_C_cost (0.5x: -28%, 1.5x: +13%).
- HSP benefit window: relative to control at 28 C it persists to ~4x HSP_prod_E_cost (+14%; break-even ~5x) and ~6x HSP_Consti_Prod (+51%; negative from 10x); at 36 C break-even ~7x cost; the 4 C penalty grows from -7% (1x) to about -15% at 3-5x and -31 to -48% at 8-10x.
- Maintenance: at 28-36 C the control loses degradation steeply with Maint_E_Rate_Const (28 C -50% at 1.5x) while the HSP arm barely changes (-4%).

**Arrhenius vs MMRT kinetics at the settled set (`kinetics_arrhenius/`, `stageA_rerun_arrhenius/`; `input_minimal_e1500_arrhenius` = same inputs with `Kinetics_Model` 0):**
- Gradient (5 seeds, day 365): Arrhenius degradation rises monotonically (12% at 4 C to 93% at 40 C; biomass 4.0 -> 29.8) with no optimum, MMRT peaks at 24 C (40% degraded; biomass 13.0). Near 15-20 C the two agree (29% vs 26-29% degraded at 15 C).
- HSP effect on biomass (Arrhenius vs MMRT): 4 C -4.2 vs -6.3%; 8 C +2.3 vs -5.4%; 24 C +5.1 vs +14.2%; 28 C +28 vs +64%; 32 C +69 vs +257%; 36 C +154 vs +584%; 40 C +318 vs +463%. Same sign pattern, benefit 2-4x smaller at 28-36 C because the Arrhenius control keeps higher biomass (11.5 vs 1.2 at 36 C); the 5-10% cold penalty is not met (4 C -4.2%, 8 C +2.3%, SD 5-8).
- Stage A collapse multiples (>= 100 cells lost) are similar at 15-28 C and higher under Arrhenius at the extremes: Maint_E_Rate_Const 4.5x at 4 C (MMRT 2.5x) and 4.0x at 36 C (3.0x); HSP_prod_E_cost 20x at 4 and 36 C (MMRT 10x, 8x); HSP_Consti_Prod no collapse to 16x (MMRT 12x at 4 and 36 C); Enz_Prod_E_cost 128x at 4 C (48x); Uptake_Prod_E_cost 48-96x (24-48x).
- HSP benefit vs HSP_prod_E_cost: break-even ~4.5x at 28 C (MMRT ~4x), ~6x at 36 C (~7x); the 4 C penalty grows later (-15% at 8x vs -23% at 5x under MMRT). The benefit peak vs Maint_E_Rate_Const shifts to higher values (28 C: 3x vs 2x; 36 C: 2x vs 1.25x).

**100x100 confirmation (`grid100/`, settled set, MMRT; wall time ~55 min total, ~11.5x the cost per run of 20x20):**
- 365-day gradient (3 seeds), 730-day gradient (2 seeds), real-climate warming (3 seeds x 6 pulses) and a focused production-rate check (Induci 0.05/0.15/0.3 x Consti 0.0003/0.001/0.0015, 7 temperatures, 3 seeds) reproduce the 20x20 results: mean |difference| in the HSP effect 2.2 points (production check; 1.7 excluding 36 C); gradient differences <= ~3 points below 36 C; control survival identical (69% at 36-40 C, 365 d).
- Seed noise is 4-5x smaller: median SD of the % change 0.9 (100x100, 3 seeds) vs 4.5 (20x20, 8 seeds); gradient SD 0.1-2 vs 2-6 points below 32 C.
- Working pair (0.05, 0.0015) at 100x100: biomass change -5.7 / -3.3 / +1.1 / +16 / +66 / +256 / +589% at 4 / 8 / 16 / 24 / 28 / 32 / 36 C, CUE 0.406. The 4 C penalty is robust (SD 0.6); the 8 C penalty is about -3% (below the 5-10% target). (0.15, 0.001): -5.4 / -3.9 / +0.5 / +11 / +58 / +231 / +509%, CUE 0.404. Induci 0.3 with Consti >= 0.001 gives a positive 8 C response (+3.6 to +4.4%, SD ~1).
- Warming with the working pair (100x100 vs 20x20): +0 / +2 / +4 / +6 C -> +2.4 / -0.9 / +8.3 / +22.1% vs +2.4 / -2.2 / +5.1 / +19.5%; per-seed means agree within 0.4 points at +2, +4, +6; the +2 C dip is essentially neutral.
- 730-day: day-365 and day-730 gaps and survival at 100x100 match 20x20 (E_resp 1500: control 48% alive at 36-40 C at day 730).

### 3.2 Temperature reference (maintenance only, original baseline)
Degradation peaks 29.2% at 24 C; HSP proportion rises 0.017 -> 0.071 (4 -> 32 C), effective maintenance Ea 27,680 -> 22,520;
at 36 / 40 C energy:C collapses (delayed), HSP -> 0 and Ea returns to 30,000, which worsens the starvation.

### 3.3 Maintenance-Ea function (slope, n, intercept), HSP cost, Maint_E_Rate_Const
- Slope (0 to -90k): no effect at 4-8 C, benefit grows with temperature; at 28 C degradation 22.8 -> 30.6%; collapse day at 40 C 73 -> never.
  HSP proportion plateau is nearly independent of slope. Smaller n = stronger, saturating mitigation (28 C: 33.6% at n=0.01 vs 23.2% at n=1.6);
  slope and n are partly interchangeable (~slope/n while HSP fraction << n).
- Intercept (the HSP-free maintenance Ea) is the strongest single lever of the Ea function: at 24 C degradation 34.9% (10k) -> 5.3% (80k);
  higher intercept moves the optimum cooler and brings collapse to cooler temperatures.
- HSP cost: monotone — free HSP is best at every temperature; at 24 C 34.4% (cost 0) -> 15.6% (2400).
- Maint_E_Rate_Const: 0.5 -> 4.0 lowers 24 C degradation 38.4 -> 16.6%; Const 4 collapses even at 4 C.
- Constitutive HSP thermal Ea: raising it **hurts** at every temperature (24 C: 32.8% at Ea 0 -> 29.2% at 30k -> 7.6% at 90k) and moves collapse cooler.
- Inducible HSP thermal Ea: null (<= ~1 point). Inducible production supplies only ~2% of HSP proportion at the original baseline.

### 3.4 HSP vs no-HSP control (stage C) — reframes stages A-B
- Earlier sweeps compared HSP settings with each other only. Against a no-HSP control, the **original baseline HSP is a net loss at every temperature <= 32 C**
  (biomass -14% at 4 C, -32% at 28 C). A benefit needs cheap HSP **and** strong mitigation (cost <= ~300 with slope <= -60k) **and** enough
  maintenance stress at warm temperatures.
- Cold penalty is small and not controlled monotonically by HSP cost; mid temperatures (12-24 C) often show the largest penalty.

### 3.5 CUE
- CUE ~= (E:C / E_max)^(1+Ferm_Cond) minus enzyme C cost; E:C settles where energy income balances the energy bill (60-77% maintenance).
  Low CUE (0.2-0.3) in the original baseline is a consequence of that bill. Levers: Maint_E_Rate_Const (0.5 -> ~0.35 at 16 C), E_max (halving: +0.045),
  E_resp (730 -> 1500: 0.29 -> 0.52 at 16 C; E_ferm weak because weighted 0.1).
- Every route to higher CUE (lower maintenance, smaller E_max, higher energy yield) removes the energy stress HSP relieves, so CUE and the warm
  benefit trade off; a higher intercept plus higher E_resp restores both (settled set).

### 3.6 Window search (cold penalty 5-10%, warm benefit, CUE >= 0.35)
Seed-replicated screens: E_max route (Const 1.5, E_max 1050): cold -5%, 28 C +7%, CUE 0.355 (near miss); E_resp 1000-1200 at Const 2-2.5: CUE 0.33-0.35;
E_resp 1100 (cost 250-350): CUE 0.36-0.37 but cold penalty only -3 to -5%; intercept 40k + E_resp 1200: cold -5 to -9%, 28 C +55%, CUE 0.34-0.35;
**intercept 45k + E_resp 1500 or 1700 meets all criteria** (settled set above). Non-monotonic 4 C / 8 C responses appear in some sets (8 C flips sign between neighbouring costs).

### 3.7 730-day gradient (settled set)
Substrate is 94-98% degraded by day 730 at 16-28 C; energy:C -> 0 in both arms at 20-32 C around day 550-670 (substrate exhaustion), so late values there are not a thermal
steady state. HSP benefit persists at 36-40 C (control collapses by day ~28-220). Cold penalty fades at 4 C by day 730.

**Re-run with the working pair (Induci 0.05, Consti 0.0015; `gradient_730_v2/`, 5 seeds):** day-365 biomass change (E_resp 1500) -6.4 / -5.4 / +2.3 / -2.3 / +5.0 / +14.2 / +64 / +257 / +584 / +462% at 4 / 8 / 12 / 16 / 20 / 24 / 28 / 32 / 36 / 40 C
(E_resp 1700: -8.6 / -3.1 / -1.3 / -0.1 / +3.4 / +7.0 / +22 / +138 / +469 / +595%); the old pair had -2.5 to -8% at 12-24 C. Day 365 CUE with HSP 0.40-0.45 up to 24 C (control 0.31-0.47) and 0.26 vs 0.07 at 40 C.
At day 730 HSP runs keep all 400 cells at 36-40 C (control 192-279) and E:C stays > 0 at 40 C; absolute biomass gain at day 730 is +13.8 / +20.4 / +20.0 / +10.0 at 28 / 32 / 36 / 40 C (1500). Substrate exhaustion at 16-32 C is unchanged.
The 365-day values agree with the 8-seed fine-grid search (`hsp_prod_fine/`).

### 3.8 Production rates at E_resp 1500 (cost 350, slope -90k, intercept 40k)
HSP_Induci_Prod 0.3 -> 10 leaves HSP content unchanged. Direct probe (16 C, day 100): inducible production is N-limited, not energy-limited — it is multiplied by
total assimilated biomass material (C+N+P allocated to biomass) and an Arrhenius factor (2.08 at 16 C), then scaled by ARROEI = min(1, N supply / N demand, E supply / E demand),
shared with inducible enzyme and osmolyte production. ARROEI was 0.43 at 0.3 and 0.043 at 3.0, so realised production is the same (~0.008 per cell per day, ~32% of assimilated mass
together with constitutive). HSP C and N are not subtracted from the cell's C/N gain; only the energy cost is charged. The shared throttle also lowers inducible enzyme production
(-19% at 0.3, -24% at 3.0 at 16 C; no change at 4 or 28 C; one seed, one day), which probably contributes to the HSP-vs-control difference at mid temperatures. HSP_Consti_Prod 0.003 raises HSP proportion at 4 C to 0.05 without a clear
cold-penalty increase and removes the 28 C benefit; 0.01 overloads (HSP 12% of biomass at 4 C, 50% at 28 C; CUE 0.26).

**Re-run at the settled set (E_resp 1500, Const 2.5, cost 300, intercept 45k, slope -90k; `hsp_prod_final/`, 5 seeds, Consti 0-0.01 x Induci 0-3):**
- HSP content: constitutive-only content is ~proportional to the rate (16 C: 0.0085 at 0.0003, 0.086 at 0.003, 0.28 at 0.01); inducible-only content saturates by Induci ~0.3 (0.034 at 16 C for 0.3, 1 and 3).
- Induci above ~0.1 only adds side effects: end-of-run enzyme pool at 16 C is 0.044 (control) for Induci <= 0.1, 0.039-0.041 at 0.3, 0.037 at 1-3 (N throttle); biomass change at 16 C worsens from about -4 to -8% (Induci 0.1) to -12 to -16% (1-3), and the 24-32 C benefit is largest at Induci 0.1 (24 C +17%, 28 C +65%, 32 C +225% at Consti 0.0003) vs 0.3 (+5.5%, +50%, +195%).
- Consti 0.0001-0.003 changes the cold penalty little (4 C -3 to -9% at Induci 0.3) but raises the 32-36 C benefit up to 0.003 (+229% / +550%); 0.01 overloads (4 C -26%, CUE 0.29).
- The settled Induci 0.3 sits above the saturation point; (Induci 0.1, Consti 0.001) gives 4/8 C -6.2/-3.6%, 24 C +15%, 28 C +64%, 32 C +256%, CUE 0.415 with no enzyme loss (not replicated beyond these 5 seeds).

**Fine-grid replicate search (`hsp_prod_fine/`, 8 seeds, Induci 0.05-0.3 x Consti 0.0003-0.002, same settled set otherwise):**
- Three sets meet the strict window + CUE >= 0.35 screen: the settled (Induci 0.3, Consti 0.0003), (0.15, 0.001) and (0.05, 0.0015). The enzyme pool relative to control at 16 C is 1.00-1.03 for Induci <= 0.15, 0.96-0.97 at 0.2 and 0.91-0.93 at 0.3 (onset of the N throttle).
- Biomass change (%) at 4 / 8 / 16 / 24 / 28 / 32 / 36 C, CUE (mean <= 32 C): settled -6.4 / -5.7 / -7.1 / +5.0 / +51 / +192 / +301, 0.397; (0.15, 0.001) -6.9 / -6.3 / -3.6 / +10 / +61 / +231 / +498, 0.413; (0.05, 0.0015) -6.2 / -5.2 / -0.7 / +15 / +64 / +257 / +581, 0.416.
- The 4 C penalty tracks HSP content (Consti 0.0003 -> 0.002 at Induci 0.05: -2.9 -> -8.6%); the 8 C response is noisy (SD 2-8) and turns positive for some sets at Consti >= 0.0015.
- Recommended region if adopted: Induci 0.05-0.15, Consti 0.001-0.0015. The warming, 730-day and real-climate runs above used the old (0.3, 0.0003) pair.

**Warming scenarios re-run with the working pair (Induci 0.05, Consti 0.0015; `climate_warming_v2/`, real daily temperature + 0/2/4/6 C, Psi -0.5, 6 pulses x 5 seeds):**
biomass change with HSP +2.4 / -2.2 / +5.1 / +19.5% (old pair -1.9 / -8.0 / -2.8 / +9.0%); seed-level SEM 0.3-0.8; CUE gain +0.03 / +0.06 / +0.08 / +0.10; energy:C gain +54 / +82 / +104 / +125; HSP degrades 0.3 / 0.3 / 1.2 / 2.9 points more substrate;
the +2 C dip persists but is smaller; at +6 C the gap builds steadily through the year (+0.2% at day 150 to +19.5% at day 365).

### 3.9 Real climate (`input/climate.csv`, 6 identical years)
With real Psi (down to -6.5 MPa) the single taxon dies back every dry season (alive < 200 by day ~195; 2-5 cells left; ~7% of substrate degraded), HSP and control
identical. Cause (from code, not directly tested): drought tolerance is normalised across taxa ((allocation - min)/(max - min)); with one taxon this is 0/0 -> 0.
With constant Psi -0.5: no mortality problem; benefit emerges only at sustained warmth (+6 C offset).

### 3.10 Drought tolerance and osmolyte production (real climate)
- Tolerance is 0 for the single taxon **and for all 100 taxa of the full-model `input/`** (identical osmolyte allocation 0.055 -> min = max -> 0/0 -> 0). Osmolyte production only affects survival through this tolerance value.
- Opt-in change in `src/microbe.py` (backup `src/microbe.py.bak_before_osmo_tol_ref`): parameter row `Osmo_tol_ref`; if > 0, tolerance = min(1, (Osmo_Consti_Prod + Osmo_Induci_Prod) x genes / Osmo_tol_ref). Default behaviour verified unchanged. Input folder `input_minimal_e1500_osmo/` (Osmo_tol_ref 0.055).
- Real-climate scan (`osmo_scan/`, 20x20, 2 pulses x 3 seeds, working HSP pair; production scaled by m, tolerance = min(m, 1)): cells alive at day 365 of 400: 1.8 / 15 / 79 / 200 / 331 / 395 / 399.7 at tolerance 0 / 0.25 / 0.5 / 0.7 / 0.8 / 0.9 / 1.0; biomass C at day 365 0.01 / 0.05 / 0.27 / 0.68 / 1.22 / 1.91 / 3.04; osmolyte share of the dry-season energy bill 0 -> 12%. m = 1.5 adds cost (15%) with no benefit.
- Even at tolerance 1 biomass is 3.0 (wet constant Psi: 7.85); with free osmolytes (cost 0) 3.64, so ~17% of the loss is osmolyte cost and the rest is moisture limitation of uptake and enzymes. Energy:C reaches 0 in dry spells when osmolytes are produced (cells survive; only C/N/P minima kill).

### 3.11 Port to the full model (`input_e1500/`)
- Built from `input/` (script `analysis/build_full_input_e1500.py`): E_resp and E_ferm x2 (nearest ten), settled maintenance/HSP block, HSP production pegged at 0.0015 / 0.05, enzyme and transporter rate ranges rescaled so community-mean per-cell totals equal the minimal model's (measured 0.00193 / 0.00019 / 0.22489 vs targets 0.002 / 0.0002 / 0.225; the original full-model totals were 0.0125 / 0.0039 / 0.092).
- Test (`full_model_e1500/`, 20x20, flat wet climate, 365 days, seed 0): with the new inputs 98 of 100 taxa persist and biomass peaks at ~58 per grid cell, but energy:C falls to 0 by day ~250 (15 C) / ~200 (28 C), HSP proportion goes to 0 with it, and only 44% / 42% of the substrate is degraded (original input: 99% / 98%, E:C 4 / 0 at day 365 after substrate exhaustion). At 15 C, 81% of the cellulose, 60% of the lignin and 28% of the hemicellulose remain at day 365 (original: 2% / 0% / 0%): degradation is limited by the 6x (inducible) and 20x (constitutive) lower enzyme investment and by energy exhaustion.
- Operational note: the base conda environment's scipy 1.15.2 no longer imports on macOS 27 (`scipy.stats` fails); the `deca_dement` environment works and reproduces earlier minimal-model results exactly.

**Enzyme-rate calibration in the full model (`full_model_e1500/enz_mult_sweep.csv`, `sub_mult_test.csv`; `input_e1500` as base, 20x20, flat wet climate, 365 d, mean of 3 random communities):**
a multiplier on `Enz_Prod` and `Constit_Prod` together (1 = matched per-cell totals). At 15 C (28 C) for multipliers 0.1 / 0.25 / 0.5 / 1 / 2 / 4 / 8 / 10: substrate degraded 25 / 28 / 34 / 48 / 72 / 90 / 98 / 99% (27 / 29 / 34 / 44 / 68 / 89 / 97 / 98%); first day E:C < 100 = 139 / 149 / 170 / 236 / 312 / 301 / 276 / 261 (110 / 115 / 125 / 155 / 221 / 243 / 218 / 209).
Minimal model reference (settled set, 16 C): degraded 6 / 12 / 29% and E:C 803 / 709 / 662 at days 100 / 200 / 365. The minimal model's degradation extent at day 365 is matched only at multipliers 0.1-0.25 (28 C: 0.5-1), where energy runs out by day ~110-195; energy viability is longest at multipliers 2-4, where 72-90% of the substrate is degraded. At every multiplier the full community degrades 24-70% by day 100 (minimal: 6%) because it starts with ~20 biomass per grid cell against 345 substrate C (minimal: 1 against 100).
Scaling the initial substrate x3 / x6 (15 C, multiplier 1 / 2): E:C at day 200 = 328 / 407 (509 / 569) vs 160 (331) at x1, first day E:C < 100 = 312 / 332 (352 / 355) vs 236 (312). Recommended multiplier 2 (`input_e1500_enz2/`, community-mean totals 0.00396 inducible / 0.00039 constitutive): longest energy-viable period without exhausting the substrate (15 C day 365: 39% of cellulose, 35% of lignin, 19% of hemicellulose remaining). With 100-day pulses E:C stays 430-520 at day 100 for multipliers 1-2.

**Adopted full-model setup (`input_e1500_enz2/`): enzyme multiplier 2 and `taxa_per_box` 0.05.**
`taxa_per_box` 0.05 vs 0.20 (multiplier 2, 20x20, flat wet, 3 communities): initial biomass 5.1 vs 20.4 per grid cell, but the community grows ~9x in the first 100 days, so the effect is modest. 15 C (28 C): degraded 33 / 49 / 67% (37 / 54 / 66%) at days 100 / 200 / 365 vs 38 / 56 / 72% (41 / 59 / 68%); energy:C at day 200 / 365 = 391 / 115 (224 / 0) vs 331 / 53 (147 / 0); first day energy:C < 100 = 324 (258) vs 313 (221); 7 of 100 taxa are not placed on a 400-cell grid at the lower density (92.7 present).
Multiplier sweep at `taxa_per_box` 0.05 (15 C; production multipliers 0.1 / 0.25 / 0.5 / 1 / 2): degraded 24 / 25 / 26 / 28 / 33% at day 100 and 25 / 27 / 32 / 46 / 67% at day 365; energy:C at day 200 = 2 / 9 / 73 / 223 / 391.
Initial enzyme pool test (`init_enzyme_test.csv`; `Enz_min`/`Enz_max` scaled, production multiplier 2): x0.1 gives 5.6 / 11.9 / 26.8% degraded and energy:C 463 / 431 / 393 at days 100 / 200 / 365; x0.02 gives 2.7 / 6.8 / 17.7% and 315 / 376 / 404 (minimal model 16 C: 6.2 / 11.6 / 29.0% and 803 / 709 / 662). The early burst and the energy exhaustion are driven mainly by the initial enzyme stock (50 enzymes at 0.04-0.05 vs 1 in the minimal model). Applied (`Enz_min` 0.004, `Enz_max` 0.005).

**Final full-model setup (`input_e1500_enz2/`, adopted):** energy yields x2, settled maintenance/HSP block, matched per-cell enzyme/transporter totals with enzyme multiplier 2, `taxa_per_box` 0.05, initial enzyme pool 0.004-0.005 (rebuild with `analysis/build_full_input_e1500.py --dst ../input_e1500_enz2 --enz-mult 2 --enz-init 0.004 0.005 --taxa-per-box 0.05`).
Flat wet test (`final_flat_test.csv`, 20x20, 3 communities), degraded at days 100 / 200 / 365 and energy:C at 365 vs the minimal model (settled set, 5 seeds): 8 C 4.4 / 9.1 / 20.1% and 406 (3.6 / 5.9 / 11.4% and 705); 15-16 C 5.6 / 11.9 / 26.8% and 393 (6.2 / 11.6 / 29.0% and 662); 22 vs 20 C 6.2 / 13.1 / 29.1% and 333 (7.1 / 13.8 / 36.7% and 634); 28 C 6.2 / 12.9 / 28.2% and 249 (7.3 / 14.5 / 39.4% and 568); 36 C 5.5 / 11.1 / 22.2% and 110 (5.7 / 10.5 / 25.6% and 466). CUE (last 30 d) 0.43 / 0.41 / 0.37 / 0.33 / 0.28 at 8 / 15 / 22 / 28 / 36 C (minimal 0.43 / 0.45 / 0.44 / 0.39 / 0.32). The full community's degradation is flatter across temperature (peak ~22 C, 29%) than the minimal model's (peak ~24 C, 40%); energy:C falls below 100 only at 36 C (day ~285).

**Full model under the real climate, drought tolerance off vs on (`full_climate_osmo/`, `input_e1500_enz2` vs `input_e1500_enz2_osmo` with `Osmo_tol_ref` 0.055 = tolerance 1 for all taxa; 20x20, 6 pulses of 365 days with the model's reseeding logic, 3 communities; runner `analysis/full_climate_run.py`, figure `full_climate_osmo.png`):**
- Without tolerance the community dies back every dry season: cell-taxa alive 1,700 -> ~90 at day 250 -> ~21 at day 365 in pulse 1 (7-10 in later pulses), taxa present 93 -> ~10 (2-5 in later pulses), biomass 0.0-0.1 at day 365, and 5-9% of the substrate degraded per pulse.
- With tolerance 1 cell-taxa alive rise through the year (3,272 at day 365 in pulse 1; ~1,100-1,300 in pulses 2-6) with no in-pulse loss of taxa; biomass 14.8 (pulse 1, from 5.1) and 6.9-10.3 (pulses 2-6, from 1.3-2.3), 7-13% degraded, energy:C 186-300 at day 365 (it reaches ~0 during the driest weeks, days ~250-295, in both settings), CUE ~0.3, osmolyte production ~8% of the dry-season energy cost. Flat wet reference for the same inputs: biomass 28.6, 26.8% degraded. Minimal model with tolerance 1: biomass 3.0 (from 1.0), 11% degraded, osmolyte share 12%.
- Taxa present at the start of each pulse falls 93 -> 83 -> 72 -> 61 -> 53 -> 47: the model's reseeding draws each taxon per grid cell with probability equal to its relative abundance (about one taxon per cell, so `taxa_per_box` only sets pulse 1), so with 400 cells rare taxa are lost at each pulse; this is probably a grid-size effect (not tested at 100x100).

**Full model under the real climate with osmolyte allocation varied across taxa: relative (original, `Osmo_tol_ref` 0) vs absolute tolerance (`full_climate_osmo_range/`; `input_e1500_enz2` with `Osmo_Induci_Prod` 0.01-0.09 uniform, `Osmo_Consti_Prod` 0.005, so total allocation 0.015-0.095 with mean 0.055 = the settled value; `Osmo_tol_ref` 0, 0.055 or 0.095; 20x20, 6 pulses, 3 communities, identical communities across the three settings; runner `analysis/full_climate_osmo_range.py`, plots `analysis/plot_full_climate_osmo_range.py`):**
- Tolerance assigned (allocation 0.015 -> 0.094): relative 0 -> 1 (mean 0.50); absolute ref 0.055 0.28 -> 1 (taxa above 0.055 capped at 1, mean 0.82); absolute ref 0.095 0.16 -> 0.99 (mean 0.58).
- Day 365 (mean of 3 communities), pulse 1 / pulses 2-6: cell-taxa alive 802 / 597 (relative), 2,233 / 1,257 (ref 0.055), 918 / 692 (ref 0.095); taxa present 60 / 27, 84 / 47, 67 / 33; biomass C 3.96 / 4.36, 10.1 / 6.6, 4.46 / 3.96; energy:C 408 / 367, 270 / 286, 399 / 407; degraded 10.1 / 6.9, 11.8 / 6.9, 10.3 / 6.7 %.
- For comparison with identical allocation (0.055 in all taxa; 3.x above): relative tolerance = 0 for all taxa gave ~21 cell-taxa and ~10 taxa at day 365; absolute tolerance 1 gave 3,272 cell-taxa, 92 taxa, biomass 14.8. Varying allocation therefore rescues the original method partly (mean tolerance 0.5), and the absolute mapping with ref 0.055 gives the best outcome of the three.
- Which taxa persist: the change in cell-taxa count from day 1 to day 365 (pulse 1) rises steeply with allocation under relative tolerance (0.01x in the lowest allocation bin, 1.9x in the highest) and under ref 0.095; with ref 0.055 low-allocation taxa are lost less (bins above ~0.06 grow about 2.2x, the 0.048 bin 1.2x, the lowest two bins 0.1-0.25x).
- Selection across pulses: the abundance-weighted mean allocation on day 1 rises from 0.054 (pulse 1) to 0.082 (relative) / 0.082 (ref 0.095) / 0.074 (ref 0.055) by pulse 6; taxa present on day 1 fall from 93 to 32 / 35 / 35 (relative / ref 0.055 / ref 0.095), lower than with identical allocation (47), presumably because low-allocation taxa are removed in addition to the reseeding effect (not separated here).
- Degradation is nearly the same in the three settings (not investigated why). Differences between relative and ref 0.095 in pulses 2-6 (e.g. biomass 4.4 vs 4.0) are within what 3 communities can resolve, because the pulses start from different communities.
- Osmolyte share of total energy cost, days 200-330 (window differs from the ~8% figure above): 4.9 % (relative), 4.3 % (ref 0.055), 5.6 % (ref 0.095).
- Caveats: allocation range, constitutive/inducible split and the reference values were chosen for this test; only 3 communities on a 20x20 grid; the allocation-by-survival result is a count ratio of cell-taxa, not a per-cell survival probability (cells also gain taxa during the wet season).

**Change to the relative tolerance method (after the runs above).** `src/microbe.py` (backup `microbe.py.bak_before_osmo_tol_lower`) now computes relative tolerance as `(allocation - Osmo_tol_lower) / (max allocation - Osmo_tol_lower)`, clipped to [0, 1], with the optional parameter `Osmo_tol_lower` defaulting to 0 (the original method used the community minimum as the lower bound). `Osmo_tol_ref` > 0 still takes precedence. Consequences: identical allocation across taxa (single taxon, or the full `input/`) now gives tolerance 1 by default, and the tolerance-0 results in 3.10, 3.11 and the `ref 0` runs of `full_climate_osmo_range/` (original method) are reproduced with `Osmo_tol_lower` set to a negative value (e.g. -1), which selects the community minimum as the lower bound (the original method; verified identical to the pre-change code, and the rerun in `full_climate_osmo_range_v2/` reproduces the `ref 0` runs of `full_climate_osmo_range/` exactly). `Osmo_tol_ref` defaults to 0 (absent row = 0). With allocation 0.015-0.095 and the new default, relative tolerance equals 0.16 -> 1, i.e. the same mapping as `Osmo_tol_ref` 0.095.

**Rerun of the varied-osmolyte comparison with the new relative method (`full_climate_osmo_range_v2/`; same inputs, 3 communities, 6 pulses; conditions: original relative `Osmo_tol_lower` -1, new relative `Osmo_tol_lower` 0, absolute `Osmo_tol_ref` 0.055, absolute 0.095; the original-method and absolute-0.055 runs reproduce `full_climate_osmo_range/` exactly):**
- Tolerance assigned (allocation 0.015 -> 0.094): original 0 -> 1 (mean 0.50); new relative 0.16 -> 1 (mean 0.58, nearly the same line as absolute 0.095); absolute 0.055 0.28 -> 1 (mean 0.82).
- Day 365, pulse 1 / pulses 2-6: cell-taxa alive 802 / 597 (original), 918 / 666 (new relative), 2,233 / 1,257 (abs 0.055), 918 / 692 (abs 0.095); taxa present 60 / 27, 69 / 31, 84 / 47, 67 / 33; biomass C 3.96 / 4.36, 4.35 / 5.92, 10.1 / 6.6, 4.46 / 3.96; degraded 10.1 / 6.9, 10.2 / 8.8, 11.8 / 6.9, 10.3 / 6.7 %.
- New relative vs original: modestly better in pulse 1 (about 14% more cell-taxa, 9 more taxa, biomass +10%) because the lowest-allocation taxa now keep some tolerance. In pulses 2-6 its biomass (5.9 vs 4.4, peak 7.5 vs 5.5) and degradation (8.8 vs 6.9%) are higher, but the pulses start from different communities (day-1 biomass 2.3 vs 1.7), so this is not attributed to tolerance (3 communities).
- Selection on allocation across pulses is about as strong with the new relative method as with the original (mean allocation on day 1 of pulse 6: 0.082 / 0.082 / 0.074 / 0.082 for original / new relative / abs 0.055 / abs 0.095); taxa present on day 1 of pulse 6: 32 / 36 / 35 / 35. Osmolyte share of energy cost days 200-330: 4.9 / 5.4 / 4.3 / 5.6 %.
- Only absolute tolerance with `Osmo_tol_ref` at or below the mean allocation (0.055) changes the outcome substantially; the three other mappings give similar results.

**Full-model HSP vs no-HSP control gradient (`full_gradient_730/`; `input_e1500_enz2`, flat climate Psi -0.5 / FM 50, 4-40 C in 4 C steps, one 730-day pulse without resupply, 20x20, 5 communities; control = `HSP_Induci_Prod` and `HSP_Consti_Prod` 0; runner `analysis/full_gradient_730.py`; minimal-model reference `gradient_730_v2/`, E_resp 1500, 5 seeds):**
- Day 365, biomass change with HSP vs control (full / minimal), %: 4 C -2.3 / -6.4; 8 C -3.1 / -5.4; 12 C -2.0 / +2.3; 16 C +3.5 / -2.3; 20 C +10.6 / +5.0; 24 C +32.7 / +14.2; 28 C +82.9 / +64; 32 C +176 / +257; 36 C +239 / +584; 40 C +244 / +462. In the full model HSP biomass exceeds the control in 5 of 5 communities from 16 C upward (day 365) and in 1 of 5 at 4-12 C, so the cold penalty is small (2-3%) and the crossover is 12-16 C.
- Day 365, substrate degraded with HSP / control, % of initial C: 4 C 16.1 / 16.2; 12 C 27.7 / 28.0; 20 C 33.5 / 31.7; 28 C 33.0 / 19.1; 36 C 26.5 / 5.4; 40 C 20.9 / 3.5 (minimal: 36 C 25.6 / 7.7). Energy:C in the control reaches 0 at 32 C and above by day 365 (HSP 134-414 across all temperatures). HSP content (HSP_prop) 0.02-0.05 (minimal 0.02-0.11).
- CUE (mean of last 30 days): the HSP gain is much smaller than in the minimal model, +0.13 at most at 28 C (day 365) against +0.29 at 32 C in the minimal model, and about 0 at day 730 (minimal +0.17 at 28 C), where energy:C is 0 at 20 C and above in both HSP and control.
- Day 730: the HSP advantage in biomass grows (+94% at 24 C, +188% at 28 C, +403% at 36 C) and degradation is 63.7 / 54.6 / 40.6 % at 12 / 28 / 36 C with HSP vs 63.2 / 19.6 / 5.4 % in the control, while the cold penalty becomes -5.5 to -6.0% at 4-8 C (HSP biomass below the control in 5 of 5 communities at 4 C, 4 of 5 at 8 C, 3 of 5 at 12 C at day 730). The minimal model has exhausted the substrate by day 730 under HSP at 20-36 C (90-98% degraded), the full model has not (40-65%), so day-730 comparisons between the two are affected by substrate depletion in the minimal model; day 365 is the cleaner comparison.
- The no-HSP control is less depressed in the warm full model than in the minimal model (biomass 7.6 vs 1.2 at 36 C, day 365; 3.5-6% degraded in both), so the HSP-to-control biomass ratio is smaller (3.4x vs 6.8x).
- Caveats: one pulse without resupply (energy:C falls to 0 at 16-20 C and above by day 730 as the labile pools are used; the experiments planned with resupply every 365 days will behave more like day 365); 5 communities on 20x20; flat wet climate (no drought); HSP and control differ only in HSP production rates, so the control also removes the HSP energy cost.

**Full-model HSP vs control gradient with 365-day substrate resupply (`full_gradient_resupply/`; same inputs and temperatures as the no-resupply gradient above; 3 consecutive pulses of 365 days at constant temperature, substrate/monomers/enzymes reset and community reseeded between pulses by the model's repopulation logic; 5 communities; runner `analysis/full_gradient_resupply.py`):**
- Pulse 1 reproduces the no-resupply day-365 values exactly (same seeds and random stream), e.g. biomass change +83 % at 28 C, +239 % at 36 C.
- Biomass change with HSP vs control on day 365 (ratio of community means), pulse 1 / 2 / 3, %: 4 C -2.3 / -1.5 / -24.7; 8 C -3.1 / -7.1 / +8.8; 12 C -2.0 / -26.1 / -23.5; 16 C +3.5 / +12.9 / +24.1; 20 C +10.6 / +15.1 / +23.6; 24 C +33 / +73 / +77; 28 C +83 / +105 / +135; 32 C +176 / +231 / +387; 36 C +239 / +279 / +435; 40 C +244 / +297 / +357. HSP biomass is higher in all 5 communities from 24 C upward in pulses 2 and 3 (from 16 C upward in pulse 1).
- Below 20 C the sign is not reliable in pulses 2-3: per-community changes at 12 C range -47 to +2 % (pulse 2) and -38 to +20 % (pulse 3), at 16 C -52 to +268 %. HSP and control start pulses 2-3 from different reseeded communities (day-1 biomass differs by 15 % in pulse 2 and 40 % in pulse 3 between paired runs), which pulse 1 does not have.
- Start of pulse: biomass 5.0 (pulse 1) -> 1.4 (pulse 2) -> 2.1-2.4 (pulse 3); cell-taxa 1,690 -> ~395 -> ~395; taxa present 92 -> ~83 -> ~71. Substrate degraded per pulse is lower in pulse 2 than in pulse 1 (e.g. 25-28 % vs 33-34 % at 20-24 C with HSP) and higher in pulse 3 at 16-28 C (36-43 % with HSP).
- Energy:C on day 365 stays 134-473 with HSP at all temperatures in all pulses; in the control it is 0 at 32 C and above (pulse 1) and 0-44 at 28 C and above in pulses 2-3. CUE advantage of HSP at 24-28 C: +0.08 to +0.14 in all three pulses (at 28 C: 0.34 vs 0.22, 0.38 vs 0.24, 0.34 vs 0.23 in pulses 1 / 2 / 3); cold CUE penalty about -0.02 or smaller.
- Cell-taxa alive on day 365 at 36-40 C: HSP about 1,900-2,200 vs control about 800 in pulses 2 and 3.
- Conclusion for the settled set: the warm-temperature benefit of HSP is robust to resupply (it grows over the pulses); the cold penalty is at most a few percent in pulse 1 and not distinguishable from reseeding noise in pulses 2-3 with 5 communities.

**Full-model MMRT vs Arrhenius (`Kinetics_Model` 1 vs 0) with 365-day resupply (`kinetics_comparison_resupply/`; Arrhenius runs in `full_gradient_resupply_arrhenius/`; `input_e1500_enz2` with only `Kinetics_Model` changed in `runtime.txt`, so each model uses its own parameters (MMRT Cp/HT0/ST0 vs Arrhenius Ea), which are not matched to each other; 3 communities (seeds 0-2), 4, 12, 20, 28, 36, 40 C, 3 pulses, HSP and control; the MMRT numbers here are the same three communities from `full_gradient_resupply/`, so they differ slightly from the 5-community values above; runner `analysis/full_gradient_resupply.py --kinetics 0`, comparison `analysis/compare_kinetics_resupply.py`):**
- Day 365, pulse 1 (MMRT / Arrhenius): biomass with HSP 20.6 / 24.4 at 4 C, 29.7 / 28.6 at 20 C, 27.6 / 30.4 at 28 C, 22.1 / 33.0 at 36 C, 18.1 / 34.7 at 40 C; control 21.0 / 25.0, 26.8 / 25.7, 13.9 / 15.4, 6.8 / 8.8, 6.0 / 7.1. At 12-20 C the two are within about 5 % in biomass, degradation and CUE.
- HSP biomass effect, pulse 1 (MMRT / Arrhenius), %: 4 C -2.0 / -2.2; 12 C -1.9 / -0.1; 20 C +10.9 / +11.2; 28 C +99 / +97; 36 C +224 / +276; 40 C +201 / +387. Pulse 2 / pulse 3 at 36 C: +213 / +238 (MMRT) vs +404 / +899 (Arrhenius); at 40 C: +270 / +275 vs +457 / +625. The three-community Arrhenius values at 36-40 C in pulses 2-3 are large and variable.
- Shape of the temperature response: MMRT HSP biomass peaks at 20-28 C and falls to 18-22 at 36-40 C (substrate degraded 17-22 %); under Arrhenius HSP biomass keeps rising with temperature (33-35 at 36-40 C, 35-38 % degraded), i.e. there is no thermal optimum for the HSP community, and the control also does a little better at 36-40 C than under MMRT but still collapses (energy:C 0).
- CUE: the HSP advantage at 28 C is similar (+0.09 MMRT / +0.11 Arrhenius in pulse 1), larger under Arrhenius at 36-40 C (0.34-0.35 vs 0.27-0.28 with HSP; control 0.24-0.25 vs 0.28-0.29). In the cold both models show a CUE penalty of about 0.02 or less.
- Cold side after resupply is noisy in both models (pulse 3, 12 C: -25 % MMRT / -10 % Arrhenius).
- As in the minimal-model comparison above, Arrhenius has no thermal optimum; unlike the minimal model (where the Arrhenius HSP benefit at 28-36 C was 2-4x smaller than under MMRT), the HSP benefit here is the same at 28 C and larger under Arrhenius at 36-40 C. The Arrhenius and MMRT parameter sets have not been matched; the agreement at 12-20 C is a property of the two parameter sets, not a calibration.

**Full model: enzyme Vmax-Km slope sweep across temperature, MMRT (`mmrt_Vmax_Km_slope`) vs Arrhenius (`Vmax_Km`) (`km_sweep/`; `input_e1500_enz2`, flat climate Psi -0.5 / FM 50, one 365-day pulse, HSP as in the input, 4-40 C in 4 C steps, slopes 7 (current) / 3 / 1 / 0.3 / 0.1, 3 communities, 20x20; runner `analysis/full_km_sweep.py`):**
- The slope scales Km only (MMRT: Km = MMRT(Km parameters) x substrate-specific factor x slope + intercept; Arrhenius: Km0 = Vmax0 x slope + intercept, then its own Arrhenius term); Vmax is unchanged. At the current slope the median Km is about 150-170 (20-28 C) against 147 (cellulose), 86 (hemicellulose), 49 (lignin) and 2-13 (the other pools) of substrate C per cell, so the large pools are near half-saturation and the small pools far below Km.
- Substrate degraded by day 365 at the MMRT optimum (24 C) / range over 12-32 C, MMRT: slope 7 29 % / 24-29 %; 3 69 % / 56-69 %; 1 88 % / 77-88 %; 0.3 91 % / 82-91 %. Arrhenius (no optimum; rises to 40 C): slope 7 19-38 % (4-40 C); 3 44-83 %; 1 65-96 %; 0.3 73-98 %. The MMRT optimum stays at 24 C for every slope.
- From slope 7 to 3 degradation scales about with 1/slope (x2.3), i.e. at the current setting the enzymes operate in the Km-limited (first-order) regime. Between slope 1 and 0.3 it barely changes (+4 percentage points at 24 C): starch, proteins, chitin and organic P are 100 % degraded, but cellulose (78-83 % at 24 C) and lignin (84-91 %) are not, which is consistent with Vmax x enzyme limiting the large recalcitrant pools once Km is small. Whether more Vmax or more enzyme would remove that limit was not tested.
- Costs of faster degradation: biomass rises from 18-30 (slope 7) to 62-77 at slopes 3 and 1 over 12-32 C (the flat wet reference was 28.6); energy:C at day 365 falls with the slope at warm temperatures (MMRT 24 C: 303 / 248 / 46 / 6 for slopes 7 / 3 / 1 / 0.3; Arrhenius 28 C: 290 / 174 / 0 / 0), and CUE falls to about 0.25 at slope 1 over 16-40 C (MMRT) so at slopes <= 1 the community is energy-starved by day 365 at warm temperatures. Slope 3 keeps energy:C at 170-473 over 4-32 C in MMRT (120-158 at 36-40 C) and CUE 0.25-0.44.
- Slope 0.1 is not usable: in 55 of the 300 runs (25 MMRT and 28 Arrhenius at slope 0.1, plus 1 at 0.3 for each model) cellulose C becomes NaN once a cell has used all its cellulose and lignin (0/0 in `Grid.degradation`, `ss7 / (ss7 + cellulose_C)`, a documented existing behaviour), and from the next day the cellulose monomer pool is NaN in every cell. Those runs were excluded; points with fewer than 2 clean communities are not plotted.
- Caveats: Km only (the sweep cannot raise Vmax separately); one pulse, flat wet climate, 3 communities; uptake Km (`Uptake_Vmax_Km`) was not varied; the earlier calibration of enzyme production, energy yields and initial enzyme pool was done at slope 7 and has not been revisited at the faster degradation.

**Finer Vmax-Km slope scan between 1 and 3 (`km_sweep_1to3/`: new runs of slopes 1.5, 2 and 2.5 in `km_sweep_fine/` merged with slopes 1 and 3 of `km_sweep/`; same inputs, climate, temperatures and 3 communities; no NaN runs in this range):**
- Substrate degraded by day 365, MMRT at 24 C (optimum) / range over 12-32 C: slope 3 69 % / 56-69; 2.5 75 % / 62-75; 2 79 % / 67-79; 1.5 84 % / 73-84; 1 88 % / 77-88. Arrhenius at 4 / 24 / 40 C: slope 3 44 / 71 / 83 %; 2.5 49 / 77 / 88; 2 55 / 83 / 91; 1.5 60 / 87 / 93; 1 65 / 90 / 96. Degradation rises smoothly as the slope falls, with the MMRT optimum at 24 C throughout.
- More than 50 % degraded across the temperature range (MMRT): slope 3 from 12 to 36 C; slope 2.5 from 12 to 40 C (49 % at 8 C, 54 % at 40 C); slope 2 from 8 to 40 C (55 % at 8 C, 59 % at 40 C); at 4 C it stays below 50 % for all slopes (32-48 %). Arrhenius is above 50 % from 12 C at slope 3, from 8 C at slope 2.5 and from 4 C at slope 2 and below.
- Energy:C on day 365 (MMRT 24 C / 28 C): slope 3 248 / 199; 2.5 195 / 153; 2 132 / 84; 1.5 77 / 51; 1 46 / 22 (Arrhenius 24 C / 28 C: 254 / 174; 178 / 106; 109 / 30; 38 / 12; 7 / 0). CUE at 24 C (MMRT): 0.30, 0.27, 0.25, 0.23, 0.25. Biomass at the optimum 67 -> 76.
- Per-substrate degradation at MMRT 24 C: cellulose 56, 63, 68, 73, 78 % and lignin 57, 65, 71, 77, 84 % for slopes 3, 2.5, 2, 1.5, 1; the other pools are 79-96 % at slope 3 (hemicellulose 79 %) and 90-99 % at slope 2 (hemicellulose 90 %).
- Trade-off: each 0.5 step down in slope adds about 4-6 percentage points of degradation at the optimum and removes about 30-60 units of energy:C at 24 C; energy:C at 28 C falls below 100 at slope 2 or lower in both models (MMRT 84, Arrhenius 30).

**Uptake Vmax-Km slope (`Uptake_mmrt_Vmax_Km_slope` MMRT / `Uptake_Vmax_Km` Arrhenius) crossed with enzyme slope 2, 3, 4 (`uptake_sweep/`; `input_e1500_enz2`, flat climate, 365 days, HSP as in the input, uptake slopes 10 / 5 (current) / 2.5 / 1, temperatures 4, 12, 20, 24, 28, 40 C, 3 communities, 432 runs, no NaN runs; runner `analysis/full_uptake_sweep.py`; the uptake-5 cells with enzyme slopes 2 and 3 reproduce the earlier `km_sweep_1to3` values exactly):**
- Uptake Km = uptake Vmax x slope, and the monomer pools (0.001-5 C per cell) are far below Km (median 117-231), so uptake per unit transporter is about pool / slope: the uptake slope acts as a rate divisor, independent of uptake Vmax. At the current slope 5 the monomer pool is 1.0-2.1 C per cell on day 365 (MMRT 24 C); it is 3.8-6.9 at slope 10, 0.3-0.5 at 2.5 and 0.06-0.09 at 1, so uptake stops limiting below about 2.5.
- Substrate degraded by day 365, MMRT 24 C (enzyme slope 2 / 3 / 4): uptake 10 74.6 / 60.9 / 47.5 %; 5 79.3 / 68.6 / 56.8; 2.5 81.1 / 72.4 / 62.6; 1 83.6 / 75.3 / 66.3. Arrhenius 28 C: 10 81.6 / 67.3 / 53.5; 5 84.8 / 75.3 / 63.1; 2.5 86.1 / 78.2 / 67.7; 1 88.0 / 82.6 / 71.7. A 10-fold change in the uptake slope is worth 9, 14 and 19 percentage points at enzyme slopes 2, 3 and 4 (MMRT 24 C), against 22.5 points for enzyme slope 4 -> 2 at uptake 5; the effect saturates below uptake 2.5 (2-4 points from 2.5 to 1).
- Biomass (MMRT 24 C) rises with faster uptake up to uptake 2.5 (enzyme 3: 57 / 67 / 73 / 73 for 10 / 5 / 2.5 / 1) and flat or lower at 1 for enzyme slope 2 (78 -> 74). Energy:C falls with faster uptake (enzyme 3: 370 / 248 / 175 / 155) and CUE falls from 0.37 to 0.27; at enzyme slope 2 energy:C at 24 C is 265 / 132 / 71 / 95 for uptake 10 / 5 / 2.5 / 1 (MMRT) and 177 / 30 / 13 / 51 at 28 C (Arrhenius).
- Combinations are not equivalent for energy: at MMRT 24 C enzyme slope 2 with uptake 10 degrades 74.6 % with energy:C 265 and CUE 0.31, versus 68.6 % with 248 and 0.30 for enzyme 3 / uptake 5 (current) and 72.4 % with 175 and 0.27 for enzyme 3 / uptake 2.5; enzyme 4 with uptake 1 gives 66.3 % with 222 and 0.29, close to enzyme 3 / uptake 5. In Arrhenius at 28 C enzyme 2 / uptake 10 gives 81.6 % with energy:C 177 against 75.3 % with 174 for enzyme 3 / uptake 5. For a given amount of degradation, a faster enzyme with slower uptake keeps more energy than a slower enzyme with faster uptake; this has been measured, the mechanism has not been tested.
- Temperature shape does not change with the uptake slope (MMRT optimum 24 C at every combination; Arrhenius rises to 40 C).
- Caveats: Km only; uptake Vmax scale not varied; 3 communities, one pulse, flat wet climate; the accumulated monomer pool at uptake 10 (4-7 C per cell, 1-2 % of the initial substrate C) is counted as degraded substrate.

**Enzyme and uptake Vmax-Km slopes both 2.5: HSP vs control gradient with 365-day resupply (`input_e1500_enz2_km25/` = `input_e1500_enz2` with `mmrt_Vmax_Km_slope`, `Vmax_Km` 7 -> 2.5 and `Uptake_mmrt_Vmax_Km_slope`, `Uptake_Vmax_Km` 5 -> 2.5, built by `analysis/build_full_input_km.py`; results in `full_gradient_resupply_km25/` (MMRT, 10 temperatures, 5 communities, 3 pulses, 200 runs), `full_gradient_resupply_km25_arrhenius/` (Arrhenius, 6 temperatures, 3 communities, 36 runs) and `kinetics_comparison_resupply_km25/`; no NaN):**
- MMRT, pulse 1, day 365 (slope 2.5/2.5 vs old 7/5): substrate degraded with HSP 42 / 59 / 72 / 78 / 81 / 81 / 81 / 77 / 72 / 62 % at 4 / 8 / 12 / 16 / 20 / 24 / 28 / 32 / 36 / 40 C (old 16-34 %); the control 42 / 59 / 72 / 77 / 79 / 74 / 62 / 42 / 21 / 10 %. Biomass with HSP 55-82 (old 21-34); energy:C with HSP 476 / 479 / 402 / 260 / 149 / 82 / 55 / 54 / 72 / 93 (old 134-452), control 0-2 from 24 C upward; CUE with HSP 0.42 at 4 C falling to 0.21-0.24 from 24 C.
- HSP effect on biomass (pulse 1 / 2 / 3, %): 4 C -3 / -11 / -9; 12 C 0 / -2 / +3; 20 C +6 / +15 / +15; 28 C +42 / +59 / +63; 36 C +224 / +276 / +273; 40 C +373 / +515 / +576 (old slopes pulse 1: 4 C -2, 20 C +11, 28 C +83, 36 C +239, 40 C +244). The relative HSP advantage at 20-32 C is smaller than at the old slopes because the control now degrades more and lives longer there (control biomass 55 vs 18 at 28 C, pulse 1; degradation 62 % vs 19 %); the advantage at 36-40 C is as large or larger (degradation difference with HSP at 36 C 51 vs 21 percentage points). Degradation with HSP is flat at 81 % from 20 to 28 C; without HSP it peaks at 79 % at 20 C.
- Later pulses reach higher degradation (HSP 88-93 % at 16-28 C in pulse 3) and lower end-of-pulse energy:C (MMRT 24 C pulse 3: 35, 28 C: 10); biomass in pulses 2-3 is 96-118 at 16-28 C. The substrate is nearly used up in the second half of each pulse at warm temperatures (substrate remaining 7-12 % at the end of pulse 3 at 20-28 C).
- Arrhenius at 2.5/2.5 (3 communities, same communities for the MMRT numbers): HSP biomass change pulse 1 +34 % (28 C), +104 % (36 C), +162 % (40 C) against MMRT +49, +261, +386 %; degradation with HSP 83 / 87 / 89 % at 28 / 36 / 40 C (pulse 1) and 87 / 97 / 99 % in pulse 3, control 66 / 47 / 37 %; the Arrhenius control survives better at warm temperatures than the MMRT control (biomass 59 / 40 / 32 vs 49 / 18 / 11 at 28 / 36 / 40 C).
- Energy:C at day 365 is low in all warm HSP runs (53-93 at 24-40 C in pulse 1, 10-83 in pulse 3) and 0 in the control from 28 C: the community is energy-limited by the end of the year at warm temperatures; the earlier calibration (enzyme multiplier, initial enzyme pool, energy yields) was done at slope 7/5 and has not been revisited.

**Real climate with slopes 2.5/2.5 and relative tolerance from 0 to the maximum allocation (`full_climate_km25_osmo0/`; `input_e1500_enz2_km25` with rows `Osmo_tol_lower,0` and `Osmo_tol_ref,0` added to `parameters.csv`; `input/climate.csv`, MMRT, 6 pulses of 365 days with the model's reseeding, 5 communities, 20x20; runner `analysis/full_climate_run.py`, figure `full_climate_km25_osmo0.png`; the previous run is `full_climate_osmo` with slopes 7/5 and `Osmo_tol_ref` 0.055, seeds 0-2):**
- All taxa have the same osmolyte allocation (0.005 + 0.05 = 0.055), so the relative scale 0 to maximum gives tolerance 1 for every taxon (verified from `TaxDroughtTol`), i.e. no added drought mortality, the same tolerance as the previous run. It is only below 1 for taxa whose allocation is below the community maximum, which needs a min-max range in `Osmo_Consti_Prod` / `Osmo_Induci_Prod`.
- Day 365 by pulse (1 / 2 / 3 / 4 / 5 / 6), new vs previous: substrate degraded 50 / 47 / 57 / 64 / 63 / 65 % (previous 12.8 / 6.9 / 7.7 / 8.3 / 9.3 / 9.1); biomass C 51 / 52 / 63 / 73 / 76 / 83 (previous 14.8 / 6.9 / 7.7 / 8.8 / 9.6 / 10.3); cell-taxa alive 7,757 / 3,689 / 2,697 / 2,274 / 1,852 / 1,571 (previous 3,272 / 1,206 / 1,284 / 1,192 / 1,204 / 1,136); energy:C 266 / 343 / 256 / 153 / 73 / 58 (previous 186 / 300 / 297 / 245 / 238 / 193); CUE 0.32 / 0.36 / 0.29 / 0.26 / 0.26 / 0.28.
- Within a pulse most degradation occurs before about day 200; it stalls during the driest weeks (about days 200-295, Psi near -4 to -6 MPa; biomass flat, no loss of cells because tolerance is 1) and resumes after day 295. Energy:C falls below 5 on 104 / 38 / 36 / 56 / 112 / 120 days in pulses 1-6 (seed 0; previous 65 / 51 / 64 / 71 / 73 / 140).
- Pulses 4-6 start from larger communities (day-1 biomass 1.8 in pulse 2, 11.5 in pulse 6, previous 1.3-2.3) and reach 63-65 % degradation, but with lower end-of-year energy:C; taxa present at the start of each pulse falls 92 / 83 / 67 / 51 / 39 / 33 (previous 93 / 83 / 72 / 61 / 53 / 47), so richness is lost faster than before with tolerance identical across taxa (reseeding follows abundance; not separated from competition here).
- Osmolyte production takes 11.8 % of total energy cost over days 200-330 (previous 4.6 %).
- Caveats: 20x20 grid, 5 communities, MMRT only; osmolyte allocation identical across taxa; the energy calibration was set at slopes 7/5.

**Real climate, slopes 2.5/2.5, osmolyte allocation range so that tolerance varies (`full_climate_km25_osmorange/`; `input_e1500_enz2_km25_osmorange` = `input_e1500_enz2_km25` with `Osmo_Induci_Prod` 0.01-0.09 (constitutive 0.005), built by `analysis/build_full_input_osmo_range.py`; `Osmo_tol_lower` 0 and `Osmo_tol_ref` 0, so tolerance = allocation / maximum allocation; 6 pulses, 5 communities, 20x20; runner `analysis/full_climate_osmo_range.py` (now overwrites existing tolerance rows), plots `analysis/plot_full_climate_km25_range.py`):**
- Allocation 0.015-0.095 (mean 0.055), tolerance 0.16-1 (mean 0.58), verified from `TaxDroughtTol`.
- Day 365 by pulse (1 / 2 / 3 / 4 / 5 / 6), range vs identical allocation (tolerance 1, `full_climate_km25_osmo0`): substrate degraded 42 / 38 / 46 / 52 / 55 / 58 % vs 50 / 47 / 57 / 64 / 63 / 65; biomass C 24 / 29 / 36 / 43 / 44 / 47 vs 51 / 52 / 63 / 73 / 76 / 83; cell-taxa alive 3,147 / 1,768 / 1,493 / 1,326 / 1,252 / 1,236 vs 7,757 / 3,689 / 2,697 / 2,274 / 1,852 / 1,571; taxa present 80 / 55 / 40 / 30 / 23 / 16 vs 92 / 81 / 66 / 50 / 38 / 32; energy:C 489 / 448 / 337 / 327 / 294 / 255 vs 266 / 343 / 256 / 153 / 73 / 58; CUE 0.47 / 0.44 / 0.39 / 0.39 / 0.39 / 0.39 vs 0.32 / 0.36 / 0.29 / 0.26 / 0.26 / 0.28.
- Cell-taxa rise to about 5,800 by day 140 in pulse 1 and fall to about 2,000 by day 295 during the dry period (biomass 32 -> 13); taxa present drops during the dry season (pulse 1 92 -> 81 by day 295; pulse 2 80 -> 55), whereas with identical tolerance 1 no taxa and no cells are lost within a pulse. Degradation still continues slowly in the dry period.
- Selection: the abundance-weighted mean allocation on day 1 rises from 0.055 (pulse 1) to 0.062 / 0.069 / 0.076 / 0.079 / 0.081 (pulses 2-6); taxa present on day 1 fall 92 / 80 / 62 / 47 / 34 / 25 (identical allocation 92 / 83 / 67 / 51 / 39 / 33). Cell-taxa in the highest allocation bins grow 3-6x through a pulse; the lowest bins end at 0-0.3x of the day-1 number.
- Energy:C is below 5 on 22 / 45 / 49 / 40 / 106 / 124 days in pulses 1-6 (seed 0; identical allocation 104 / 38 / 36 / 56 / 112 / 120). Osmolyte production is 10.2 % of energy cost over days 200-330 (identical 11.8 %).
- Interpretation: with the same total osmolyte spend, varying allocation lowers the community tolerance from 1 to 0.58 on average, so the dry season kills cells and taxa again, the community is about half the size, and the survivors (higher allocation) have more energy per unit biomass.
- Caveats: allocation range (0.01-0.09 inducible, constitutive fixed) chosen for this test; relative scale to the community maximum; 20x20, 5 communities, MMRT.

**Real climate, slopes 2.5/2.5, ranges for both constitutive and inducible osmolyte production (`full_climate_km25_osmorange2/`; `input_e1500_enz2_km25_osmorange2` = `input_e1500_enz2_km25` with `Osmo_Consti_Prod` 0.002-0.005 and `Osmo_Induci_Prod` 0.04-0.08; `Osmo_tol_lower` 0, `Osmo_tol_ref` 0; 6 pulses, 5 communities, 20x20; figure `km25_range_community.png` with the identical-allocation and the inducible-range runs side by side):**
- Total allocation 0.043-0.084 (mean 0.0635), tolerance 0.52-1 (mean 0.76; verified = allocation / maximum).
- Day 365 by pulse (1 / 2 / 3 / 4 / 5 / 6): substrate degraded 45 / 42 / 49 / 54 / 61 / 63 % (inducible-only range 0.015-0.095, tolerance mean 0.58: 42 / 38 / 46 / 52 / 55 / 58; identical allocation: 50 / 47 / 57 / 64 / 63 / 65); biomass C 32 / 36 / 40 / 46 / 52 / 52 (24 / 29 / 36 / 43 / 44 / 47; 51 / 52 / 63 / 73 / 76 / 83); cell-taxa alive 4,332 / 2,066 / 1,607 / 1,362 / 1,231 / 1,089 (3,147 / 1,768 / 1,493 / 1,326 / 1,252 / 1,236; 7,757 / 3,689 / 2,697 / 2,274 / 1,852 / 1,571); taxa present 90 / 68 / 50 / 40 / 33 / 24 (80 / 55 / 40 / 30 / 23 / 16; 92 / 81 / 66 / 50 / 38 / 32); energy:C 449 / 441 / 349 / 319 / 280 / 249 (489 / 448 / 337 / 327 / 294 / 255; 266 / 343 / 256 / 153 / 73 / 58); CUE 0.45 / 0.43 / 0.39 / 0.39 / 0.36 / 0.36 (0.47-0.39; 0.32-0.26).
- Pulse 1 still rises to about 6,100 cell-taxa by day 170 and falls to about 2,900 by day 295 (inducible-only range: 5,800 -> 2,000); biomass 35 -> 19; the dry-season loss of taxa is smaller (pulse 1 92 -> 90, pulse 2 81 -> 68).
- Selection on allocation is weaker: the abundance-weighted mean allocation on day 1 rises from 0.0634 (pulse 1) to 0.073 (pulse 6) (all-taxa mean 0.0635; inducible-only range 0.055 -> 0.081). Taxa present on day 1: 92 / 81 / 61 / 50 / 38 / 30.
- Energy:C below 5 on 2 / 28 / 21 / 38 / 53 / 110 days in pulses 1-6 (seed 0; identical allocation 104 / 38 / 36 / 56 / 112 / 120). Osmolyte production is 9.2 % of energy cost over days 200-330 (identical allocation 11.8 %, inducible-only range 10.2 %).
- Across the three settings biomass, taxa present and degradation follow the mean tolerance (1, 0.76, 0.58) in the order identical > this range > inducible-only range, and energy:C follows the reverse order; the count of cell-taxa alive in pulses 5-6 does not follow it (identical 1,852 / 1,571, this range 1,231 / 1,089, inducible-only range 1,252 / 1,236).

**Repeat of the osmolyte-range run (`input_e1500_enz2_km25_osmorange2`) on a 50x50 grid with 12 pulses (`full_climate_km25_osmorange2_50x50_12p/`; 5 communities, MMRT, `Osmo_tol_lower` 0, `Osmo_tol_ref` 0; wall time 47 min; the climate was extended to 12 years by repeating year 1 (`analysis_output/climate_12_years.csv`; `input/climate.csv` holds 6 identical years, and a first attempt stopped with a KeyError at day 2190); a 100x100 grid was not run: about 31x the per-day cost of 20x20 and 3.7 GB per run, estimated 3.5-4.5 h for 12 pulses with 2-3 communities in parallel):**
- Day 365 by pulse 1-12 (50x50): substrate degraded 44 / 43 / 57 / 64 / 68 / 69 / 70 / 70 / 70 / 71 / 71 / 71 %; biomass C 32 / 37 / 47 / 53 / 57 / 58 / 61 / 60 / 60 / 61 / 61 / 60; taxa present 99 / 91 / 78 / 63 / 53 / 42 / 33 / 28 / 25 / 22 / 18 / 17; energy:C 435 / 423 / 348 / 212 / 135 / 104 / 97 / 91 / 84 / 79 / 69 / 59; CUE 0.46 / 0.44 / 0.40 / 0.35 / 0.35 / 0.36 / 0.37 / 0.38 / 0.38 / 0.38 / 0.38 / 0.38. The 20x20 run, pulses 1-6: degraded 45 / 42 / 49 / 54 / 61 / 63; biomass 32 / 36 / 40 / 46 / 52 / 52; taxa present 90 / 68 / 50 / 40 / 33 / 24; energy:C 449 / 441 / 349 / 319 / 280 / 249.
- Degradation and biomass in pulses 1-2 agree between the grids (taxa present at day 365 already differ: 99 / 91 vs 90 / 68); from pulse 3 the larger grid degrades more (57 / 64 / 68 / 69 vs 49 / 54 / 61 / 63 % in pulses 3-6) and has more biomass and taxa present, and from pulse 4 lower end-of-year energy:C. Taxa present on day 1 of pulses 1-6: 100 / 95 / 88 / 77 / 65 / 52 (20x20: 92 / 81 / 61 / 50 / 38 / 30); at 50x50 all 100 taxa are placed at the start (2,500 cells), at 20x20 about 92.
- After about pulse 7 the community approaches a steady state in degradation (70-71 %) and biomass (60-61), while taxa present keep falling (33 -> 17) and end-of-year energy:C keeps falling (97 -> 59). The abundance-weighted mean allocation on day 1 plateaus at 0.074 (all-taxa mean 0.0635; 20x20 pulse 6 0.073). Osmolyte production is 8.8 % of energy cost over days 200-330 (20x20 9.2 %).
- Energy: the minimum over the year of the mean energy:C is 92 / 96 / 70 / 29 / 9 / 4 / 2 / 0.6 / 0.3 / 0.2 / 0.1 / 0.1 across pulses 1-12 (20x20, pulses 1-6: 86 / 79 / 47 / 42 / 40 / 38), so the dry-season energy collapse appears from pulse 5-6 at the larger grid but not in the first six pulses at 20x20.
- The five communities behave differently: the SD across communities of the degradation at day 365 is 28 percentage points from pulse 4 (20x20: 16-19). One community (seed 1) stays at 19-25 % degraded with biomass 14-20 and energy:C near 200; the other four reach 76-89 % degraded by pulse 6-12 and three of them have energy:C 0 at day 365 from pulse 6 onward. The 5-community means mix the two behaviours.
- Caveats: 5 communities (one of them an outlier); climate repeated identically; 50x50 vs 20x20 comparison covers pulses 1-6 only; no 100x100 check.

**Why energy:C falls with each pulse (`full_climate_km25_osmorange2_50x50_12p/`: `fungal_share.csv`, `energy_decline.png`; `start_biomass_test/`; scripts `analysis/fungal_share.py`, `analysis/start_biomass_test.py`, `analysis/plot_energy_decline.py`):**
- Reseeding (`Grid.repopulation`) gives every new pulse fresh substrate, monomers and enzymes, and places each taxon in each cell with probability equal to its share of the cell-taxa counts of the previous pulse (about 1 taxon per cell). It carries no biomass over, but each placed microbe takes its initial cell mass from `Microbes_init`: C = 1 for bacteria and C = 25 for fungi (half of `max_size_b` 2 and `max_size_f` 50).
- The fungal share of the reseeded cells rises over the pulses (1 / 2 / 10 / 24 / 37 / 47 / 54 / 59 / 62 / 65 / 67 / 69 % of cells in pulses 1-12), and fungi gain share within each pulse (4 / 15 / 28 / 43 / 50 / 56 % of cell-taxa on day 365 in pulses 1-6), so the next pulse is seeded with more fungi. Starting biomass per cell therefore rises from 1.6 (pulse 2) to 17-24 (pulse 12; 5.0 in pulse 1, where about 4 bacterial rows and 0.03 fungal rows per cell are placed). The only community that stays bacterial (seed 1) keeps a starting biomass of 1.3-2.3 and a flat energy:C of 170-220.
- Across the 55 community-pulses of pulses 2-12 the energy:C at day 365 correlates with starting biomass at Spearman -0.81 (fungal share -0.81, peak biomass -0.77, substrate remaining on day 200 +0.76). Mean substrate remaining on day 200 falls from 70 % (pulse 1) to 36 % (pulse 12) and the biomass peak rises from 35 to 99, while the day on which energy:C first falls below 50 stays about the same (day 233-249, set by the dry season).
- Controlled single-pulse test (20x20, 5 communities, real climate, only the start changed; energy:C on day 365 / substrate degraded on day 200): one bacterial cell per cell, C 1: 567 / 19 %; normal pulse-1 start (5 C): 449 / 31 %; same start x3 (15 C): 324 / 46 %; x10 (50 C): 11 / 70 %; one bacterial cell per cell scaled to C 25: 270 / 56 %; one fungal cell per cell, C 25: 140 / 56 %. With equal starting biomass the fungal start ends lower than the bacterial one (140 vs 270, two communities at 0 vs one), so the fungal identity or the smaller number of starting rows adds to the effect but was not separated.
- Interpretation: a larger starting biomass uses up the labile substrate earlier in the year (and with a fixed substrate supply per pulse), so by the dry season (about days 200-295) and after rewetting little substrate is left to restore energy; the end-of-year energy:C then falls. The decline is a consequence of the reseeding scheme (fungal cells start 25x larger and fungi accumulate share), not of the passage of pulses as such.
- Caveat: the mechanism was tested at 20x20 with 5 communities and one climate year; the 50x50 12-pulse statistics come from 5 communities.

**Bacteria only (`fb` 0 in `runtime.txt`; `input_e1500_enz2_km25_osmorange2_bac/`; results in `full_climate_km25_osmorange2_bac_20x20_12p/`; 20x20, 12 pulses, 5 communities (seeds 0-4; `fb` changes the random draws, so the communities are not the same ones as with fungi), real climate extended to 12 years, `Osmo_tol_lower` 0, `Osmo_tol_ref` 0, osmolyte ranges as in `osmorange2`):**
- Energy:C on day 365 does not decline: 502 / 543 / 533 / 523 / 499 / 496 / 531 / 535 / 519 / 509 / 493 / 501 over pulses 1-12 (with fungi at 20x20: 449 / 441 / 349 / 319 / 280 / 249 over pulses 1-6; at 50x50: 435 -> 59 over 12 pulses). The minimum over the year of the mean energy:C stays at 127-141 (with fungi 86 -> 37 at 20x20, 92 -> 0.05 at 50x50).
- Starting biomass stays at 1.0 per cell from pulse 2 (5.1 in pulse 1), biomass at day 365 stays at 25-29 (with fungi 32 -> 52 at 20x20 over 6 pulses, up to 60 at 50x50), substrate remaining on day 200 stays at 75-78 % (with fungi 69 -> 49 %), and degradation per pulse stays at 36-42 % (with fungi 42 -> 63 % over 6 pulses at 20x20, up to 71 % at 50x50).
- Taxa present at day 365 still falls (99 / 80 / 56 / 42 / 34 / 29 / 26 / 22 / 19 / 18 / 17 / 16), so the loss of taxa across pulses is not caused by fungi (with fungi at 20x20: 90 / 68 / 50 / 40 / 33 / 24).
- Per community: three of the five (seeds 2-4) have energy:C of 590-690 and 37-59 % degraded in every pulse; seed 0 has 480-550 and 33-42 %; seed 1 (the low-degradation community) degrades 17-25 % and its energy:C varies between 0 and 211 across pulses.
- Conclusion: the decline of energy:C with each pulse requires the fungal cells (C = 25 at reseeding and increasing share); without them the pulses repeat each other, but the community also degrades only about 40 % of the substrate per year.
- Caveats: 20x20 only, 5 communities, different communities from the fungal runs.

**More communities at 50x50: 15 communities (`full_climate_km25_osmorange2_50x50_12p_15c/` = the 5 original communities (seeds 0-4) plus seeds 5-9 (`..._b2`) and 10-14 (`..._b3`), same input, 12 pulses, 50x50; plots `km25_range_community.png`, `energy_decline.png`; table `fungal_share.csv`):**
- Day 365, mean (SE) over 15 communities, pulses 1 / 2 / 3 / 4 / 6 / 9 / 12: substrate degraded 45 (3) / 44 (4) / 59 (6) / 66 (6) / 70 (7) / 72 (7) / 72 (7) %; biomass C 31 (2) / 36 (4) / 48 (6) / 53 (6) / 56 (6) / 59 (6) / 59 (6); taxa present 99 / 90 / 78 / 64 / 44 / 26 / 18 (SE about 1); energy:C 456 (27) / 443 (27) / 354 (37) / 231 (46) / 121 (43) / 84 (34) / 69 (31); CUE 0.46 / 0.44 / 0.39 / 0.36 then 0.36. The 5-community means were close (energy:C 435 -> 59, degradation 44 -> 71 %); seeds 5-14 alone gave energy:C 467 -> 74 and degradation 45 -> 73 %.
- The mean energy:C falls until about pulse 10 and then levels off at 69-73 (pulses 10-12; the 5-community means still looked to be falling, 79 -> 59), while degradation and biomass flatten from about pulse 6-7; the minimum over the year of the mean energy:C is 94 / 96 / 84 / 49 / 21 / 12 / 3 / 2 / 0.9 / 0.8 / 0.7 / 0.5 across pulses 1-12.
- The community-to-community spread remains large (SE of the degradation 6-7 percentage points, of energy:C 31-46): by pulse 12 twelve of the 15 communities are fungus-dominated (starting biomass 14-25 C per cell, fungal share 55-100 % of cells) and 10 of the 15 have energy:C of 0 at day 365 (7 of 15 in pulse 6); three communities (seeds 1, 5, 7) stay mostly bacterial (starting biomass 2-3 C per cell, fungal share 5-9 %), degrade only 21-37 % in pulse 12 and keep energy:C of 100-376.
- The explanation of the energy:C decline holds with 15 communities: the fungal share of the reseeded cells rises from 1 % (pulse 1) to 67 % (pulse 12), the starting biomass from 1.6 (pulse 2) to 17 per cell, and the Spearman correlation of the energy:C at day 365 with the starting biomass is -0.77 (fungal share -0.76; pulses 2-12 x 15 communities, n = 165).

**How fungi gain share and come to dominate (`fungal_tracking/`: `baseline` (12 pulses) and three variants (6 pulses), 20x20, 5 communities (seeds 0-4), `input_e1500_enz2_km25_osmorange2`, real climate extended to 12 years; recorder `analysis/fungal_tracking.py` (reads the model state, reproduces the earlier run exactly), figures from `analysis/plot_fungal_dominance.py`: `fungal_within_pulse.png`, `fungal_ratchet.png`, `fungal_variants.png`):**
- Model differences between fungi and bacteria (death rates are identical, 0.001 with beta 10): (1) stoichiometry (`Cfrac/Nfrac/Pfrac` fungi 0.9 / 0.09 / 0.01, bacteria 0.825 / 0.16 / 0.015, target N:C 0.100 vs 0.194); (2) cell mass at reseeding and at division (C 25, division above `max_size_f` 50, vs C 1 and `max_size_b` 2); (3) in `Grid.reproduction` the biomass of every fungal cell of a taxon is set to the taxon's grid average each day (translocation) and daughters always spread in the +x direction, while bacterial daughters disperse at random within one cell.
- Within a pulse (pulse 1, mean of 5 communities, day 1 -> day 365): fungal cells 13 -> 172 and biomass 0.84 -> 17.6 C per cell (x21), bacterial cells 1,675 -> 4,160 and biomass 4.2 -> 14.8 (x3.5); fungal share of biomass 17 % -> 53 %. Births per cell present over the year are 3.1 (fungi) vs 2.1 (bacteria) while deaths are similar (1.6 vs 1.5, drought/basal); there are no starvation deaths in any group in any pulse.
- Nitrogen: bacterial biomass N:C falls from the target 0.194 to 0.150-0.152 within about 30 days of each pulse and stays there in all 12 pulses (77 % of target); fungal N:C stays at 0.098-0.108 (target 0.100) until about day 100 and falls to 0.08-0.09 later in the year.
- Reseeding turns the within-pulse gain into a ratchet: the fungal share of cell-taxa on day 1 / day 365 / in the reseeding frequency is 0.8 / 4.2 / 2.5 % (pulse 1), 2.6 / 14.1 / 9.6 % (pulse 2) ... 66 / 68 / 68 % (pulse 12); the fungal share of biomass C on day 1 rises 17 / 38 / 65 / 78 / 87 / 91 % over pulses 1-6 (cells are 25x larger) and reaches 97 % by pulse 12. Median growth factor of the number of cells per fungal taxon is 6 / 15 / 7 / 4 / 3 / 3 over pulses 1-6 (bacterial taxa 1.5 / 1.9 / 1.5 / 1 / 1 / 1). The gain per pulse ends (day 1 and day 365 shares converge, about 62-68 %) once fungal cells fill the grid.
- The community narrows: taxa present on day 365 fall from 8.4 fungal and 82 bacterial (pulse 1) to 2.8 fungal and 6 bacterial (pulse 12); in community 0 in pulse 12, four fungal taxa of 11 present hold 84 % of the cumulative cell counts (pulse 1: 3 % fungal among 91 taxa).
- Variants (fungal share of biomass C on day 1, pulses 1-6): baseline 17 / 38 / 65 / 78 / 87 / 91 %; fungal C:N:P set equal to bacterial 17 / 26 / 27 / 34 / 35 / 39 %; fungal size set equal to bacterial (`max_size_f` 2) 17 / 27 / 40 / 46 / 48 / 50 %; both 17 / 19 / 19 / 19 / 20 / 20 % (about the 19 % that fungal taxa make up of the taxa). Total starting biomass in pulse 6 is 10.9 (baseline), 2.2, 1.0 and 1.0 C per cell, and end-of-pulse biomass 52 / 26 / 45 / 30. With both equalised, fungi gain no share, so translocation and directed spread give no advantage by themselves here (they cannot be switched off by a parameter and stay in all variants).
- Interpretation: the lower N and P requirement per unit C gives fungi the competitive edge in this N-limited system (largest effect: without it the fungal share of cells stays at 1-8 %), the 25x larger cell mass multiplies it (without it the share of biomass levels at about 50 %), and reseeding in proportion to cell-taxa counts carries each pulse's gain into the next start. Not tested: a direct manipulation of the N supply; one grid size (20x20) and 5 communities.

**Test of the N explanation: raising the N supply (`fungal_tracking/subN_x0.5`, `subN_x2`, `subN_x4`, `subN_x10` vs `baseline`; `fungal_tracking.py --sub-n-mult K` multiplies the N of every initial substrate (initial substrate N 8.3 per cell against 346 C, C:N about 42; `Input_NH4` is 0 and unused); 20x20, 5 communities (same communities in every run), 6 pulses; figure `n_supply.png`, script `analysis/plot_n_supply.py`):**
- Bacterial N:C (days 100-365, mean of pulses 1-6; target 0.194): 0.150 (x0.5), 0.150 (x1), 0.196 (x2), 0.260 (x4), 0.273 (x10); fungal N:C (target 0.100): 0.059, 0.090, 0.169, 0.170, 0.159. Pulse-1 biomass growth of bacteria (day 365 / day 1): x1.1 (x0.5), x3.5 (x1), x6.6 (x2), x6.4 (x4), x5.2 (x10); of fungi x14, x21, x15, x12, x9.
- Fungal share of biomass C on day 1, pulses 1-6: x0.5 17 / 46 / 84 / 94 / 97 / 98 %; x1 17 / 39 / 65 / 78 / 87 / 91; x2 17 / 31 / 53 / 70 / 78 / 82; x4 17 / 20 / 33 / 48 / 57 / 62; x10 17 / 21 / 34 / 40 / 49 / 54. Fungal share of cell-taxa on day 365 in pulse 6: 79 / 53 / 22 / 11.5 / 8 % (x0.5 / x1 / x2 / x4 / x10).
- Total starting biomass C per cell in pulse 6: 18.8 / 10.9 / 5.1 / 2.8 / 2.6; energy:C on day 365, pulse 1 -> 6: x0.5 509 -> 187; x1 449 -> 249; x2 389 -> 199; x4 346 -> 250; x10 264 -> 258 (the decline over the pulses is -322, -200, -190, -96 and -6). Substrate degraded on day 365, pulse 1 -> 6: x0.5 31 -> 63 %; x1 45 -> 63; x2 52 -> 55; x4 50 -> 54; x10 42 -> 37.
- Conclusion: the fungal takeover and the starting-biomass ratchet are controlled by the N supply. Lowering N makes it stronger and faster, raising N weakens it, and at x10 the decline of energy:C across pulses is gone. The fungal share still rises at x10 (17 -> 54 % of biomass), so N is the main driver but not the only one (cell mass at reseeding is). Raising N also changes the first pulse itself (bacteria grow x1.1 to x6.6, energy:C in pulse 1 falls from 509 to 264) and lowers the degradation at x10.
- Caveats: only organic N in the substrates was changed (NH4 starts at 0); 20x20, 5 communities.

**N tracker (`analysis/n_tracker.py`, `analysis/plot_n_budget.py`):** records, per grid cell and day, the N pools (litter substrates, DeadMic + DeadEnz substrates, organic monomers, NH4, microbes by group, enzymes (implicit: enzyme pool C x N_cost)) and fluxes (decay, uptake by group, N into enzymes, enzyme turnover into DeadEnz, N of killed cells, N released to NH4 by living cells (MicLoss), leaching, and the C released in the same mortality step by cells with too little N (overflow C) by group), and a daily closure check. It reads the model state and does not change it. Closure test (seed 0, 2 pulses): largest daily residual 5e-7 (float32 rounding), and the change in total N over pulse 1 (-0.0636) equals the cumulative leaching (0.0636), so N is conserved apart from leaching. HSP and osmolyte N costs are not deducted from any pool in the model; they only enter the N limitation of inducible production. N is added to the system at each pulse start because reseeded cells bring their initial N.

**N budgets from the N tracker (`n_tracking/x1` = baseline and `n_tracking/x4` = N x4; 20x20, 5 communities, 6 pulses, real climate; figures `n_pools.png`, `n_fluxes.png`; script `analysis/plot_n_budget.py`; per grid cell and year, mean of communities):**
- Closure: largest daily residual 1.3e-5 (x1) and 1.8e-5 (x4); N is conserved apart from leaching, which is 0.05 N per year (0.5 % of the 9.3 N present) at x1 and 12.7 N (38 % of 33.7) at x4.
- x1, pulse 1 / 3 / 6: N released by decay 8.9 / 8.5 / 10.2; N uptake by bacteria 5.5 / 1.6 / 0.6 and by fungi 3.4 / 6.9 / 9.6 (fungi 38 / 81 / 94 % of the uptake); N in cells killed 3.2 / 0.8 / 0.3 (bacteria) and 1.3 / 2.6 / 4.3 (fungi), i.e. 58 % of the N taken up by bacteria in pulse 1 and 37-45 % of that taken up by fungi; N into enzymes 1.3 / 1.8 / 2.4, all of which returns through enzyme turnover (1.3 / 1.8 / 2.4); N released to NH4 by living cells 0.09 (bacteria) and 0.02 (fungi) in pulse 1 and about 0 later.
- Where the N is at the end of the year, x1 pulse 1 / 3 / 6: litter substrates 1.9 / 1.9 / 1.0; dead microbes and enzymes (substrate) 3.3 / 3.0 / 4.1 (36-44 % of the 9.3-9.5 N); bacteria 2.2 / 0.8 / 0.3; fungi 1.6 / 3.0 / 4.0. The dead pool fills mainly in days 150-295 (dry season) and holds the N of killed cells until it decays.
- C released to respiration by cells with too little N (overflow), % of the C they take up, pulse 1 / 3 / 6: bacteria 17.5 / 15.5 / 16.0 % (about 20-25 % between days 50 and 200) against fungi 0.3 / 1.7 / 0.5 % at x1; at x4 bacteria 4.3 / 2.0 / 2.9 % and fungi 0.3 / 0 / 0 %. This is the stoichiometric C loss (`MicLoss` in `Grid.mortality`) that lowers bacterial CUE when N is short.
- x4, pulse 1: N uptake by bacteria 46 and fungi 14, N released to NH4 by living cells 31 (bacteria) and 11 (fungi), leached 12.7; N:C of the bacterial biomass 0.26 (above the 0.194 target) so the surplus N is released and leached.
- N cycling: the N of the substrates (8.3) is recycled several times in a year at x1 (N released by decay 8.9, uptake 8.9) because killed cells return their N to the dead pool; fungi take about all of the decay N from pulse 3 on (uptake = decay in the pulse-6 curve).
- Not in the budget: HSP and osmolyte N costs (not deducted anywhere in the model); P was not tracked.

**Why necromass (DeadMic, DeadEnz) accumulates (`necromass_check/substrate_decay.csv`, `necromass_decay.png`; scripts `analysis/necromass_check.py`, `analysis/plot_necromass.py`; 20x20, 5 communities, 2 pulses, `osmorange2` input, real climate):**
- Both necromass pools are degradable in the model (DeadMic needs 3 enzymes, DeadEnz 6, as for the other substrates, and each has a monomer with its own transporter) and nothing in the code treats them specially.
- First-order rate (decay / mass, per day), days 1-100 / 150-250 / 300-365 in pulse 1: DeadMic 0.0063 / 0.0035 / 0.0035 (turnover 160 / 280 / 280 days); DeadEnz 0.0071 / 0.0067 / 0.0090; for comparison starch 0.0035-0.0043, chitin 0.0050-0.0058, proteins 0.0035-0.0058, cellulose 0.0007-0.0016, hemicellulose 0.0016-0.0028, lignin 0.0005-0.0009. Every pool turns over only 0.1-0.9 % per day because the pool masses (0.4-30 per cell for DeadMic) are far below Km (about 58 per cell at slope 2.5, estimated from the earlier probe at slope 7), so the rate is proportional to the pool mass (first-order regime).
- DeadMic grows from 1.4 (day 100) to 30 (day 295) per cell (C+N+P) because the dead cells arrive in the dry season (days 150-295, Psi about -5 MPa), when the enzyme Vmax is multiplied by about 0.35 and the uptake by about 0.12; the decay stays at 0.03-0.05 per day then and rises to 0.1 per day after rewetting (day 295), which is not long enough to clear the pool (25.8 per cell on day 365 in pulse 1, 23.5 in pulse 2). `Grid.repopulation` resets the substrates, so this necromass, and its N (36-44 % of all N at the end of a pulse), is discarded at each pulse boundary.
- Monomers: the pool left after uptake is small in the wet season (Mon_DeadMic 0.018 per cell, Mon_DeadEnz 0.013) and large in the dry season (Mon_DeadMic 0.67 on average over days 150-250, up to about 1.6; inflow 0.046 per day), i.e. uptake limits the necromass monomers only in the dry season; after rewetting the pool falls to 0.4 and declines.
- Not done: no change to the decay of necromass was tested (options: a lower Km for the necromass rows, a separate first-order necromass decay, carrying necromass over between pulses).

**Minimal model (one taxon, one substrate), N limitation with a C-rich substrate (`input_minimal_e1500_crich` = `input_minimal_e1500` with `IdealSub` C 200, N 5, P 0.5 instead of 100, 5, 0.5; `input_minimal_e1500_km25`, `input_minimal_e1500_crich_km25` = the same with the four Vmax-Km slopes at 2.5 (built by `analysis/build_full_input_km.py`); flat 15 C, Psi -0.5, one 365-day pulse, 20x20 (400 identical cells), 3 seeds, `n_tracker.py --keep-osmo`; results in `minimal_n/` (slopes 7/5) and `minimal_n_km25/` (2.5/2.5), figure `minimal_n.png`, script `analysis/plot_minimal_n.py`):**
- The settled minimal model is already N-limited: with the 100:5:0.5 substrate the biomass N:C falls from 0.194 (target) to 0.143 by day 365 (below 97 % of the target from day 17) and 11.9 % of the C taken up is released to respiration by cells short of N; the N:C of the uptake is 0.051 against the 0.194 needed. Cells level off at N:C 0.143-0.147 in all four runs (what this level is: see the next paragraph on the minimum quota).
- Slopes 7/5, reference (C:N 20) / C-rich (C:N 40): biomass C on day 365 8.7 / 3.7; substrate C degraded 26.7 (26.7 %) / 22.3 (11.2 %); N uptake 1.34 / 0.53; N:C of the uptake 0.051 / 0.026; overflow 11.9 / 41.0 % of the C uptake, starting on day 99 / 51; CUE of the model (uptake-weighted) 0.47 / 0.59, effective CUE (model CUE minus the overflow share) 0.35 / 0.18.
- Slopes 2.5/2.5, reference / C-rich: biomass C 18.5 / 6.0; substrate C degraded 53.1 (53.1 %) / 32.7 (16.4 %); N uptake 2.97 / 0.92; N:C of the uptake 0.053 / 0.028; overflow 13.7 / 42.0 %, starting on day 51 / 31; CUE of the model 0.50 / 0.61, effective CUE 0.36 / 0.19.
- Faster kinetics (2.5/2.5) double the use of the ideal substrate (degradation 26.7 -> 53.1 %) but leave the N limitation unchanged (N:C 0.143-0.145, overflow share 12-14 % and 41-42 %, effective CUE 0.35-0.36 and 0.18-0.19); the C-rich substrate still degrades only 16.4 % of its C (32.7 C, against 53.1 C of the 100-C substrate) and reaches a third of the biomass. Lower biomass means less enzyme production, which probably limits the decay further; not isolated.
- The limit is the N supplied per unit of C (set by the substrate stoichiometry), not the amount of N: 4.2-4.4 of the 5 N are still in the substrate on day 365 in the C-rich runs, NH4 is 0 and nothing leaches.
- The CUE reported by the model (`CUE_system`) is computed in `metabolism` before the overflow, which happens in `mortality`; the overflow is not in it. The N:C and the overflow oscillate with a period of 8-16 days (shorter at 2.5/2.5), probably from synchronised division of the 400 identical cells (not confirmed).
- Caveats: one temperature, one year, 3 seeds of 400 identical cells.

**The same two minimal-model comparisons with a fungus instead of a bacterium (`fb` 1 in `runtime.txt`: `input_minimal_e1500_fb1`, `input_minimal_e1500_crich_fb1`, `input_minimal_e1500_km25_fb1`, `input_minimal_e1500_crich_km25_fb1`, each the matching `fb` 0 input with only `fb` changed; results in `minimal_n_fb1/` (slopes 7/5) and `minimal_n_km25_fb1/` (2.5/2.5), figure `minimal_n.png` made with `analysis/plot_minimal_n.py ... fun`; flat 15 C, 365 days, 3 seeds, `n_tracker.py --keep-osmo`; N closes, residual below 1e-6):**
- A fungal taxon is placed in a cell with probability `taxa_per_box` x `max_size_b` / `max_size_f` = 0.04, so it starts in about 14 of the 400 cells with C = 25 each (grid-mean 1.0 C per cell, as for the bacterium) and has to spread (occupied cells day 1 / 100 / 365: 14 / 26 / 66 and 14 / 26 / 245 for C 100 and C 200 at slopes 7/5; 14 / 48 / 101 and 14 / 48 / 312 at 2.5/2.5).
- Fungus, slopes 7/5, C 100 (C:N 20) / C 200 (C:N 40): biomass C on day 365 7.0 / 20.5; substrate C degraded 19.1 (19.1 %) / 48.3 (24.2 %); N uptake 1.04 / 1.34; N:C of the uptake 0.053 / 0.027; biomass N:C 0.119 / 0.054; overflow 0 / 0 % of the C uptake; effective CUE 0.40 / 0.45. Slopes 2.5/2.5: biomass 12.2 / 41.9; degraded 31.7 (31.7 %) / 95.1 (47.6 %); N uptake 1.85 / 2.79; biomass N:C 0.121 / 0.054; overflow 0 / 0.5 %; effective CUE 0.41 / 0.45.
- The same comparisons for the bacterium (above): biomass 8.7 / 3.7 (7/5) and 18.5 / 6.0 (2.5/2.5); degraded C 26.7 / 22.3 and 53.1 / 32.7; overflow 11.9 / 41.0 % and 13.7 / 42.0 %; effective CUE 0.35 / 0.18 and 0.36 / 0.19.
- Minimum N:C: `Microbe.minimum_cell_quota` sets `MinRatios` = optimal mass fraction - tolerance (Cfrac 0.825 / 0.90 minus 0.09, Nfrac 0.16 / 0.09 minus 0.04, Pfrac 0.015 / 0.010 minus 0.005, for bacteria / fungi), i.e. C 0.735, N 0.120, P 0.010 (bacteria) and C 0.810, N 0.050, P 0.005 (fungi). In `Grid.mortality` a living cell is rebalanced when N / (C+N+P) < MinRatios N, which corresponds to a trigger N:C of about 0.138 (bacteria) and 0.053 (fungi) (with the measured P:C of 0.0136 and 0.0119); the rebalancing then sets the cell to N:C = MinRatios N / MinRatios C = 0.163 (bacteria) and 0.062 (fungi) by releasing the excess C (targets 0.194 and 0.100). The observed N:C lies between trigger and reset: a bacterium in the C 100 run oscillates between 0.138 and 0.162 (seed 0, days 250-365), the population mean in the C 200 runs sits at 0.143-0.148, and the fungus is pinned at 0.0535-0.054 just above its trigger; this is the level at which the C-rich fungal run sits from about day 140. The fungus is N-replete with the C:N 20 substrate (biomass N:C 0.12 above the 0.10 target) and C-limited, and with C:N 40 it reaches its quota: new biomass N:C = N:C of the uptake / CUE = 0.027 / 0.45 = 0.06, which is between the trigger (0.053) and the reset value (0.062), so almost no overflow is needed (overflow only flickers on days 140-175). A substrate richer in C than C:N 40 would be expected to push the fungus into overflow (not tested).
- So the same substrate has opposite effects: more C helps the fungus (2.2x and 3x more C degraded than with the C 100 substrate) and penalises the bacterium (and the fungal runs show no sawtooth of N:C). This is the stoichiometric reason for the fungal advantage found in the multi-taxon runs.
- Caveats: the fungus has not colonised the whole grid after a year (16-78 % of the cells), so the year-end values mix stoichiometry with spatial spread; one temperature; 3 seeds.

**Energy balance and C allocation of the minimal model (`minimal_e_km25_fb1/` = fungus and `minimal_e_km25/` = bacterium, each with `ref/` (IdealSub C 100) and `crich/` (C 200); Km slopes 2.5/2.5, flat 15 C, 365 days, 3 seeds; recorder `analysis/e_tracker.py` (reads model state, closes the daily energy budget to 0.003 on fluxes of about 500), figures `minimal_e.png` from `analysis/plot_minimal_e.py`):**
- No cell is ever short of energy in any of the four cases: the E pool is never floored at zero (0 cell-days) and the cap at E_max (1400) is never reached. The model's E:C of the biomass is 34-68 % of E_max.
- The cells route about half of the C they take up to energy in every case: 47-63 % (fungus C 100: 52 / 48 / 54 / 63 % in days 1-30 / 31-100 / 101-200 / 201-365; fungus C 200: 51 / 48 / 51 / 54; bacterium C 100: 53 / 44 / 47 / 51; bacterium C 200: 52 / 35 / 36 / 40). The fraction is set by the E:C (fraction routed = 1 - (E:C / E_max)^1.1).
- Energy spending (per unit biomass): maintenance 47-77 % of the energy spent, heat shock proteins 14-16 %, transport 7-36 %, enzymes and osmolytes about 1.5 %.
- Fate of the C taken up, days 31-100 / 201-365 (energy / enzymes / released by cells short of N or P / retained as biomass, % of uptake): bacterium C 100: 44 / 0.7 / 23 / 32 and 51 / 1.1 / 12 / 36; bacterium C 200: 35 / 0.6 / 54 / 11 and 40 / 0.9 / 41 / 18; fungus C 100: 48 / 0.7 / 0 / 51 and 63 / 2.4 / 0 / 34; fungus C 200: 48 / 0.7 / 0 / 51 and 54 / 1.2 / 0 / 45 (released 3.2 % in days 101-200).
- Energy produced / spent (days 1-30 / 31-100 / 101-200 / 201-365): fungus C 100 2.32 / 2.16 / 1.47 / 0.92; fungus C 200 2.46 / 2.30 / 1.73 / 1.62; bacterium C 100 2.14 / 1.95 / 1.53 / 1.55; bacterium C 200 2.27 / 1.71 / 1.37 / 1.41. Only the fungus on the C 100 substrate goes below 1, after about day 270, when the C taken up falls with the substrate: its E:C falls to 34 % of E_max and the fraction of uptake C routed to energy rises to 0.85 (63 % of the C uptake) and the retained share falls to 34 %.
- Energy limitation and N limitation do not occur together in these runs: the N-limited bacterium has an energy surplus (produced / spent 1.4-1.9) and a higher E:C (up to 950 = 68 % of E_max on the C 200 substrate); the fungus on the C 100 substrate is N-replete (biomass N:C 0.12 above its 0.10 target) and the only case in which energy produced falls below energy spent. In the bacterium each release of C by cells short of N lowers the C but not the E of the cell, so E:C jumps and the fraction routed to energy falls in a sawtooth with the synchronous divisions of the identical cells.
- Caveats: the three seeds and the 400 identical cells; flat 15 C; the energy shortfall indicator is the floor of the E pool at zero (the throttle of the production of enzymes, HSP and osmolytes by the available energy, `ARREO*`, is internal to `Grid.metabolism` and not recorded).

**Ongoing NH4 and PO4 inputs (code change) and the minimal-model scan (`nh4_scan_minimal/`; `src/grid.py` changed (backup `src/grid.py.bak_before_nh4_po4_input`): `Grid.__init__` reads the NH4 and PO4 entries of `MonInput` and `Grid.uptake` adds `Input_NH4` as N and `Input_PO4` as P to the NH4 and PO4 pools of every grid cell each day before uptake; the organic-monomer and FermProd inputs of `sub_mon_inputs.csv` stay unused; before the change `Input_NH4` and `Input_PO4` had no effect (test: identical output); scripts `analysis/nh4_scan.py` (stages `nh4` and `po4`; `--values`, `--tag`) and `analysis/n_tracker.py` (now with P pools and fluxes, the inputs, the limiting element of each C release, cells below the N and P triggers, `--set NAME=VALUE`)):**
- Verification: with both inputs 0 the new code reproduces the earlier `minimal_n_km25/crich` run exactly (difference 0.0 over 32 columns); with inputs the N and P budgets close to 5e-7 when the inputs are counted.
- Setup: minimal bacterium, `input_minimal_e1500_km25` (C 100) and `input_minimal_e1500_crich_km25` (C 200), Km slopes 2.5/2.5, flat 15 C, 365 days, 3 seeds; Input_NH4 0 and 0.0005-0.05 (log scale) with Input_PO4 = 0.09 x Input_NH4; relieved = C released by cells short of N below 2 % of the C uptake (year total) and daily mean biomass N:C above the trigger 0.138 on all days 31-365 (the N:C part never binds: the daily mean stays at 0.138-0.139 even with no input; the share of cell-days below the trigger is reported as well).
- C released by cells short of N, % of the C uptake, Input_NH4 0 / 0.0005 / 0.001 / 0.002 / 0.003 / 0.005 / 0.008 / 0.0125 / 0.02: C 100: 13.7 / 11.7 / 9.4 / 7.1 / 4.2 / 0 / 0 / 0 / 0; C 200: 42.0 / 40.1 / 37.2 / 34.6 / 31.4 / 25.6 / 18.9 / 12.8 / 1.9. Refinement: C 100 2.9 % at 0.0035, 1.6 % at 0.004, 0.06 % at 0.0045; C 200 10.1 % at 0.015, 5.8 % at 0.0175, 1.9 % at 0.02, 0 at 0.03 and 0.05.
- First input meeting the criterion: C 100 (C:N 20): Input_NH4 0.004 N per cell per day (with PO4 0.00036); C 200 (C:N 40): 0.02 (with PO4 0.0018), 5x more. Inputs per year: 1.46 N and 0.13 P (0.29 of the 5 N initially in the substrate) for C 100, 7.3 N and 0.66 P (1.46 x the initial substrate N) for C 200. No P-limited release at these levels with the tied PO4 (0 %); the share of cell-days below the N trigger 0.37 % and 0.16 %.
- Effect at those inputs (without input): biomass C on day 365 31.9 (18.5) and 51.7 (6.0); substrate C degraded 74.6 % (53.1 %) and 59.9 % (16.4 %); effective CUE 0.43 (0.36) and 0.44 (0.19); N uptake 5.06 (2.97) and 8.36 (0.92). More N than needed does not help the C 100 substrate (biomass 34.3 at 0.005 and 31.9 at 0.02: C-limited; biomass N:C rises to 0.24); for C 200 biomass is 51.3 at 0.03 and 47.6 at 0.05.
- Leaching (N leached because of the input as % of the N input, mean of 3 seeds; zero-input leaching subtracted): 45 % (C 100 at 0.004) and 38 % (C 200 at 0.02); P leached because of the input at the tied PO4: 30 % and 37 %. The leaching of N falls from about 68 % (0.0005) to 38 % (0.02) as the uptake increases.
- PO4 stage (PO4 as a multiple of 0.09 x Input_NH4, at 0.004 and 0.02): for C 200, 0 x / 0.25 x / 0.5 x / 1 x / 2 x / 4 x: C released by cells short of P 34.7 / 19.8 / 11.2 / 0 / 0 / 0 % (short of N 0 / 0 / 0 / 1.9 / 2.2 / 2.3 %), biomass 8.8 / 20.7 / 36.5 / 51.7 / 50.9 / 49.4: once N is supplied the P becomes limiting and at least the tied PO4 is needed; for C 100: 0 x gives 1.1 % short of P (and 0.7 % short of N), from 0.25 x on 0 % short of P; biomass 30.3 / 31.5 / 32.1 / 31.9 / 31.9 / 31.5. More P than the tied amount does not help.
- Caveats: the PO4 leaching attributable to the input is not meaningful below an input of about 0.001 (the PO4 pool also receives P released by cells); 3 seeds of 400 identical cells; flat 15 C, one year; criterion decided by the C-release part only.

**Minimal model with Input_NH4 0.02 and Input_PO4 0.002 (inputs set in four new copies of the Km 2.5/2.5 folders: `input_minimal_e1500_km25_in`, `input_minimal_e1500_crich_km25_in`, `input_minimal_e1500_km25_fb1_in`, `input_minimal_e1500_crich_km25_fb1_in`; only the two parameter rows differ; flat 15 C, 365 days, 3 seeds; results and figures in `minimal_n_km25_in/` (bacterium), `minimal_n_km25_fb1_in/` (fungus), `minimal_e_km25_in/`, `minimal_e_km25_fb1_in/`; `plot_minimal_n.py` now also draws the N input and the N leached):**
- Bacterium, C 100 / C 200, without input -> with input: biomass C on day 365 18.5 -> 32.3 / 6.0 -> 51.5; substrate C degraded 53 -> 77 % / 16 -> 60 %; C released by cells short of N 13.7 -> 0 % / 42.0 -> 1.8 % of the uptake; effective CUE 0.36 -> 0.43 / 0.19 -> 0.44; biomass N:C 0.143 -> 0.243 / 0.145 -> 0.142; N uptake 3.0 -> 8.9 / 0.9 -> 8.3 per cell per year; N input 7.3, N leached 3.0 (C 100) and 2.8 (C 200).
- With the C 200 substrate the constant input stops being enough at the end of the year: the N uptake rises to 0.038 / 0.046 / 0.052 per day on days 300 / 340 / 365 (input 0.02, the rest from decay), the biomass N:C falls from 0.19 (day 175) to 0.142 and the release of C starts on day 349 (20 % of the daily uptake on day 355, 7.8 % over the last 30 days), all N-limited (0 % P-limited). The year total (1.8 %) hides this.
- Fungus, C 100 / C 200, without input -> with input: biomass 12.2 -> 11.5 / 41.9 -> 40.8; degraded 32 -> 32 % / 48 -> 48 %; release 0 / 0.5 -> 0; biomass N:C 0.121 -> 0.181 / 0.054 -> 0.158; N uptake 1.8 -> 8.5 / 2.8 -> 7.5; N released by the living cells (C 100) 5.8; N leached 6.5 / 2.9. The input changes nothing for the fungus except N storage and N loss: it is C-limited, not N-limited.
- Energy balance with the input: no energy shortfall in any of the four cases; energy produced / spent 1.5 / 1.7 (bacterium, days 201-365), 0.91 / 1.53 (fungus C 100 / C 200); C routed to energy 57 / 54 % (bacterium), 67 / 56 % (fungus); bacterium C retained as biomass 42 / 43 % (previously 36 / 18 %); E:C smooth, without the sawtooth of the N-limited runs. The late-year energy deficit of the fungus on the C 100 substrate is unchanged (E:C 31 % of E_max, 31 % of the C retained).
- Caveats: constant input per cell and day; 3 seeds of 400 identical cells; the NH4 pool averages 0.05-0.12 and the PO4 pool 0.005-0.012 per cell.

**Inputs set in the main Km 2.5/2.5 folders, and the minimal model under the real climate (year 1 of `input/climate.csv`, one 365-day pulse; `Input_NH4` 0.02 and `Input_PO4` 0.002 now set in the eight folders with Km slopes 2.5/2.5: `input_minimal_e1500_km25`, `_crich_km25`, `_km25_fb1`, `_crich_km25_fb1` and `input_e1500_enz2_km25`, `_km25_osmorange`, `_km25_osmorange2`, `_km25_osmorange2_bac` (the `_in` copies are now duplicates of the first four); each README notes the change, and earlier results from these folders used inputs of 0; results and figures in `minimal_n_clim_km25/`, `minimal_n_clim_km25_fb1/`, `minimal_e_clim_km25/`, `minimal_e_clim_km25_fb1/` (the plot scripts take `--climate CSV` and `--label TEXT`); 3 seeds):**
- The climate year: mean temperature 15.7 C; Psi below -1 MPa on 223 days (days 44-362) and below -4 MPa on 103 days; Psi -0.5 on day 100, -5.2 on day 200, -5.3 on day 250, -1.3 on day 295. A single taxon has drought tolerance 1 in the default relative scale (`Osmo_tol_lower` 0), so there is no drought mortality (all 400 bacterial cells alive on day 365).
- Bacterium, C 100 / C 200 (flat 15 C with the same inputs in brackets): biomass C on day 365 7.1 (32.3) / 10.8 (51.5); substrate C degraded 22.7 % (76.8 %) / 15.9 % (59.7 %). Fungus, C 100 / C 200: 4.7 (11.5) / 10.2 (40.8); degraded 15.6 % (31.6 %) / 15.0 % (47.9 %); the fungus occupied 64 and 133 of the 400 cells.
- No C is released by cells short of N or P in any of the four cases (0 % of the uptake); biomass N:C 0.16-0.26 on average; N input 7.3 per cell, N leached 4.8-7.1.
- Instead the cells become energy-limited in the dry season: energy produced / spent is below 1 on 124-163 days (minimum 0.26 / 0.29 for the bacterium, 0.13 / 0.26 for the fungus), E:C falls to 0 (minimum of the mean), the fraction of uptake C routed to energy reaches 96 % (days 201-365: 77 / 69 % bacterium, 89 / 70 % fungus) and the retained share of the C taken up falls to about 4 % at the peak; energy shortfalls (E pool floored at zero) start on day 245-266 and peak on day 246-268: 6.9 % / 0.7 % (bacterium) and 14.1 % / 1.1 % (fungus) of the cell-days; osmolytes take 10-14 % of the energy spent on days 101-365. The deficit ends at the rewetting on day 295, when the energy produced spikes (up to 4.6 per unit of energy spent).
- The C 100 substrate (less C per cell) is more affected than C 200, and the fungus more than the bacterium; the N input has no role in this energy limitation.
- The cells shed the surplus N (N released by living cells per cell and year, C 100 / C 200: bacterium 3.7 / 2.0, fungus 4.3 / 3.6), which then leaches; the fungus's biomass N:C resets to 0.062 whenever it reaches about 0.21.
- Caveats: one climate year, one pulse; 3 seeds of 400 identical cells; no drought mortality because tolerance is 1.

**What drives the burst of energy production at day 295-300 under the real climate (`spike_diagnosis/spike.csv`, `spike_diagnosis.png`; scripts `analysis/diagnose_spike.py`, `analysis/plot_spike.py`; minimal bacterium on the C 100 substrate, Km 2.5/2.5 with the inputs, 3 seeds; counterfactual climates run on the same communities):**
- It is the rewetting. Psi rises from -4.1 (day 294) to -1.3 (day 295) and -0.5 MPa (day 296); the moisture factor on uptake (`exp(0.5 (Psi - wp_fc))`) rises from 0.13-0.22 to 0.85 and 1.0 and the factor on enzyme Vmax from 0.35-0.47 to 0.92 and 1.0. The C taken up per unit biomass rises from 0.008-0.011 (days 288-294) to 0.039-0.043 (days 295-296) and the energy produced per unit biomass from 10-14 to 50-54 (5.8 times the mean of days 260-280, peak on day 296).
- It is a rate effect on a stock: during the dry season the monomer C pool seen by uptake grows from 1.2-1.6 (days 100-200) to 2.4-2.9 per cell (days 250-295) because uptake is suppressed more (about 0.15) than decay (about 0.4); the share of the pool taken per day is 1.3 % in days 201-293 and 8-11 % after rewetting. Uptake is first-order in the pool, so the uptake scales with the product of the moisture factor and the pool: 0.18 x 2.51 on day 288 and 1.0 x 2.28 on day 296 give a 5.2-fold increase against 4.7-fold observed in the C uptake. The pool is drawn down from 2.4 to 1.2 by day 310, which ends the burst; the C released by decay rises from 0.022 to 0.075 per day over days 294-300.
- Counterfactuals (peak energy per unit biomass in days 285-315 relative to days 260-280): real climate 5.8x; Psi held at its day-285 value from day 285 (no rewetting) 1.4x; Psi -0.5 all year (real temperature) 1.0x (no stockpile either: pool 1.1); temperature set to 15 C from day 285 (real Psi) 6.0x; temperature 15 C all year 7.4x. The warm days after rewetting (up to 28 C on day 302) are not the trigger.
- The dry-season energy limitation is therefore an uptake limitation by moisture (the monomer pool is twice as large as in the wet season while the cells are energy-limited), not a shortage of substrate.
- Caveats: one substrate, bacterium, 3 seeds; the fungal runs show an N-uptake peak at the same time (`minimal_n_clim_km25_fb1/minimal_n.png`) but were not diagnosed separately.

**Lowering the moisture sensitivity of the uptake (`uptake_psi_scan/`; code change: the hardcoded 0.5 in the uptake moisture factor `f_psi = exp(0.5 (Psi - wp_fc))` of `Grid.uptake` is now the optional parameter row `Uptake_psi_slope` (default 0.5; files `src/grid.py`, `src/initialization.py`, backups `*.bak_before_uptake_psi_slope`; with the row absent or set to 0.5 a stored real-climate run is reproduced exactly); scripts `analysis/uptake_psi_scan.py`, `analysis/e_tracker.py` (now with `--set` and the substrate C); minimal model, real climate year 1, Km 2.5/2.5 with Input_NH4 0.02 / Input_PO4 0.002, slopes 0.5 / 0.35 / 0.25 / 0.1 / 0 (factor at Psi -5 MPa 0.14 / 0.25 / 0.37 / 0.67 / 1.0), bacterium and fungus on the C 100 and C 200 substrates, 3 seeds):**
- Bacterium, C 100, slope 0.5 / 0.35 / 0.25 / 0.1 / 0: % of cell-days with an energy shortfall 6.9 / 0 / 0 / 0 / 2.9; biomass C on day 365 7.1 / 8.3 / 8.3 / 8.0 / 7.2; substrate C degraded 22.7 / 26.1 / 26.6 / 26.8 / 25.5 %; days with energy produced < spent 144 / 124 / 127 / 128 / 124 and the minimum of produced / spent 0.26 / 0.27 / 0.27 / 0.27 / 0.27 (unchanged); the burst at the rewetting (peak energy per biomass relative to days 260-280) 5.8 / 3.8 / 2.8 / 2.3 / 1.7; C routed to energy on days 201-293 91 / 88 / 85 / 84 / 85 %; the C released for lack of N or P stays 0.
- Bacterium, C 200: shortfall 0.7 / 0 / 0 / 0 / 0 %; biomass 10.8 / 11.2 / 11.2 / 11.2 / 11.0; degraded 15.9 / 16.7 / 16.8 / 17.2 / 17.5 %; days produced < spent 124 / 116 / 112 / 118 / 114.
- Fungus, C 100: shortfall 14.1 / 18.6 / 14.1 / 37.3 / 39.4 % (worse at slopes <= 0.1: at slope 0 the shortfall lasts from day 254 to day 365 and the energy produced / spent falls to about 0.1 by day 365, instead of days 245-294 at slope 0.5); biomass 4.7 / 5.1 / 5.2 / 4.9 / 5.0; degraded 15.5 / 17.2 / 18.1 / 15.4 / 15.2 %; days produced < spent 163 / 151 / 146 / 208 / 206. Fungus, C 200: shortfall 1.1 / 0 / 0 / 0 / 0 %, degraded 15.0 / 15.7 / 15.9 / 16.2 / 16.1 %. Why the fungus does worse at the lowest slopes was not examined.
- The energy deficit of the dry season does not depend mainly on the uptake moisture factor: with no moisture effect on uptake (slope 0) the bacterium on the C 100 substrate still has energy produced < spent on days 186-294 (produced / spent 0.69 on days 201-293 against 0.60 at slope 0.5, i.e. about a quarter of the deficit removed). Per unit biomass (days 1-100 / 201-293 at slope 0): C uptake 0.028 / 0.009 (biomass 2 -> 5, the decay flux limited by the enzyme moisture factor), maintenance 5.8 / 7.5 (mean temperature 12 -> 20 C), osmolytes 0.1 / 2.5 (10-20 % of the energy spent), energy produced 19.9 / 10.6.
- Caveats: one climate year; 3 seeds of identical cells; the enzyme moisture factor `exp(0.25 (Psi - wp_fc))` and the osmolyte production were not changed.

**Moisture sensitivities 0.25 for enzymes and uptake, C 200 substrate, real climate (`minimal_n_clim_psi25(_fb1)/`, `minimal_e_clim_psi25(_fb1)/` single year; `minimal_clim6_psi25/` six pulses; scripts `analysis/plot_minimal_n.py`, `plot_minimal_e.py`, `plot_minimal_pulses.py`; input folders `input_minimal_e1500_crich_km25(_fb1)`, Input_NH4 0.02, Input_PO4 0.002, Km slopes 2.5/2.5):**
- Code: the enzyme moisture slope is now the optional parameter `Enzyme_psi_slope` (default 0.25, which was the hardcoded value, so earlier runs are unchanged); `Uptake_psi_slope` is 0.25 in the eight Km 2.5/2.5 input folders.
- Year 1, bacterium: biomass C 11.2 on day 365, 16.8% of substrate C degraded, E:C 25% of E_max in days 201-365, energy produced / spent 1.12 (days 101-200) and 1.27 (days 201-365), 0 cell-days with an energy shortfall, no C released by cells short of N or P. Fungus: biomass 10.5, 15.9% degraded, 0 shortfall, same energy pattern.
- Six pulses (reseeded at day 365, 3 seeds): the bacterium repeats itself (biomass on day 365 11.2-11.6, ~17% degraded, 108-112 days per pulse with energy produced < spent, no overflow, no shortfall).
- Fungus: only pulse 1 is comparable (start biomass 0.9 per cell, 10.5 on day 365). Pulses 2-6 start at C = 25 per cell in all 400 cells (reseeding places a single taxon in every cell at Microbes_init mass), reach biomass ~82-83, degrade ~100% of the substrate by day ~320, and show 13% of cell-days with an energy shortfall and energy produced / spent ~0.01 after the substrate runs out. This is a reseeding artifact, not a property of the fungus.
- Caveats: climate repeats exactly each pulse; 3 seeds of 400 identical cells; fungal pulses >= 2 need a fixed-biomass reseeding rule before they mean anything.

## 4. What to re-analyze with the simplified model before porting

| Priority | Item | Why | Suggested test |
|---|---|---|---|
| Done | Stage-A sweeps (maintenance, enzyme/uptake/HSP production and energy costs) — re-run at the settled set (`stageA_rerun/`) | Done at flat 15 C, single seed, E_resp 730. The settled energy yield (~2x) is equivalent to halving all energy costs and E_max, so collapse thresholds and rankings will shift. | Re-run vs temperature (4-40 C), HSP and control, 5 seeds, at the settled set |
| High | HSP_Consti_Ea (and HSP_Induci_Ea) | Run at the original constitutive-dominated mix, single seed, no control. HSP content is now set by the constitutive rate at 0.0003. | Re-run at the settled set; confirm inducible Ea remains null |
| Done | HSP production-rate test (3.8) | Re-run at the settled set (`hsp_prod_final/`) | Possible follow-up: replicate around Induci 0.05-0.2, Consti 0.0005-0.002 |
| High | HSP N cost (`HSP_prod_N_cost`), substrate N supply, and the shared N/E throttle (ARROEI) | Inducible HSP production is capped by N supply and the same factor throttles inducible enzyme/osmolyte production, so HSP-vs-control differences include an enzyme effect; the full model's N supply will differ | Sweep `HSP_prod_N_cost` and substrate N with enzyme production tracked; check the constitutive throttle (ARROEHC) the same way |
| High | Enz_Prod, Constit_Prod, Uptake_C_cost, transporter costs vs temperature | Never varied with temperature or against a control; Uptake_C_cost was fragile at 15 C | One-parameter x temperature sweeps at the settled set |
| Medium | n (Ea saturation) and slope -60k vs -90k at the final set | n fixed at 0.2 since stage B; slope -90k is near the Ea = 0 limit | Small grid with seeds |
| Done (settled set) | Arrhenius vs MMRT (`Kinetics_Model` 0 vs 1) | Gradient and stage-A sweeps repeated under Arrhenius (see 3.8) | Full model: enzymes have different Ea / MMRT parameters, so repeat key checks there |
| Medium | Multi-substrate / N and P limitation | Minimal scenario has one substrate with ample N and P, so only energy limits growth | Needs the full model or a 2-3 substrate minimal variant |
| Medium | Competition among taxa and relative drought tolerance | Tolerance is relative across taxa (0 for a single taxon); HSP vs non-HSP competition cannot be seen in one taxon | Two-taxon minimal test, or test in the full model |
| Medium | Horizon | Biomass still growing at day 365 (67-77% substrate left at 15.7 C); gaps are year-1 snapshots | Check key conclusions with 2-3 pulses of real-climate runs |
| Low | Cooling scenarios (-2, -4 C) | Cold penalty not yet tested under real climate | Real-climate runs with negative offsets |
| Not examined | Osmolyte parameters, `death_rate_bac`, moisture/drought parameters, `Maint_E_Rate_Tref` | Excluded by design | Revisit in the full model |

Findings unlikely to change (mechanistic, held across baselines): direction of the maintenance-Ea parameters (steeper slope, smaller n, lower intercept all help; intercept strongest);
HSP cost monotone harmful; constitutive HSP Ea higher is harmful; inducible Ea null at low inducible rates; CUE set by the energy balance; collapse = energy exhaustion.

## 5. Caveats that apply to everything above

- 400 identical clones synchronise: step-like events (E:C and CUE jumps, biomass dips) and +/- 3-5 point jitter in biomass-change values; 8 C and 4 C responses are least reliable.
- Percent change vs a collapsing control is inflated (e.g. +305% at 36 C is +3.6 absolute biomass units); prefer absolute differences for warm-side claims.
- Early sweeps (stages A, B) were single-seed; final screens used 5-8 seeds.
- Corrected during the work: the early N-limitation diagnosis (real cause: energy exhaustion); the early 100-day death attributed to nutrients (drought); the claim that inducible HSP is the
  larger pathway (not at baseline); "control collapses at 28-36 C" (true at 36 C, depleted but not collapsed at 28 C).
- Code changes that affect the full model: `monomer.py` deterministic identity monomer-transporter mapping when `n_uptake == n_monomers` and `Uptake_per_monomer == 1`;
  `grid.py` publishes `HSP_prop`, `maint_E_cost_Ea`, `tax_maint_E_cost` (diagnostic only).

## 6. Where the results are

`analysis_output/`: `maint_sweep`, `enzprod_sweep(_fine)`, `hspprod_sweep`, `hspprodEcost_sweep`, `enzprodEcost_sweep`, `uptakeprod_sweep`, `uptakeprodEcost_sweep` (stage A);
`maint_thermal`, `maint_hsp_slope_thermal`, `maint_hsp_n_thermal`, `hsp_consti_ea_thermal`, `hsp_induci_ea_thermal` (B1); `hsp_inducible_*_thermal` (B2);
`hsp_window_screen1/2`, `hsp_window_replicates`, `cue_levers`, `monomer_energy`, `cue_window_replicates`, `eresp_window_replicates(2,3)`, `slope_intercept_replicates`,
`eresp_int40_replicates`, `eresp_int_high_replicates`, `hsp_prod_e1500` (C); `cue_timeseries`, `gradient_730`, `climate_run_*`, `climate_warming` (D).
Scripts in `analysis/`: `sweep_minimal.py`, `sweep_temperature.py`, `sweep_hsp_effect.py`, `plot_*.py`, `cue_timeseries.py`, `gradient_730.py`, `climate_run.py`, `climate_warming_summary.py`.

## 4. Maintenance sensitivity with endogenous respiration off vs on (2026-10-15)

Terminology: **endogenous respiration (R_endo)** = respiring free biomass C (at the respiration energy yield of the necromass monomer) to pay energy shortfalls; switch `R_endo` in `parameters.csv` (0 = off: shortfalls are floored away and the energy is not accounted for; named `E_deficit_erosion` until 2026-10-16).

**Setup (`maint_endo_sweep/`; runner `analysis/maint_endo_sweep.py`, figures and table `analysis/plot_maint_endo.py`):** model at commit 22cbb3c (intracellular HSP/osmolyte pools, C+N+P-only rebalancing, Uptake/Enzyme_psi_slope 0.25, CUE after rebalancing losses and floored at 0; re-run of the 61b611f sweep, all other outputs identical); `input_minimal_e1500_crich_km25_in` (bacterium, C-rich substrate, Km slopes 2.5, Input_NH4 0.02 / Input_PO4 0.002); `Maint_E_Rate_Const` 0.2-5 x 2.5 (15 values, stage-A multiples); flat 4 / 15 / 28 / 36 C, Psi -0.5; 20x20, 365 days, seeds 0-2; `R_endo` 0 and 1; HSP production as in the folder (no no-HSP control). 360 runs. Not directly comparable with `stageA_rerun/` (different input folder and model version).

- **Identical until the energy pool runs out.** Below the multiple at which cells first have an energy shortfall the two settings give the same results (no shortfalls occur). Shortfalls start at 2x (4 C), 4x (15 C), 3-3.5x (28 C) and 2x (36 C) baseline.
- **Beyond it the outcome is binary.** R_endo off: cells stay alive (about 400 of 400 at 15-28 C; at 4 C and at 36 C >= 4x about 277, from basal death without reproduction) at E = 0, with no enzyme production and no growth, on unpaid energy (per cell per year: about 20,000-25,000 at 15 C, 56,000-62,000 at 28 C, 25,000-51,000 at 36 C, 1,700-5,200 at 4 C). R_endo on: the shortfall is paid from biomass (1-4 C per cell per year eroded) and the population dies out: day 286 (15 C, 4x), 218 (28 C, 3.5x), 305-325 (36 C, 2x); at 4 C, 3x it is still declining at day 365 (biomass 0.12) and extinct from 3.5x (days 256-322).
- **Degradation and energy:C are the same in both settings at every multiple**: in both, enzyme production stops once E = 0 (the unpaid energy covers only maintenance and transport), so degradation is set by the enzyme pool built before the shortfall.
- **Collapse multiple** (>= 100 of 400 cells lost, seed mean): R_endo on 3x / 4x / 3.5x / 2x at 4 / 15 / 28 / 36 C; R_endo off 3x (4 C) and 4x (36 C) from basal death of non-dividing cells, none at 15 and 28 C.
- **Baseline (1x)**: degraded 4.0 / 51 / 94 / 33 %, biomass C 2.9 / 41 / 68 / 20, CUE 0.24 / 0.41 / 0.32 / 0.29, energy:C 342 / 595 / 457 / 420 at 4 / 15 / 28 / 36 C. Degradation falls steeply with maintenance at 28 and 36 C (28 C: 94% at 1x, 55% at 2x, 27% at 3x).
- **Below baseline (0.2-0.8x)**: degradation and biomass rise steadily as maintenance falls (15 C: 51 -> 68% degraded, biomass 41 -> 57; 28 C: 94 -> 99%; 36 C: 33 -> 59%; 4 C: 4.0 -> 4.6%); energy:C rises to 560-790.
- **CUE including rebalancing losses**: fast growth at low maintenance makes cells N-limited on the C-rich substrate (biomass N:C at the 0.138 minimum), and the C respired by rebalancing lowers CUE: 15 C 0.32 at 0.2x and 0.35 at 0.4x against 0.44 at 0.6x (CUE peaks at 0.6x); 28 C about 12% lower than without the rebalancing loss up to 1x (0.32 at 1x); 36 C 7% lower at 0.2x; none at 4 C or at >= 1.5x. With R_endo on, CUE is 0 at the extinction threshold and missing once no cells are left in the last 30 days.

### 4b. Same sweep without the HSP effect on maintenance Ea (2026-10-16)

**Setup (`maint_endo_sweep_slope0/`; comparison figures and table in `maint_endo_sweep_slope0/comparison/` from `analysis/plot_hsp_slope_comparison.py`):** as section 4 but `Maint_E_HSP_Ea_Slope` 0 (runner option `--set Maint_E_HSP_Ea_Slope=0`): cells still produce HSP and pay its costs (HSP share of biomass unchanged: 0.017 / 0.028 / 0.043 / 0.053 at 4 / 15 / 28 / 36 C, 1x), but maintenance Ea stays at the 45,000 intercept instead of about 38,000 / 34,100 / 29,100 / 26,100.

- **At 1x (R_endo off), with vs without the HSP effect**: degraded 4.0 / 4.0, 51 / 45, 94 / 57, 33 / 5.8 %; biomass C 2.9 / 2.8, 41 / 33, 68 / 30, 20 / 2.3; CUE 0.24 / 0.24, 0.41 / 0.38, 0.32 / 0.25, 0.29 / 0.11 at 4 / 15 / 28 / 36 C. Without the effect, 36 C cells run out of energy by about day 130 at 1x.
- **Energy shortfalls begin** at 2 / 4 / 3 / 2x with the HSP effect and 2 / 3 / 1.75 / 0.8x without it; with R_endo on the population is lost from 3 / 4 / 3.5 / 2x vs 3 / 3 / 1.75 / 0.8x (at 36 C it dies out even at the 2.5 baseline).
- **No effect at 4 C** (maintenance Arrhenius reference 273 K, so Ea hardly matters near 0 C; HSP content is also lowest there). The effect grows with temperature, so the small penalty of `Maint_E_Rate_Const` 2.5 at 15-36 C in section 4 is largely the HSP-lowered Ea.

## 5. Minimal model under the climate scenarios, endogenous respiration off vs on (2026-10-16)

**Setup (`climate_endo/`; runner `analysis/climate_endo_run.py`, figure and table `analysis/plot_climate_endo.py`):** model at commit 22cbb3c; `input_minimal_e1500_crich_km25_in` (bacterium, baseline maintenance and HSP); first year of `input/climate_data/climate_ambient_temp.csv`, `climate_moderate_warm.csv` (+5 C) and `climate_extreme_warm.csv` (+10 C), all with the same Psi record (below -1 MPa on 223 days, -4 to -6.5 MPa around days 200-295); one 365-day pulse, 20x20, seeds 0-2; R_endo off and on.

- **R_endo makes no difference**: no cell had an energy shortfall in any run (energy:C of living biomass stays at or above about 190, lowest in the driest weeks around day 270), so the two settings are identical.
- **Drought controls the year.** Biomass peaks around day 150-170 and falls through the dry season; drought mortality cuts the living cells to 84 / 148 / 187 of 400 (ambient / moderate / extreme) around day 292, and they recover after rewetting (day ~295). Drought tolerance is about 0.38 throughout (osmolyte share 0.027 of biomass with Osmo_tol_K 0.05).
- **Warming increases decomposition**: day-365 substrate degraded 8.0 / 13.9 / 17.9 %, biomass C 0.93 / 2.3 / 3.1, cells alive 234 / 395 / 399 (seed SD of degradation <= 0.2 points). Warmer runs grow more before the dry season and regrow faster after it. CUE (last 30 days) is 0.48 in all three; daily CUE falls from about 0.48 to 0.2-0.25 in the dry season.

### 5b. Fungus under the climate scenarios (2026-10-16)

**Setup (`climate_endo_fungus/`):** as section 5 with `input_minimal_e1500_crich_km25_fb1_in` (identical parameters, `fb` 1); model at commit ed80ea0 (R_endo renaming only).

- The fungus starts in 14 of 400 cells (0.9 C per grid cell, the same total as the bacterium); placement scales fungal occupancy by max_size_b / max_size_f (see OPEN_ISSUES, fungal distribution).
- **R_endo makes no difference** (no shortfalls; energy:C of living biomass at or above about 140).
- Day 365: substrate degraded 7.9 / 13.6 / 17.1 % (bacterium 8.0 / 13.9 / 17.9), biomass C 1.13 / 2.6 / 3.2, occupied cells 12 / 29 / 34 (3-9 at the end of the dry season), CUE 0.49 / 0.49 / 0.48 (ambient / +5 C / +10 C). Per grid cell the fungus decomposes about as much as the bacterium.

### 5c. Six pulses (six identical climate years) (2026-10-16)

**Setup (`climate_endo_6p/` bacterium, `climate_endo_fungus_6p/` fungus):** as 5 / 5b with 6 pulses of 365 days (`climate_endo_run.py --pulses 6`; `sweep_minimal.run_one` now runs several pulses with `Grid.repopulation` between them, as `dementpy.py`); per-pulse table `climate_endo_pulses.csv`.

- **With one taxon, nothing carries over between pulses.** Repopulation resets substrates, monomers, enzymes and the HSP/osmolyte pools and re-places the community from the taxa's cell counts; with a single taxon every cell is refilled with the initial cell. The climate years are identical, so the bacterial pulses are replicates: day-365 degradation 8.0-8.2 / 13.8-14.0 / 17.6-18.0 % across pulses, R_endo off = on; pulse 1 reproduces section 5.
- **Fungus, after the repopulation fix (commit a42b5f4; `climate_endo_fungus_6p/`)**: pulses 2-6 start in 12-20 occupied cells (0.8-1.3 C per grid cell; pulse 1: 14 cells, 0.9) and are replicates of pulse 1: day-365 degradation 7.9-8.7 / 13.6-15.7 / 16.1-19.4 % (ambient / +5 C / +10 C; more pulse-to-pulse spread than the bacterium because few cells are occupied), no energy shortfalls, R_endo off = on.
- **Before the fix (`climate_endo_fungus_6p_before_repop_fix/`)** repopulation refilled all 400 cells at the initial fungal cell mass (25 C per grid cell, 28 times the pulse-1 start): degradation 79-80 / 97 / 99 %, the warm runs ran out of energy after the substrate was gone (shortfalls on 29-68 days per pulse at +5 C, 139-140 at +10 C), and R_endo lowered day-365 biomass C (9.3-11.1 vs 12.9-14.5 at +5 C; 1.6-1.9 vs 3.4-4.1 at +10 C). Those results are an artifact and are kept only for reference.

## 6. Why decomposition is low under the real climate: moisture test (2026-10-16)

**Setup (`moisture_test/`; figure and table `analysis/plot_moisture_test.py`):** bacterial minimal model as section 5 (ambient climate, first year, 20x20, 3 seeds, R_endo on); variants: real climate; Psi held at -0.5 MPa with the real temperature (`moisture_test/climate_ambient_wetpsi.csv`); real climate with `Enzyme_psi_slope` 0, with `Uptake_psi_slope` 0, with both 0, and with both 0 plus `Osmo_tol_K` 1e-9 (drought tolerance about 1, so no drought mortality); fungus with Psi -0.5.

| Variant (bacterium) | Degraded, day 365 | Max biomass C | Cells killed |
|---|---|---|---|
| Real climate | 8.0 % | 4.1 | 1,026 |
| No moisture effect on enzymes | 14.5 % | 6.6 | 1,317 |
| No moisture effect on uptake | 8.3 % | 4.2 | 1,121 |
| No moisture effect on enzymes or uptake | 16.7 % | 6.0 | 1,442 |
| ... and no drought mortality | 42.6 % | 33.6 | 145 |
| Psi held at -0.5 MPa | 43.7 % | 35.2 | 143 |

- **Moisture is the limit (8.0 -> 43.7 % without dry soil; fungus 7.9 -> 34.4 %).** The dry season (Psi -4 to -6.5 MPa, days ~200-295) falls in the warm months, when growth and decomposition would be fastest. Two moisture effects act together: drought mortality (tolerance about 0.38, death probability about 0.026 per cell per day; the population falls to 84 of 400 cells) and the moisture factor on enzyme Vmax. Removing either alone gains about 7 points (14.5 % without the enzyme factor; 15.3 % without drought mortality, `beta_bac` = 0, section 7); removing both gives 42.6 %. (Corrected 2026-10-16: the first version attributed ~26 points to drought mortality, which holds only with the moisture factors already off.)
- The moisture factor on uptake adds almost nothing on its own (8.3 %), because the monomer supply is set by decay, not by uptake capacity.
- Up to about day 150 all variants are alike (about 5-6 % degraded): the first months are slow in every case (cool spring, small starting biomass). The remaining 1 point between "no drought mortality" and Psi -0.5 is the inducible osmolyte response and its costs.
- The size of the drought effect is therefore set mainly by drought tolerance (`Osmo_tol_K`, osmolyte production and the mortality parameters `death_rate_bac`/`beta_bac`), which are uncalibrated (OPEN_ISSUES, calibration).

## 7. Drought-mortality parameters under the scenario climates (2026-10-16)

**Setup (`drought_sweep/`; `analysis/drought_sweep.py`, `analysis/plot_drought_sweep.py`):** bacterial minimal model, first year, 20x20, 2 seeds, R_endo on, ambient / +5 C / +10 C climates; one at a time: `beta_bac` (= `beta_fun`) 10 (baseline), 5, 2.5, 1, 0 and `Osmo_tol_K` 0.05 (baseline), 0.02, 0.01, 0.005, 0.0025. Drought mortality per cell per day is `death_rate * (1 - (1 - tol) * (Psi - wp_fc) * beta)` with `tol = s / (s + Osmo_tol_K)` (s = osmolyte share of biomass, about 0.028 here), so the excess mortality is set by the product `(1 - tol) * beta` (baseline 0.62 x 10 = 6.2). `wp_fc` and `death_rate` were not varied: `wp_fc` also sets the enzyme, uptake, osmolyte-induction and leaching responses, and `death_rate` the basal mortality.

| Setting | (1 - tol) x beta | Ambient | +5 C | +10 C | Cells killed (ambient) |
|---|---|---|---|---|---|
| Baseline (beta 10, K 0.05) | 6.2 | 8.0 % | 13.9 % | 18.0 % | 1,026 |
| beta 5 | 3.1 | 10.2 | 17.5 | 22.8 | 776 |
| beta 2.5 | 1.6 | 12.0 | 21.3 | 27.4 | 468 |
| beta 1 | 0.6 | 13.7 | 25.3 | 30.6 | 284 |
| beta 0 (no drought mortality) | 0 | 15.3 | 27.6 | 24.6 | 161 |
| K 0.02 (tol 0.60) | 4.0 | 9.6 | 16.0 | 21.1 | 914 |
| K 0.01 (tol 0.75) | 2.5 | 10.8 | 19.3 | 25.0 | 656 |
| K 0.005 (tol 0.86) | 1.4 | 12.2 | 22.3 | 27.6 | 452 |
| K 0.0025 (tol 0.93) | 0.7 | 13.9 | 25.2 | 30.7 | 296 |

(Day-365 % of the substrate degraded, mean of 2 seeds; seed spread 0.1-0.5 points except beta 0 at +10 C, 22.3 and 26.9.)

- **Both parameters act through `(1 - tol) * beta`, and settings with the same product give about the same result** (beta 2.5 and K 0.005: 12.0/21.3/27.4 and 12.2/22.3/27.6 %). They differ in what tolerance responds to: with a low `Osmo_tol_K` tolerance saturates (0.86-0.93), so differences in osmolyte content between cells, and extra osmolyte production, buy little tolerance; with a lower `beta` tolerance keeps its baseline sensitivity to osmolyte content.
- **Relieving drought mortality alone recovers only part of the gap**: at most about 15 % at ambient (vs 43.7 % with Psi held at -0.5), because the enzyme moisture factor still slows decay through the dry season (the curves separate only after about day 190 and all flatten during the dry season).
- **The warming response is kept, and gets larger in absolute terms**: +10 C vs ambient is 2.2-2.3x at every setting from beta 10 to beta 1 (18.0 vs 8.0; 30.6 vs 13.7 %).
- **Some drought mortality is needed at +10 C.** With beta 0 the population fills all 400 cells and outgrows its monomer supply after the dry season: biomass becomes N-rich (N:C 0.33), enzyme production and CUE fall to zero for 43 days, and degradation (24.6 %) ends below beta 1 (30.6 %). Mortality supplies turnover and necromass that keep the community productive.
- No energy shortfalls in any run (R_endo has no effect).

## 8. Moisture slopes of enzyme activity and uptake (2026-10-16)

**Setup (`psi_slope_sweep/`; `analysis/param_sweep.py`, `analysis/plot_param_sweep.py`):** bacterial minimal model, first year, 20x20, 2 seeds, R_endo on, beta 1 (new baseline), ambient / +5 C / +10 C; one at a time `Enzyme_psi_slope` and `Uptake_psi_slope` 0, 0.1, 0.25 (baseline), 0.5, 1. Below wp_fc each process is multiplied by exp(slope x (Psi - wp_fc)); at Psi -5 MPa that is 1 / 0.67 / 0.37 / 0.14 / 0.02. Mortality is reported as dead biomass C (all deaths) and drought kills; the new step diagnostics `Grid.Starved` and `Grid.Dead_C` (no effect on results; baseline runs identical to section 7, beta 1) also count threshold (starvation) deaths.

| Setting | Degraded day 365 (amb / +5 / +10, %) | Mean biomass C (amb / +10) | Dead C, year (amb / +10) |
|---|---|---|---|
| Enzyme 0 | 33.5 / 65.4 / 80.1 | 8.2 / 17.7 | 5.5 / 12.2 |
| Enzyme 0.1 | 21.2 / 40.9 / 50.6 | 5.7 / 11.9 | 4.2 / 9.1 |
| Baseline (0.25 / 0.25) | 13.7 / 25.3 / 30.6 | 4.3 / 8.3 | 3.6 / 6.6 |
| Enzyme 0.5 | 6.8 / 10.9 / 13.7 | 3.0 / 4.7 | 2.5 / 4.0 |
| Enzyme 1 | 5.3 / 8.4 / 10.6 | 2.4 / 3.8 | 1.9 / 2.9 |
| Uptake 0 | 14.2 / 25.5 / 26.1 | 4.5 / 7.9 | 3.6 / 6.9 |
| Uptake 0.1 | 14.2 / 25.6 / 25.3 | 4.5 / 7.9 | 3.6 / 6.1 |
| Uptake 0.5 | 12.2 / 22.0 / 23.7 | 4.0 / 6.8 | 2.9 / 5.7 |
| Uptake 1 | 9.1 / 15.4 / 18.8 | 2.9 / 5.0 | 2.2 / 3.9 |

- **`Enzyme_psi_slope` is the dominant moisture control on decomposition.** Degradation falls from 80 to 11 % (+10 C) across the range; the curves separate only after about day 150, as the soil dries. Without a moisture effect on enzymes (slope 0) the +10 C run degrades 80 % in the first year.
- **`Uptake_psi_slope` matters little up to the baseline**; at 0.5-1 it lowers degradation by 2-12 points, with biomass falling through the dry season and recovering after rewetting.
- **Energy trap at the end of the dry season.** With a steep enzyme slope (0.5, 1) the community's energy pool reaches zero during the dry season (E:C 450 at day 180 -> 0 by about day 220 at slope 1), enzyme production stops, the enzyme pool decays, and the community does not recover after rewetting: no enzymes, so no monomers, so no energy (biomass shrinks to 0.2-0.9 per grid cell through R_endo; at +10 C slope 1 also 276 starvation deaths). The same trap appears at +10 C with uptake slope 0-0.1 in one of two seeds (22 vs 30 % degraded): without the uptake restriction the large dry-season community depletes the monomer pool (0.4 vs 1.4), its energy pool hits zero at about day 290, and it never restarts (as with beta 0, section 7). Recovery after drought therefore depends on whether the energy reserve lasts through the dry season.
- **CUE** drops through the dry season at every setting (to about 0.1 at baseline), the more the steeper the slopes; at enzyme slope >= 0.5 or uptake slope 1 it reaches 0. Annual uptake-weighted CUE is not a useful summary here: when dry-season uptake is negligible the annual value reflects the wet season only (it is lowest at the baseline).
- **Mortality** with beta 1 is mostly basal: drought kills are 270-320 cells per year at every setting (no trend), and dead biomass C follows biomass.

## 9. What limits degradation (2026-10-16)

**Setup (`degradation_limits/`; `analysis/param_sweep.py`, `analysis/plot_limit_sweep.py`):** bacterial minimal model, first year, 20x20, 2 seeds, R_endo on, new baseline `Enzyme_psi_slope` = `Uptake_psi_slope` = 0.1 and beta 1; ambient / +5 C / +10 C. New step diagnostics in `Grid.Diag` (no effect on results; checked identical): enzyme saturation S / (Km + S) on the largest substrate; uptake capacity of the community relative to the organic monomer pool (>= 1: the whole pool can be taken up in a day); transporter saturation M / (Km + M); mineral N and P pools. Test runs, one at a time: `Uptake_Vmax_Km` (uptake Km = slope x transporter Vmax) 2.5 (baseline), 1, 0.5, 0.25, 0.1, 0.025, 0.0025; `Enz_Prod_min`/`max` 0.002 (baseline), 0.004, 0.008; `Vmax_Km` (enzyme Km slope) 2.5, 0.25; and `Uptake_Vmax_Km` 0.25 with `Enz_Prod` 0.004 and 0.008.

**Diagnosis at the baseline:**
- **Degradation is enzyme-limited.** Enzymes are 92-96 % saturated by substrate all year (S about 200 vs Km 10-18 per grid cell), so decay is close to enzyme pool x Vmax, and the enzyme pool (0.03-0.07 per grid cell) follows biomass.
- **Growth is limited by uptake kinetics, not by monomer supply.** Transporters are about 1 % saturated (monomer pool about 1.2 vs uptake Km 35-60 per grid cell, because Km = 2.5 x Vmax): the community can take up only 2-30 % of the monomer pool per day (capacity ratio 0.10 / 0.17 / 0.19 annual mean), and organic monomers accumulate (2.5-3.5 per grid cell in the first 50 days, about 1.2 later). At low M the uptake rate is about Vmax / Km x M = M / 2.5 per unit transporter, so raising transporter Vmax does not help; only the Km slope (or intercept) does.
- Not N- or P-limited: NH4 and PO4 accumulate, and biomass N:C rises above its initial value (0.19 -> 0.27-0.33; C-limited cells).
- With slopes 0.1 biomass already grows during the dry season (+2.1 / +4.4 / +4.1 C per grid cell from the start to the end of Psi < -4 MPa).

| Run | Degraded day 365 (amb / +5 / +10, %) | Capacity ratio (amb) | Dry-season biomass change (amb / +10) | Min E:C (amb) |
|---|---|---|---|---|
| Baseline (uptake Km slope 2.5) | 22.9 / 42.6 / 52.2 | 0.10 | +2.1 / +4.1 | 424 |
| Uptake Km slope 1 | 26.1 / 52.2 / 63.9 | 0.31 | +2.5 / +5.6 | 377 |
| Uptake Km slope 0.5 | 26.9 / 54.8 / 68.3 | 0.64 | +2.0 / +5.7 | 314 |
| Uptake Km slope 0.25 | 25.3 / 55.6 / 69.8 | 1.18 | +1.3 / +6.0 | 179 |
| Uptake Km slope 0.1 | 20.4 / 54.8 / 70.7 | 2.3 | +0.3 / +6.1 | 1-21 |
| Uptake Km slope 0.025 | 17.4 / 56.4 / 70.9 | 6.4 | -3.8 / +5.9 | 0 (119 starved) |
| Enz_Prod 0.004 | 29.1 / 58.0 / 73.6 | 0.13 | +3.8 / +7.7 | 461 |
| Enz_Prod 0.008 | 51.4 / 98.0 / 99.6 | 0.18 | +8.3 / +20.3 | 527 |
| Enzyme Km slope 0.25 | 29.2 / 59.6 / 76.7 | 0.13 | +3.9 / +8.0 | 582 |
| Uptake Km slope 0.25 + Enz_Prod 0.004 | 35.3 / 79.8 / 93.1 | 1.56 | +3.6 / +10.4 | 380 |
| Uptake Km slope 0.25 + Enz_Prod 0.008 | 74.7 / 100 / 100 | 2.59 | +13.1 / +14.5 | 510 (warm runs exhaust the substrate, then E -> 0) |

- **Lowering the uptake Km slope removes the uptake limitation**: at 0.25 the capacity ratio is about 1 at ambient (1.1-1.2 in the dry season) and 2-3 when warm, and the monomer pool is drawn down to about 0.06. The warm runs gain 13-18 points; ambient only 2-4 points, because biomass is then set by the monomer supply, i.e. by decay, i.e. by the enzyme pool.
- **Too fast uptake destabilises the ambient run**: at slope 0.1 and below the supply-limited community has no monomer buffer, its energy reserve runs out in the dry season (one of two seeds at 0.1; both at 0.025 and 0.0025, with 117-225 starvation deaths) and it does not recover (the energy trap of section 8). Slope 0.25-0.5 is the range where uptake is (nearly) non-limiting and the ambient run keeps an energy margin.
- **Once uptake is not limiting, enzyme production sets degradation**: doubling `Enz_Prod` with uptake slope 0.25 gives 35 / 80 / 93 %; the enzyme Km slope helps less (saturation is already above 0.9; its effect comes through the growth feedback).
- Ambient stays low in every case mainly because enzyme activity is slow at the cool ambient temperatures (mean enzyme Vmax about 4-7 per unit enzyme vs higher in the warm runs) and the community starts small (1 C per grid cell; the first ~150 days are slow in every run).

## 10. Enzyme Km slope (`Vmax_Km`), ambient climate (2026-10-16)

**Setup (`enzyme_km_sweep/`; `analysis/param_sweep.py`, `analysis/plot_param_sweep.py`):** bacterial minimal model, first year, ambient climate only, 20x20, 3 seeds, R_endo on; baseline slopes 0.1, beta 1, `Uptake_Vmax_Km` 0.25 (set in the km25 inputs, commit 284d31d); `Vmax_Km` 5, 2.5 (baseline), 1, 0.5, 0.25, 0.1, 0.025, 0.0025. Enzyme Km = `Vmax_Km` x enzyme Vmax (MMRT, before the moisture factor; floor `Km_min` 0.01), so Km also rises with temperature.

| Vmax_Km | Degraded, day 365 | Mean biomass C | Mean enzyme pool | Enzyme saturation | Dry-season biomass change | Min E:C |
|---|---|---|---|---|---|---|
| 5 | 18.5 % | 5.3 | 0.028 | 0.88 | -0.6 | 11 |
| 2.5 (baseline) | 25.6 % | 7.0 | 0.036 | 0.94 | +1.3 | 189 |
| 1 | 31.2 % | 8.3 | 0.042 | 0.97 | +2.7 | 307 |
| 0.5 | 33.4 % | 8.8 | 0.044 | 0.99 | +3.3 | 338 |
| 0.25 | 34.5 % | 9.1 | 0.045 | 0.99 | +3.5 | 359 |
| 0.1 | 35.8 % | 9.4 | 0.047 | 1.00 | +3.8 | 371 |
| 0.025 | 36.4 % | 9.6 | 0.048 | 1.00 | +4.0 | 376 |
| 0.0025 | 36.7 % | 9.6 | 0.048 | 1.00 | +4.1 | 377 |

(Seed spread about 1 point.)

- **Small changes in saturation have a large effect through growth.** Saturation goes only from 0.94 to 1.00 between slope 2.5 and 0.1, but degradation rises from 25.6 to 35.8 %: with uptake no longer limiting, growth follows the monomer supply, so faster decay means more biomass, more enzyme and faster decay again (the curves separate from about day 100 and diverge most after the dry season).
- **The response saturates below a slope of about 0.25-0.5** (34-37 %): enzymes are then fully saturated and decay is enzyme pool x Vmax, independent of Km and of the amount of substrate left.
- **Saturation is lowest in the warm months** (dry season), because Km follows the MMRT Vmax (Enz_Km_temp_sens = 1) and rises with temperature, while the moisture factor lowers Vmax but not Km. With a high slope this offsets part of the temperature response of decay; with a low slope it vanishes.
- **Slope 5 approaches the dry-season energy trap** (biomass falls in the dry season, min E:C 11, CUE near 0 at the end of the dry season).
- With S about 200 per grid cell, even the baseline is close to zero-order kinetics; at low slopes decay no longer slows as the litter is used up until very little is left.

## 11. Enzyme production, ambient climate (2026-10-16)

**Setup (`enzyme_production_sweep/`; `analysis/param_sweep.py`, `analysis/plot_param_sweep.py`):** bacterial minimal model, first year, ambient climate, 20x20, 3 seeds, R_endo on; baseline slopes 0.1, beta 1, `Uptake_Vmax_Km` 0.25, `Vmax_Km` 1 (set in the km25 inputs, commit d421c72). One at a time: constitutive production `Constit_Prod_min`/`max` (enzyme C per unit biomass C per day) 0.0001-0.005 (baseline 0.0002); inducible production `Enz_Prod_min`/`max` (enzyme C per unit of C taken up into biomass) 0.002 (baseline)-0.1; `Enz_loss_rate` (enzyme turnover per day) 0.04 (baseline), 0.03, 0.02, 0.01. New diagnostics in `Grid.Diag`: production by route and the joint C/N/E limitation factors of constitutive and inducible production (no effect on results; checked identical).

**Baseline:** enzyme production is never resource-limited (limitation factors 1.0 all year), so it is set by the rate parameters alone. Constitutive production supplies about 94 % of it (0.0002 x biomass C per day vs 0.002 x biomass-C uptake). Enzyme C is 0.9 % of C uptake; the enzyme pool (0.04-0.09 per grid cell) is about 0.5 % of biomass; day-365 degradation 31.2 %.

| Setting | Degraded day 365 | Day 50 % degraded | Max biomass C | Enzyme C / C uptake |
|---|---|---|---|---|
| Constit_Prod 0.0001 (x0.5) | 5.7 % (community collapses, 280 starvation deaths) | - | 4.1 | 0.6 % |
| Constit_Prod 0.00015 (x0.75) | 10.9 % (collapses) | - | 5.8 | 0.7 % |
| **Baseline** | **31.2 %** | - | 18.9 | 0.9 % |
| Constit_Prod 0.0003 (x1.5) | 99.1 % | 283 | 79.5 | 1.3 % |
| Constit_Prod 0.0005 (x2.5) | 100 % | 183 | 64.0 | 2.0 % |
| Constit_Prod 0.002 (x10) | 100 % | 87 | 48.7 | 5.1 % |
| Enz_Prod 0.005 (x2.5) | 53.5 % | 355 | 39.4 | 1.0 % |
| Enz_Prod 0.01 (x5) | 100 % | 257 | 78.6 | 1.4 % |
| Enz_Prod 0.05 (x25) | 100 % | 52 | 37.9 | 4.1 % |
| Enz_loss_rate 0.03 | 82.6 % | 308 | 64.7 | 0.8 % |
| Enz_loss_rate 0.02 | 100 % | 226 | 73.8 | 1.1 % |

- **The response is a threshold, not a gradient.** A factor of 2 in constitutive production separates collapse (x0.5, 5.7 %) from complete degradation within the year (x1.5, 99 %); the same holds for inducible production and enzyme turnover.
- **Why:** with enzymes saturated by substrate and uptake not limiting, nothing in the model depends on population density until the substrate runs out. Enzyme production, decay and uptake all scale with biomass, so the per-capita growth rate is fixed by enzyme per biomass, enzyme Vmax (temperature), yield and maintenance, and the community grows or declines exponentially (per-capita growth days 60-180: 0.006 per day at baseline, 0.0105 at x1.5, 0.003 at x0.75 before the collapse). Small changes in that rate compound over the year. Below the threshold the community cannot cover maintenance, its energy pool runs out (E:C 0 by about day 200 at x0.5) and it does not recover (the energy trap).
- The lower uptake Km (section 9) and enzyme Km (section 10) slopes removed two density-dependent feedbacks (monomer competition at the transporters, substrate saturation of enzymes), which is part of why the system is now this sensitive. The same compounding will apply to temperature, through the MMRT Vmax: warming responses are likely to be threshold-like as well.
- After the substrate is exhausted the community starves (closed litter system): biomass, CUE and enzyme saturation fall to zero, and the runs with the fastest degradation have starvation deaths late in the year.
- Inducible production is somewhat self-limiting (tied to uptake): at Enz_Prod 0.1 decay outpaces uptake and organic monomers briefly accumulate to about 80 per grid cell.

## 12. Enzyme Km slope with the higher enzyme production; shape of the degradation curve (2026-10-16)

**Setup (`enzyme_km_sweep_prod/`):** as section 10 but with `Constit_Prod` 0.00025 and `Enz_Prod` 0.005 (minimal km25 inputs, commit e72b8d9); ambient climate, 3 seeds; `Vmax_Km` 0.25, 0.5, 1 (baseline), 2.5, 5, 10, 20. Shape is summarised as the % of the substrate degraded in each quarter of the year (days 1-91, 92-182, 183-273, 274-365). Diagnostic (`analysis/initial_biomass_shape_test.py`, not a model option): starting biomass x5 and x20 at `Vmax_Km` 2.5, 5 and 10, 1 seed.

| Vmax_Km | Degraded day 365 | Day 25 / 50 / 75 % | Degraded per quarter (%) | Max biomass C |
|---|---|---|---|---|
| 0.25 | 99.9 % | 206 / 277 / 314 | 6 / 13 / 29 / 51 | 80 |
| 0.5 | 99.8 % | 208 / 280 / 317 | 6 / 13 / 28 / 53 | 81 |
| **1 (baseline)** | **98.2 %** | 213 / 288 / 325 | 6 / 12 / 26 / 54 | 78 |
| 2.5 | 83.7 % | 228 / 307 / 351 | 5 / 11 / 21 / 47 | 66 |
| 5 | 61.0 % | 253 / 339 / - | 5 / 9 / 15 / 32 | 46 |
| 10 | 34.4 % | 313 / - / - | 4 / 6 / 9 / 15 | 21 |
| 20 | 11.6 % (collapse in the dry season) | - | 3 / 4 / 4 / 1 | 6 |

| Starting biomass | Vmax_Km | Degraded day 365 | Degraded per quarter (%) |
|---|---|---|---|
| x1 | 2.5 / 5 / 10 | 83.5 / 61.0 / 35.0 % | 5/11/21/47; 5/9/15/32; 4/6/9/15 |
| x5 | 2.5 / 5 / 10 | 100 / 98.1 / 74.4 % | 10/23/42/24; 9/20/32/37; 8/15/21/31 |
| x20 | 2.5 / 5 / 10 | 100 / 99.9 / 98.7 % | 26/62/12/0; 24/53/22/0; 21/41/30/6 |

- **The Km slope sets how much is degraded, not the shape.** `Vmax_Km` 2.5 gives about 84 % in the first year, but at every slope degradation accelerates through the year (about half of it in the last quarter): biomass grows exponentially from the 1 C per grid cell inoculum until the substrate runs out (no plateau; biomass at day 365 is the maximum).
- A higher Km slope damps decay as the substrate is used up (saturation 0.98 -> 0.8-0.9), but because decay stays proportional to the enzyme pool, and the enzyme pool to biomass, it does not create a plateau in biomass. Very high slopes (20) push the community into the dry-season energy trap.
- **A larger starting biomass shortens the slow start but does not produce a steady decline**: the curve is still sigmoidal, with the fast phase earlier (x5) or the substrate exhausted in the first half-year (x20).
- A steady decline needs density dependence that the model now lacks: something that lowers the per-capita return as biomass or the enzyme pool grows (for example decay saturating in enzyme, i.e. reverse Michaelis-Menten or ECA kinetics; uptake competition that depends on biomass; or density-dependent mortality).

## 13. Decay with enzyme saturation (runtime switch `Decay_MM_Form`) (2026-10-16)

**New option (`src/grid.py`, `Grid.degradation`):** runtime.txt row `Decay_MM_Form` (provisional name; absent = 0). 0 = forward Michaelis-Menten, Decay = Vmax E S / (Km + S), saturating in substrate and linear in enzyme (the original model; checked bit-identical). 1 = the same times an enzyme-saturation factor K_E S / (K_E S + E), with K_E = `Enz_Sat_K` (parameters.csv; enzyme C per unit substrate C+N+P; E = enzyme acting on the substrate). With scarce enzyme it equals form 0; with enzyme in excess decay approaches Vmax K_E S^2 / (Km + S), about Vmax K_E S for S >> Km, i.e. first order in substrate and no longer increasing with enzyme (the reverse Michaelis-Menten limit), so Vmax K_E is the maximum relative decay rate per day and keeps the MMRT temperature and moisture responses. All input folders now carry `Decay_MM_Form 0` and an inactive `Enz_Sat_K` 0.002. Uptake is unchanged (forward Michaelis-Menten). New step diagnostic `d_enz_sat_factor` (the factor on the largest substrate).

**Runs (`enzyme_saturation/`; ambient unless stated, minimal model, first year, R_endo on, Enz_Prod 0.005, Vmax_Km 2.5, Uptake_Vmax_Km 0.25):**

| Run | Degraded day 365 | Day 50 % | Degraded per quarter (%) | Max biomass C |
|---|---|---|---|---|
| Form 0, defaults (Constit_Prod 0.00025) | 83.7 % | 307 | 5 / 11 / 21 / 47 | 66 |
| Form 1, Constit 0.00025, K_E 0.0005 / 0.001 / 0.002 / 0.005 / 0.01 / 0.05 | 14 / 27 / 38 / 52 / 62 / 76 % | - (to 0.005: 358) | e.g. K_E 0.005: 5 / 9 / 14 / 24 | 7-58 |
| Form 1, Constit 0.001, K_E 0.001 | 61.5 % | 295 | 10 / 20 / 17 / 15 | 31 |
| Form 1, Constit 0.001, K_E 0.002 | 82.1 % | 200 | 14 / 31 / 22 / 15 | 39 |
| Form 1, Constit 0.002, K_E 0.001 | 67.1 % | 250 | 15 / 23 / 16 / 14 | 27 |
| Form 1, Constit 0.002, K_E 0.002 | 86.2 % | 163 | 22 / 34 / 19 / 11 | 33 |
| Form 1, Constit 0.005, K_E 0.002 | 88.6 % | 136 | 32 / 32 / 15 / 9 | 23 |

Climate response (2 seeds), day-365 degraded (day 50 %): form 0 defaults 83 % (308) / 100 % (223) / 100 % (188) for ambient / +5 C / +10 C; form 1 Constit 0.001, K_E 0.001: 62 % (295) / 70 % (237) / 72 % (214); K_E 0.002: 82 % (200) / 87 % (161) / 88 % (146).

- **With the production rates of form 0, form 1 only slows decay** (the curve still accelerates; enzyme never far exceeds K_E S). **With higher enzyme production it gives the steady decline**: biomass and the enzyme pool level off by about day 200 (density dependence), the saturation factor falls to 0.1-0.2, and decay settles near Vmax K_E S, an exponential decline modulated by temperature and moisture.
- **The warming response becomes graded instead of threshold-like**: +10 C shortens the time to 50 % by 50-80 days instead of exhausting the substrate by day 250-300 (form 0).
- Cost: at saturation most of the extra enzyme is wasted (factor 0.1-0.2), CUE falls (0.34 vs 0.43), and the warm runs with K_E 0.002 come close to the energy trap at the end of the dry season (min E:C 0.1-2 at +5/+10 C); lower enzyme production that still saturates (for example Constit 0.0005) would cost less.

## 14. Enzyme Km slope x enzyme production options, forward Michaelis-Menten (2026-10-16)

**Setup (`km_production/`, runs in `runs.tar.gz`; `analysis/param_sweep.py`, `analysis/plot_km_production.py`):** bacterial minimal model, first year, 20x20, 2 seeds, Decay_MM_Form 0, R_endo on, Uptake_Vmax_Km 0.25, slopes 0.1, beta 1. Enzyme production scaled together (`Constit_Prod` / `Enz_Prod` = x1 0.00025 / 0.005 (defaults), x1.1, x1.2, x1.4) crossed with `Vmax_Km` 2.5-13 at ambient; five options near 80 % then run under the three climates.

Ambient, day-365 % degraded:

| Production \ Vmax_Km | 2.5 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 11 |
|---|---|---|---|---|---|---|---|---|---|
| x1 (defaults) | **83.2** | 77.8 | 68.2 | 61.1 | | | | | |
| x1.1 | 99.4 | 97.6 | 91.2 | **82.5** | | | | | |
| x1.2 | 100 | 100 | 99.6 | 97.2 | 93.2 | 86.9 | **80.0** | | |
| x1.4 | 100 | 100 | 100 | 100 | | 99.7 | | 97.3 | 91.0 |

Options near 80 % (degraded per quarter at ambient; day-365 degraded and day of 50 % at +5 / +10 C):

| Option | Ambient | Quarters (%) | +5 C | +10 C | Biomass C at day 365, +10 C |
|---|---|---|---|---|---|
| x1, Vmax_Km 2.5 (defaults) | 83.2 % (day 308) | 5 / 11 / 20 / 46 | 100 % (223) | 100 % (188) | 21 |
| x1, Vmax_Km 3 | 77.8 % (314) | 5 / 10 / 19 / 43 | 100 % (230) | 100 % (194) | 33 |
| x1.1, Vmax_Km 5 | 82.5 % (306) | 5 / 11 / 21 / 45 | 99.9 % (224) | 99.9 % (192) | 34 |
| x1.2, Vmax_Km 8 | 80.0 % (306) | 5 / 11 / 21 / 42 | 98.7 % (232) | 99.4 % (202) | 49 |
| x1.4, Vmax_Km 13 | 80.1 % (303) | 5 / 12 / 22 / 41 | 97.1 % (232) | 98.2 % (208) | 57 |

- **Production and the Km slope trade off along a ridge**: about +10 % production offsets +2.5 to +3 in `Vmax_Km` near 80 %. Options on the ridge have almost the same ambient curve (the shape is set by growth from the small inoculum, section 12).
- **Higher Km slope with higher production damps the warming response a little**: enzyme saturation is lower (0.69-0.78 vs 0.91 at ambient) and falls further as Km rises with temperature, so +10 C reaches 50 % 20 days later (208 vs 188) and the litter is not completely used up (98 vs 100 %); the community crashes less after the substrate runs out (biomass at day 365 at +10 C 57 vs 21). The ambient curve is slightly less back-loaded (last quarter 41 vs 46 %).
- The warm runs still exhaust the substrate within the year, and min E:C reaches about 0 at +10 C in every option (end-of-year starvation once the substrate is gone).

## 15. Multi-year runs of the chosen parameters; MMRT vs Arrhenius (2026-10-16)

**Setup (`multiyear/`, runs in `runs.tar.gz`; `analysis/climate_endo_run.py` now takes `--runtime`, `analysis/plot_multiyear.py`):** bacterial minimal model with the chosen parameters (production x1.4: `Constit_Prod` 0.00035, `Enz_Prod` 0.007; `Vmax_Km` 13; minimal km25 inputs, commit 5fbae43), Decay_MM_Form 0, R_endo on, three climates. (a) Six yearly pulses (litter reset and community re-placed each year; 2 seeds); (b) one continuous 3-year pulse (no reset; 1 seed). Each with `Kinetics_Model` 1 (MMRT) and 0 (Arrhenius; same Vmax0, Ea and Km slopes).

| Kinetics | Climate | Degraded day 365 (pulses 1-6) | Day 50 % | Biomass C at day 365 |
|---|---|---|---|---|
| MMRT | Ambient | 80.1-81.8 % | 300-303 | 61-62 |
| MMRT | +5 C | 96.9-97.1 % | 232-234 | 65-66 |
| MMRT | +10 C | 98.0-98.2 % | 205-208 | 54-57 |
| Arrhenius | Ambient | 94.2-95.0 % | 258-261 | 69-72 |
| Arrhenius | +5 C | 99.0-99.1 % | 199-201 | 52-54 |
| Arrhenius | +10 C | 98.2 % | 163-165 | 4-5 |

- **The six pulses are replicates** with both kinetics (repopulation restores the t0 community; identical climate years): no drift between years.
- **Arrhenius decomposes faster at the same parameters**: 95 vs 81 % at ambient, 50 % reached 43 days earlier. Both rate constants are similar in the cool season (enzyme Vmax 4.6 vs 4.2 at 11 C) but Arrhenius keeps rising with temperature (mean enzyme Vmax 12.0 vs 6.8 above 30 C), whereas MMRT flattens above about 20-25 C and falls at the highest temperatures; transporter Vmax is also about 40 % higher under Arrhenius (22 vs 16 annual mean at ambient).
- **The warming shift in timing is similar** (day 50 % 303 -> 207 MMRT, 259 -> 164 Arrhenius), but under Arrhenius the +10 C community exhausts the substrate by about day 250 and then crashes (biomass 4-5 at day 365, E:C 0); under MMRT it keeps 54-57.
- Arrhenius would need its own calibration to reach about 80 % at ambient (for example a higher `Vmax_Km` or lower production).
- **Continuous 3-year run (no litter input):** the substrate falls to 20 % remaining at the end of year 1 (MMRT ambient), then to 10 % by day 392 and 5 % by day 411; the community then starves (biomass 0 during year 2, 330-365 starvation deaths per climate). There is no levelling off at 20 %; with a single substrate and no inputs, the remaining litter is used up within weeks of the second spring.

## 16. MMRT vs Arrhenius across a constant-temperature gradient (2026-10-16)

**Setup (`kinetics_gradient/`; `analysis/kinetics_gradient.py`, `analysis/plot_kinetics_gradient.py`):** bacterial minimal model with the chosen parameters (section 15), flat climate at 0-40 C in 5 C steps, Psi -0.5 MPa (no moisture effects), first year, 20x20, 2 seeds, R_endo on, `Kinetics_Model` 1 (MMRT) and 0 (Arrhenius). The kinetics switch changes only the enzyme and transporter rate constants (degradation and uptake); Km follows Vmax in both (Vmax_Km 13, Uptake_Vmax_Km 0.25).

| T (C) | Enzyme Vmax MMRT / Arrh | Transporter Vmax MMRT / Arrh | Degraded day 365 MMRT / Arrh | Day 50 % MMRT / Arrh |
|---|---|---|---|---|
| 0 | 1.1 / 2.5 | 3.8 / 11.1 | 2.7 / 23.0 % | - / - |
| 5 | 2.3 / 3.4 | 7.7 / 14.6 | 12.6 / 55.3 % | - / 353 |
| 10 | 3.9 / 4.4 | 13.0 / 19.1 | 67.5 / 95.2 % | 328 / 269 |
| 15 | 5.7 / 5.7 | 18.9 / 24.7 | 99.4 / 99.7 % | 231 / 214 |
| 20 | 7.1 / 7.4 | 23.4 / 31.7 | 99.7 / 99.8 % | 193 / 176 |
| 25 | 7.6 / 9.5 | 25.2 / 40.3 | 99.6 / 99.6 % | 184 / 148 |
| 30 | 7.2 / 12.0 | 23.8 / 50.9 | 99.4 / 99.5 % | 197 / 128 |
| 35 | 6.0 / 15.1 | 19.7 / 63.8 | 98.1 / 99.5 % | 240 / 114 |
| 40 | 4.4 / 18.9 | 14.6 / 79.3 | 53.8 / 99.4 % | 355 / 102 |

- **The enzyme rate constants agree at 10-20 C** (the reference range of the two parameter sets) and diverge outside it: MMRT is lower in the cold (1.1 vs 2.5 at 0 C), peaks at about 25 C and declines above it (4.4 at 40 C), while Arrhenius rises monotonically (18.9 at 40 C). The transporter rate constant is higher under Arrhenius at every temperature (about 1.3-1.4 x at 10-20 C, 3 x at 0 C, 5 x at 40 C).
- **Decomposition follows the rate constants with amplification by growth**: under MMRT the community decomposes the litter fastest at 25 C (50 % at day 184) and slows on both sides, failing to reach 50 % within the year at <= 5 C and at 40 C (54 % degraded); under Arrhenius time to 50 % falls monotonically (353 -> 102 days from 5 to 40 C).
- **Cold:** at 0-5 C the MMRT community barely grows (maximum biomass 1.7-8.5 vs 19-48, CUE 0.18-0.34, min E:C 37 at 0 C), so differences of about 2x in the rate constants become 4-8x differences in decomposition.
- **Warm:** above 15-20 C the substrate is exhausted within the year under both (except MMRT at 40 C) and the community then starves; under Arrhenius at 30-40 C almost all cells die after the substrate runs out (13-385 starvation deaths; biomass at day 365 near 0), under MMRT biomass at day 365 stays higher because decomposition is slower.
- Annual CUE declines with temperature above 10 C under both (0.43 -> 0.29-0.32 at 40 C), more steeply under Arrhenius.

## 17. Arrhenius Vmax matched to MMRT at 15-20 C (2026-10-16)

**Change:** in the minimal km25 inputs (MMRT parameters unchanged) `Vmax0_min`/`max` 55,000,000 -> 50,370,000 (x0.916) and `Uptake_Vmax0_min`/`max` 55,000,000 -> 41,745,000 (x0.759). Both are used only by the Arrhenius branch (and set the Arrhenius Km through the Km slopes). The factors are the inverse geometric-mean Arrhenius/MMRT ratios at 15-20 C (1 C steps) of the enzyme Vmax on the main substrate (IdealSub; new diagnostic `d_deg_Vmax_main`, plus `d_deg_Km_main`) and of the transporter Vmax.

- **The first matching attempt used the wrong measure.** Section 16's "enzyme Vmax" is the mean over all three substrates; Arrhenius gives each substrate its own Ea (IdealSub 36,000, DeadMic 37,000, DeadEnz 35,000 J/mol), so the mean matched MMRT (ratio 1.01-1.05 at 15-20 C) while the rate on the litter substrate was 8-12 % higher. Scaling by the mean (x0.982) and the transporter alone left Arrhenius 13-15 days ahead at 15-20 C; the transporter scaling alone had no effect (uptake is not limiting at Uptake_Vmax_Km 0.25).
- **Matched on IdealSub, the curves overlap at 12.5-20 C** (`kinetics_gradient_matched/`, 2 seeds; day of 25 / 50 / 75 % degraded):

| T (C) | MMRT | Arrhenius, original | Arrhenius, matched |
|---|---|---|---|
| 10 | 252 / 328 / - | 206 / 269 / 315 | 227 / 296 / 346 |
| 12.5 | 205 / 268 / 313 | - | 200 / 261 / 306 |
| 15 | 176 / 231 / 272 | 162 / 214 / 254 | 176 / 232 / 274 |
| 17.5 | 156 / 206 / 245 | - | 157 / 208 / 247 |
| 20 | 145 / 193 / 231 | 131 / 176 / 213 | 141 / 188 / 226 |
| 22.5 | 139 / 186 / 223 | - | 128 / 172 / 208 |
| 25 | 138 / 184 / 221 | 109 / 148 / 183 | 116 / 157 / 192 |

- Outside the match range the two still diverge as expected from their temperature functions: Arrhenius is faster in the cold (0-10 C) and keeps accelerating above 20 C, while MMRT slows above about 25 C (day 50 % at 40 C: 107 vs 355 days).

## 18. Rounded Arrhenius Vmax (5e7 / 4.2e7): gradient and climate runs (2026-10-16)

**Change:** minimal km25 inputs `Vmax0` 5.037e7 -> 5.0e7, `Uptake_Vmax0` 4.1745e7 -> 4.2e7 (commit 2eff07f). The MMRT runs are unaffected (these parameters are Arrhenius-only; same inputs and seeds) and are reused from sections 15-17. Outputs: `kinetics_gradient_rounded/`, `multiyear_arrh_matched/`.

- **Gradient (constant T, Psi -0.5):** day of 50 % degraded MMRT / Arrhenius: 10 C 328 / 298.5; 12.5 C 267.5 / 262.5; 15 C 230.5 / 234.5; 17.5 C 206 / 209; 20 C 193 / 189.5; 22.5 C 185.5 / 172; 25 C 184 / 158; 40 C 354.5 / 107. Rounding moved the Arrhenius timing by 0.5-2.5 days.
- **Six yearly pulses, real climate (2 seeds):**

| Climate | Degraded day 365, MMRT / Arrhenius | Day 50 %, MMRT / Arrhenius | Biomass C at day 365, MMRT / Arrhenius |
|---|---|---|---|
| Ambient | 80.1-81.8 / 85.5-86.5 % | 300-303 / 288-291 | 61-62 / 64-65 |
| +5 C | 96.9-97.1 / 98.5-98.6 % | 232-234 / 220-222 | 65-66 / 66-68 |
| +10 C | 98.0-98.2 / 98.4-98.5 % | 205-208 / 176-178 | 54-57 / 14-16 |

  (Before matching, Arrhenius at ambient was 94-95 % with day 50 % at 258-261.) The pulses replicate under both kinetics. The remaining gap comes from the parts of the year outside 12.5-20 C: Arrhenius is faster both in the cold season and in the warm season above about 20 C, so it stays 12 days ahead at ambient, 11 at +5 C and 30 at +10 C, and at +10 C its community again exhausts the substrate and crashes (biomass 14-16 vs 54-57 at day 365).
- **Continuous 3-year run (Arrhenius):** 20 % remaining at day 349 (ambient) / 290 (+5 C) / 237 (+10 C) and 5 % at day 396 / 335 / 299; the community starves in year 2, as under MMRT.

## 19. Is energy limiting with the chosen parameters? (2026-10-17)

From the runs of sections 15-18 (six yearly pulses on the climate data, MMRT and matched Arrhenius; constant-temperature gradient), split into the time while substrate remains (> 5 % left) and after it is exhausted. Energy shortfall = a day on which endogenous respiration (R_endo) had to pay for an energy deficit; production limitation = a day on which constitutive or inducible production (enzymes, osmolytes, HSP) was capped by the C/N/E limit.

- **While substrate remains, energy is essentially never limiting.** No shortfalls and no capped production in any climate run (except about 2 days just before exhaustion at +10 C under Arrhenius) or at 5-40 C; the lowest E:C is 165-540 (E_max ratio 1,400). The exception is MMRT at 0 C, where rates are so slow that the community runs short on 30 days and production is energy-capped on 211 days (min E:C 37).
- **Energy shapes allocation all the time:** 53-60 % of the C taken up is respired to make energy (the demand curve keeps E:C well below E_max), more where rates are slow (MMRT: 78 % at 0 C, 62-67 % at 35-40 C) and less where they are fast (Arrhenius: 45 % at 40 C). This is a cost on CUE rather than a hard limit.
- **After the substrate is exhausted, energy becomes limiting** and the community lives on R_endo until it starves: real climate, days per year with shortfalls MMRT / Arrhenius: ambient 0 / 0 (substrate never exhausted), +5 C 0 / 3, +10 C 19 / 66 (R_endo C 4 / 41 per grid cell); gradient: from 17.5 C upward (20-140 days). No energy was left unpaid in any run.

## 20. Making ambient nearly energy-limited and +5 C energy-limited (2026-10-17)

**Setup (`energy_limitation/`; `analysis/param_sweep.py`, `analysis/energy_limitation_summary.py`, `analysis/plot_energy_limitation.py`):** bacterial minimal model, MMRT, first year, 2 seeds, R_endo on, real climate (ambient and +5 C; +10 C for the last set). Energy limitation is scored only while substrate remains (> 5 %): days with a shortfall paid by R_endo, days with production capped by the C/N/E limit, minimum E:C. Sets: (a) `Enz_Prod_E_cost` 600-9,600 x production x1 / x1.5 / x2 (relative to Constit_Prod 0.00035 / Enz_Prod 0.007); (b) `Maint_E_Rate_Const` 2.5-15; (c) production x1.25 / x1.5 / x2 x maintenance 5-8; (d) production x1.25 with `Vmax_Km` 10-13 and maintenance 7-8.

- **Enzyme energy cost does not produce energy limitation.** Even at 16 x the cost (9,600) there are no shortfalls and no capped production; min E:C stays above 140. Cells cover any cost that scales with their activity by routing more of their uptake C to energy, and enzyme spending pays for itself through decay; the cost only slows growth (ambient 80 -> 66.5 % at x1 production).
- **Maintenance does, but with a threshold**: costs that scale with biomass rather than activity can outrun the supply. At current production, maintenance 5 gives min E:C 295 (ambient) / 205 (+5 C) and no limitation; 6 gives ambient barely limited (2.5 shortfall days, 80 capped days, min E:C 80) and +5 C limited (55 days), but only 39 / 64 % degraded; 7.5 and above collapse both.
- **Higher production with higher maintenance gives the requested contrast at about 80 % degraded:**

| Production x1.25 (Constit 0.0004375, Enz_Prod 0.00875) | Ambient: degraded / min E:C / limited days | +5 C: degraded / limited days | +10 C: degraded / limited days |
|---|---|---|---|
| Vmax_Km 13, maintenance 7 | 78 % / 221 / 0 | 81 % / 111 | 77 % / 134 |
| Vmax_Km 12, maintenance 7.5 | 79 % / 192 / 0 | 79 % / 125 | - |
| Vmax_Km 11, maintenance 7.75 | 82 % / 184 / 0 | 82 % / 125 | - |
| Vmax_Km 10, maintenance 7.5 | 91 % / 147 / 1.5 | 88 % / 125 | 87 % / 125 |
| Vmax_Km 13, maintenance 8 | 65 % / 114 / 31.5 | 67 % / 134 | 65 % / 104 |

- **"Energy-limited" at +5 C is a collapse, not a slowdown.** The warmer community grows faster and larger early in the year; late in the dry season its maintenance exceeds the supply, its energy pool runs out (from about day 230-250 at +5 C, 210 at +10 C), it shrinks through R_endo to about 1-2 C per grid cell and does not recover, and decomposition stalls with 20-30 % of the substrate left. Ambient stays above the threshold and keeps decomposing through the autumn, so all three climates end the year at about the same % degraded (the warming response of annual decomposition disappears).
- The model has no graded energy-limited state: a community either covers its maintenance or enters the energy trap (sections 8, 11). A persistent low-activity energy-limited state would need a mechanism such as maintenance that falls when energy is scarce (dormancy).

## 21. Lower enzyme Km slope with higher maintenance (2026-10-17)

**Setup:** as section 20 (production x1.25), `Vmax_Km` 10 / 9 / 8 x `Maint_E_Rate_Const` 7.5-9, three climates, 2 seeds (`energy_limitation/runs_ecost5.tar.gz`, `lower_km_x_maintenance_x1.25.txt`). Limited days counted while substrate remains (> 5 %).

| Vmax_Km, maintenance | Ambient: degraded / limited days (first day) / min E:C | +5 C: degraded / limited days (from day) / biomass day 365 | +10 C: degraded / limited days |
|---|---|---|---|
| 12, 7.5 (section 20) | 79 % / 0 / 192 | 79 % / 125 (239) / 1.1 | - |
| 10, 8 | 86 % / 0 / 185 | 85 % / 132 (234) / 1.6 | 84 % / 112 |
| 10, 8.5 | 80 % / 4.5 (267) / 131 | 80 % / 135 (231) / 1.9 | 79 % / 106 |
| 9, 8.5 | 87 % / 0 / 159 | 86 % / 137 (230) / 1.5 | 85 % / 118 |
| 9, 9 | 81 % / 23 (260) / 98 | 82 % / 142 (224) / 1.9 | 80 % / 146 |
| 8, 9 | 89 % / 8.5 (309) / 130 | 88 % / 141 (225) / 1.7 | 87 % / 130 |

- A lower Km slope raises decay per unit biomass, so a higher maintenance is needed for the same degradation; near 80 % this brings ambient to the edge of energy limitation in the late dry season (Vmax_Km 10 / maintenance 8.5: one seed limited on days 267-281, min E:C 131; Vmax_Km 9 / maintenance 9: 23 days, min E:C 98).
- **It does not prevent the warm-climate collapse**: at +5 and +10 C the community still runs out of energy from about day 225-240 (+5 C) and 200-210 (+10 C), shrinks to 1-2 C per grid cell and does not recover; annual degradation is the same in all three climates. In runs that degrade more than 90 % at ambient, the ambient limited days come at the end of the year, when the substrate is nearly gone (days 309-361).

## 22. Energy limitation when near-complete degradation is allowed (2026-10-17)

**Setup (`energy_limitation/runs_ecost6.tar.gz`, `production_x_km_x_maintenance_wide.txt`; `analysis/energy_limitation_grid_summary.py`):** production x1.25 / x1.5 / x2 x `Vmax_Km` 4 / 6 / 8 x `Maint_E_Rate_Const` 8 / 10 / 12, ambient and +5 C (and +10 C for the chosen option), 2 seeds, MMRT, R_endo on. Degradation at ambient may now approach 100 %. Limited days counted while substrate remains (> 5 %).

- **With high production the substrate is used up before the dry season** (all x2 runs, x1.5 at +5 C): no energy limitation while substrate remains, only starvation afterwards. Warming then only advances the exhaustion date.
- **x1.5 / Vmax_Km 6-8 / maintenance 8-10** make ambient barely limited in the dry season (1.5-29 days, min E:C 1-122) but +5 C exhausts the substrate first and is not limited before exhaustion: the opposite of the target.
- **Production x1.25, Vmax_Km 6, maintenance 8** gives the requested pattern:

| Climate | Degraded day 365 | Day 95 % degraded | Limited days before / after 95 % (first day) | Min E:C before 95 % | Max biomass C |
|---|---|---|---|---|---|
| Ambient | 99.0 % | 316-319 | 12-17 / 47-50 (302-304) | 57-60 | 44-45 |
| +5 C | 97.9-98.1 % | 239-241 | 18-23 / 93-101 (218-221) | 1-2 | 53 |
| +10 C | 98.1-98.2 % | 212-213 | 15-16 / 101-105 (196-198) | 3 | 54-55 |

  At ambient the energy pool stays above about 170 through the dry season and the community is limited only as the last of the substrate runs out (from about day 303). At +5 and +10 C the larger, warmer community runs out of energy in the dry season (from about day 220 / 197) while 10-20 % of the substrate is left; the remaining substrate is still used up within about 20 days, after which the community starves. Warming advances the 95 % date by 78 / 105 days.

## 23. Rounded production values with Vmax_Km 6, maintenance 8 (2026-10-17)

`Constit_Prod` 0.0004 with `Enz_Prod` 0.007 / 0.008 / 0.009 (vs 0.0004375 / 0.00875 in section 22); three climates, 2 seeds (`energy_limitation/runs_ecost7.tar.gz`). Limited days while substrate > 5 % (first day); after = once below 5 %.

| Constit / Enz_Prod | Ambient: degraded, day 95 %, limited (first), min E:C | +5 C: degraded, day 95 %, limited before / after (first) | +10 C: degraded, day 95 %, limited before / after |
|---|---|---|---|
| 0.0004 / 0.007 | 89 %, -, 0, 207 | 92 %, -, 123 / 0 (243) | 90 %, -, 154 / 0 |
| 0.0004 / 0.008 | 95 %, -, 11.5 (354), 109 | 94 %, -, 130 / 0 (236) | 94 %, -, 120 / 0 |
| 0.0004 / 0.009 | 98 %, 349, 6 (343), 96 | 96 %, 299, 70 / 67 (230) | 96 %, 253, 50 / 71 |
| 0.0004375 / 0.00875 | 99 %, 318, 14.5 (303), 59 | 98 %, 240, 21 / 97 (220) | 98 %, 213, 16 / 103 |

- Rounding `Constit_Prod` down to 0.0004 (-9 %) slows the community more than rounding `Enz_Prod` up to 0.009 (+3 %) speeds it, because constitutive production supplies most of the enzyme. With 0.0004 / 0.009 the pattern is kept: ambient decomposes through the dry season (min E:C about 220-250 then) and is limited only near the end of the litter (from about day 343); +5 and +10 C run out of energy in the dry season (from day 230 / 203) and collapse with more substrate left than in section 22 (about 10 % at the crash vs 5-10 %), reaching 96 % by day 365.
- With `Enz_Prod` 0.008 or 0.007 the warm communities collapse earlier relative to their substrate and end at 92-94 %, with ambient similar or slightly ahead.

## 24. New minimal-model defaults: Constit_Prod 0.00045, Enz_Prod 0.008, Vmax_Km 6, Maint_E_Rate_Const 8 (2026-10-17)

Set in the 8 minimal km25 input folders (commit b14b846; Uptake_Vmax_Km 0.25, slopes 0.1, beta 1, Arrhenius Vmax0 5e7 / Uptake_Vmax0 4.2e7 unchanged). Check run: MMRT, real climate, first year, 2 seeds (`energy_limitation/runs_defaults_check.tar.gz`, figure `energy_limitation_defaults_c0.00045_e0.008_km6_maint8.png`).

| Climate | Degraded day 365 | Day 50 % / 95 % | Energy-limited days while > 5 % left (first day; % left then) | Min E:C while > 5 % left | Limited days after |
|---|---|---|---|---|---|
| Ambient | 99.2 % | 223-224 / 315-318 | 18-23 (291-299; 12-15 %) | 68-85 | 48-51 |
| +5 C | 98.0 % | 165 / 242 | 24 (218; 12 %) | 2.3-2.4 | 90-102 |
| +10 C | 98.2 % | 147 / 214 | 15-17 (197-199; 9-10 %) | 2.3-2.6 | 99-104 |

Practically the same as production x1.25 (section 22). Ambient becomes energy-limited only at the end of the dry season with 12-15 % of the litter left; the warm communities run out of energy in the dry season and collapse with 10-12 % left, after which the remaining litter is used up and the community starves. The Arrhenius match of section 17-18 is unaffected (it compares rate constants, which do not depend on production, maintenance or the Km slope).

## 25. HSP and osmolyte effects and production on/off, current defaults (2026-10-17)

**Setup (`hsp_osmo_test/`; `analysis/hsp_osmo_test.py`, `analysis/plot_hsp_osmo_test.py`):** current minimal defaults (section 24), MMRT, real climate, first year, 2 seeds, R_endo on. Effect off (production and its C/N/E costs kept): HSP `Maint_E_HSP_Ea_Slope` 0 (maintenance Ea stays 45,000); osmolytes `Osmo_tol_K` 1e9 (drought tolerance about 0). Production off: both production rates 0 (no cost, no effect).

| Variant | Degraded day 365 (amb / +5 / +10) | Energy-limited days while > 5 % left, ambient (first day) | Drought kills, ambient | Annual CUE, ambient |
|---|---|---|---|---|
| Baseline | 99.2 / 98.0 / 98.2 % | 20.5 (295) | 253 | 0.26 |
| HSP effect off | 58.0 / 69.2 / 53.1 % | 143 (223) | 281 | 0.18 |
| HSP production off | 67.1 / 74.6 / 60.4 % | 131.5 (233) | 276 | 0.19 |
| Osmolyte effect off | 99.2 / 98.0 / 98.1 % | 7 (312) | 357 | 0.28 |
| Osmolyte production off | 99.5 / 98.9 / 98.9 % | 13 (290) | 361 | 0.30 |
| Both effects off | 57.4 / 69.0 / 53.1 % | 142.5 (220) | 380 | 0.19 |
| Both production off | 95.5 / 83.1 / 68.7 % | 104 (256) | 377 | 0.22 |

- **HSP is what keeps the community viable at Maint_E_Rate_Const 8.** HSP holds the maintenance Ea at about 30,000-34,000 instead of 45,000 (about 30 % less maintenance at 20 C); without its effect the community runs out of energy in the dry season at every temperature (from day 220 at ambient, 135-140 at +10 C) and decomposition stalls at 53-69 %. Producing HSP without the effect costs a little more than not producing it (58 vs 67 % at ambient), so the benefit outweighs the cost by far.
- **Feedback in the warm runs:** once energy runs short in the dry season, HSP production is capped, the HSP share falls (turnover 0.1/day) and the maintenance Ea climbs back towards 40,000-44,000, which raises maintenance further and accelerates the collapse.
- **Osmolytes do not pay for themselves with beta 1.** Their protection lowers drought kills (253 vs 357 at ambient) but hardly changes decomposition; their production costs energy in the dry season, and removing it relieves the warm runs (min E:C at +5 / +10 C 40 / 29 instead of 2; energy-limited days 9 / 4.5 instead of 24 / 16) and speeds decomposition slightly (day 95 % at ambient 304 vs 317).
- Without both, the community does better than without HSP alone at warm temperatures, because it no longer pays the osmolyte energy cost.

## 26. Osmolyte sweeps: beta, Osmo_prod_E_cost, Osmo_Induci_Prod (2026-10-17)

**Setup (`osmolyte_sweeps/`; `analysis/param_sweep.py`, `analysis/plot_osmo_sweeps.py`):** current minimal defaults, MMRT, real climate (three scenarios), first year, 2 seeds, R_endo on. `beta_bac` = `beta_fun` 1 (default), 2, 3, 5, 10, each with osmolytes on and with osmolyte production off; `Osmo_prod_E_cost` 600 (default), 300, 150, 60; `Osmo_Induci_Prod` 0.05 (default), 0.025, 0.01, 0.005. Limited days and min E:C while > 5 % of the substrate is left.

| Setting | Day 95 % degraded (amb / +5 / +10) | Limited days (amb / +5 / +10) | Min E:C (amb / +5 / +10) | Drought kills, ambient |
|---|---|---|---|---|
| beta 1, osmolytes on (default) | 317 / 242 / 214 | 20.5 / 24 / 16 | 77 / 2 / 2 | 253 |
| beta 1, osmolyte production off | 304 / 231 / 207 | 13 / 9 / 4.5 | 101 / 40 / 29 | 361 |
| beta 3, on / off | 326 / 246 / 215; 318 / 235 / 210 | 5.5 / 24 / 15; 0 / 11.5 / 7.5 | 115 / 3 / 5; 191 / 48 / 27 | 530; 801 |
| beta 10, on / off | 361 / 259 / 223; - (92 %) / 258 / 219 | 0.5 / 26.5 / 19.5; 0 / 3.5 / 9.5 | 162 / 31 / 19; 293 / 139 / 93 | 1,424; 2,291 |
| Osmo_prod_E_cost 300 | 308 / 236 / 211 | 16.5 / 15 / 11.5 | 73 / 5 / 11 | 259 |
| Osmo_prod_E_cost 150 | 304 / 234 / 208 | 12.5 / 13 / 6.5 | 69 / 11 / 16 | 278 |
| Osmo_prod_E_cost 60 | 302 / 232 / 208 | 12.5 / 14 / 9 | 82 / 22 / 28 | 264 |
| Osmo_Induci_Prod 0.005-0.025 | 317 / 242 / 214 | 14.5-17 / 22 / 14-15.5 | 75-80 / 2-3 / 2 | 261-278 |

- **A steeper beta raises drought kills about linearly (osmolytes cut them by about 35-40 % at every beta) but changes decomposition little up to beta 5** (day 95 % at ambient 317 -> 337). Only at beta 10 do osmolytes pay for themselves at ambient (96 vs 92 % degraded at day 365); in the warm runs osmolyte production still costs more than it saves (production off decomposes slightly faster at every beta). Higher beta also thins the community, which relieves the ambient energy limitation (limited days 20.5 -> 0.5 at beta 10) while the warm runs stay limited.
- **`Osmo_Induci_Prod` has almost no effect**: the osmolyte pool (2-2.5 % of biomass) comes from constitutive production (0.005 x biomass C per day, against inducible 0.05 x biomass-C uptake x the drought factor, about ten times smaller); cost and protection are set by `Osmo_Consti_Prod`.
- **Lowering `Osmo_prod_E_cost` keeps the protection and removes most of the energy penalty**: at 60-150 the warm runs approach osmolyte-production-off (min E:C 11-28 vs 2; limited days 6.5-14 vs 16-24; day 95 % 3-10 days earlier) with drought kills unchanged (about 260-280). This also weakens the +5 C energy limitation that the current defaults were tuned for (section 24).

## 27. Osmo_prod_E_cost 150 across beta (2026-10-17)

**Decision:** `Osmo_prod_E_cost` lowered from 600 to 150 in the eight minimal km25 input folders, so osmolytes give a net benefit. The full-model km25 folders are unchanged (600).

**Setup (`osmolyte_sweeps/osmo_cost150_beta.*`; `analysis/plot_osmo_cost_beta.py`):** `beta_bac` = `beta_fun` 1, 2, 3, 5, 10 at cost 150, compared with the cost-600 and osmolyte-production-off runs from section 26 (same design: MMRT, real climate, first year, 2 seeds, R_endo on).

| beta | Day 95 % degraded, amb / +5 / +10: cost 150 (cost 600; off) | Limited days, amb / +5 / +10: cost 150 (cost 600; off) | Min E:C, +5: cost 150 (600; off) | Drought kills, ambient: cost 150 (600; off) |
|---|---|---|---|---|
| 1 | 304 / 234 / 208 (317 / 242 / 214; 304 / 231 / 207) | 12.5 / 13 / 6.5 (20.5 / 24 / 16; 13 / 9 / 4.5) | 11 (2; 40) | 278 (253; 361) |
| 2 | 307 / 235 / 209 (320 / 244 / 215; 310 / 234 / 208) | 5.5 / 12 / 8 (8.5 / 22.5 / 15.5; 3 / 11 / 5.5) | 14 (4; 29) | 377 (403; 579) |
| 3 | 312 / 236 / 210 (326 / 246 / 215; 318 / 235 / 210) | 9 / 14.5 / 9 (5.5 / 24 / 15; 0 / 11.5 / 7.5) | 19 (3; 48) | 518 (530; 801) |
| 5 | 320 / 240 / 212 (337 / 248 / 217; 334 / 241 / 211) | 0 / 16.5 / 12 (3.5 / 24 / 15; 0 / 11 / 8) | 21 (11; 66) | 767 (809; 1,221) |
| 10 | 345 / 246 / 214 (361 / 259 / 223; not reached, 92 % / 258 / 219) | 0 / 9 / 11 (0.5 / 26.5 / 19.5; 0 / 3.5 / 9.5) | 73 (31; 139) | 1,418 (1,424; 2,291) |

- **At cost 150, osmolytes are neutral at beta 1 and a net benefit from beta 2 upward.** At ambient, day 95 % is 3 / 6 / 14 days earlier than with osmolyte production off at beta 2 / 3 / 5, and at beta 10 degradation reaches 98.4 % vs 92.1 %. The warm runs are within 1-3 days of osmolyte production off at beta 1-5 and gain at beta 10 (+5: 246 vs 258; +10: 214 vs 219). Cost 150 is faster than cost 600 everywhere (12-17 days at ambient, 5-13 days in the warm runs).
- **Protection is unchanged by the cost:** drought kills at cost 150 match cost 600 (osmolytes cut kills by 23-38 % vs off, more at higher beta), and dry-season drought tolerance stays about 0.38 at ambient and 0.24-0.31 in the warm runs.
- **Energy limitation:** the warm runs remain energy-limited at every beta (+5: 9-16.5 days, min E:C 11-21 for beta 1-5; +10: 6.5-12 days), but much less than at cost 600 (min E:C about 2-11). Ambient is limited for 12.5 days at beta 1 (min E:C 69, about the same as with osmolyte production off) and not at all from beta 5.

## 28. HSP cold cost at the final set (2026-10-19)

**Decisions recorded the same day** (see `SENSITIVITY_PLAN.md`): beta 1; endogenous respiration on by default (`R_endo,1` in every input folder, code default 1); HSP should carry a modest cold cost; Morris design for the final sensitivity analysis; no fungi; CUE 0.26-0.35 accepted, higher CUE a secondary benchmark.

**Setup (`hsp_cold_cost/`; `analysis/kinetics_gradient.py`, `analysis/plot_hsp_cold_cost.py`):** final minimal set (Osmo_prod_E_cost 150), MMRT, constant temperature 0-40 C with wet Psi (-0.5 MPa), one 365-day pulse, 2 seeds, R_endo on. Each setting is compared with its control (HSP_Consti_Prod = HSP_Induci_Prod = 0). The control does not depend on HSP cost, production or the Ea slope; for the Tref variant, Maint_E_Rate_Const was rescaled (8 -> 16, i.e. 8 x exp(45000/R x (1/273 - 1/283)) = 16.1) so that the control maintenance curve is unchanged (its results match the original control).

Mean biomass C over the year, HSP vs control (%):

| Setting | 0 C | 5 C | 7.5 C | 10 C | 12.5 C | 15 C | 20 C | 25 C | Crossover |
|---|---|---|---|---|---|---|---|---|---|
| Current (HSP_prod_E_cost 300) | -18 | +1 | +9 | +12 | +10 | +9 | +16 | +44 | ~5 C |
| HSP_prod_E_cost 450 | -24 | -6 | +5 | +7 | +8 | +8 | +15 | +44 | ~6 C |
| HSP_prod_E_cost 600 | -28 | -10 | -2 | +2 | +7 | +7 | +15 | +43 | ~9 C |
| HSP_prod_E_cost 1200 | -40 | -28 | -23 | -14 | 0 | +4 | +12 | +40 | ~12.5 C |
| HSP_Consti_Prod 0.003 | -24 | 0 | +13 | +14 | +13 | +12 | +22 | +55 | ~5 C |
| Ea slope -60000 | -18 | -3 | +4 | +5 | +6 | +5 | +11 | +36 | ~6 C |
| Tref 283 K, Maint_E_Rate_Const 16 | -27 | -13 | -12 | -6 | 0 | +2 | +7 | +25 | ~12.5 C |

- **At the final set HSP costs biomass only below about 5 C**, so it is a net benefit in every climate scenario (ambient daily temperature 5.3-30.7 C, median 15.7). Controls collapse at 30 C and above in every setting (degradation 14 % vs 93-100 % with HSP).
- **Production-cost route** (matches the rationale that HSP diverts resources from enzymes and growth): HSP_prod_E_cost 600 (the original value, before the settled set lowered it to 300) moves the crossover to about 9 C with a modest winter cost (-10 % at 5 C, -2 % at 7.5 C) and leaves the warm benefit unchanged; 1200 gives an extreme cold penalty. More constitutive HSP strengthens the benefit without a cold cost; a shallower Ea slope shrinks the benefit everywhere.
- **Tref route** (Tref is where HSP's Ea reduction has no effect on maintenance; below it HSP raises maintenance): Tref 283 K with the rescaled rate constant gives a crossover at about 12.5 C, a cold cost of 6-13 % at 5-10 C and a smaller warm benefit (+7 % at 20 C, +25 % at 25 C). It is a thermal cost, not a resource cost, and it raises the maintenance of HSP-carrying cells by roughly 10-20 % at every temperature.
- Degradation by day 365 is nearly unaffected between 10 and 25 C (99-100 %); the HSP effect shows in biomass and timing (day 50 % 8-35 days earlier with HSP at 10-25 C in the current setting).
- Either change will shift the climate-scenario energy calibration (section 24-27), so it needs the climate re-baseline (plan step 1.1) before the Morris runs.

**28b. CUE and energy limitation for the HSP options (2026-10-19).** Annual CUE is uptake-weighted; limited days and minimum E:C are counted while > 5 % substrate is left.

*Constant-temperature gradient (wet, one pulse):*
- CUE: HSP raises CUE above the control from 7.5 C (current setting +0.01-0.03 at 7.5-20 C, +0.12 at 25 C, +0.17 at 30 C). HSP_prod_E_cost 600 lowers CUE vs the current setting by 0.024 at 5 C and by 0.005-0.01 at 7.5-15 C (0.32-0.34 at 7.5-15 C); 1200 by 0.07 at 5 C. The Tref 283 K option lowers CUE by 0.01-0.03 at 5-20 C and by 0.03 / 0.12 at 25 / 30 C.
- Energy limitation: no HSP setting is energy-limited during decomposition between 5 and 20 C; all are limited for most of the year at 0 C (kinetic limb; HSP runs 200-235 days vs 140 for the control). Controls become limited from 20 C (8.5 days; 165 at 25 C). Minimum E:C at 5-25 C: current 340-470; cost 600 285-440; cost 1200 140-375; Tref option 176 at 20 C, 28 at 25 C, 1.5 at 30 C. Endogenous respiration after the substrate is exhausted (about 60 C units per cell at 12.5-25 C) is the same in every setting.

*Real climate (first year, 2 seeds):*

| Setting | CUE (amb / +5 / +10) | Limited days | Min E:C | Day 95 % degraded | % degraded, day 365 |
|---|---|---|---|---|---|
| Current (HSP_prod_E_cost 300) | 0.287 / 0.303 / 0.296 | 12.5 / 13 / 6.5 | 69 / 11 / 16 | 304 / 234 / 208 | 99.4 / 98.6 / 98.7 |
| HSP_prod_E_cost 450 | 0.280 / 0.300 / 0.295 | 8.5 / 16.5 / 10.5 | 68 / 7 / 8 | 308 / 237 / 211 | 99.3 / 98.5 / 98.4 |
| HSP_prod_E_cost 600 | 0.269 / 0.297 / 0.292 | 12 / 20 / 19 | 64 / 3 / 2 | 314 / 241 / 217 | 99.2 / 98.2 / 98.1 |
| Tref 283 K, rate constant 16, cost 300 | 0.220 / 0.248 / 0.242 | 62 / 155 / 168 | 36 / 0 / 0 | 340 / - / - | 97.7 / 92.9 / 90.0 |
| Tref 283 K, rate constant 16, cost 600 | 0.212 / 0.239 / 0.230 | 63 / 156 / 170 | 29 / 0 / 0 | 348 / - / - | 97.2 / 90.8 / 87.4 |

- HSP_prod_E_cost 600 leaves ambient as before (12 limited days, min E:C 64) and deepens the warm-run energy limitation (+5: 20 days, min E:C 3; +10: 19 days, min E:C 2), restoring the +5 C tipping that the osmolyte cost cut (section 27) had weakened. CUE falls by 0.02 at ambient and by less than 0.01 in the warm runs.
- The Tref option breaks the calibration: CUE 0.21-0.25, ambient limited for about 60 days, warm runs limited for 155-170 days with degradation 87-93 %. It would need maintenance re-tuned.

## 29. N limitation check (2026-10-19)

**NH4 inputs are on.** `Input_NH4` 0.02 (N per grid cell per day) and `Input_PO4` 0.002 are added to the monomer pools every day before uptake (`grid.py` `uptake()`); NH4 and PO4 leach at `Leaching` 0.1 x exp(0.5 (Psi - wp_fc)) per day. `Init_NH4` is 0. The substrate in the `crich` folders is IdealSub C:N:P 200:5:0.5 (C:N 40).

**N limitation does occur.** Bacterial minimum N share is `Nfrac_b - Nrange` = 0.16 - 0.04 = 0.12 of C+N+P (biomass N:C about 0.138). When a cell falls below it, rebalancing (`mortality()`) cuts C to the minimum ratio and respires the excess, which counts against CUE. Days at the N floor while > 5 % substrate is left (current defaults, 2 seeds):
- Real climate: ambient 0; +5 C 45 (from about day 160); +10 C 67 (from about day 125).
- Constant temperature (wet, one pulse): 0 at 0-10 C and 35-40 C; 31-75 days at 12.5-30 C, starting at day 110-235.
- This is the step down in CUE (0.39 -> 0.26) at 12.5-30 C in `hsp_cold_cost/hsp_cue_timeseries.png`: N demand of peak growth exceeds the organic N arriving with C (C:N 40) plus 0.02 NH4, and the NH4 pool is drawn down to about the daily input. Part of the warm-scenario response in the current runs is therefore N limitation.

**More N removes the floor but exposes an uptake artefact** (`nh4/`: Input_NH4 0.05 / 0.1 / 0.2, and IdealSub 200:10:1 with Input_NH4 0 and 0.02):

| Case | Days at N floor (amb / +5 / +10) | Max biomass N:C | CUE | Energy-limited days | Min E:C | Day 95 % degraded |
|---|---|---|---|---|---|---|
| C:N 40, NH4 0.02 (current) | 0 / 45 / 67 | 0.29 / 0.22 / 0.20 | 0.287 / 0.303 / 0.296 | 12.5 / 13 / 6.5 | 69 / 11 / 16 | 304 / 234 / 208 |
| C:N 40, NH4 0.05 | 0 / 0 / 0 | 0.34 | 0.263 / 0.288 / 0.300 | 10 / 32.5 / 21 | 80 / 0.07 / 0.2 | 329 / 252 / 220 |
| C:N 40, NH4 0.1 | 0 / 0 / 0 | 0.35 | 0.190 / 0.254 / 0.271 | 278 / 146 / 123 | 0 / 0 / 0 | not reached / not reached / 332 |
| C:N 20, NH4 0 | 0 / 5 / 10.5 | 0.19 / 0.30 / 0.32 | 0.294 / 0.279 / 0.289 | 1.5 / 28 / 28 | 178 / 0.9 / 0.2 | 326 / 265 / 237 |
| C:N 20, NH4 0.02 | 0 / 0 / 0 | 0.33 | 0.267 / 0.293 / 0.302 | 14.5 / 29.5 / 20.5 | 75 / 0.1 / 0.04 | 322 / 251 / 223 |

- NH4 uptake is not regulated by N demand: cells take up NH4 until their C share reaches its minimum (Cfrac_b - Crange = 0.735, biomass N:C about 0.34), rebalancing then returns the excess N to the NH4 pool, and it is taken up again. Every unit of monomer mass taken up, NH4 included, is charged `Active_Uptake_E_cost` and the transporter-production energy cost, so this cycle drains energy. At 0.1-0.2 NH4 the communities collapse; on the gradient the cost falls hardest where growth is slow (10 C: 64 % degraded at NH4 0.1, 6 % at 0.2; corrected, see section 30).
- Luxury N uptake already happens with the current inputs during slow growth (biomass N:C 0.26-0.33 at ambient and at 10 C).
- Removing the floor by either route makes the warm runs more energy-limited (min E:C about 0) and 15-25 days slower to 95 % degraded, and removes the CUE step on the gradient (CUE at 15-25 C 0.33-0.35 with NH4 0.05, vs 0.30-0.33).
- Not fixed yet; options are in OPEN_ISSUES.md. The energy calibration (sections 24-28) will need re-checking after the fix.

## 30. Passive mineral uptake (2026-10-19)

**Decision and code change.** NH4 and PO4 enter cells by passive diffusion: their uptake is no longer charged the active-transport energy cost (`Active_Uptake_E_cost`) or the transporter-production energy cost (`Uptake_Prod_E_cost`, from which the transporter C cost output is derived). Optional switch `Mineral_Uptake_Passive` (parameters.csv, default 1; row added to every input folder); 0 reproduces the earlier behaviour exactly (checked: identical daily output for Input_NH4 0.02 and 0.05, three climates, 2 seeds). `grid.py`: `Taxon_Transport_Uptake` (biomass C + N + P taken up, without NH4 and PO4 when the switch is on) replaces `Taxon_Total_Uptake` in the two transport energy costs; `Taxon_Total_Uptake` is unchanged elsewhere (it still scales inducible HSP production).

**Correction to section 29.** A test folder for the C:N 20 runs was created as a link to the C-rich input folder in the analysis workspace, so that folder's `initial_substrates.csv` was overwritten with 200:10:1 from 11:53 to 12:20 on 2026-10-19. Affected and re-run: the gradient runs at Input_NH4 0.1 and 0.2 (10 C at 0.2: 6 % degraded, not 5 %; other conclusions unchanged). Unaffected: the climate NH4 runs, the gradient at NH4 0.05, the C:N 20 runs and all HSP runs. The repository input folders were never touched.

**Results** (final minimal set, osmolyte cost 150, HSP cost 300, R_endo on, 2 seeds; real climate, first year):

| Case | Days at N floor (amb / +5 / +10) | Max biomass N:C | CUE | Days with endogenous respiration | Min E:C | Day 95 % degraded | % degraded, day 365 |
|---|---|---|---|---|---|---|---|
| Charged, NH4 0.02 (previous) | 0 / 45 / 67 | 0.29 / 0.22 / 0.20 | 0.287 / 0.303 / 0.296 | 14 / 12 / 7 | 69 / 11 / 16 | 304 / 234 / 208 | 99.4 / 98.6 / 98.7 |
| Passive, NH4 0.02 | 0 / 48 / 69 | 0.29 / 0.22 / 0.20 | 0.289 / 0.300 / 0.292 | 14 / 12 / 7 | 74 / 26 / 16 | 303 / 230 / 207 | 99.4 / 98.8 / 98.7 |
| Passive, NH4 0.05 | 0 / 0 / 0 | 0.34 | 0.265 / 0.300 / 0.304 | | 24 / 2 / 2 | 309 / 239 / 214 | 98.9 / 98.0 / 98.0 |
| Passive, NH4 0.1 | 0 / 0 / 0 | 0.35 | 0.278 / 0.279 / 0.283 | | 59 / 1 / 2 | 306 / 244 / 215 | 99.3 / 97.4 / 97.7 |
| Passive, NH4 0.2 | 0 / 0 / 0 | 0.35 | 0.296 / 0.286 / 0.288 | 11 / 16 / | 74 / 10 / 11 | 299 / 234 / 206 | 99.4 / 98.5 / 98.6 |
| Charged, NH4 0.1 / 0.2 (for comparison) | 0 | 0.35 | 0.16-0.27 | | 0 | not reached in most | 8-95 |

- Passive uptake removes the energy drain of the N hoarding cycle: no collapse at any NH4 input (degradation 97-99 %, biomass and timing close to the current runs). Cells still take up NH4 to the C-share minimum (biomass N:C about 0.35, C:N about 2.9) when N is plentiful; this no longer costs energy.
- With Input_NH4 0.02 the warm runs are still N-limited (48 / 69 days at the floor); 0.05 or more removes the floor in every climate and on the gradient.
- With more NH4 the flagged energy-limited days rise (ambient 14.5 -> 21 / 54 / 127 at 0.05 / 0.1 / 0.2) through inducible production running below potential (`r_induci` < 1) while E:C stays high, not through energy shortfalls (days with endogenous respiration unchanged). From here on the limited-days metric should count only endogenous-respiration days, or be split by limiting resource.
- Constant-temperature gradient (wet, one pulse): passive uptake also changes the cold limb. At 0 C endogenous respiration falls from 176 days to 0 and degradation rises from 2.0 to 3.9 %; at 5 C CUE is 0.31 vs 0.27 and degradation 22.8 vs 18.3 % (the NH4 transport cost was a large share of the cold energy budget). With NH4 0.05: no N floor at any temperature, CUE 0.33 / 0.35 / 0.34 / 0.33 at 10 / 15 / 20 / 25 C (the step down at 12.5-30 C is gone), day 50 % within 3 days of the current runs.
- The HSP cold-cost screen (section 28) and the energy calibration (sections 24-27) were made with charged mineral uptake and N limitation in the warm runs; both need re-running once Input_NH4 is chosen.

## 31. Energy credit for C respired in stoichiometric rebalancing, R_stoich (2026-10-19)

**Code change.** When a cell falls below its minimum N or P share, rebalancing (`mortality()`) cuts its C to the minimum ratio and respires the excess. With `R_stoich` on, that C yields energy at the necromass respiration yield (`Mon_DeadMic` E_resp, 700 per unit C, the same yield as R_endo), credited to the cell's energy pool just before the E:C cap (energy above E_max x C is dissipated there). The C is still respired and still counts against CUE. Optional switch `R_stoich` (parameters.csv; default 1; row added to every input folder), independent of `R_endo`. New outputs: `Taxon_R_stoich_E` (grid), `R_stoich_E_Taxon` (output.py, per taxon per day), `R_stoich_E` (sweep series), `E_stoich_g` (e_tracker, which now separates the credit from `E_cap_g` and includes it in `resid_E`).

**Checks** (final minimal set, passive mineral uptake, R_endo on, 2 seeds, real climate first year):
- `R_stoich` 0 reproduces the previous code exactly (identical daily output, Input_NH4 0.02 and 0.05, three climates).
- Energy balance with `R_stoich` 1 (e_tracker, +10 C, seed 0): daily residual at most 2e-4 (6e-7 of the day's uptake energy); credit = 700 x rebalancing C on every day it occurs; no credit dissipated at the cap (E:C 500-650 during N limitation, cap 1400).

**Effect under N limitation** (Input_NH4 0.02; ambient / +5 / +10 C; ambient has no rebalancing, so it is unchanged):

| R_stoich | Days at N floor | Energy credited (per cell, year; share of uptake energy) | Days with endogenous respiration | Min E:C | CUE | Day 95 % degraded | % degraded, day 365 |
|---|---|---|---|---|---|---|---|
| 0 | 0 / 48 / 69 | 0 | 14.5 / 12 / 7 | 74 / 26 / 16 | 0.289 / 0.300 / 0.292 | 303 / 230 / 207 | 99.4 / 98.8 / 98.7 |
| 1 | 0 / 52 / 72 | - / 7,260 (4.0 %) / 12,300 (7.1 %) | 14.5 / 7 / 3 | 74 / 51 / 83 | 0.289 / 0.304 / 0.297 | 303 / 229 / 205 | 99.4 / 99.0 / 99.0 |

- The credit arrives on 90-100 days in the warm runs (peak about 240-310 energy units per cell per day) and halves the endogenous-respiration days; the minimum E:C rises from 26 / 16 to 51 / 83. Cells route less uptake C to energy (uptake energy -2 to -5 %), so CUE rises slightly (+0.004-0.005).
- With Input_NH4 0.05 (no N floor) rebalancing is rare: the credit is 0.01 % of uptake energy and the outputs are unchanged to three digits.

## 32. Re-run with Input_NH4 0.05, passive mineral uptake and R_stoich (2026-10-19)

**Change:** `Input_NH4` 0.02 -> 0.05 in the eight minimal km25 input folders (full-model folders unchanged at 0.02 until the port). All runs: final minimal set (osmolyte cost 150, beta 1, maintenance 8), `Mineral_Uptake_Passive` 1, `R_stoich` 1, `R_endo` 1, MMRT, 2 seeds (`hsp_cold_cost_nh4_0.05/`). Energy limitation is counted as days with endogenous respiration while > 5 % substrate is left (section 30).

**HSP cold cost, constant-temperature gradient** (wet, one pulse). Mean biomass, HSP vs control (%):

| HSP_prod_E_cost | 0 C | 5 C | 7.5 C | 10 C | 12.5 C | 15 C | 20 C | 25 C | Crossover |
|---|---|---|---|---|---|---|---|---|---|
| 300 (current) | -6 | -2 | +5 | +11 | +11 | +16 | +29 | +60 | ~6 C |
| 450 | -11 | -7 | -0.5 | +7 | +9 | +14 | +28 | +56 | ~7.5 C |
| 600 | -15 | -14 | -7 | +2 | +7 | +11 | +25 | +54 | ~10 C |

- CUE: at 5 C 0.308 / 0.297 / 0.284 for cost 300 / 450 / 600 (control 0.309); at 10-20 C 0.33-0.35 (300) and 0.32-0.34 (600), control 0.25-0.31; at 0 C HSP lowers CUE (0.20 / 0.19 / 0.18 vs 0.24). The N-limitation step is gone; a small dip remains at 17.5-25 C with cost 300 (P floor reached for 4-31 days at peak growth; PO4 drawn down to the daily input; R_stoich credits energy then).
- No HSP run is energy-limited below 30 C; the control is limited from 20 C (15 days; 173 at 25 C) and collapses at 30 C and above.

**Real climate, HSP cost** (ambient / +5 / +10 C):

| HSP_prod_E_cost | CUE | Endogenous-respiration days | Min E:C | Day 95 % degraded | % degraded, day 365 |
|---|---|---|---|---|---|
| 300 | 0.265 / 0.300 / 0.304 | 21 / 22.5 / 17 | 24 / 1.9 / 1.9 | 309 / 239 / 214 | 98.9 / 98.0 / 98.0 |
| 450 | 0.252 / 0.290 / 0.294 | 34 / 28.5 / 25 | 40 / 1.3 / 0.5 | 316 / 247 / 221 | 98.9 / 97.5 / 97.5 |
| 600 | 0.249 / 0.280 / 0.284 | 21 / 38 / 32 | 40 / 0.05 / 0.2 | 321 / 257 / 231 | 98.8 / 96.9 / 96.9 |

**HSP/osmolyte on-off at cost 300** (ambient / +5 / +10 C): baseline 98.9 / 98.0 / 98.0 % degraded; HSP effect off 79 / 68 / 52 %; HSP production off 89 / 75 / 59 %; osmolyte production off 99.2 / 98.4 / 98.4 % (day 95 % 2-4 days earlier, CUE 0.286 / 0.309 / 0.310, drought kills 361 / 367 / 349 vs 265 / 305 / 295); osmolyte effect off within 0.2 points of baseline. Endogenous-respiration days 21 / 22.5 / 17 (baseline) vs 105 / 168 / 156 with the HSP effect off.

- Compared with the earlier calibration (NH4 0.02, charged mineral uptake; sections 26-28): ambient is closer to the energy limit (21 vs about 14 endogenous-respiration days, min E:C 24 vs 69-77) and +5 C is clearly tipped (min E:C 2 vs 11); ambient CUE 0.265 (was 0.287), warm-run CUE unchanged (0.30). The calibration targets still hold: ambient close to energy-limited (late in the season), +5 C energy-limited, about 98-99 % degraded.
- HSP stays essential at every climate; osmolytes still cost slightly more than they save at beta 1 (as in section 27).
- P is close to limiting in the warm runs (biomass P:C minimum 0.013 vs a floor of about 0.012) and limits briefly on the wet constant-temperature gradient. `Input_PO4` 0.002 gives an N:P input ratio of 25 (substrate 10).

## 33. Defaults set: HSP_prod_E_cost 600, Input_PO4 0.005 (2026-10-19)

**Decision:** `HSP_prod_E_cost` 300 -> 600 (modest HSP cold cost: crossover about 10 C, -7 to -14 % biomass at 5-7.5 C; section 32) and `Input_PO4` 0.002 -> 0.005 (N:P input 10, as in the substrate) in the eight minimal km25 input folders. With `Input_NH4` 0.05, passive mineral uptake, `R_stoich` 1 and `R_endo` 1, this is the current minimal baseline.

**Check** (real climate first year, 2 seeds; ambient / +5 / +10 C; the PO4 0.002 column is the same run with the old P input):

| Input_PO4 | Min biomass P:C | Rebalancing energy credit | CUE | Endogenous-respiration days | Min E:C | Day 95 % degraded | % degraded, day 365 |
|---|---|---|---|---|---|---|---|
| 0.002 | 0.017 / 0.014 / 0.014 | 0 / 0.8 / 4.7 | 0.249 / 0.280 / 0.284 | 21 / 38 / 32 | 40 / 0.05 / 0.2 | 321 / 257 / 231 | 98.8 / 96.9 / 96.9 |
| 0.005 | 0.019 / 0.019 / 0.019 | 0 / 0 / 0 | 0.252 / 0.278 / 0.283 | 22.5 / 40 / 33 | 49 / 0.16 / 0.22 | 321 / 259 / 230 | 98.9 / 96.8 / 96.9 |

- No nutrient limitation remains: biomass P:C stays at about 0.019 (floor about 0.012) and N:C at about 0.20 (floor 0.138); on the constant-temperature gradient at 15-25 C the P floor is no longer reached (P:C 0.019; it was 0.012-0.013 with PO4 0.002) and the rebalancing credit is about 0.
- Baseline targets: ambient close to energy-limited (22.5 endogenous-respiration days late in the season, min E:C 49), +5 and +10 C energy-limited (33-40 days, min E:C about 0.2), 96.8-98.9 % degraded, CUE 0.25-0.28.

## 34. Step 1.1: re-baseline, 5 pulses, R_endo on/off, grid size (2026-10-19)

**Setup (`baseline_1p1/`; `analysis/climate_endo_run.py`):** current minimal baseline (section 33), MMRT, real climate (pulse p uses climate year p), 20x20, 5 pulses x 5 seeds, R_endo 0 and 1 (R_stoich 1 in both); grid-size check at 50x50 (2 seeds) and 100x100 (1 seed), one pulse, R_endo 1. Run time here: 25 s per 5-pulse run at 20x20; 50x50 and 100x100 one-pulse runs about 20 and 35 s.

**Pulses and seeds (20x20, R_endo 1; ambient / +5 / +10 C).** Every pulse reproduces the first within seed noise: % degraded 98.8-98.9 / 96.8-96.9 / 96.6-96.9; day 95 % 319-321 / 255-258 / 229-236; CUE 0.247-0.250 / 0.277-0.279 / 0.281-0.285; pulses start at the same biomass (1.01-1.02). Across-seed SD: degradation <= 0.16 points, day 95 % <= 3.5 days, CUE <= 0.006.

**R_endo off vs on** (pulse means): degradation within 0.2 points; day 95 % within 1-5 days; CUE higher with R_endo off (0.264 / 0.300 / 0.300 vs 0.249 / 0.278 / 0.284), because with R_endo off energy shortfalls are floored away instead of paid by respiring biomass; min E:C 44 / 0 / 0 (off) vs 46 / 0.1 / 0.2 (on).

**Grid size** (pulse 1, R_endo 1):

| Side | % degraded | Day 50 % / 95 % | CUE | Max biomass | Drought kills per cell | Endogenous respiration C (share of respiration) | Days with any cell in deficit |
|---|---|---|---|---|---|---|---|
| 20 | 98.87 / 96.76 / 96.93 | 223, 319 / 167, 258 / 150, 229 | 0.249 / 0.279 / 0.283 | 41.3 / 50.6 / 53.4 | 0.70 / 0.78 / 0.80 | 1.7 (1.2 %) / 34.6 (20.3 %) / 35.9 (21.1 %) | 21 / 40 / 32 |
| 50 | 98.87 / 96.78 / 96.90 | 223, 320 / 167, 258 / 150, 230 | 0.247 / 0.277 / 0.284 | 40.5 / 50.4 / 53.2 | 0.70 / 0.78 / 0.79 | 2.3 (1.6 %) / 35.0 (20.5 %) / 36.7 (21.5 %) | 39 / 47 / 34 |
| 100 | 98.87 / 96.77 / 96.85 | 223, 320 / 167, 258 / 150, 231 | 0.249 / 0.279 / 0.284 | 40.9 / 50.2 / 52.8 | 0.70 / 0.80 / 0.80 | 2.2 (1.5 %) / 34.4 (20.2 %) / 37.7 (21.9 %) | 48 / 46 / 36 |

- Every amount-based output is grid-size invariant (differences within seed noise), so 20x20 with 2 seeds and one pulse is adequate for the Morris screening.
- The count of days with endogenous respiration counts a day if any cell is in deficit, so it grows with grid size (ambient 21 -> 39 -> 48) while the amount respired does not. Use the endogenous-respiration C, or its share of respiration, as the energy-limitation metric from here on (ambient about 1.5 %, +5 / +10 C about 20 %). The day counts in sections 30-33 are 20x20 values and comparable with each other only.

