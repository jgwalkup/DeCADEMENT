#!/usr/bin/env python3
"""
export_summary_data.py
======================
Export community-averaged summary data for all metrics used by the DeCA-DEMENT
figure scripts, without generating any figures.  One set of CSVs is written per
climate scenario.

WHAT IS EXPORTED
----------------
For every scenario the following are computed (mean ± SD across communities)
and written to <out-dir>/<scenario_label>/:

System-level time series  (index = day)
  system_resp.csv               – total respiration
  system_substrate_C.csv        – total substrate C (summed across types)
  system_monomer_C.csv          – total monomer C (summed across types)
  system_microbial_biomass_C.csv– total microbial biomass C (summed across taxa)
  system_microbial_energy.csv   – total microbial energy (summed across taxa)
  system_enzyme_C.csv           – total enzyme C (summed across types)
  system_CUE.csv                – emergent system CUE
  system_CAE.csv                – emergent system CAE

Substrate / monomer type breakdowns  (index = type, columns = day)
  substrates_by_type.csv        – substrate C per type × day
  monomers_by_type.csv          – monomer C per type × day

Carbon allocation  (index = taxon, columns = day;  per alive microbe)
  c_uptake_taxon.csv            – total C uptake
  c_biosynthesis_taxon.csv      – C to biosynthesis (BC_Uptake)
  c_energy_taxon.csv            – C to energy metabolism (uptake − biosynthesis)
  c_energy_fraction_taxon.csv   – energy C / total C  (0–1)
  c_biosyn_fraction_taxon.csv   – biosynthesis C / total C  (0–1)
  c_enzyme_taxon.csv            – C to enzyme production
  c_hsp_taxon.csv               – C to HSP production
  c_osmolyte_taxon.csv          – C to osmolyte production
  c_transport_taxon.csv         – C to transporter production

Carbon allocation fractions within biosynthesis  (index = taxon, columns = day)
  c_enzyme_biosyn_frac.csv      – enzyme C / biosynthesis C
  c_hsp_biosyn_frac.csv         – HSP C / biosynthesis C
  c_osmolyte_biosyn_frac.csv    – osmolyte C / biosynthesis C
  c_transport_biosyn_frac.csv   – transporter C / biosynthesis C

Energy allocation  (index = taxon, columns = day;  per alive microbe)
  e_production_taxon.csv        – energy produced
  e_cost_taxon.csv              – total energy cost
  e_gain_taxon.csv              – net energy gain
  e_transport_taxon.csv         – energy to transport
  e_maintenance_taxon.csv       – energy to maintenance
  e_hsp_taxon.csv               – energy to HSP production
  e_enzyme_taxon.csv            – energy to enzyme production
  e_osmolyte_taxon.csv          – energy to osmolyte production
  e_cost_fraction_taxon.csv     – total cost / production  (0–1)
  e_transport_frac.csv          – transport / total cost
  e_maintenance_frac.csv        – maintenance / total cost
  e_hsp_frac.csv                – HSP / total cost
  e_enzyme_frac.csv             – enzyme / total cost
  e_osmolyte_frac.csv           – osmolyte / total cost

Efficiency  (index = taxon, columns = day)
  cue_taxon.csv                 – per-taxon CUE
  cae_taxon.csv                 – per-taxon CAE

HSP  (index = taxon, columns = day)
  hsp_pool_taxon.csv            – total HSP pool C per taxon
  hsp_biomass_ratio_taxon.csv   – HSP pool C / biomass C per taxon

Each file written contains the community mean.  A paired *_sd.csv is written
alongside with the cross-community standard deviation.

PULSE SUMMARY
-------------
pulse_summary_mean.csv / pulse_summary_sd.csv  (written to the scenario root)
  One row per pulse for each system-level variable.

USAGE
-----
  python export_summary_data.py --base /path/to/output \\
      --treatment-map ambient_dir=Ambient mod_dir="Moderate Warming" ext_dir="Extreme Warming"

  python export_summary_data.py \\
      --scenario-dirs /out/ambient /out/mod /out/ext \\
      --scenario-labels Ambient "Moderate Warming" "Extreme Warming" \\
      --out-dir /path/to/summary
"""

from __future__ import annotations

import argparse
import gc
import os
import sys
import warnings
from concurrent.futures import ProcessPoolExecutor, wait, FIRST_COMPLETED
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

# ── Hardcoded fallback ────────────────────────────────────────────────────────
BASE = Path("src/output")

_HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(_HERE))

try:
    from scenario_utils import (
        load_series,
        pick_runs,
        discover_scenarios,
        add_scenario_args,
        read_runtime,
        read_gridsize,
        year_range,
        filter_cols,
    )
except ImportError as e:
    sys.exit(f"ERROR: could not import scenario_utils: {e}")


# ── Arithmetic helpers ────────────────────────────────────────────────────────

def safe_div(a: pd.DataFrame | pd.Series,
             b: pd.DataFrame | pd.Series) -> pd.DataFrame | pd.Series:
    with np.errstate(divide="ignore", invalid="ignore"):
        r = a / b
    if isinstance(r, pd.DataFrame):
        r.replace([np.inf, -np.inf], np.nan, inplace=True)
    else:
        r = r.replace([np.inf, -np.inf], np.nan)
    return r


def _squeeze(df: pd.DataFrame | None) -> pd.Series | None:
    """Collapse a 1-row DataFrame to a Series."""
    if df is None:
        return None
    if isinstance(df, pd.Series):
        return df
    if df.shape[0] == 1:
        return df.iloc[0]
    if df.shape[1] == 1:
        return df.iloc[:, 0]
    return df.sum(axis=0)


# ── Per-run loader ────────────────────────────────────────────────────────────

def compute_run(run_dir: Path) -> dict | None:
    """Load all series for one replicate and compute derived quantities."""

    def L(name):
        return load_series(run_dir, name)

    # ── Raw loads ──────────────────────────────────────────────────────────────
    SUB  = L("SubstratesSeries")      # substrates × days
    MON  = L("MonomersSeries")        # monomers  × days
    MIC  = L("MicrobesSeries")        # taxa      × days  (biomass C)
    ENE  = L("EnergySeries")          # taxa      × days  (energy)
    ENZ  = L("EnzymesSeries")         # enzymes   × days
    RESP = L("RespSeries")            # 1         × days
    CUE_sys = L("CUE_system")         # 1         × days
    CAE_sys = L("CAE_system")         # 1         × days
    CUE_tax = L("CUE_TaxonSeries")    # taxa      × days
    CAE_tax = L("CAE_TaxonSeries")    # taxa      × days

    # C allocation
    CU   = L("C_Uptake_Taxon")        # taxa × days
    BC   = L("BC_Uptake_Taxon")       # taxa × days
    EC   = L("Enzyme_C_Taxon")        # taxa × days
    HC   = L("HSP_C_Taxon")           # taxa × days
    OC   = L("Osmo_C_Taxon")          # taxa × days
    TC   = L("Transport_C_Taxon")     # taxa × days

    # E allocation
    EP   = L("E_Prod_Taxon")          # taxa × days
    EGN  = L("Microbe_E_Gain")        # taxa × days
    ECT  = L("Microbe_E_Cost")        # taxa × days
    ET   = L("Transport_E_Taxon")     # taxa × days
    EM   = L("Maint_E_Taxon")         # taxa × days
    EH   = L("HSP_E_Taxon")           # taxa × days
    EE   = L("Enzyme_E_Taxon")        # taxa × days
    EO   = L("Osmo_E_Taxon")          # taxa × days

    # HSP
    HSP  = L("Taxon_HSP_Series")      # (taxa × HSP_type) × days

    # ── Derived: C allocation ──────────────────────────────────────────────────
    energy_C      = (CU - BC).clip(lower=0) if (CU is not None and BC is not None) else None
    energy_frac   = safe_div(energy_C, CU).clip(0, 1) if energy_C is not None else None
    biosyn_frac   = safe_div(BC, CU).clip(0, 1)       if (BC is not None and CU is not None) else None
    enz_bf        = safe_div(EC, BC).clip(0, 1)        if (EC is not None and BC is not None) else None
    hsp_bf        = safe_div(HC, BC).clip(0, 1)        if (HC is not None and BC is not None) else None
    osm_bf        = safe_div(OC, BC).clip(0, 1)        if (OC is not None and BC is not None) else None
    trans_bf      = safe_div(TC, BC).clip(0, 1)        if (TC is not None and BC is not None) else None

    # ── Derived: E allocation ──────────────────────────────────────────────────
    e_total_cost = None
    if any(x is not None for x in [ET, EM, EH, EE, EO]):
        parts = [x for x in [ET, EM, EH, EE, EO] if x is not None]
        idx0, cols0 = parts[0].index, parts[0].columns
        e_total_cost = parts[0].copy()
        for p in parts[1:]:
            cc = sorted(set(e_total_cost.columns) & set(p.columns))
            e_total_cost[cc] = e_total_cost[cc].add(p[cc], fill_value=0)

    def _efrac(x):
        if x is None or e_total_cost is None:
            return None
        cc = sorted(set(x.columns) & set(e_total_cost.columns))
        return safe_div(x[cc], e_total_cost[cc]).clip(0, 1)

    cost_frac  = safe_div(e_total_cost, EP).clip(0, 1) if (e_total_cost is not None and EP is not None) else None

    # net-energy fraction = 1 - cost fraction  (matches energy_allocation_figures)
    e_net_frac = (1 - cost_frac).clip(0, 1) if cost_frac is not None else None

    # HSP production fraction = HSP E cost / total E cost  (matches hsp_annotated_e_heatmap)
    hsp_prod_frac = None
    if EH is not None and e_total_cost is not None:
        cc = sorted(set(EH.columns) & set(e_total_cost.columns))
        hsp_prod_frac = safe_div(EH[cc], e_total_cost[cc]).clip(0, 1)

    # ── Derived: combined C cost ────────────────────────────────────────────────
    # Sum of carbon allocation costs (enzyme + HSP + osmolyte), per taxon.
    # Matches absolute_c_pathway_figures' total_cost (community-mean is taken
    # downstream; summing per-taxon here is identical since the mean is linear).
    c_total_cost = None
    c_parts = [x for x in [EC, HC, OC] if x is not None]
    if c_parts:
        c_total_cost = c_parts[0].copy()
        for p in c_parts[1:]:
            cc = sorted(set(c_total_cost.columns) & set(p.columns))
            c_total_cost[cc] = c_total_cost[cc].add(p[cc], fill_value=0)

    # ── Derived: HSP pool ──────────────────────────────────────────────────────
    hsp_pool = None
    hsp_biomass_ratio = None
    if HSP is not None:
        hsp_pool = HSP.groupby(HSP.index).sum()   # collapse HSP types
        if MIC is not None:
            common_taxa = sorted(set(hsp_pool.index) & set(MIC.index))
            common_days = sorted(set(hsp_pool.columns) & set(MIC.columns))
            if common_taxa and common_days:
                h = hsp_pool.loc[common_taxa, common_days].astype(float)
                m = MIC.loc[common_taxa, common_days].astype(float)
                hsp_biomass_ratio = safe_div(h, m).clip(lower=0)

    # ── Derived: mass-specific respiration ──────────────────────────────────────
    # Community respiration per unit total microbial biomass C (1-D, over days).
    resp_per_biomass = None
    _resp_s = _squeeze(RESP)
    if _resp_s is not None and MIC is not None:
        _bio_tot = MIC.sum(axis=0)            # total biomass C per day
        cc = sorted(set(_resp_s.index) & set(_bio_tot.index))
        if cc:
            resp_per_biomass = safe_div(_resp_s[cc], _bio_tot[cc])

    # ── Derived: monomer-pool breakdown + necromass substrate pools ──────────────
    # Non-plant monomers are tracked individually; plant-litter monomer = the rest.
    # Candidate label lists guard against naming variants across model versions.
    def _pick(df, candidates):
        if df is None:
            return None
        for c in candidates:
            if c in df.index:
                return df.loc[c]
        return None

    _NONPLANT_MON = {
        "nh4":      ["NH4"],
        "po4":      ["PO4"],
        "fermprod": ["FermProd"],
        "deadmic":  ["Mon_DeadMic", "DeadMic"],
        "deadenz":  ["Mon_DeadEnz", "DeadEnz"],
    }
    nh4_mon = po4_mon = ferm_mon = deadmic_mon = deadenz_mon = None
    plant_litter_monomer = None
    if MON is not None:
        _matched = []
        _resolved = {}
        for nm, cands in _NONPLANT_MON.items():
            lbl = next((c for c in cands if c in MON.index), None)
            _resolved[nm] = MON.loc[lbl] if lbl is not None else None
            if lbl is not None:
                _matched.append(lbl)
        nh4_mon, po4_mon, ferm_mon = _resolved["nh4"], _resolved["po4"], _resolved["fermprod"]
        deadmic_mon, deadenz_mon  = _resolved["deadmic"], _resolved["deadenz"]
        _plant_rows = [r for r in MON.index if r not in set(_matched)]
        plant_litter_monomer = MON.loc[_plant_rows].sum(axis=0) if _plant_rows else None

    deadmic_sub = _pick(SUB, ["DeadMic"])
    deadenz_sub = _pick(SUB, ["DeadEnz"])

    # ── Community-mean-per-alive-microbe pathway series (1-D over days) ──────────
    # Mean across taxa per day, mirroring absolute_c/absolute_e compute_run's
    # _comm().  Aggregated across communities as 1-D series, the resulting mean
    # and SD match raw-mode rep_mean/rep_std of the per-community taxon-mean
    # exactly (day-indexed, so no taxon-reference anchoring is involved).
    def _comm_mean(df):
        return df.mean(axis=0, skipna=True) if df is not None else None

    return dict(
        # pools
        SUB=SUB, MON=MON, MIC=MIC, ENE=ENE, ENZ=ENZ,
        RESP=_squeeze(RESP),
        resp_per_biomass=resp_per_biomass,
        # monomer breakdown + necromass pools (1-D series)
        plant_litter_monomer=plant_litter_monomer,
        nh4_monomer=nh4_mon, po4_monomer=po4_mon, fermprod_monomer=ferm_mon,
        deadmic_monomer=deadmic_mon, deadenz_monomer=deadenz_mon,
        deadmic_substrate=deadmic_sub, deadenz_substrate=deadenz_sub,
        # community-mean-per-alive-microbe pathway series (1-D)
        comm_uptake_C_per_microbe=_comm_mean(CU),
        comm_biosynthesis_C_per_microbe=_comm_mean(BC),
        comm_enzyme_C_per_microbe=_comm_mean(EC),
        comm_hsp_C_per_microbe=_comm_mean(HC),
        comm_osmolyte_C_per_microbe=_comm_mean(OC),
        comm_production_E_per_microbe=_comm_mean(EP),
        comm_transport_E_per_microbe=_comm_mean(ET),
        comm_maintenance_E_per_microbe=_comm_mean(EM),
        comm_hsp_E_per_microbe=_comm_mean(EH),
        comm_enzyme_E_per_microbe=_comm_mean(EE),
        comm_osmolyte_E_per_microbe=_comm_mean(EO),
        CUE_sys=_squeeze(CUE_sys), CAE_sys=_squeeze(CAE_sys),
        CUE_tax=CUE_tax, CAE_tax=CAE_tax,
        # C
        CU=CU, BC=BC, energy_C=energy_C,
        energy_frac=energy_frac, biosyn_frac=biosyn_frac,
        EC=EC, HC=HC, OC=OC, TC=TC,
        enz_bf=enz_bf, hsp_bf=hsp_bf, osm_bf=osm_bf, trans_bf=trans_bf,
        # E
        EP=EP, EGN=EGN, ECT=ECT,
        ET=ET, EM=EM, EH=EH, EE=EE, EO=EO,
        e_total_cost=e_total_cost, cost_frac=cost_frac,
        e_net_frac=e_net_frac, hsp_prod_frac=hsp_prod_frac,
        c_total_cost=c_total_cost,
        trans_ef=_efrac(ET), maint_ef=_efrac(EM),
        hsp_ef=_efrac(EH), enz_ef=_efrac(EE), osm_ef=_efrac(EO),
        # HSP
        hsp_pool=hsp_pool, hsp_biomass_ratio=hsp_biomass_ratio,
    )


# ── Worker count ──────────────────────────────────────────────────────────────

def _resolve_jobs(requested: int) -> int:
    """Resolve the worker count.

    requested == 0  → auto: use $SLURM_CPUS_PER_TASK if set, else os.cpu_count().
    requested >= 1  → use that many (clamped to the number of available CPUs).
    """
    avail = int(os.environ.get("SLURM_CPUS_PER_TASK", 0)) or (os.cpu_count() or 1)
    if requested <= 0:
        return max(1, avail)
    return max(1, min(requested, avail))


# ── Streaming aggregation (memory independent of community count) ───────────────

# Pool keys whose per-run frames are also summed across rows → system totals.
_POOL_KEYS = ("SUB", "MON", "MIC", "ENE", "ENZ")


def _fold(store: dict, key: str, val):
    """Fold one run's value into running sum / sumsq / count (numpy, in-place).

    The FIRST run to carry a key fixes the reference index (and columns, for
    2-D), and every later run is aligned onto it — exactly matching the legacy
    _agg_df / _agg_series, which aggregate over frames[0]'s shape rather than a
    union.  Cells missing from a run contribute NaN (counted out), so the
    finalised mean/SD ignore them, ddof=1.

    Accumulation is done on numpy arrays with in-place adds rather than pandas
    DataFrame arithmetic: the latter allocates a fresh block-managed frame on
    every fold, which serialises the whole pipeline at scale.
    """
    is_df = isinstance(val, pd.DataFrame)
    rec = store.get(key)
    if rec is None:
        if is_df:
            idx, cols = val.index, val.columns
            shape = (len(idx), len(cols))
        else:
            idx, cols = val.index, None
            shape = (len(idx),)
        rec = store[key] = {
            "is_df": is_df, "index": idx, "columns": cols,
            "sum": np.zeros(shape, dtype="float64"),
            "sumsq": np.zeros(shape, dtype="float64"),
            "count": np.zeros(shape, dtype="float64"),
        }

    # Align to the reference shape, then drop to a numpy array. Skip the (costly)
    # reindex when the run already matches the reference exactly.
    if rec["is_df"]:
        if val.index.equals(rec["index"]) and val.columns.equals(rec["columns"]):
            arr = val.to_numpy(dtype="float64", copy=False)
        else:
            arr = val.reindex(index=rec["index"],
                              columns=rec["columns"]).to_numpy(dtype="float64")
    else:
        if val.index.equals(rec["index"]):
            arr = val.to_numpy(dtype="float64", copy=False)
        else:
            arr = val.reindex(rec["index"]).to_numpy(dtype="float64")

    mask = ~np.isnan(arr)
    arr0 = np.where(mask, arr, 0.0)
    rec["sum"]   += arr0
    rec["sumsq"] += arr0 * arr0
    rec["count"] += mask        # bool adds 1.0 where present


def _finalize_store(store: dict) -> dict:
    """Turn sum/sumsq/count accumulators into (mean, std) with ddof=1.

    var = (Σx² − (Σx)²/n) / (n−1) over non-NaN samples — algebraically the
    same as np.nanstd(ddof=1); clipped at 0 to absorb float cancellation.
    NaN where n < 2 (std) or n < 1 (mean).  Accumulators are numpy arrays;
    results are rebuilt as pandas using the stored reference index/columns.
    """
    out = {}
    for key, rec in store.items():
        s, ss, c = rec["sum"], rec["sumsq"], rec["count"]
        with np.errstate(invalid="ignore", divide="ignore"):
            cpos  = np.where(c > 0, c, np.nan)
            mean  = np.where(c >= 1, s / cpos, np.nan)
            denom = np.where(c >= 2, c - 1.0, np.nan)
            var   = np.clip((ss - (s * s) / cpos) / denom, 0.0, None)
            std   = np.sqrt(var)
        if rec["is_df"]:
            mean = pd.DataFrame(mean, index=rec["index"], columns=rec["columns"])
            std  = pd.DataFrame(std,  index=rec["index"], columns=rec["columns"])
        else:
            mean = pd.Series(mean, index=rec["index"])
            std  = pd.Series(std,  index=rec["index"])
        out[key] = (mean, std)
    return out


class _StreamAgg:
    """Folds runs one at a time into mean/SD accumulators, then serves results.

    Memory is O(number of distinct keys) — a few hundred MB of accumulators —
    independent of how many communities are processed, so it scales to
    thousands of runs.  Provides the same three reductions the writers need:
      get_series(key)   1-D series aggregation        (RESP, CUE_sys, …)
      get_df(key)       full taxa/types × days agg    (all _taxon / _by_type)
      get_sumrows(key)  sum-across-rows → series agg  (system_*_C totals)
    """

    def __init__(self):
        self._df, self._sr, self._ser = {}, {}, {}
        self.n_loaded = 0
        self._fdf = self._fsr = self._fser = None

    def fold(self, d):
        if d is None:
            return
        self.n_loaded += 1
        for key, val in d.items():
            if val is None:
                continue
            if isinstance(val, pd.Series):
                _fold(self._ser, key, val)
            elif isinstance(val, pd.DataFrame):
                _fold(self._df, key, val)
                if key in _POOL_KEYS:
                    _fold(self._sr, key, val.sum(axis=0))

    def finalize(self):
        self._fdf = _finalize_store(self._df)
        self._fsr = _finalize_store(self._sr)
        self._fser = _finalize_store(self._ser)
        self._df = self._sr = self._ser = None   # release accumulators

    def get_df(self, key):
        return self._fdf.get(key, (None, None))

    def get_series(self, key):
        return self._fser.get(key, (None, None))

    def get_sumrows(self, key):
        return self._fsr.get(key, (None, None))


def _stream_load_and_fold(run_dirs: list[Path], n_jobs: int, acc: _StreamAgg):
    """Load runs (bounded parallel) and fold each into *acc*, discarding as we go.

    Outstanding futures are capped at n_jobs + buffer so completed run dicts are
    folded and freed immediately instead of piling up in the result queue — that
    cap, not the community count, sets peak memory.
    """
    n = len(run_dirs)

    if n_jobs <= 1 or n <= 1:
        for i, rd in enumerate(run_dirs):
            print(f"  [{i+1}/{n}] {rd.name}", end=" ", flush=True)
            d = compute_run(rd)
            ok = d is not None
            acc.fold(d)
            del d
            print("ok" if ok else "FAILED")
        return

    max_outstanding = n_jobs + 4
    print(f"  Streaming {n} runs across {n_jobs} workers "
          f"(≤{max_outstanding} in flight)...", flush=True)
    done = submitted = 0
    with ProcessPoolExecutor(max_workers=n_jobs) as ex:
        futures = {}
        while submitted < min(max_outstanding, n):
            rd = run_dirs[submitted]
            futures[ex.submit(compute_run, rd)] = rd.name
            submitted += 1
        while futures:
            finished, _pending = wait(set(futures), return_when=FIRST_COMPLETED)
            for fut in finished:
                name = futures.pop(fut)
                try:
                    d = fut.result()
                    ok = d is not None
                except Exception as e:
                    d, ok = None, False
                    print(f"  [run {name}] FAILED: {e}", flush=True)
                acc.fold(d)
                del d
                done += 1
                print(f"  [{done}/{n}] {name} {'ok' if ok else 'FAILED'}", flush=True)
                if submitted < n:
                    rd = run_dirs[submitted]
                    futures[ex.submit(compute_run, rd)] = rd.name
                    submitted += 1


# ── CSV writer ────────────────────────────────────────────────────────────────

def write_pair(out_dir: Path, stem: str,
               mean: pd.DataFrame | pd.Series | None,
               std:  pd.DataFrame | pd.Series | None) -> None:
    """Write mean CSV and SD CSV for one variable."""
    if mean is None:
        return
    mean_path = out_dir / f"{stem}_mean.csv"
    sd_path   = out_dir / f"{stem}_sd.csv"
    if isinstance(mean, pd.Series):
        mean.to_csv(mean_path, header=["mean"])
        if std is not None:
            std.to_csv(sd_path, header=["sd"])
    else:
        mean.to_csv(mean_path)
        if std is not None:
            std.to_csv(sd_path)


# ── Pulse aggregation ─────────────────────────────────────────────────────────

def pulse_agg(s: pd.Series | None, n_pulses: int,
              days_per_pulse: int, mode: str = "sum") -> pd.Series:
    if s is None or s.empty:
        return pd.Series([np.nan] * n_pulses,
                         index=[f"Pulse {p}" for p in range(1, n_pulses + 1)])
    values = []
    for p in range(1, n_pulses + 1):
        d0, d1 = (p - 1) * days_per_pulse + 1, p * days_per_pulse
        win = s.loc[d0:d1]
        if win.empty:
            values.append(np.nan)
        elif mode == "sum":
            values.append(win.sum())
        else:
            values.append(win.iloc[-1])
    return pd.Series(values, index=[f"Pulse {p}" for p in range(1, n_pulses + 1)])


# ── Scenario-list resolver ────────────────────────────────────────────────────

def _parse_treatment_map(entries):
    if not entries:
        return {}
    result = {}
    for entry in entries:
        if "=" not in entry:
            print(f"  WARNING: --treatment-map entry has no '=': {entry!r}  (ignored)")
            continue
        stem, label = entry.split("=", 1)
        result[stem.strip()] = label.strip()
    return result


def resolve_scenario_list(args) -> list[tuple[str, list[Path]]]:
    tmap = _parse_treatment_map(getattr(args, "treatment_map", None))

    if args.scenario_dirs:
        paths = [Path(p) for p in args.scenario_dirs]
        raw_labels = (args.scenario_labels
                      if args.scenario_labels and len(args.scenario_labels) == len(paths)
                      else [p.stem for p in paths])
        result = []
        for raw_label, sdir in zip(raw_labels, paths):
            label = tmap.get(sdir.stem, raw_label)
            runs  = pick_runs(sdir, min_cols=100)
            if not runs:
                print(f"  WARNING: no valid runs found in {sdir}")
            result.append((label, runs))
        return result

    base = Path(args.base) if args.base else BASE
    found = discover_scenarios(base)
    if not found:
        sys.exit(f"ERROR: no scenario subdirectories found under {base}")

    pos_labels = args.scenario_labels or []
    result = []
    for i, (stem, sdir) in enumerate(found):
        if stem in tmap:
            label = tmap[stem]
        elif i < len(pos_labels):
            label = pos_labels[i]
        else:
            label = stem
        runs = pick_runs(sdir, min_cols=100)
        result.append((label, runs))
    return result


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        # Some versions of scenario_utils.add_scenario_args() also define
        # --treatment-map.  'resolve' lets this script's definition override the
        # shared one instead of raising a conflicting-option error.
        conflict_handler="resolve",
    )
    add_scenario_args(parser)
    parser.add_argument(
        "--treatment-map", nargs="+", default=None, metavar="DIR=LABEL",
        help="Map directory stems to display labels (DIR_STEM=LABEL). "
             "Takes precedence over --scenario-labels.",
    )
    parser.add_argument(
        "--out-dir", default=None, metavar="DIR",
        help="Root directory to write summary CSVs "
             "(default: <base>/export_summary_data).",
    )
    parser.add_argument(
        "--jobs", "-j", type=int, default=1, metavar="N",
        help="Number of parallel worker processes for per-run loading. "
             "1 (default) = serial; 0 = auto ($SLURM_CPUS_PER_TASK or all cores); "
             "N = use N workers (clamped to available CPUs).",
    )

    args = parser.parse_args()

    scenario_list = resolve_scenario_list(args)
    if not scenario_list:
        sys.exit("ERROR: no scenarios resolved.")

    n_jobs = _resolve_jobs(args.jobs)
    if n_jobs > 1:
        print(f"Parallel per-run loading: {n_jobs} workers\n")

    base = Path(args.base) if args.base else BASE
    root_out = Path(args.out_dir) if args.out_dir else base / "export_summary_data"
    root_out.mkdir(parents=True, exist_ok=True)

    # Runtime
    first_run = next((r for _, runs in scenario_list for r in runs), None)
    rt = read_runtime(first_run)
    n_pulses       = rt["pulse"]
    days_per_pulse = rt["end_time"]
    total_days     = n_pulses * days_per_pulse

    col_start, col_end = year_range(args.years, n_pulses, days_per_pulse)
    print(f"\nRuntime: {n_pulses} pulses × {days_per_pulse} days = {total_days} total days")
    if args.years:
        print(f"Filtering to days {col_start}–{col_end}")
    print(f"Output root: {root_out}\n")

    pulse_rows_mean = {}
    pulse_rows_sd   = {}

    acc = None
    for label, run_dirs in scenario_list:
        n = len(run_dirs)
        print(f"── {label}  ({n} communities) ──────────────────────────────")

        # Stream runs into mean/SD accumulators (memory independent of count).
        acc = None
        gc.collect()
        acc = _StreamAgg()
        _stream_load_and_fold(run_dirs, n_jobs, acc)
        acc.finalize()

        loaded = acc.n_loaded
        print(f"  Loaded {loaded}/{n}\n")

        # ── Output directory for this scenario ────────────────────────────────
        safe_label = label.replace(" ", "_").replace("/", "-")
        out_dir = root_out / safe_label
        out_dir.mkdir(parents=True, exist_ok=True)

        # ── Self-contained metadata ───────────────────────────────────────────
        # Record gridsize and runtime so summary-mode figure scripts can apply
        # gridsize normalisation and draw year labels without the raw runs.
        _gs = read_gridsize(run_dirs[0]) if run_dirs else 1
        pd.Series(
            {"gridsize": _gs, "pulse": n_pulses, "end_time": days_per_pulse,
             "col_start": col_start, "col_end": col_end},
            name="value",
        ).to_csv(out_dir / "summary_meta.csv", header=["value"])

        def write(stem: str, key: str, is_series: bool = False,
                  col_filter: bool = True):
            """Aggregate and write mean+SD CSVs for one key (from the stream)."""
            if is_series:
                mean, std = acc.get_series(key)
                if mean is not None and col_filter:
                    mean = mean.loc[col_start:col_end]
                    if std is not None:
                        std = std.loc[col_start:col_end]
            else:
                mean, std = acc.get_df(key)
                if mean is not None and col_filter:
                    cols = [c for c in mean.columns if col_start <= c <= col_end]
                    mean = mean[cols]
                    if std is not None:
                        std = std[cols]
            write_pair(out_dir, stem, mean, std)

        # ── System 1-D series ──────────────────────────────────────────────────
        print("  Writing system-level series...", flush=True)

        # Respiration (1-D)
        write("system_resp",               "RESP",    is_series=True)
        write("system_resp_per_biomass",   "resp_per_biomass", is_series=True)
        write("system_CUE",                "CUE_sys", is_series=True)
        write("system_CAE",                "CAE_sys", is_series=True)

        # Pool totals: sum across types/taxa → 1-D  (streamed sum-across-rows)
        for stem, key in [("system_substrate_C",      "SUB"),
                           ("system_monomer_C",        "MON"),
                           ("system_microbial_biomass_C", "MIC"),
                           ("system_microbial_energy", "ENE"),
                           ("system_enzyme_C",         "ENZ")]:
            m, s = acc.get_sumrows(key)
            if m is not None:
                m = m.loc[col_start:col_end]
                if s is not None:
                    s = s.loc[col_start:col_end]
            write_pair(out_dir, stem, m, s)

        # Monomer-pool breakdown: plant-litter total + each non-plant pool (1-D).
        # system_monomer_C above remains the full total across all monomer types.
        print("  Writing monomer-pool breakdown...", flush=True)
        write("system_plant_litter_monomer_C", "plant_litter_monomer", is_series=True)
        write("system_nh4_monomer_C",          "nh4_monomer",      is_series=True)
        write("system_po4_monomer_C",          "po4_monomer",      is_series=True)
        write("system_fermprod_monomer_C",     "fermprod_monomer", is_series=True)
        write("system_deadmic_monomer_C",      "deadmic_monomer",  is_series=True)
        write("system_deadenz_monomer_C",      "deadenz_monomer",  is_series=True)

        # Necromass substrate pools, tracked separately (1-D).
        # system_substrate_C above remains the full total across all substrate types.
        print("  Writing necromass substrate pools...", flush=True)
        write("system_deadmic_substrate_C", "deadmic_substrate", is_series=True)
        write("system_deadenz_substrate_C", "deadenz_substrate", is_series=True)

        # ── Type-level breakdowns ──────────────────────────────────────────────
        print("  Writing type-level breakdowns...", flush=True)
        write("substrates_by_type", "SUB")
        write("monomers_by_type",   "MON")
        write("microbes_by_taxon",  "MIC")   # per-taxon biomass (taxa × days)

        # ── Community-mean-per-alive-microbe pathway series (1-D) ───────────────
        # Mean across taxa, aggregated across communities; matches raw-mode
        # absolute_c/absolute_e exactly (mean and SD).  total_cost / net_gain are
        # derived in the figure scripts from these components.
        print("  Writing community-mean pathway series...", flush=True)
        write("comm_uptake_C_per_microbe",       "comm_uptake_C_per_microbe",       is_series=True)
        write("comm_biosynthesis_C_per_microbe", "comm_biosynthesis_C_per_microbe", is_series=True)
        write("comm_enzyme_C_per_microbe",       "comm_enzyme_C_per_microbe",       is_series=True)
        write("comm_hsp_C_per_microbe",          "comm_hsp_C_per_microbe",          is_series=True)
        write("comm_osmolyte_C_per_microbe",     "comm_osmolyte_C_per_microbe",     is_series=True)
        write("comm_production_E_per_microbe",   "comm_production_E_per_microbe",   is_series=True)
        write("comm_transport_E_per_microbe",    "comm_transport_E_per_microbe",    is_series=True)
        write("comm_maintenance_E_per_microbe",  "comm_maintenance_E_per_microbe",  is_series=True)
        write("comm_hsp_E_per_microbe",          "comm_hsp_E_per_microbe",          is_series=True)
        write("comm_enzyme_E_per_microbe",       "comm_enzyme_E_per_microbe",       is_series=True)
        write("comm_osmolyte_E_per_microbe",     "comm_osmolyte_E_per_microbe",     is_series=True)

        # ── Carbon allocation ──────────────────────────────────────────────────
        print("  Writing C allocation...", flush=True)
        for stem, key in [
            ("c_uptake_taxon",          "CU"),
            ("c_biosynthesis_taxon",    "BC"),
            ("c_energy_taxon",          "energy_C"),
            ("c_energy_fraction_taxon", "energy_frac"),
            ("c_biosyn_fraction_taxon", "biosyn_frac"),
            ("c_enzyme_taxon",          "EC"),
            ("c_hsp_taxon",             "HC"),
            ("c_osmolyte_taxon",        "OC"),
            ("c_transport_taxon",       "TC"),
            ("c_total_cost",            "c_total_cost"),
            ("c_enzyme_biosyn_frac",    "enz_bf"),
            ("c_hsp_biosyn_frac",       "hsp_bf"),
            ("c_osmolyte_biosyn_frac",  "osm_bf"),
            ("c_transport_biosyn_frac", "trans_bf"),
        ]:
            write(stem, key)

        # ── Energy allocation ──────────────────────────────────────────────────
        print("  Writing E allocation...", flush=True)
        for stem, key in [
            ("e_production_taxon",   "EP"),
            ("e_gain_taxon",         "EGN"),
            ("e_cost_taxon",         "ECT"),
            ("e_transport_taxon",    "ET"),
            ("e_maintenance_taxon",  "EM"),
            ("e_hsp_taxon",          "EH"),
            ("e_enzyme_taxon",       "EE"),
            ("e_osmolyte_taxon",     "EO"),
            ("e_total_cost_taxon",   "e_total_cost"),
            ("e_cost_fraction_taxon","cost_frac"),
            ("e_net_fraction",       "e_net_frac"),
            ("e_transport_frac",     "trans_ef"),
            ("e_maintenance_frac",   "maint_ef"),
            ("e_hsp_frac",           "hsp_ef"),
            ("e_enzyme_frac",        "enz_ef"),
            ("e_osmolyte_frac",      "osm_ef"),
        ]:
            write(stem, key)

        # ── Efficiency ────────────────────────────────────────────────────────
        print("  Writing efficiency...", flush=True)
        write("cue_taxon", "CUE_tax")
        write("cae_taxon", "CAE_tax")

        # ── HSP ───────────────────────────────────────────────────────────────
        print("  Writing HSP...", flush=True)
        write("hsp_pool_taxon",          "hsp_pool")
        write("hsp_biomass_ratio_taxon", "hsp_biomass_ratio")
        write("hsp_production_fraction",  "hsp_prod_frac")

        # ── Pulse summary ─────────────────────────────────────────────────────
        print("  Computing pulse summaries...", flush=True)
        resp_m, resp_s = acc.get_series("RESP")
        sub_m,  sub_s  = acc.get_sumrows("SUB")
        mic_m,  mic_s  = acc.get_sumrows("MIC")
        ene_m,  ene_s  = acc.get_sumrows("ENE")
        cue_m,  cue_s  = acc.get_series("CUE_sys")
        cae_m,  cae_s  = acc.get_series("CAE_sys")

        pulse_mean_row = {}
        pulse_sd_row   = {}
        for vname, series_m, series_s, mode in [
            ("Resp",    resp_m,  resp_s,  "sum"),
            ("SubC",    sub_m,   sub_s,   "end"),
            ("MicC",    mic_m,   mic_s,   "end"),
            ("MicE",    ene_m,   ene_s,   "end"),
            ("CUE_sys", cue_m,   cue_s,   "end"),
            ("CAE_sys", cae_m,   cae_s,   "end"),
        ]:
            pa_m = pulse_agg(series_m, n_pulses, days_per_pulse, mode)
            pa_s = pulse_agg(series_s, n_pulses, days_per_pulse, mode)
            for pidx in pa_m.index:
                pulse_mean_row.setdefault(pidx, {})[f"{vname}_{label}"] = pa_m[pidx]
                pulse_sd_row.setdefault(pidx, {})[f"{vname}_{label}"]   = pa_s[pidx]

        pulse_rows_mean[label] = pulse_mean_row
        pulse_rows_sd[label]   = pulse_sd_row

        files = sorted(out_dir.glob("*.csv"))
        print(f"  Wrote {len(files)} files to {out_dir}\n")

        # Drop this scenario's finalised aggregates before the next scenario.
        acc = None
        gc.collect()

    # ── Cross-scenario pulse summary ──────────────────────────────────────────
    print("Writing cross-scenario pulse summary...")
    all_pulse_idx = [f"Pulse {p}" for p in range(1, n_pulses + 1)]

    # Build flat DataFrames: rows = pulses, columns = variable_scenario
    mean_records = {p: {} for p in all_pulse_idx}
    sd_records   = {p: {} for p in all_pulse_idx}
    for label in [l for l, _ in scenario_list]:
        for pidx in all_pulse_idx:
            mean_records[pidx].update(pulse_rows_mean.get(label, {}).get(pidx, {}))
            sd_records[pidx].update(pulse_rows_sd.get(label, {}).get(pidx, {}))

    pulse_mean_df = pd.DataFrame(mean_records).T
    pulse_sd_df   = pd.DataFrame(sd_records).T
    pulse_mean_df.index.name = "pulse"
    pulse_sd_df.index.name   = "pulse"
    pulse_mean_df.to_csv(root_out / "pulse_summary_mean.csv")
    pulse_sd_df.to_csv(root_out   / "pulse_summary_sd.csv")
    print(f"  pulse_summary_mean.csv")
    print(f"  pulse_summary_sd.csv")
    print(f"\nDone.  All output written to {root_out}\n")


if __name__ == "__main__":
    main()
