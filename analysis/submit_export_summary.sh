#!/bin/bash
# =============================================================================
# submit_export_summary.sh
# SLURM job for export_summary_data.py (streaming per-run aggregation) on HPC3.
#
# Project / base / output are ARGUMENTS so one script serves all scenario
# directories (mmrtTempSens, mmrtNoTemp, arrTempSens, arrNoTemp, …).
#
# USAGE
# -----
#   sbatch submit_export_summary.sh --project DIR --base DIR --out-dir DIR [OPTIONS]
#
# Required:
#   --project  DIR   Project root (contains analysis/, src/, input/); the job cd's here
#   --base     DIR   Raw simulation root to summarise (the temp_gradient dir)
#   --out-dir  DIR   Directory to write summary CSVs
#
# Optional:
#   --conda-env NAME Conda environment            (default: decadementenv)
#   --jobs      INT  Worker processes; 0 = use $SLURM_CPUS_PER_TASK (default: 0)
#   --help           Show this help and exit
#
# SBATCH resources (cpus/mem/time/account/partition) are set as directives below.
# Override them at submit time with sbatch flags, e.g.:
#   sbatch --cpus-per-task=16 --time=03:00:00 submit_export_summary.sh ...
#
# EXAMPLES
# --------
#   sbatch submit_export_summary.sh \
#       --project /pub/$USER/DeCADEMENT20260824mmrtTempSens \
#       --base    /pub/$USER/DeCADEMENT20260824mmrtTempSens/simulations_20260826mmrtTempSens/temp_gradient \
#       --out-dir /pub/$USER/DeCADEMENT20260824mmrtTempSens/simulations_20260826mmrtTempSens/summary
#
# NOTES
# -----
#   • Memory is independent of community count (streaming fold): cores are the
#     lever, not --mem.  I/O on shared storage tends to be the ceiling past
#     ~16-32 workers, so cpus-per-task=32 is a sensible default.
#   • --mem peak ≈ (workers + 4) runs in flight + ~200 MB accumulators; 48G is
#     conservative headroom — trim it after checking `seff <jobid>`.
#   • Communities with corrupt taxon identities (duplicate labels) are skipped
#     and named in the log; the run still completes over the good communities.
# =============================================================================

#SBATCH --job-name=export_summary
#SBATCH --account=allisons_lab
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=48G
#SBATCH --time=06:00:00
#SBATCH --output=logs/export_summary_%j.out
#SBATCH --error=logs/export_summary_%j.err

set -euo pipefail

# ── Defaults ─────────────────────────────────────────────────────────────────
PROJECT=""
BASE=""
OUTDIR=""
CONDA_ENV="decadementenv"
JOBS=0          # 0 → export_summary_data.py auto-detects SLURM_CPUS_PER_TASK

# ── Argument parsing ─────────────────────────────────────────────────────────
while [[ $# -gt 0 ]]; do
    case "$1" in
        --project)   PROJECT="$2";   shift 2 ;;
        --base)      BASE="$2";      shift 2 ;;
        --out-dir)   OUTDIR="$2";    shift 2 ;;
        --conda-env) CONDA_ENV="$2"; shift 2 ;;
        --jobs|-j)   JOBS="$2";      shift 2 ;;
        --help|-h)
            sed -n '2,50p' "$0"; exit 0 ;;
        *)
            echo "ERROR: unknown or unexpected argument '$1'" >&2
            echo "       run with --help for usage." >&2
            exit 1 ;;
    esac
done

# ── Validate required args ───────────────────────────────────────────────────
missing=0
for pair in "PROJECT:--project" "BASE:--base" "OUTDIR:--out-dir"; do
    var="${pair%%:*}"; flag="${pair##*:}"
    if [[ -z "${!var}" ]]; then
        echo "ERROR: ${flag} is required." >&2; missing=1
    fi
done
[[ "${missing}" -eq 1 ]] && { echo "Run with --help for usage." >&2; exit 1; }

if [[ ! -d "${PROJECT}" ]]; then
    echo "ERROR: --project directory does not exist: ${PROJECT}" >&2; exit 1
fi
if [[ ! -d "${PROJECT}/analysis" ]]; then
    echo "ERROR: no analysis/ under --project: ${PROJECT}" >&2
    echo "       (is this the project root?)" >&2; exit 1
fi
if [[ ! -d "${BASE}" ]]; then
    echo "ERROR: --base directory does not exist: ${BASE}" >&2; exit 1
fi

# Resolve the JOBS value: blank/"auto" → 0
if [[ "${JOBS}" == "auto" || -z "${JOBS}" ]]; then JOBS=0; fi

# ── Environment ──────────────────────────────────────────────────────────────
cd "${PROJECT}"
mkdir -p logs

module load miniconda3
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "${CONDA_ENV}"

echo "================================================================"
echo "  export_summary_data.py"
echo "================================================================"
echo "  Host        : $(hostname)"
echo "  Start       : $(date)"
echo "  Project     : ${PROJECT}"
echo "  Base        : ${BASE}"
echo "  Out         : ${OUTDIR}"
echo "  Conda env   : ${CONDA_ENV}"
echo "  CPUs alloc  : ${SLURM_CPUS_PER_TASK:-<not set>}"
echo "  --jobs      : ${JOBS}  (0 = auto-detect SLURM_CPUS_PER_TASK)"
echo "================================================================"
echo

python analysis/export_summary_data.py \
    --base    "${BASE}" \
    --out-dir "${OUTDIR}" \
    --jobs    "${JOBS}"

echo
echo "Finished    : $(date)"
