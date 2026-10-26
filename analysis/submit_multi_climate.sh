#!/bin/bash
# =============================================================================
# submit_multi_climate.sh
# SLURM array job for run_multi_climate_simulations.py on UCI HPC3 / Green Planet
#
# Each array task simulates ONE community across ALL climate scenarios.
# Array size is auto-detected from the community directory at submit time.
#
# USAGE
# -----
#   bash submit_multi_climate.sh [OPTIONS]
#
# Required:
#   --community-dir DIR    Directory of community subdirectories (must exist)
#   --output-dir    DIR    Base output directory for simulation results
#
# Climate source (mutually exclusive — same as run_multi_climate_simulations.py):
#   --climate FILE|DIR     Climate CSV or directory of CSVs inside --input-dir
#   --static-temp T [T …]  One or more constant temperatures (°C)
#   --static-psi  P [P …]  One or more constant water potentials (must match --static-temp)
#
# Branching (when runtime.txt sets branch_pulse > 0):
#   --spinup-climate FILE  Climate CSV for the shared spin-up (required with branch_pulse)
#   --reuse-spinup         Load matching saved spin-up states from <output-dir>/spinup/
#   --no-save-spinup       Do not write spin-up states to <output-dir>/spinup/
#   --full-branch-output   Write full series in every scenario folder (default: spin-up days
#                          once per community in <output-dir>/spinup/<community>_s<seed>/)
#
# Optional:
#   --src-dir    DIR       DEMENTpy source directory      (default: src)
#   --input-dir  DIR       Input CSVs / runtime.txt dir   (default: input_full_port_bac)
#   --seed       INT       Base random seed                (default: 0)
#   --conda-env  NAME      Conda environment name          (default: decadementenv)
#   --python     PATH      Python interpreter for the tasks (default: the conda env's
#                          python, found at submit time; see NOTES)
#   --tasks      LIST      Explicit array task list instead of --range, e.g. 422-424,431,502
#   --exclude    NODES     Nodes to keep tasks off, e.g. hpc3-19-03 or hpc3-19-03,hpc3-14-01
#                          (a node that fails every task within seconds drains the whole array)
#                          (task N = community N in the sorted pool), for resubmitting failures
#   --partition  NAME      SLURM partition                 (default: standard)
#   --account    NAME      SLURM account / PI group        (default: allisons_lab)
#   --cpus       INT       CPUs per task                   (default: 1)
#   --mem        STR       Memory per task                 (default: 8G)
#   --time       HH:MM:SS  Wall-clock limit per task       (default: 12:00:00)
#   --max-array  INT       Max simultaneously running tasks (default: 64)
#   --log-dir    DIR       Directory for SLURM log files   (default: logs/slurm)
#   --mail-user  ADDRESS   Email address for SLURM notifications (default: none, no email)
#   --mail-type  TYPES     Events to email about, e.g. END,FAIL or ALL (default: END,FAIL;
#                          only used with --mail-user). For an array job SLURM sends one
#                          message for the whole array, not one per task.
#   --dry-run              Print the sbatch command without submitting
#
# EXAMPLES
# --------
#   # Three static temperature treatments, free partition
#   bash submit_multi_climate.sh \
#       --community-dir community_pool/communities \
#       --output-dir    /pub/$USER/deca_dement/simulations \
#       --static-temp   10 15 20 \
#       --static-psi    -0.5 -0.5 -0.5 \
#       --conda-env     decadementenv \
#       --partition     free \
#       --mem           8G \
#       --time          12:00:00
#
#   # Seasonal heatwaves branched from a 5-year ambient spin-up
#   # (runtime.txt: pulse 6, branch_pulse 5)
#   bash submit_multi_climate.sh \
#       --community-dir  community_pool_port/communities \
#       --output-dir     /pub/$USER/deca_dement/sim_heatwave \
#       --climate        input_full_port_bac/climate_data \
#       --spinup-climate input_full_port_bac/climate_data/climate_ambient_temp.csv \
#       --time           06:00:00
#
#   # Climate directory
#   bash submit_multi_climate.sh \
#       --community-dir community_pool/communities \
#       --output-dir    /pub/$USER/deca_dement/simulations \
#       --climate       climate_scenarios \
#       --conda-env     decadementenv \
#       --partition     standard \
#       --account       allisons_lab \
#       --time          24:00:00
#
# NOTES
# -----
#   • Run this script from the project root (same directory as src/, input_full_port_bac/, etc.)
#   • Communities must have been assembled from the same input folder (assemble_communities.py
#     --input-dir input_full_port_bac): the per-taxon trait matrices are stored in the community
#     directories; everything else (Enz_loss_rate, Enzymes_per_sub, mineral inputs, energy yields,
#     runtime.txt) is read from --input-dir at run time.
#   • Each task runs every scenario in sequence: at 100x100 a 5-pulse run takes about 2-3 h per
#     scenario on one core (4 GB), at 40x40 about 20-30 min (1 GB).
#   • With branch_pulse B in runtime.txt each task simulates pulses 1..B once with the spin-up
#     climate, then pulses B+1..pulse for every scenario from that state, so the wall time is about
#     (B + n_scenarios x (pulse - B)) pulses instead of n_scenarios x pulse. The spin-up state is
#     saved to <output-dir>/spinup/ (about 0.5 GB per community at 40x40 after 5 pulses, several
#     times that at 100x100; --no-save-spinup skips it). Measured on a Mac at 40x40, 6 pulses with
#     branch_pulse 5 and 3 scenarios: about 13 min per community (straight: about 28 min).
#     A branched task holds the model plus a serialized copy of its state in memory: at 100x100
#     allow about 16G.
#   • Branched output (csv) is split: <output-dir>/spinup/<community>_s<seed>/ holds the days before
#     the branch and the complete Initialization folder once; each scenario folder holds the branch
#     days only (about 150 MB per community at 40x40 with 3 scenarios, instead of about 340 MB).
#     Read full series with analysis/branch_output.py: read_table(run_dir, 'MicrobesSeries').
#   • Seeds are --seed + the community's position in the sorted pool, so a task's seed (and its
#     output) is the same as in a run over the whole pool.
#   • Each task handles exactly one community; all scenarios for that community
#     run sequentially within the task.  Wall-time limit must cover all scenarios.
#   • Array throttle (--max-array) protects scheduler and filesystem from floods.
#   • Logs go to <log-dir>/deca_%A_%a.{out,err}  where %A = job ID, %a = task ID.
#   • Tasks do not run conda: the conda env's python is located once at submit time (on the login
#     node, where numpy/pandas imports are checked) and every task calls it directly. Hundreds of
#     simultaneous 'conda activate' calls crashed conda at task start (2026-10-26).
#   • To resubmit every failed or cancelled task of an array:
#       LIST=$(sacct -j <jobID> -X -n -P --format=JobID,State | awk -F'|' '$2 !~ /COMPLETED|RUNNING|PENDING/ {split($1,a,"_"); gsub(/[][]|%.*/,"",a[2]); print a[2]}' | paste -sd, -)
#       bash submit_multi_climate.sh ... --tasks "$LIST"
#   • If a task fails, resubmit only that community with:
#       bash submit_multi_climate.sh ... --range FAILED_IDX FAILED_IDX
# =============================================================================

set -euo pipefail

# ── Defaults ─────────────────────────────────────────────────────────────────
SRC_DIR="src"
INPUT_DIR="input_full_port_bac"
COMMUNITY_DIR=""
OUTPUT_DIR=""
CLIMATE=""
SPINUP_CLIMATE=""
REUSE_SPINUP=0
NO_SAVE_SPINUP=0
FULL_BRANCH_OUTPUT=0
STATIC_TEMP=""         # space-separated values collected during arg parse
STATIC_PSI=""
SEED=0
CONDA_ENV="decadementenv"
PARTITION="standard"
ACCOUNT="allisons_lab"
CPUS=1
MEM="8G"
TIME="12:00:00"
MAX_ARRAY=64
LOG_DIR="logs/slurm"
MAIL_USER=""
MAIL_TYPE="END,FAIL"
DRY_RUN=0
RANGE_START=""
TASKS=""
PYTHON_BIN=""
EXCLUDE=""
RANGE_END=""

# ── Argument parsing ──────────────────────────────────────────────────────────
# We collect --static-temp and --static-psi as multi-value lists manually
# because bash doesn't support nargs='+' natively.
COLLECTING_TEMP=0
COLLECTING_PSI=0

while [[ $# -gt 0 ]]; do
    case "$1" in
        --src-dir)       SRC_DIR="$2";       shift 2; COLLECTING_TEMP=0; COLLECTING_PSI=0 ;;
        --input-dir)     INPUT_DIR="$2";     shift 2; COLLECTING_TEMP=0; COLLECTING_PSI=0 ;;
        --community-dir) COMMUNITY_DIR="$2"; shift 2; COLLECTING_TEMP=0; COLLECTING_PSI=0 ;;
        --output-dir)    OUTPUT_DIR="$2";    shift 2; COLLECTING_TEMP=0; COLLECTING_PSI=0 ;;
        --climate)       CLIMATE="$2";       shift 2; COLLECTING_TEMP=0; COLLECTING_PSI=0 ;;
        --spinup-climate) SPINUP_CLIMATE="$2"; shift 2; COLLECTING_TEMP=0; COLLECTING_PSI=0 ;;
        --reuse-spinup)  REUSE_SPINUP=1;     shift;   COLLECTING_TEMP=0; COLLECTING_PSI=0 ;;
        --no-save-spinup) NO_SAVE_SPINUP=1;  shift;   COLLECTING_TEMP=0; COLLECTING_PSI=0 ;;
        --full-branch-output) FULL_BRANCH_OUTPUT=1; shift; COLLECTING_TEMP=0; COLLECTING_PSI=0 ;;
        --seed)          SEED="$2";          shift 2; COLLECTING_TEMP=0; COLLECTING_PSI=0 ;;
        --conda-env)     CONDA_ENV="$2";     shift 2; COLLECTING_TEMP=0; COLLECTING_PSI=0 ;;
        --python)        PYTHON_BIN="$2";    shift 2; COLLECTING_TEMP=0; COLLECTING_PSI=0 ;;
        --tasks)         TASKS="$2";         shift 2; COLLECTING_TEMP=0; COLLECTING_PSI=0 ;;
        --exclude)       EXCLUDE="$2";       shift 2; COLLECTING_TEMP=0; COLLECTING_PSI=0 ;;
        --partition)     PARTITION="$2";     shift 2; COLLECTING_TEMP=0; COLLECTING_PSI=0 ;;
        --account)       ACCOUNT="$2";       shift 2; COLLECTING_TEMP=0; COLLECTING_PSI=0 ;;
        --cpus)          CPUS="$2";          shift 2; COLLECTING_TEMP=0; COLLECTING_PSI=0 ;;
        --mem)           MEM="$2";           shift 2; COLLECTING_TEMP=0; COLLECTING_PSI=0 ;;
        --time)          TIME="$2";          shift 2; COLLECTING_TEMP=0; COLLECTING_PSI=0 ;;
        --max-array)     MAX_ARRAY="$2";     shift 2; COLLECTING_TEMP=0; COLLECTING_PSI=0 ;;
        --log-dir)       LOG_DIR="$2";       shift 2; COLLECTING_TEMP=0; COLLECTING_PSI=0 ;;
        --mail-user)     MAIL_USER="$2";     shift 2; COLLECTING_TEMP=0; COLLECTING_PSI=0 ;;
        --mail-type)     MAIL_TYPE="$2";     shift 2; COLLECTING_TEMP=0; COLLECTING_PSI=0 ;;
        --range)         RANGE_START="$2"; RANGE_END="$3"; shift 3;
                         COLLECTING_TEMP=0; COLLECTING_PSI=0 ;;
        --dry-run)       DRY_RUN=1;          shift;   COLLECTING_TEMP=0; COLLECTING_PSI=0 ;;
        --static-temp)
            COLLECTING_TEMP=1; COLLECTING_PSI=0; shift
            ;;
        --static-psi)
            COLLECTING_PSI=1; COLLECTING_TEMP=0; shift
            ;;
        -[0-9]*|--*=[0-9-]*)
            # Negative number or --flag=value form — treat as a bare value
            if [[ "${COLLECTING_TEMP}" -eq 1 ]]; then
                STATIC_TEMP="${STATIC_TEMP} $1"; shift
            elif [[ "${COLLECTING_PSI}" -eq 1 ]]; then
                STATIC_PSI="${STATIC_PSI} $1"; shift
            else
                echo "ERROR: unexpected argument '$1'" >&2
                exit 1
            fi
            ;;
        -*)
            # Unknown flag — stop collecting temp/psi
            COLLECTING_TEMP=0; COLLECTING_PSI=0
            echo "WARNING: unknown flag '$1' (ignored)" >&2
            shift
            ;;
        *)
            # Bare value — could be a multi-value for --static-temp or --static-psi
            if [[ "${COLLECTING_TEMP}" -eq 1 ]]; then
                STATIC_TEMP="${STATIC_TEMP} $1"; shift
            elif [[ "${COLLECTING_PSI}" -eq 1 ]]; then
                STATIC_PSI="${STATIC_PSI} $1"; shift
            else
                echo "ERROR: unexpected positional argument '$1'" >&2
                exit 1
            fi
            ;;
    esac
done

# Trim leading spaces from multi-value lists
STATIC_TEMP="${STATIC_TEMP# }"
STATIC_PSI="${STATIC_PSI# }"

# ── Validate required args ────────────────────────────────────────────────────
if [[ -z "${COMMUNITY_DIR}" ]]; then
    echo "ERROR: --community-dir is required." >&2; exit 1
fi
if [[ -z "${OUTPUT_DIR}" ]]; then
    echo "ERROR: --output-dir is required." >&2; exit 1
fi
if [[ -n "${STATIC_TEMP}" && -z "${STATIC_PSI}" ]] || \
   [[ -z "${STATIC_TEMP}" && -n "${STATIC_PSI}" ]]; then
    echo "ERROR: --static-temp and --static-psi must be supplied together." >&2; exit 1
fi
if [[ -n "${STATIC_TEMP}" && -n "${CLIMATE}" ]]; then
    echo "ERROR: --static-temp/--static-psi cannot be combined with --climate." >&2; exit 1
fi

# Resolve to absolute paths from CWD (where the user runs this script)
PROJECT_ROOT="$(pwd)"
SRC_DIR="$(realpath "${SRC_DIR}")"
INPUT_DIR="$(realpath "${INPUT_DIR}")"
COMMUNITY_DIR="$(realpath "${COMMUNITY_DIR}")"
# OUTPUT_DIR may not exist yet — don't realpath it
OUTPUT_DIR="${OUTPUT_DIR}"   # passed as-is; Python script creates it

# ── Count communities to size the array ──────────────────────────────────────
# Count subdirs that contain EnzGenes.csv (same criterion as get_community_dirs)
N_COMMUNITIES=$(find "${COMMUNITY_DIR}" -maxdepth 2 -name "EnzGenes.csv" \
                | awk -F/ '{print $(NF-1)}' | sort -u | wc -l)

if [[ "${N_COMMUNITIES}" -eq 0 ]]; then
    echo "ERROR: no community directories (containing EnzGenes.csv) found in:" >&2
    echo "       ${COMMUNITY_DIR}" >&2
    exit 1
fi

# Apply --tasks or --range if given, otherwise use full 1..N span
if [[ -n "${TASKS}" && -n "${RANGE_START}" ]]; then
    echo "ERROR: --tasks and --range cannot be combined." >&2; exit 1
fi
if [[ -n "${TASKS}" ]]; then
    if [[ ! "${TASKS}" =~ ^[0-9]+(-[0-9]+)?(,[0-9]+(-[0-9]+)?)*$ ]]; then
        echo "ERROR: --tasks must look like 422-424,431,502 (got '${TASKS}')." >&2; exit 1
    fi
    ARRAY_START="${TASKS%%[-,]*}"
    ARRAY_SPEC="${TASKS}%${MAX_ARRAY}"
else
    if [[ -n "${RANGE_START}" && -n "${RANGE_END}" ]]; then
        ARRAY_START="${RANGE_START}"
        ARRAY_END="${RANGE_END}"
    else
        ARRAY_START=1
        ARRAY_END="${N_COMMUNITIES}"
    fi
    ARRAY_SPEC="${ARRAY_START}-${ARRAY_END}%${MAX_ARRAY}"
fi

# ── Python interpreter for the tasks (no conda calls inside the jobs) ─────────
if [[ -z "${PYTHON_BIN}" ]]; then
    CAND="${HOME}/.conda/envs/${CONDA_ENV}/bin/python"
    if [[ -x "${CAND}" ]]; then
        PYTHON_BIN="${CAND}"
    elif command -v conda >/dev/null 2>&1; then
        ENV_PREFIX=$(conda env list 2>/dev/null | awk -v e="${CONDA_ENV}" '$1==e {print $NF}')
        [[ -n "${ENV_PREFIX}" ]] && PYTHON_BIN="${ENV_PREFIX}/bin/python"
    fi
fi
if [[ -z "${PYTHON_BIN}" || ! -x "${PYTHON_BIN}" ]]; then
    echo "ERROR: could not find the python of conda env '${CONDA_ENV}'." >&2
    echo "       Looked for ${HOME}/.conda/envs/${CONDA_ENV}/bin/python and in 'conda env list'." >&2
    echo "       Pass it with --python /path/to/env/bin/python." >&2
    exit 1
fi
if ! "${PYTHON_BIN}" -c "import numpy, pandas" >/dev/null 2>&1; then
    echo "ERROR: ${PYTHON_BIN} cannot import numpy and pandas." >&2; exit 1
fi

# ── Build the python command args that each task will run ─────────────────────
PYTHON_ARGS=(
    "analysis/run_multi_climate_simulations.py"
    "--src-dir"       "${SRC_DIR}"
    "--input-dir"     "${INPUT_DIR}"
    "--community-dir" "${COMMUNITY_DIR}"
    "--output-dir"    "${OUTPUT_DIR}"
    "--seed"          "${SEED}"
    "--range"         "\${SLURM_ARRAY_TASK_ID}" "\${SLURM_ARRAY_TASK_ID}"
    "--quiet"
)

if [[ -n "${CLIMATE}" ]]; then
    PYTHON_ARGS+=("--climate" "${CLIMATE}")
fi

if [[ -n "${SPINUP_CLIMATE}" ]]; then
    PYTHON_ARGS+=("--spinup-climate" "${SPINUP_CLIMATE}")
fi

if [[ "${REUSE_SPINUP}" -eq 1 ]]; then
    PYTHON_ARGS+=("--reuse-spinup")
fi

if [[ "${NO_SAVE_SPINUP}" -eq 1 ]]; then
    PYTHON_ARGS+=("--no-save-spinup")
fi

if [[ "${FULL_BRANCH_OUTPUT}" -eq 1 ]]; then
    PYTHON_ARGS+=("--full-branch-output")
fi

if [[ -n "${STATIC_TEMP}" ]]; then
    PYTHON_ARGS+=("--static-temp" ${STATIC_TEMP})
    PYTHON_ARGS+=("--static-psi"  ${STATIC_PSI})
fi

# Convert array to a single string for embedding in the heredoc
PYTHON_CMD="${PYTHON_BIN} ${PYTHON_ARGS[*]}"

# ── Build optional SBATCH directives ─────────────────────────────────────────
ACCOUNT_LINE=""
if [[ -n "${ACCOUNT}" ]]; then
    ACCOUNT_LINE="#SBATCH --account=${ACCOUNT}"
fi

EXCLUDE_LINE=""
if [[ -n "${EXCLUDE}" ]]; then
    EXCLUDE_LINE="#SBATCH --exclude=${EXCLUDE}"
fi

MAIL_LINES=""
if [[ -n "${MAIL_USER}" ]]; then
    MAIL_LINES="#SBATCH --mail-user=${MAIL_USER}
#SBATCH --mail-type=${MAIL_TYPE}"
fi

# ── Create log directory ──────────────────────────────────────────────────────
mkdir -p "${LOG_DIR}"

# ── Print submit summary ──────────────────────────────────────────────────────
echo "================================================================"
echo "  DeCA-DEMENT multi-climate SLURM array job"
echo "================================================================"
echo "  Project root  : ${PROJECT_ROOT}"
echo "  Community dir : ${COMMUNITY_DIR}"
echo "  Output dir    : ${OUTPUT_DIR}"
echo "  Communities   : ${N_COMMUNITIES}"
echo "  Array range   : ${ARRAY_SPEC}  (${MAX_ARRAY} concurrent max)"
echo "  Partition     : ${PARTITION}"
echo "  CPUs/task     : ${CPUS}"
echo "  Mem/task      : ${MEM}"
echo "  Wall time     : ${TIME}"
echo "  Conda env     : ${CONDA_ENV}"
echo "  Python        : ${PYTHON_BIN}"
echo "  Log dir       : ${LOG_DIR}"
if [[ -n "${ACCOUNT}" ]];   then echo "  Account       : ${ACCOUNT}"; fi
if [[ -n "${MAIL_USER}" ]]; then echo "  Email         : ${MAIL_USER} (${MAIL_TYPE})"; fi
if [[ -n "${EXCLUDE}" ]]; then echo "  Exclude nodes : ${EXCLUDE}"; fi
if [[ -n "${CLIMATE}" ]];   then echo "  Climate       : ${CLIMATE}"; fi
if [[ -n "${SPINUP_CLIMATE}" ]]; then echo "  Spin-up clim. : ${SPINUP_CLIMATE}"; fi
if [[ "${REUSE_SPINUP}" -eq 1 ]]; then echo "  Reuse spin-up : yes"; fi
if [[ "${NO_SAVE_SPINUP}" -eq 1 ]]; then echo "  Save spin-up  : no"; fi
if [[ -n "${STATIC_TEMP}" ]]; then
    echo "  Static temp   : ${STATIC_TEMP}"
    echo "  Static psi    : ${STATIC_PSI}"
fi
echo "================================================================"

if [[ "${DRY_RUN}" -eq 1 ]]; then
    echo "  DRY RUN — job script shown below, not submitted."
    echo ""
fi

# ── Generate and submit the job script ───────────────────────────────────────
JOB_SCRIPT=$(cat <<SLURM_SCRIPT
#!/bin/bash
#SBATCH --job-name=deca_dement
#SBATCH --array=${ARRAY_SPEC}
#SBATCH --partition=${PARTITION}
${ACCOUNT_LINE}
${MAIL_LINES}
${EXCLUDE_LINE}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=${CPUS}
#SBATCH --mem=${MEM}
#SBATCH --time=${TIME}
#SBATCH --output=${LOG_DIR}/deca_%A_%a.out
#SBATCH --error=${LOG_DIR}/deca_%A_%a.err

# ── Environment ------------------------------------------------------------------
# The conda env's python is called directly (located at submit time); no conda
# commands run here, because hundreds of simultaneous 'conda activate' calls
# crashed conda at task start.
if ! "${PYTHON_BIN}" -c "import numpy, pandas"; then
    echo "ERROR: ${PYTHON_BIN} cannot import numpy/pandas on \$(hostname)." >&2
    exit 1
fi

# ── Working directory ------------------------------------------------------------
# Run from the project root so relative paths (src/, input/, analysis/) resolve.
cd "${PROJECT_ROOT}"

# ── Diagnostic info --------------------------------------------------------------
echo "======================================================================"
echo "  Job     : \${SLURM_JOB_ID}   Task: \${SLURM_ARRAY_TASK_ID}"
echo "  Node    : \$(hostname)"
echo "  Started : \$(date)"
echo "  Community index: \${SLURM_ARRAY_TASK_ID} / ${N_COMMUNITIES}"
echo "======================================================================"

# ── Run -------------------------------------------------------------------------
${PYTHON_CMD}

EXIT_CODE=\$?

echo ""
echo "======================================================================"
echo "  Finished : \$(date)"
echo "  Exit code: \${EXIT_CODE}"
echo "======================================================================"

exit \${EXIT_CODE}
SLURM_SCRIPT
)

if [[ "${DRY_RUN}" -eq 1 ]]; then
    echo "${JOB_SCRIPT}"
    echo ""
    echo "  (dry run — use without --dry-run to submit)"
else
    JOB_ID=$(echo "${JOB_SCRIPT}" | sbatch --parsable)
    echo ""
    echo "  Submitted job array: ${JOB_ID}"
    echo "  Monitor with:"
    echo "    squeue -u \$USER"
    echo "    squeue -j ${JOB_ID}"
    echo "  Cancel with:"
    echo "    scancel ${JOB_ID}"
    echo "  Check a task log:"
    echo "    tail -f ${LOG_DIR}/deca_${JOB_ID}_${ARRAY_START}.out"
fi
