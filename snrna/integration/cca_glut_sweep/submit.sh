#!/bin/bash
# Submit transfer array -> aggregate for one normalization method.
# Data prep (prep.R) must be run locally first; the following files must
# exist on the cluster before submitting:
#   snrna/integration/{SCRIPT_NAME}/{NORM}/obj_bf_common.qs
#   snrna/integration/{SCRIPT_NAME}/{NORM}/obj_chk_common.qs
#   snrna/integration/{SCRIPT_NAME}/{NORM}/sct_features.rds  (SCT only)
#
# Usage: ./slurm/submit.sh [lognorm|sct] [script-name]
# Can be run from any directory.
#
# Example:
#   ./slurm/submit.sh lognorm snrna_zaremba-chk_integration-cca_glut
#   ./slurm/submit.sh sct     snrna_zaremba-chk_integration-cca_glut_dach2

set -euo pipefail

# Absolute paths derived from script location (slurm/ -> integration/ -> project root)
SLURM_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
INTEGRATION_DIR="$(dirname "${SLURM_DIR}")"
PROJECT_ROOT="$(dirname "$(dirname "${INTEGRATION_DIR}")")"

CONDA_SH="${HOME}/miniforge3/etc/profile.d/conda.sh"
CONDA_ENV="r44"
ACTIVATE=". ${CONDA_SH} && conda activate ${CONDA_ENV}"

NORM=${1:-lognorm}
SCRIPT_NAME=${2:-snrna_zaremba-chk_integration-cca_glut}
NORM_DIR="${INTEGRATION_DIR}/${SCRIPT_NAME}/${NORM}"
LOG_DIR="${SLURM_DIR}/logs/${SCRIPT_NAME}/${NORM}"
mkdir -p "$LOG_DIR"

# Verify prep files are present before submitting
for f in obj_bf_common.qs obj_chk_common.qs; do
  if [[ ! -f "${NORM_DIR}/${f}" ]]; then
    echo "ERROR: missing prep file ${NORM_DIR}/${f}" >&2
    echo "Run prep.R locally first." >&2
    exit 1
  fi
done
if [[ "${NORM}" == "sct" && ! -f "${NORM_DIR}/sct_features.rds" ]]; then
  echo "ERROR: missing ${NORM_DIR}/sct_features.rds" >&2
  exit 1
fi

# shellcheck disable=SC1090
. "${CONDA_SH}" && conda activate "${CONDA_ENV}"

NPARAMS=$(Rscript --vanilla -e \
  "source('${SLURM_DIR}/config.R'); cat(nrow(transfer_params))")
echo "Submitting ${NPARAMS} transfer tasks for norm=${NORM} script=${SCRIPT_NAME}"

# --- 1. Transfer array ---
XFER_JID=$(sbatch --parsable \
  --job-name="xfer_${NORM}" \
  --partition=medium \
  --array="1-${NPARAMS}" \
  --time=8:00:00 \
  --mem=64G \
  --cpus-per-task=48 \
  --chdir="${PROJECT_ROOT}" \
  --output="${LOG_DIR}/transfer_%A_%a.out" \
  --wrap="${ACTIVATE} && Rscript ${SLURM_DIR}/transfer.R \
            --norm ${NORM} \
            --script-name ${SCRIPT_NAME}")
echo "Transfer array job: ${XFER_JID}"

# --- 2. Integrate array (parallel to transfer) ---
INT_JID=$(sbatch --parsable \
  --job-name="int_${NORM}" \
  --partition=medium \
  --array="1-${NPARAMS}" \
  --time=8:00:00 \
  --mem=128G \
  --cpus-per-task=48 \
  --chdir="${PROJECT_ROOT}" \
  --output="${LOG_DIR}/integrate_%A_%a.out" \
  --wrap="${ACTIVATE} && Rscript ${SLURM_DIR}/integrate.R \
            --norm ${NORM} \
            --script-name ${SCRIPT_NAME}")
echo "Integrate array job: ${INT_JID}"

# --- 3. Aggregate transfer ---
AGG_JID=$(sbatch --parsable \
  --job-name="agg_${NORM}" \
  --partition=short \
  --dependency=afterok:${XFER_JID} \
  --time=1:00:00 \
  --mem=32G \
  --cpus-per-task=1 \
  --chdir="${PROJECT_ROOT}" \
  --output="${LOG_DIR}/aggregate_%j.out" \
  --wrap="${ACTIVATE} && Rscript ${SLURM_DIR}/aggregate.R \
            --norm ${NORM} \
            --script-name ${SCRIPT_NAME}")
echo "Aggregate transfer job: ${AGG_JID}"

# --- 4. Aggregate integration ---
AGG_INT_JID=$(sbatch --parsable \
  --job-name="agg_int_${NORM}" \
  --partition=short \
  --dependency=afterok:${INT_JID} \
  --time=1:00:00 \
  --mem=32G \
  --cpus-per-task=1 \
  --chdir="${PROJECT_ROOT}" \
  --output="${LOG_DIR}/aggregate_integrate_%j.out" \
  --wrap="${ACTIVATE} && Rscript ${SLURM_DIR}/aggregate_integrate.R \
            --norm ${NORM} \
            --script-name ${SCRIPT_NAME}")
echo "Aggregate integration job: ${AGG_INT_JID}"
