#!/bin/bash
# SAMap on SCTransform-normalized input for the plusSATB2_noPre glut-only finch x chicken
# comparison (finch = hybrid labels + Glut-SATB2-1; chicken = Ex_* minus Ex_Pre_*).
# One job, four steps, each skipped if its output already exists:
#   1. export raw counts for exactly these cell sets      (samap_sct_export_inputs.py)
#   2. SCTransform each species independently, r44 env    (toolkit R/sct_transform_export.R)
#   3. reassemble SCT h5ads with the original obs         (samap_sct_assemble_h5ad.py)
#   4. SAMap with its own normalization disabled          (run_samap_sct.py)
#SBATCH --job-name=samap_sct_plusSATB2_noPre
#SBATCH --partition=long
#SBATCH --cpus-per-task=32
#SBATCH --mem=500G
#SBATCH --time=24:00:00
#SBATCH --output=/private/groups/colquittlab/song-system-grn/snrna/integration/gg_glutonly_hybrid/samap/logs/samap_sct_plusSATB2_noPre_%j.out
#SBATCH --error=/private/groups/colquittlab/song-system-grn/snrna/integration/gg_glutonly_hybrid/samap/logs/samap_sct_plusSATB2_noPre_%j.err
set -euo pipefail
export OPENBLAS_NUM_THREADS=16 OMP_NUM_THREADS=16 MKL_NUM_THREADS=16 NUMEXPR_NUM_THREADS=16
SAMAP_DIR=/private/groups/colquittlab/song-system-grn/snrna/integration/gg_glutonly_hybrid/samap
DATA=$SAMAP_DIR/data
PY=$HOME/miniforge3/envs/SAMap/bin/python
RS=$HOME/miniforge3/envs/r44/bin/Rscript
SCT_R=/private/groups/colquittlab/finch-integration-toolkit/R/sct_transform_export.R
cd "$SAMAP_DIR"
echo "host=$(hostname) start=$(date)"

echo "== step 1: export raw counts"
[ -s "$DATA/chicken_noPre_matrix.mtx" ] && [ -s "$DATA/finch_plusSATB2_matrix.mtx" ] \
  || $PY samap_sct_export_inputs.py --subset plusSATB2_noPre

for sp in finch_plusSATB2 chicken_noPre; do
  echo "== step 2: SCTransform $sp  $(date)"
  [ -s "$DATA/${sp}_sct_matrix.mtx" ] || $RS "$SCT_R" "$DATA/$sp" "$DATA/${sp}_sct"
done

echo "== step 3: assemble h5ads  $(date)"
[ -s "$DATA/finch_plusSATB2_sct.h5ad" ] && [ -s "$DATA/chicken_noPre_sct.h5ad" ] \
  || $PY samap_sct_assemble_h5ad.py --subset plusSATB2_noPre

echo "== step 4: SAMap  $(date)"
$PY run_samap_sct.py --ncpus 32 \
    --finch_h5ad "$DATA/finch_plusSATB2_sct.h5ad" \
    --chicken_h5ad "$DATA/chicken_noPre_sct.h5ad" \
    --out_dir "$SAMAP_DIR/results/SCT_plusSATB2_noPre"
echo "end=$(date)"
