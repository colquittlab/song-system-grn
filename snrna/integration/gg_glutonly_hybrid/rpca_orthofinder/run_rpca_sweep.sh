#!/bin/bash
# RPCA label-transfer sweep on the OrthoFinder-ortholog inputs (prepare_inputs.R), with the
# same settings as the original RBH-ortholog run in ../rpca/results/cca: SCTransform per
# species, reduction rpca, dims 1:40, k.filter 200 (nominal; Seurat ignores k.filter under
# SCT), max.features 200, k.score 30, k.anchor 2-100. Runs the unchanged
# finch-integration-toolkit/R/cca_transfer.R; writes only to results/cca/ here.
#
#   bash snrna/integration/gg_glutonly_hybrid/rpca_orthofinder/run_rpca_sweep.sh
set -euo pipefail
cd "$(dirname "$0")"
TOOLKIT=${FINCH_INTEGRATION_TOOLKIT:-/ssd/brad/rstudio/multiome/finch-integration-toolkit}
RENV_LIB=${RENV_LIB:-/ssd/brad/rstudio/multiome/song-system-grn/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu}
TAG=gg_glutonly_hybrid_plusSATB2_noPre_orthofinder
mkdir -p results/cca logs
export R_LIBS="$RENV_LIB" OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4
for KA in 2 5 10 20 30 50 75 100; do
  Rscript --no-init-file "$TOOLKIT/R/cca_transfer.R" "$TAG" data/seurat results/cca 40 "$KA" rpca 200 SCT 200 30 \
    > "logs/rpca_ka${KA}.log" 2>&1 &
done
wait
echo "done: $(ls results/cca/*_matrix.csv | wc -l) matrices"
