#!/bin/bash
# All-cell cross-species integrations (finch x chicken = gg, finch x mouse = yao) with
# SCTransform normalization, matching the glut-only pipeline. OMA 1:1 orthologs for GSI and
# label transfer; SAMap uses its own BLAST maps. Outputs sit beside the log-normalized ones:
#   GSI    composite_scoring/results/gsi_corr_<pair>_adult_hybrid_SCT.csv
#   CCA    rpca_sweep/results/cca/<pair>_adult_hybrid_cca_SCT_*_ka5_kfNA_d30_*.csv
#   SAMap  datasets/samap_bf-adult_gg-adult/results_hybrid_SCT, datasets/samap_bf-mature_yao/results_adult_hybrid_SCT
# Run from anywhere; submits SLURM jobs with dependencies and prints their IDs.
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
LIVE=/private/groups/colquittlab/song-system-grn/snrna/integration
TK=/private/groups/colquittlab/finch-integration-toolkit
PY=$HOME/miniforge3/envs/saturn/bin/python
SPY=$HOME/miniforge3/envs/SAMap/bin/python
RS=$HOME/miniforge3/envs/r44/bin/Rscript
SCT=$LIVE/composite_scoring/data/sct
LOG=$LIVE/logs
sb() { sbatch --parsable --output=$LOG/%x_%j.out --error=$LOG/%x_%j.err "$@"; }
dep() { local IFS=:; echo "afterok:$*"; }

# ---- GSI on SCT-corrected counts (subsampled, ortholog genes) ----
for pair in gg yao; do
  E=$(sb -J gsisct_export_${pair} -p medium -c 8 --mem=300G -t 6:00:00 --wrap "$PY $HERE/allcell_sct.py gsi-export --pair $pair")
  S1=$(sb -J sct_gsi_${pair}_finch -p long -c 8 --mem=300G -t 24:00:00 --dependency=afterok:$E --wrap "$RS $TK/R/sct_transform_export.R $SCT/gsi_${pair}/finch $SCT/gsi_${pair}/finch_sct")
  S2=$(sb -J sct_gsi_${pair}_ref   -p long -c 8 --mem=500G -t 48:00:00 --dependency=afterok:$E --wrap "$RS $TK/R/sct_transform_export.R $SCT/gsi_${pair}/ref $SCT/gsi_${pair}/ref_sct")
  sb -J gsisct_run_${pair} -p medium -c 8 --mem=300G -t 6:00:00 --dependency=$(dep $S1 $S2) --wrap "$PY $HERE/allcell_sct.py gsi-run --pair $pair"
done

# ---- CCA label transfer with SCT (same anchors settings: cca, ka5, kfNA, d30) ----
for t in gg_adult_hybrid yao_adult_hybrid; do
  BASE_DIR=$LIVE/rpca_sweep TAG=$t NORM=SCT sb -J cca_sct_${t%%_*} --export=ALL $TK/run_cca.sbatch
done

# ---- SAMap on SCT-corrected counts (complete transcriptomes, SAM normalization off) ----
declare -A ASM
for sp in bf gg mm; do
  X=$(sb -J samapsct_export_$sp -p medium -c 4 --mem=300G -t 6:00:00 --wrap "$PY $HERE/allcell_sct.py samap-export --species $sp")
  S=$(sb -J sct_samap_$sp -p long -c 8 --mem=500G -t 48:00:00 --dependency=afterok:$X --wrap "$RS $TK/R/sct_transform_export.R $SCT/samap/$sp $SCT/samap/${sp}_sct")
  ASM[$sp]=$(sb -J samapsct_assemble_$sp -p medium -c 4 --mem=300G -t 6:00:00 --dependency=afterok:$S --wrap "$PY $HERE/allcell_sct.py samap-assemble --species $sp")
done
RUN=$HERE/gg_glutonly_hybrid/samap/run_samap_sct.py
sb -J samapsct_gg -p long -c 32 --mem=500G -t 48:00:00 --dependency=$(dep ${ASM[bf]} ${ASM[gg]}) \
   --wrap "export OPENBLAS_NUM_THREADS=16 OMP_NUM_THREADS=16; $SPY $RUN --ncpus 32 --finch_h5ad $SCT/samap/bf_sct.h5ad --chicken_h5ad $SCT/samap/gg_sct.h5ad --ref_id gg --ref_name chicken --maps $LIVE/datasets/samap_bf-adult_gg-adult/maps/ --out_dir $LIVE/datasets/samap_bf-adult_gg-adult/results_hybrid_SCT"
sb -J samapsct_yao -p long -c 32 --mem=500G -t 48:00:00 --dependency=$(dep ${ASM[bf]} ${ASM[mm]}) \
   --wrap "export OPENBLAS_NUM_THREADS=16 OMP_NUM_THREADS=16; $SPY $RUN --ncpus 32 --finch_h5ad $SCT/samap/bf_sct.h5ad --chicken_h5ad $SCT/samap/mm_sct.h5ad --ref_id mm --ref_name mouse --maps /private/groups/colquittlab/saturn/samap_bf-dev_lamanno/maps_full/ --out_dir $LIVE/datasets/samap_bf-mature_yao/results_adult_hybrid_SCT"
echo submitted
