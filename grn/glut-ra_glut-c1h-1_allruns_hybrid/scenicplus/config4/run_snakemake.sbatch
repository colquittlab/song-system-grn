#!/bin/bash
#SBATCH --job-name=scenicplus_ra_c1h1_all_hybrid_config4
#SBATCH --partition=medium
#SBATCH --cpus-per-task=40
#SBATCH --mem=600G
#SBATCH --time=12:00:00
#SBATCH --output=logs/scenicplus_config4_%j.out
#SBATCH --error=logs/scenicplus_config4_%j.err
# No --account: prism doesn't need one here. --mem/--time are set by hand, not measured; tf_to_gene and
# region_to_gene (the GBM steps) dominate wall time, so check `sacct` after the first run.
# The whole workflow runs inside this one allocation (--cores 40); it is not a per-rule submit.
#
# Submit from configN/Snakemake (so relative config/ and workflow/ resolve), after `mkdir -p logs`:
#   cd /private/groups/colquittlab/scenicplus/motor-pathway_multiome/motor-pathway_multiome_seurat_cellbender.0.05_preprocess_cr/glut-ra_glut-c1h-1_allruns_hybrid/scenicplus/config4/Snakemake && mkdir -p logs && sbatch ../run_snakemake.sbatch
set -euo pipefail
set +u   # conda's activate scripts read unset variables
source /private/home/${USER}/miniforge3/etc/profile.d/conda.sh
conda activate scenicplus
set -u
# ACC_GEX.h5mu comes from config1 via prepare_shared_upstream.sh (the search space and region-to-gene fit are this config's own);
# without it Snakemake would silently rebuild the accessibility + expression object.
test -s /private/groups/colquittlab/scenicplus/motor-pathway_multiome/motor-pathway_multiome_seurat_cellbender.0.05_preprocess_cr/glut-ra_glut-c1h-1_allruns_hybrid/scenicplus/shared/outs/ACC_GEX.h5mu || { echo "missing shared input /private/groups/colquittlab/scenicplus/motor-pathway_multiome/motor-pathway_multiome_seurat_cellbender.0.05_preprocess_cr/glut-ra_glut-c1h-1_allruns_hybrid/scenicplus/shared/outs/ACC_GEX.h5mu -- run prepare_shared_upstream.sh first" >&2; exit 1; }
mkdir -p /private/groups/colquittlab/scenicplus/motor-pathway_multiome/motor-pathway_multiome_seurat_cellbender.0.05_preprocess_cr/glut-ra_glut-c1h-1_allruns_hybrid/scenicplus/config4/outs /private/groups/colquittlab/scenicplus/motor-pathway_multiome/motor-pathway_multiome_seurat_cellbender.0.05_preprocess_cr/glut-ra_glut-c1h-1_allruns_hybrid/scenicplus/config4/tmp
snakemake --cores 40 --rerun-incomplete --printshellcmds
