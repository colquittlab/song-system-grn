#!/bin/bash
#SBATCH --job-name=scenicplus_hybrid_config24
#SBATCH --partition=medium
#SBATCH --cpus-per-task=40
#SBATCH --mem=600G
#SBATCH --time=08:00:00
#SBATCH --output=logs/scenicplus_config24_%j.out
#SBATCH --error=logs/scenicplus_config24_%j.err
# No --account: prism doesn't need one here. --mem/--time are set by hand, not measured; tf_to_gene and
# region_to_gene (the GBM steps) dominate wall time, so check `sacct` after the first run.
# The whole workflow runs inside this one allocation (--cores 40); it is not a per-rule submit.
#
# Submit from configN/Snakemake (so relative config/ and workflow/ resolve), after `mkdir -p logs`:
#   cd /private/groups/colquittlab/scenicplus/motor-pathway_multiome/motor-pathway_multiome_seurat_cellbender.0.05_preprocess_cr/ra-arco-hvc-nc_hybrid/scenicplus/config24/Snakemake && mkdir -p logs && sbatch ../run_snakemake.sbatch
set -euo pipefail
set +u   # conda's activate scripts read unset variables
source /private/home/${USER}/miniforge3/etc/profile.d/conda.sh
conda activate scenicplus
set -u
# These come from config1 via prepare_shared_upstream.sh; without them Snakemake would silently redo the 3.5 h
# region-to-gene fit and the 35 GB object, per config.
for f in ACC_GEX.h5mu search_space.tsv region_to_gene_adj.tsv; do
    test -s /private/groups/colquittlab/scenicplus/motor-pathway_multiome/motor-pathway_multiome_seurat_cellbender.0.05_preprocess_cr/ra-arco-hvc-nc_hybrid/scenicplus/shared/outs/$f || { echo "missing shared input /private/groups/colquittlab/scenicplus/motor-pathway_multiome/motor-pathway_multiome_seurat_cellbender.0.05_preprocess_cr/ra-arco-hvc-nc_hybrid/scenicplus/shared/outs/$f -- run prepare_shared_upstream.sh first" >&2; exit 1; }
done
mkdir -p /private/groups/colquittlab/scenicplus/motor-pathway_multiome/motor-pathway_multiome_seurat_cellbender.0.05_preprocess_cr/ra-arco-hvc-nc_hybrid/scenicplus/config24/outs /private/groups/colquittlab/scenicplus/motor-pathway_multiome/motor-pathway_multiome_seurat_cellbender.0.05_preprocess_cr/ra-arco-hvc-nc_hybrid/scenicplus/config24/tmp
snakemake --cores 40 --rerun-incomplete --printshellcmds
