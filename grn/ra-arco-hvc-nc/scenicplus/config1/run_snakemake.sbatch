#!/bin/bash
#SBATCH --job-name=scenicplus_hybrid_config1
#SBATCH --cpus-per-task=40
#SBATCH --mem=400G
#SBATCH --time=72:00:00
#SBATCH --output=logs/scenicplus_config1_%j.out
#SBATCH --error=logs/scenicplus_config1_%j.err
# EDIT for prism before submitting: --partition / --account, --mem and --time (unmeasured guesses;
# tf_to_gene and region_to_gene, the GBM steps, dominate wall time).
# The whole workflow runs inside this one allocation (--cores 40); it is not a per-rule submit.
#
# Submit from configN/Snakemake (so relative config/ and workflow/ resolve), after `mkdir -p logs`:
#   cd /private/groups/colquittlab/scenicplus/motor-pathway_multiome/motor-pathway_multiome_seurat_cellbender.0.05_preprocess_cr/ra-arco-hvc-nc_hybrid/scenicplus/config1/Snakemake && mkdir -p logs && sbatch ../run_snakemake.sbatch
set -euo pipefail
set +u   # conda's activate scripts read unset variables
source /private/home/${USER}/miniforge3/etc/profile.d/conda.sh
conda activate scenicplus
set -u
mkdir -p /private/groups/colquittlab/scenicplus/motor-pathway_multiome/motor-pathway_multiome_seurat_cellbender.0.05_preprocess_cr/ra-arco-hvc-nc_hybrid/scenicplus/config1/outs /private/groups/colquittlab/scenicplus/motor-pathway_multiome/motor-pathway_multiome_seurat_cellbender.0.05_preprocess_cr/ra-arco-hvc-nc_hybrid/scenicplus/config1/tmp
snakemake --cores 40 --rerun-incomplete --printshellcmds
