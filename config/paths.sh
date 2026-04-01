#!/usr/bin/env bash
# config/paths.sh
#
# Central configuration for all external data paths used by shell scripts.
# Source this file at the top of each shell script:
#   source "$(dirname "$0")/../../config/paths.sh"   # adjust depth as needed
#
# Edit to match your local environment. See data/README.md for download
# instructions for reference files.

# ---------------------------------------------------------------------------
# Zebra finch (lonStrDom2 / GCF_005870125.1) reference assembly
# Source: NCBI RefSeq GCF_005870125.1
#   https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_005870125.1/
# ---------------------------------------------------------------------------

GENOME_FA_LONSTR="/mnt/nest/assembly/lonStrDom2/ncbi/GCF_005870125.1_lonStrDom2_genomic_ucsc_only.fna"
CHROM_SIZES_LONSTR="/mnt/nest/assembly/lonStrDom2/ucsc/chrom.sizes.ucsc"

# ---------------------------------------------------------------------------
# cisTarget motif databases
# Source: SCENIC+ v10 — https://resources.aertslab.org/cistarget/
# ---------------------------------------------------------------------------

CISTARGET_SINGLETONS_DIR="/home/brad/nest/cistarget/v10nr_clust_public/singletons"

# Path to the create_cisTarget_databases script repository
# Source: https://github.com/aertslab/create_cisTarget_databases
CREATE_CISTARGET_SCRIPT_DIR="/home/brad/repos/create_cisTarget_databases"
