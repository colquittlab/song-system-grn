#!/usr/bin/env bash
# Build a cisTarget database on the consensus regions of the hybrid-labeled run (ra-arco-hvc-nc_hybrid).
#
# The existing database (ra-arco-hvc-nc_seurat-clustering, 517,176 regions) was built on the Nov 2024 consensus and is left
# untouched. The consensus peaks were regenerated from the hybrid labels on 2026-10-05 (499,348 regions; only 25,023 are
# identical to the old ones), so SCENIC+ runs on the current cisTopic object need a database on these regions: they are
# otherwise matched to the old regions by the 40 % overlap rule and 14 % of them get no motif scores.
#
# Same recipe as create_cistarget_db.sh (1 kb padded background, v10nr_clust singletons, Cluster-Buster), nothing else.
# Runs locally: ~20.5k motifs x 499k regions. Needs ~60 GB of disk for the three feather files.
#
#   nohup ./create_cistarget_db_hybrid.sh > build.log 2>&1 &
#
# The consensus BED is copied into the database directory first, so the database records the regions it was built on.
set -euo pipefail

BASE_DIR="/hdd/jupyter/brad/scenicplus/motor-pathway_multiome/motor-pathway_multiome_seurat_cellbender.0.05_preprocess_cr/ra-arco-hvc-nc_hybrid"
SRC_BED="${BASE_DIR}/pycisTopic/scATAC/consensus_peak_calling/consensus_regions.bed"

# paths.sh still points at the pre-reorganization nest layout and at ~/repos; these are the current locations
GENOME_FASTA="/mnt/nest/common/assembly/lonStrDom2/ncbi/GCF_005870125.1_lonStrDom2_genomic_ucsc_only.fna"
CHROMSIZES="/mnt/nest/common/assembly/lonStrDom2/ucsc/chrom.sizes.ucsc"
SCRIPT_DIR="/home/brad/ssd/repos/create_cisTarget_databases"
CBDIR="/home/brad/nest/cistarget/v10nr_clust_public/singletons"
export PATH="/home/brad/micromamba/envs/scenicplus/bin:/opt/cactus-bin-v2.6.12/bin:${PATH}"   # cbust, bedtools, python deps

DATABASE_PREFIX="ra-arco-hvc-nc_hybrid"
OUT_DIR="/hdd/jupyter/brad/scenicplus/cistarget/${DATABASE_PREFIX}"
REGION_BED="${OUT_DIR}/consensus_regions.bed"
FASTA_FILE="${OUT_DIR}/lonStrDom2_1kb_bg_padding.fa"
MOTIF_LIST="${OUT_DIR}/motifs.txt"
THREADS="${THREADS:-40}"

mkdir -p "${OUT_DIR}"
if ls "${OUT_DIR}"/*.feather >/dev/null 2>&1; then
    echo "${OUT_DIR} already has feather files; remove them to rebuild" >&2
    exit 1
fi

cp -p "${SRC_BED}" "${REGION_BED}"
echo "$(date +%T) regions: $(wc -l < "${REGION_BED}")"

"${SCRIPT_DIR}/create_fasta_with_padded_bg_from_bed.sh" \
    "${GENOME_FASTA}" \
    "${CHROMSIZES}" \
    "${REGION_BED}" \
    "${FASTA_FILE}" \
    1000 \
    yes

ls "${CBDIR}" | grep -v '\.meme$' > "${MOTIF_LIST}"   # the old motifs.txt lists the cluster-buster (.cb) files
echo "$(date +%T) motifs: $(wc -l < "${MOTIF_LIST}")"

"${SCRIPT_DIR}/create_cistarget_motif_databases.py" \
    -f "${FASTA_FILE}" \
    -M "${CBDIR}" \
    -m "${MOTIF_LIST}" \
    -o "${OUT_DIR}/${DATABASE_PREFIX}" \
    --bgpadding 1000 \
    -t "${THREADS}"

echo "$(date +%T) done"
ls -lh "${OUT_DIR}"
