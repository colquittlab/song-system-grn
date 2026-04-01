

source "$(dirname "$0")/../../config/paths.sh"  # load shared reference paths

BASE_DIR="/hdd/brad/jupyter/multiome/motor-pathway/motor-pathway_multiome_seurat_no-sct_preprocess/"
REGION_BED=$BASE_DIR/"pycisTopic/pycisTopic_hvc-nc_glut/scATAC/consensus_peak_calling/consensus_regions.bed"
GENOME_FASTA="${GENOME_FA_LONSTR}"
CHROMSIZES="${CHROM_SIZES_LONSTR}"

SCRIPT_DIR="${CREATE_CISTARGET_SCRIPT_DIR}"

${SCRIPT_DIR}/create_fasta_with_padded_bg_from_bed.sh \
        ${GENOME_FASTA} \
        ${CHROMSIZES} \
        ${REGION_BED} \
        lonStrDom2_hvc-nc_glut_seurat-clustering_1kb_bg_padding.fa \
        1000 \
        yes
