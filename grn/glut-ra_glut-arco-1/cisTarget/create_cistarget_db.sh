

source "$(dirname "$0")/../../../config/paths.sh"  # load shared reference paths

BASE_DIR="/hdd/jupyter/brad/scenicplus/motor-pathway_multiome/ra-arco-hvc-nc_seurat-clustering_241123/"
REGION_BED=$BASE_DIR/"pycisTopic/scATAC/consensus_peak_calling/consensus_regions.bed"
GENOME_FASTA="${GENOME_FA_LONSTR}"
CHROMSIZES="${CHROM_SIZES_LONSTR}"
SCRIPT_DIR="${CREATE_CISTARGET_SCRIPT_DIR}"

DATABASE_PREFIX="ra-arco-hvc-nc_seurat-clustering"
OUT_DIR=$DATABASE_PREFIX
CBDIR="${CISTARGET_SINGLETONS_DIR}"
FASTA_FILE="${OUT_DIR}/lonStrDom2_1kb_bg_padding.fa"
MOTIF_LIST="${OUT_DIR}/motifs.txt"

if [ ! -d $OUT_DIR ]
then
    mkdir $OUT_DIR
fi

${SCRIPT_DIR}/create_fasta_with_padded_bg_from_bed.sh \
        ${GENOME_FASTA} \
        ${CHROMSIZES} \
        ${REGION_BED} \
        ${FASTA_FILE} \
        1000 \
        yes

ls ${CBDIR} > $OUT_DIR/motifs.txt

"${SCRIPT_DIR}/create_cistarget_motif_databases.py" \
    -f ${FASTA_FILE} \
    -M ${CBDIR} \
    -m ${MOTIF_LIST} \
    -o ${OUT_DIR}/${DATABASE_PREFIX} \
    --bgpadding 1000 \
    -t 40
