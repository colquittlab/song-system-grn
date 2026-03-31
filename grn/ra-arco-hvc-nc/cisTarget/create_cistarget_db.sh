

BASE_DIR="/hdd/jupyter/brad/scenicplus/motor-pathway_multiome/ra-arco-hvc-nc_seurat-clustering_241123/"
REGION_BED=$BASE_DIR/"pycisTopic/scATAC/consensus_peak_calling/consensus_regions.bed"
GENOME_FASTA="/mnt/nest/assembly/lonStrDom2/ncbi/GCF_005870125.1_lonStrDom2_genomic_ucsc_only.fna"
CHROMSIZES="/mnt/nest/assembly/lonStrDom2/ucsc/chrom.sizes.ucsc"
SCRIPT_DIR="/home/brad/repos/create_cisTarget_databases"

DATABASE_PREFIX="ra-arco-hvc-nc_seurat-clustering"
OUT_DIR=$DATABASE_PREFIX
CBDIR="/home/brad/nest/cistarget/v10nr_clust_public/singletons"
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
