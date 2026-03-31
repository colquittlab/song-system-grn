

BASE_DIR="/hdd/brad/jupyter/multiome/motor-pathway/motor-pathway_multiome_seurat_no-sct_preprocess/"
REGION_BED=$BASE_DIR/"pycisTopic/pycisTopic_hvc-nc_glut/scATAC/consensus_peak_calling/consensus_regions.bed"
GENOME_FASTA="/mnt/nest/assembly/lonStrDom2/ncbi/GCF_005870125.1_lonStrDom2_genomic_ucsc_only.fna"
CHROMSIZES="/mnt/nest/assembly/lonStrDom2/ucsc/chrom.sizes.ucsc"

SCRIPT_DIR="/home/brad/repos/create_cisTarget_databases"

${SCRIPT_DIR}/create_fasta_with_padded_bg_from_bed.sh \
        ${GENOME_FASTA} \
        ${CHROMSIZES} \
        ${REGION_BED} \
        lonStrDom2_hvc-nc_glut_seurat-clustering_1kb_bg_padding.fa \
        1000 \
        yes
