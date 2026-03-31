#! /usr/bin/env nextflow

chrom_pattern = /_[0-9]+$/

process chunkMAFs {

    executor 'slurm'
    queue 'medium'
    time '6:00:00'
    cpus 1
    memory '200G'

    input:
        tuple val(chrom), path(maf_file), path(refseq_file)

    output:
        tuple val(chrom), path("chunk_${root}_*.maf")

    script:
        root = maf_file.baseName

    """
    # Get chromosome size from refseq fasta index
    samtools faidx ${refseq_file}

    # Generate 1Mb window BED
    bedtools makewindows -g ${refseq_file}.fai -w 1000000 > windows.bed

    # Split MAF by windows
    mafSplit windows.bed chunk_${root}_ ${maf_file}
    """
}

process runMafFilter {

    publishDir(
        path: "${params.results_dir}",
        mode: 'copy',
        saveAs: { filename ->
            filename.replace('.maf', '.maf')
        }
    )

    memory { 50.GB * task.attempt }
    time { 24.h * task.attempt }
    errorStrategy { task.exitStatus in [143,137,104,134,139,9,140] ? 'retry' : 'terminate' }
    maxRetries 3
    executor 'slurm'
    queue 'long'
    cpus 1

    input:
        tuple val(chrom), val(chunk_name), path(maf_file)
        path config_template

    output:
        tuple val(chrom), path("${chunk_name}.filter.maf")
        path "${chunk_name}.trash_aln.maf", optional: true

    script:
    """
    # Create per-chunk config by substituting input/output file names
    sed -e 's|INPUT_MAF|${maf_file}|g' \
        -e 's|OUTPUT_MAF|${chunk_name}.filter.maf|g' \
        -e 's|TRASH_MAF|${chunk_name}.trash_aln.maf|g' \
        ${config_template} > maffilter_${chunk_name}.txt

    maffilter param=maffilter_${chunk_name}.txt
    """
}

workflow {

    // Load MAF files split by chromosome
    mafs_ch = channel.fromPath("${params.maf_dir}/*.maf")
        .map { file ->
            def chrom = file.baseName.replaceAll(".maf", "")
            tuple(chrom, file)
        }
        .filter { chrom, _file ->
            params.chroms ? chrom in params.chroms : true
        }

    // Combine MAFs with refseq for chunking
    maf_refseq_ch = mafs_ch.map { chrom, maf_file ->
        tuple(chrom, maf_file, file(params.refseq_file))
    }

    // Split MAF files into smaller chunks
    chunkMAFs(maf_refseq_ch)

    // Map chunks back to their chromosome name
    chunks_ch = chunkMAFs.out.flatMap { chrom, chunks ->
            chunks instanceof List
                ? chunks.collect { chunk -> tuple(chrom, chunk.baseName, chunk) }
                : [tuple(chrom, chunks.baseName, chunks)]
        }

    config_template = file(params.config_template)

    runMafFilter(chunks_ch, config_template)
}
