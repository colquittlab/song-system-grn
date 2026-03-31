#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
 * Extract and split MAF files from a HAL alignment
 *
 * This pipeline:
 * 1. Extracts MAF for the entire genome using cactus-hal2maf
 * //2. Gets chromosome sizes from HAL using halStats
 * //3. Splits into per-chromosome MAF files using mafSplit
 */

process cactusHal2Maf {

    publishDir './results', mode: 'copy'

    memory { 400.GB * task.attempt }
    time { 24.h * task.attempt }
    errorStrategy { task.exitStatus in [143,137,104,134,139,9,140] ? 'retry' : 'terminate' }
    maxRetries 3
    executor 'slurm'
    queue 'long'
    cpus 1

    //executor 'local'

    input:
        path hal_file
        val species

    output:
        path "genome.maf"

    script:
    """
    cactus-hal2maf \\
        --chunkSize 100000 \\
        --batchCores 96 \\
        --batchCount 10 \\
        --noAncestors \\
        --filterGapCausingDupes \\
        --outType norm single \\
        --batchSystem slurm \\
        --maxLocalJobs 800 \\
        --slurmTime 200:00:00 \\
        --slurmPartition long \\
        --refGenome ${species} \\
        --logFile avian.maf.gz.log \\
        ./js \\
        ${hal_file} \\
        genome.maf
    """
}

process halStats {

    executor 'local'

    input:
        path hal_file
        val species

    output:
        path "genome.bed"

    script:
    """
    halStats --chromSizes ${species} ${hal_file} | awk 'BEGIN{OFS="\t"}{print \$1, 0, \$2}' > genome.bed
    """
}

process mafSplit {

    publishDir './results/chrom_mafs', mode: 'copy'

    memory { 50.GB * task.attempt }
    time { 4.h * task.attempt }
    errorStrategy { task.exitStatus in [143,137,104,134,139,9,140] ? 'retry' : 'terminate' }
    maxRetries 3
    executor 'slurm'
    queue 'medium'
    cpus 1

    input:
        path maf_file

    output:
        path "chrom_mafs/*.maf"

    script:
    """
    mkdir -p chrom_mafs
    mafSplit -byTarget -useFullSequenceName dummy.bed chrom_mafs/ ${maf_file}
    """
}

workflow {

    // Validate required parameters
    if (!params.hal_file) {
        error "Please provide --hal_file parameter"
    }
    if (!params.species) {
        error "Please provide --species parameter"
    }

    hal_ch = channel.fromPath(params.hal_file, checkIfExists: true)

    // Extract MAF and get chromosome sizes in parallel
    cactusHal2Maf(hal_ch, params.species)
    //halStats(hal_ch, params.species)

    // Split into per-chromosome MAF files
    //mafSplit(cactusHal2Maf.out)
}
