#! /usr/bin/env nextflow

process extractChromSizes {

    executor 'local'

    input:
        path hal_file
        val species

    output:
        path "chrom.sizes"

    script:
    """
    halStats --chromSizes ${species} ${hal_file} > chrom.sizes
    """
}

process runMafFilter {

    memory { 50.GB * task.attempt }
    time { 2.h * task.attempt }
    errorStrategy { task.exitStatus in [143,137,104,134,139,9,140] ? 'retry' : 'terminate' }
    maxRetries 3
    executor 'slurm'
    queue 'medium'
    cpus 1

    input:
        tuple val(chrom), val(chunk_name), path(maf_file)
        path config_template

    output:
        tuple val(chrom), path("${chunk_name}.filter.maf")

    script:
    """
    sed -e 's|INPUT_MAF|${maf_file}|g' \
        -e 's|OUTPUT_MAF|${chunk_name}.filter.maf|g' \
        -e 's|TRASH_MAF|${chunk_name}.trash_aln.maf|g' \
        ${config_template} > maffilter_${chunk_name}.txt

    maffilter param=maffilter_${chunk_name}.txt
    """
}

process runPhyloPSubtree {

    memory { 80.GB * task.attempt }
    time { 12.h * task.attempt }
    errorStrategy { task.exitStatus == 1 ? 'ignore' : (task.exitStatus in [143,137,104,134,139,9,140] ? 'retry' : 'terminate') }
    maxRetries 3
    executor 'slurm'
    queue 'long'
    cpus 1

    input:
        tuple val(set_name), val(chromname), path(infile), path(neutral_model)
        val mode
        val method
        val subtree

    output:
        tuple val(set_name), path("scores.wig")

    script:
    """
    phyloP --wig-scores --msa-format MAF --method ${method} --mode ${mode} --subtree ${subtree} \
            --chrom ${chromname} ${neutral_model} ${infile} > scores.wig
    """
}

process runPhyloPAll {

    memory { 80.GB * task.attempt }
    time { 6.h * task.attempt }
    errorStrategy { task.exitStatus == 1 ? 'ignore' : (task.exitStatus in [143,137,104,134,139,9, 140] ? 'retry' : 'terminate') }
    maxRetries 5
    executor 'slurm'
    queue 'long'
    cpus 1

    input:
        tuple val(set_name), val(chromname), path(infile), path(neutral_model)
        val mode
        val method

    output:
        tuple val(set_name), path("scores.wig")

    script:
    """
    phyloP --wig-scores --msa-format MAF --method ${method} --mode ${mode} \
            --chrom ${chromname} ${neutral_model} ${infile} > scores.wig
    """
}

process combineWigsBySet {

    publishDir(
        path: "${params.results_dir}",
        mode: 'copy',
        saveAs: { filename ->
            "${subdir}/${filename}"
        }
    )

    executor 'local'

    input:
        tuple val(set_name), path("scores.wig")
        path subdir
    output:
        path "phyloP_${set_name}.wig"
    script:
    """
    cat scores.wig* > phyloP_${set_name}.wig
    """
}

process combineAllWigs {

    publishDir(
        path: "${params.results_dir}",
        mode: 'copy',
        saveAs: { filename ->
            "${subdir}/${filename}"
        }
    )

    executor 'local'

    input:
        path "phyloP_*.wig"
        path subdir
    output:
        path "phyloP.wig"
    script:
    """
    cat phyloP_*.wig > phyloP.wig
    """
}

process convertWigToBigWig {

    publishDir(
        path: "${params.results_dir}",
        mode: 'copy',
        saveAs: { filename ->
            "${subdir}/${filename}"
        }
    )

    executor 'local'

    input:
        path infile
        path chrom_sizes
        path subdir
    output:
        path "phyloP.bw"
    script:
    """
    wigToBigWig ${infile} ${chrom_sizes} phyloP.bw
    """
}

workflow {

    // Extract chromosome sizes for target genome
    extractChromSizes(params.hal_file, params.target_species)

    // Define chromosome to set mapping
    def chrom_to_set = [:]
    params.chrom_sets.each { set_name, chrom_str ->
        chrom_str.split(',').each { chrom ->
            chrom_to_set[chrom.trim()] = set_name
        }
    }

    // Load pre-chunked MAF files from chunked_maf_dir (output of chunkMAFs pipeline)
    // Directory structure: chunked_maf_dir/<chrom>/chunk_<chrom>_*.maf
    chunks_ch = channel.fromPath("${params.chunked_maf_dir}/*/chunk_*.maf")
        .map { file ->
            def parent = file.parent.name
            // Filename format from mafSplit: chunk_<parent>_<chrom>.<num>.maf
            def chrom = file.baseName.replaceAll("^chunk_${parent}_", '').replaceAll(/\.\d+$/, '')
            tuple(chrom, file)
        }
        .filter { chrom, _file -> chrom_to_set.containsKey(chrom) }
        .filter { chrom, chunk -> chunk.readLines().size() >= 3 }

    // Get chromosome set models
    rescaled_models_ch = channel.fromPath("${params.rescaled_models_dir}/rescaled_*.mod")
        .map { file ->
            def set_name = file.baseName.replaceAll("rescaled_", "").replaceAll(".mod", "")
            tuple(set_name, file)
        }
        .filter { set_name, _file -> set_name != 'unassigned' }

    subdir = file("./results/${params.mode}_${params.method}_${params.subtree}")

    // Optionally filter chunks with maffilter
    if (params.maffilter_config) {
        maffilter_config = file(params.maffilter_config)

        maffilter_input = chunks_ch
            .map { chrom, chunk -> tuple(chrom, chunk.baseName, chunk) }

        runMafFilter(maffilter_input, maffilter_config)

        filtered_chunks = runMafFilter.out
            .filter { chrom, maf -> maf.readLines().size() >= 3 }
    } else {
        filtered_chunks = chunks_ch
    }

    // Map chunks to their chromosome set and appropriate model
    chunks_with_models = filtered_chunks
        .map { chromname, chunk ->
            def set_name = chrom_to_set.get(chromname, 'unassigned')
            tuple(set_name, chromname, chunk)
        }
        .combine(rescaled_models_ch, by: 0)
        .map { set_name, chromname, chunk, model ->
            tuple(set_name, chromname, chunk, model)
        }

    if (params.subtree == "all") {
        runPhyloPAll(chunks_with_models, params.mode, params.method)
        ch_wigs_by_set = runPhyloPAll.out.groupTuple(by: 0)
    } else {
        runPhyloPSubtree(chunks_with_models, params.mode, params.method, params.subtree)
        ch_wigs_by_set = runPhyloPSubtree.out.groupTuple(by: 0)
    }

    // Combine wigs within each chromosome set
    combineWigsBySet(ch_wigs_by_set, subdir)

    // Combine all sets into final wig
    combineAllWigs(combineWigsBySet.out.collect(), subdir)

    // Convert to BigWig
    convertWigToBigWig(combineAllWigs.out, extractChromSizes.out, subdir)
}
