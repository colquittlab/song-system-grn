#! /usr/bin/env nextflow

process preprocessBed {

    executor 'local'

    input:
        path features_bed

    output:
        path "features_by_chrom/*.bed"

    script:
    """
    mkdir features_by_chrom

    # Trim to 3 columns and split by chromosome
    awk -F'\\t' 'BEGIN{OFS="\\t"} {
        print \$1, \$2, \$3 > "features_by_chrom/"\$1".bed"
    }' ${features_bed}
    """
}

process runPhyloPFeaturesSubtree {

    memory { 80.GB * task.attempt }
    time { 12.h * task.attempt }
    errorStrategy { task.exitStatus == 1 ? 'ignore' : (task.exitStatus in [143,137,104,134,139,9,140] ? 'retry' : 'terminate') }
    maxRetries 3
    executor 'slurm'
    queue 'long'
    cpus 1

    input:
        tuple val(set_name), val(chromname), path(infile), path(neutral_model), path(features_bed)
        val mode
        val method
        val subtree

    output:
        tuple val(set_name), path("scores.txt")

    script:
    """
    phyloP --features ${features_bed} --msa-format MAF --method ${method} --mode ${mode} --gff-scores --subtree ${subtree} \
            --chrom ${chromname} ${neutral_model} ${infile} > scores.txt
    """
}

process runPhyloPFeaturesAll {

    memory { 80.GB * task.attempt }
    time { 12.h * task.attempt }
    errorStrategy { task.exitStatus == 1 ? 'ignore' : (task.exitStatus in [143,137,104,134,139,9, 140] ? 'retry' : 'terminate') }
    maxRetries 5
    executor 'slurm'
    queue 'long'
    cpus 1

    input:
        tuple val(set_name), val(chromname), path(infile), path(neutral_model), path(features_bed)
        val mode
        val method

    output:
        tuple val(set_name), path("scores.txt")

    script:
    """
    phyloP --features ${features_bed} --msa-format MAF --method ${method} --mode ${mode} --gff-scores \
            --chrom ${chromname} ${neutral_model} ${infile} > scores.txt
    """
}

process combineResultsBySet {

    publishDir(
        path: "${params.results_dir}",
        mode: 'copy',
        saveAs: { filename ->
            "${subdir}/${filename}"
        }
    )

    executor 'local'

    input:
        tuple val(set_name), path("scores.txt")
        path subdir
    output:
        path "phyloP_features_${set_name}.txt"
    script:
    """
    cat scores.txt* > phyloP_features_${set_name}.txt
    """
}

process combineAllResults {

    publishDir(
        path: "${params.results_dir}",
        mode: 'copy',
        saveAs: { filename ->
            "${subdir}/${filename}"
        }
    )

    executor 'local'

    input:
        path "phyloP_features_*.txt"
        path subdir
    output:
        path "phyloP_features.txt"
    script:
    """
    cat phyloP_features_*.txt > phyloP_features.txt
    """
}

process convertToBedGraph {

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
        path subdir
    output:
        path "phyloP_features.bedGraph"
    script:
    """
    # GFF cols: 1=chr, 4=start(1-based), 5=end, 6=score; convert start to 0-based for bedGraph
    awk -F'\\t' 'BEGIN{OFS="\\t"} !/^#/{print \$1, \$4-1, \$5, \$6}' ${infile} \
        | sort -k1,1 -k2,2n > phyloP_features.bedGraph
    """
}

process convertToBigWig {

    publishDir(
        path: "${params.results_dir}",
        mode: 'copy',
        saveAs: { filename ->
            "${subdir}/${filename}"
        }
    )

    executor 'local'

    input:
        path bedgraph
        path chrom_sizes
        path subdir
    output:
        path "phyloP_features.bw"
    script:
    """
    bedGraphToBigWig ${bedgraph} ${chrom_sizes} phyloP_features.bw
    """
}

workflow {

    // Define chromosome to set mapping
    def chrom_to_set = [:]
    params.chrom_sets.each { set_name, chrom_str ->
        chrom_str.split(',').each { chrom ->
            chrom_to_set[chrom.trim()] = set_name
        }
    }

    // Preprocess BED: fix score column to integer and split by chromosome
    preprocessBed(file(params.features_bed))

    // Create channel of per-chromosome BED files
    features_ch = preprocessBed.out.flatten()
        .map { file ->
            def chrom = file.baseName
            tuple(chrom, file)
        }

    // Load pre-chunked MAF files from chunked_maf_dir (output of chunkMAFs pipeline)
    // Directory structure: chunked_maf_dir/<chrom>/chunk_<chrom>_*.maf
    chunks_ch = channel.fromPath("${params.chunked_maf_dir}/*/chunk_*.maf")
        .map { file ->
            def chrom = file.parent.name
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

    subdir = file("./results/features_${params.mode}_${params.method}_${params.subtree}")

    // Join chunks with their chromosome's BED, then map to set and model
    chunks_with_features = chunks_ch
        .combine(features_ch, by: 0)
        .map { chromname, chunk, bed ->
            def set_name = chrom_to_set.get(chromname, 'unassigned')
            tuple(set_name, chromname, chunk, bed)
        }
        .combine(rescaled_models_ch, by: 0)
        .map { set_name, chromname, chunk, bed, model ->
            tuple(set_name, chromname, chunk, model, bed)
        }

    if (params.subtree == "all") {
        runPhyloPFeaturesAll(chunks_with_features, params.mode, params.method)
        ch_results_by_set = runPhyloPFeaturesAll.out.groupTuple(by: 0)
    } else {
        runPhyloPFeaturesSubtree(chunks_with_features, params.mode, params.method, params.subtree)
        ch_results_by_set = runPhyloPFeaturesSubtree.out.groupTuple(by: 0)
    }

    // Combine results within each chromosome set
    combineResultsBySet(ch_results_by_set, subdir)

    // Combine all sets into final results
    combineAllResults(combineResultsBySet.out.collect(), subdir)

    // Convert to bedGraph then bigWig for UCSC genome browser
    convertToBedGraph(combineAllResults.out, subdir)
    convertToBigWig(convertToBedGraph.out, file(params.chrom_sizes), subdir)
}
