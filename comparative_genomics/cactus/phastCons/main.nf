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

// Use chooseLines to randomly select rho_sample_size chunk paths per set.
// Runs once per chromosome set; finds MAF files via shell find so no
// Nextflow channel-close dependency is needed.
process sampleForRho {

    executor 'local'

    input:
        tuple val(set_name), val(chrom_list)

    output:
        tuple val(set_name), path("sample_list.txt")

    script:
    def chroms = chrom_list.split(',').collect { it.trim() }.join(' ')
    """
    for chrom in ${chroms}; do
        for f in ${params.chunked_maf_dir}/genome/chunk_genome_\${chrom}.*.maf; do
            [ -f "\$f" ] && [ \$(wc -l < "\$f") -ge 3 ] && echo "\$f"
        done
    done | shuf -n ${params.rho_sample_size} > sample_list.txt
    """
}

// Estimate rho (scaling factor) for conserved/non-conserved models per chunk
process estimateRho {

    errorStrategy 'ignore' // some chunks fail to converge
    executor 'slurm'
    queue 'long'
    time '24:00:00'
    cpus 1

    input:
        tuple val(set_name), val(chromname), path(maf_file)

    output:
        tuple val(set_name), path("newtree.cons.mod"), path("newtree.noncons.mod")

    script:
    """
    phastCons \
        --target-coverage ${params.target_coverage} \
        --expected-length ${params.expected_length} \
        --msa-format MAF \
        --seqname ${chromname} \
        --gc 0.4 --no-post-probs \
        --estimate-rho newtree \
        ${maf_file} ${params.rescaled_models_dir}/rescaled_${set_name}.mod
    """
}

// Average cons and noncons models separately per chromosome set.
// Files staged with the path pattern so phyloBoot can glob them.
process averageModelsBySet {

    publishDir(
        path: "${params.results_dir}",
        mode: 'copy',
        saveAs: { filename -> "${subdir}/models/${filename}" }
    )

    executor 'local'

    input:
        tuple val(set_name), path("newtree.cons.*.mod"), path("newtree.noncons.*.mod"), path(subdir)

    output:
        tuple val(set_name), path("ave.cons.mod"), path("ave.noncons.mod")

    script:
    """
    phyloBoot --read-mods newtree.cons.*.mod* --output-average ave.cons.mod
    phyloBoot --read-mods newtree.noncons.*.mod* --output-average ave.noncons.mod
    """
}

// Run phastCons on each chunk using the averaged models for its chromosome set
process runPhastCons {

    memory { 40.GB * task.attempt }
    time { 6.h * task.attempt }
    errorStrategy { task.exitStatus == 1 ? 'ignore' : (task.exitStatus in [143,137,104,134,139,9,140] ? 'retry' : 'terminate') }
    maxRetries 3
    executor 'slurm'
    queue 'long'
    cpus 1

    input:
        tuple val(set_name), val(chromname), path(maf_file), path(cons_model), path(noncons_model)

    output:
        tuple val(set_name), path("chunk.bed"), path("chunk.wig")

    script:
    """
    phastCons \
        --target-coverage ${params.target_coverage} \
        --expected-length ${params.expected_length} \
        --msa-format MAF \
        --most-conserved chunk.bed --score \
        --seqname ${chromname} \
        ${maf_file} ${cons_model},${noncons_model} \
        > chunk.wig
    """
}

process combineBySet {

    publishDir(
        path: "${params.results_dir}",
        mode: 'copy',
        saveAs: { filename -> "${subdir}/${filename}" }
    )

    executor 'local'

    input:
        tuple val(set_name), path("chunk.bed"), path("chunk.wig"), path(subdir)

    output:
        tuple val(set_name), path("phastCons_${set_name}.bed"), path("phastCons_${set_name}.wig")

    script:
    """
    cat chunk.bed* | sort -k1,1 -k2,2n > phastCons_${set_name}.bed
    cat chunk.wig* > phastCons_${set_name}.wig
    """
}

process combineAllSets {

    publishDir(
        path: "${params.results_dir}",
        mode: 'copy',
        saveAs: { filename -> "${subdir}/${filename}" }
    )

    executor 'local'

    input:
        path "phastCons_*.bed"
        path "phastCons_*.wig"
        path subdir

    output:
        tuple path("phastCons.bed"), path("phastCons.wig")

    script:
    """
    cat phastCons_*.bed | sort -k1,1 -k2,2n > phastCons.bed
    cat phastCons_*.wig > phastCons.wig
    """
}

process convertWigToBigWig {

    publishDir(
        path: "${params.results_dir}",
        mode: 'copy',
        saveAs: { filename -> "${subdir}/${filename}" }
    )

    executor 'local'

    input:
        tuple path(wig_file), path(chrom_sizes)
        path subdir

    output:
        path "phastCons.bw"

    script:
    """
    wigToBigWig ${wig_file} ${chrom_sizes} phastCons.bw
    """
}

workflow {

    // Extract chromosome sizes for target genome
    extractChromSizes(params.hal_file, params.target_species)

    // Build chrom -> set_name map from chrom_sets parameter
    def chrom_to_set = [:]
    params.chrom_sets.each { set_name, chrom_str ->
        chrom_str.split(',').each { chrom ->
            chrom_to_set[chrom.trim()] = set_name
        }
    }

    subdir = file("./results/tc${params.target_coverage}_el${params.expected_length}")

    // Load pre-chunked MAF files (output of chunkMAFs pipeline)
    // Directory structure: chunked_maf_dir/<chrom>/chunk_<chrom>_*.maf
    chunks_ch = channel.fromPath("${params.chunked_maf_dir}/*/chunk_*.maf")
        .map { file ->
            def parent = file.parent.name
            def chromname = file.baseName.replaceAll("^chunk_${parent}_", '').replaceAll(/\.\d+$/, '')
            tuple(chromname, file)
        }
        .filter { chromname, _file -> chrom_to_set.containsKey(chromname) }
        .filter { _chromname, chunk -> chunk.readLines().size() >= 3 }
        .map { chromname, file ->
            def set_name = chrom_to_set[chromname]
            tuple(set_name, chromname, file)
        }

    // Run sampleForRho once per chromosome set using a simple value channel —
    // no file-scanning operators needed, so no channel-close issues
    sets_ch = channel.from(
        params.chrom_sets.collect { set_name, chrom_str -> tuple(set_name, chrom_str) }
    )
    sampleForRho(sets_ch)

    // Read selected paths and reconstruct (set_name, chromname, maf_file) tuples
    sample_ch = sampleForRho.out
        .flatMap { set_name, sample_list ->
            sample_list.readLines()
                .findAll { it.trim() }
                .collect { path_str ->
                    def maf = file(path_str.trim())
                    def parent = maf.parent.name
                    def chrom = maf.baseName.replaceAll("^chunk_${parent}_", '').replaceAll(/\.\d+$/, '')
                    tuple(set_name, chrom, maf)
                }
        }

    // Estimate rho on sampled chunks only
    estimateRho(sample_ch)
    // output: tuple(set_name, cons_mod, noncons_mod)

    // Group models by chromosome set for separate averaging
    cons_by_set_ch = estimateRho.out
        .map { set_name, cons_mod, noncons_mod -> tuple(set_name, cons_mod) }
        .groupTuple()

    noncons_by_set_ch = estimateRho.out
        .map { set_name, cons_mod, noncons_mod -> tuple(set_name, noncons_mod) }
        .groupTuple()

    // Average models separately per chromosome set
    averageModelsBySet(
        cons_by_set_ch
            .join(noncons_by_set_ch)
            .map { set_name, cons_mods, noncons_mods -> tuple(set_name, cons_mods, noncons_mods, subdir) }
    )
    // output: tuple(set_name, ave.cons.mod, ave.noncons.mod)

    // Pair each chunk with the averaged models for its chromosome set
    run_ch = chunks_ch
        .combine(averageModelsBySet.out, by: 0)
        .map { set_name, chromname, maf_file, ave_cons, ave_noncons ->
            tuple(set_name, chromname, maf_file, ave_cons, ave_noncons)
        }

    // Run phastCons with the appropriate set-specific averaged models
    runPhastCons(run_ch)
    // output: tuple(set_name, bed, wig)

    // Combine beds and wigs within each chromosome set
    by_set_ch = runPhastCons.out
        .groupTuple(by: 0)
        .map { set_name, beds, wigs -> tuple(set_name, beds, wigs, subdir) }

    combineBySet(by_set_ch)

    // Combine all sets into final outputs
    all_beds = combineBySet.out.map { it[1] }.collect()
    all_wigs = combineBySet.out.map { it[2] }.collect()

    combineAllSets(all_beds, all_wigs, subdir)

    // Convert to BigWig
    convertWigToBigWig(
        combineAllSets.out
            .map { bed, wig -> wig }
            .combine(extractChromSizes.out)
            .map { wig, chrom_sizes -> tuple(wig, chrom_sizes) },
        subdir
    )
}
