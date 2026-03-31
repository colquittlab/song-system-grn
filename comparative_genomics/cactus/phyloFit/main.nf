
process extractCodons {

    memory { 10.GB * task.attempt }
    time { 10.m * task.attempt }
    errorStrategy { task.exitStatus in [143,137,104,134,139,9, 140] ? 'retry' : 'terminate' }
    maxRetries 3
    executor 'slurm'
    queue 'short'
    cpus 1

    input:
        tuple path(maf_file), path(gff_file)    

    output:
        path "4d-codons.ss"

    script:
    """
    msa_view ${maf_file} --4d --features ${gff_file} > 4d-codons.ss
    """

}

process extractGenomeFasta { 
    
    executor 'local'

    input:
        path hal_file
        val species

    output:
        path "ancestral.fasta"

    script:
    """
    hal2fasta ${hal_file} ${species} > ancestral.fasta
    """
}


process faSplit {

    //memory { 10.GB * task.attempt }
    //time { 10.m * task.attempt }
    //errorStrategy { task.exitStatus in [143,137,104,134,139,9, 140] ? 'retry' : 'terminate' }
    //maxRetries 3
    executor 'local'
    //queue 'short'
    //cpus 1

    input:
        path fasta_file

    output:
        path "ancestral_fasta_split/*"

    script:
    """
    mkdir ancestral_fasta_split
    faSplit byname ${fasta_file} ancestral_fasta_split/
    """
} 

process buildRepeatMaskerLib {

    executor 'local'

    output:
        val true

    script:
    """
    echo ">dummy" > dummy.fa
    echo "ACGT" >> dummy.fa
    RepeatMasker -species ${params.repeatmasker_species} dummy.fa -xsmall -pa 1 || true
    """
}

process getAncestralRepeats {

    memory { 50.GB * task.attempt }
    time { 3.h * task.attempt }
    errorStrategy { task.exitStatus in [143,137,104,134,139,9, 140] ? 'retry' : 'terminate' }
    maxRetries 3
    executor 'slurm'
    queue 'medium'
    cpus 4

    input:
        path fasta_file
        val lib_ready

    output:
        path "*.out"

    script:
    """
    RepeatMasker -species ${params.repeatmasker_species} ${fasta_file} -xsmall -pa ${task.cpus}
    """
}   

process generateAncestralRepeatsBed {

    memory { 10.GB * task.attempt }
    time { 10.m * task.attempt }
    errorStrategy { task.exitStatus in [143,137,104,134,139,9, 140] ? 'retry' : 'terminate' }
    maxRetries 3
    executor 'slurm'
    queue 'short'
    cpus 1

    input:
        path repeats_file

    output:
        path "ancestral_repeats.bed"

    script:
    """
    grep -Ev "tRNA|SINE/tRNA.*|Low_complexity" ${repeats_file} | awk '{OFS="\t"}{if (NR>3) print \$5, \$6, \$7, \$11}' > ancestral_repeats.bed
    """
}

process combineAncestralRepeatsBed {

    publishDir './results', mode: 'copy'

    executor 'local'

    input:
        path "ancestral_repeats.bed*"

    output:
        path "ancestral_repeats_combined.bed"

    script:
    """
    cat ancestral_repeats.bed* | sort -k1,1 -k2,2n > ancestral_repeats_combined.bed
    """
}

process subSampleBed {

    executor 'local'

    input:
        path ancestral_repeats_bed_file
        val bed_lines

    output:
        path "subsampled-ancestral-repeats.bed"

    script:
    """
    shuf -n ${bed_lines} ${ancestral_repeats_bed_file} > subsampled-ancestral-repeats.bed
    """
}

process getAncestralChromSizes {

    executor 'local'

    input:
        path hal_file
        val species

    output:
        path "ancestral.chrom.sizes"

    script:
    """
    halStats --chromSizes ${species} ${hal_file} > ancestral.chrom.sizes
    """
}

process randomSelectBedPositions {

    executor 'local'

    input:
            path bed_file
	        val n_positions

    output:
            path "random_selected.bed"

    script:
    """
    # Expand bed intervals to individual positions
    awk 'BEGIN{OFS="\t"}{
        for(i=\$2; i<\$3; i++) {
           print \$1, i, i+1
        }
    }' ${bed_file} | shuf -n ${n_positions} | sort -k1,1 -k2,2n > random_selected.bed
    """
}

process filterBedByChromSizes {

    executor 'local'

    input:
        path bed_file
        path chrom_sizes

    output:
        path "filtered.bed"

    script:
    """
    # Filter bed to only include regions within valid chromosome ranges
    # Also clamp any regions that extend beyond chromosome ends
    awk 'BEGIN{OFS="\t"}
    NR==FNR {sizes[\$1]=\$2; next}
    {
        if (\$1 in sizes) {
            start = \$2
            end = \$3
            max_size = sizes[\$1]
            # Skip if start is beyond chromosome end
            if (start >= max_size) next
            # Clamp end to chromosome size
            if (end > max_size) end = max_size
            # Ensure valid interval
            if (start >= 0 && end > start) {
                print \$1, start, end
            }
        }
    }
    ' ${chrom_sizes} ${bed_file} > filtered.bed
    """
}


process liftOverRepeatsBed {

    memory { 100.GB * task.attempt }
    time { 10.m * task.attempt }
    errorStrategy { task.exitStatus in [143,137,104,139,134,9,140] ? 'retry' : 'terminate' }
    maxRetries 3
    executor 'slurm'
    queue 'short'
    cpus 1

    input:
        path ancestral_repeats_bed_file

    output:
        path "ancestral_repeats_lifted.bed"

    script:
    """
    halLiftover ${params.hal_file} ${params.ancestral_species} ${ancestral_repeats_bed_file} ${params.target_species} ancestral_repeats_lifted.bed
    """
}
process splitBedByChromosome {

    executor 'local'

    input:
        path ancestral_repeats_bed_file

    output:
        path "ancestral_repeats_by_chromosome/*"

    script:
    """
    splitFileByColumn ${ancestral_repeats_bed_file} ancestral_repeats_by_chromosome/ 
    """
}

process splitBedByChromosomeSet {

    executor 'local'

    input:
            path bed_file
	            val chrom_sets_json

    output:
            path "bed_by_set/*"

    script:
        """
    #!/usr/bin/env python3
    import json
    import os
    from collections import defaultdict
    
    # Parse chromosome sets
    chrom_sets = json.loads('${chrom_sets_json}')
    
    # Create mapping from chromosome to set
    chrom_to_set = {}
    for set_name, chrom_str in chrom_sets.items():
        for chrom in chrom_str.split(','):
            chrom_to_set[chrom.strip()] = set_name
    
    # Create output directory
    os.makedirs("bed_by_set", exist_ok=True)
    
    # Group bed entries by chromosome set
    set_entries = defaultdict(list)
    
    with open("${bed_file}", 'r') as f:
        for line in f:
            if line.strip():
                chrom = line.split('\t')[0]
                set_name = chrom_to_set.get(chrom, 'unassigned')
                set_entries[set_name].append(line)
    
    # Write out bed files for each set
    for set_name, entries in set_entries.items():
        with open(f"bed_by_set/{set_name}.bed", 'w') as out:
            for entry in sorted(entries):
                out.write(entry)
    """
}
    
process extractMafFromHal {

    publishDir './results/mafs', mode: 'copy'

    memory { 200.GB * task.attempt }
    time { 12.h * task.attempt }
    errorStrategy { task.exitStatus in [143,137,104,134,139,9,140] ? 'retry' : 'terminate' }
    maxRetries 3
    executor 'slurm'
    queue 'long'
    cpus 48

    input:
        tuple path(hal_file), val(ref_species), val(prefix), path(bed_file)

    output:
        tuple val(prefix), path("${prefix}_*.maf")

    script:
    def chrom = bed_file.getBaseName().toString().split(".bed")[0]
    """
    cactus-hal2maf --batchCores 48  --refGenome ${ref_species} --onlyOrthologs --bedRanges ${bed_file} ./jobStore ${hal_file} ${prefix}_${chrom}.maf
    """
}

process splitBedByChromosome_lifted {

    executor 'local'

    input:
        path ancestral_repeats_bed_file

    output:
        path "ancestral_repeats_lifted_by_chromosome/*"

    script:
    """
    splitFileByColumn ${ancestral_repeats_bed_file} ancestral_repeats_lifted_by_chromosome/ 
    """
}

process extractFeatures {

    memory { 50.GB * task.attempt }
    time { 1.h * task.attempt }
    errorStrategy { task.exitStatus in [143,137,104,134,139,9, 140] ? 'retry' : 'terminate' }
    maxRetries 3
    executor 'slurm'
    queue 'medium'
    cpus 1

    input:
        tuple path(maf_file), path(bed_file)

    output:
        path "features.ss"

    script:
    """
    msa_view ${maf_file} --features ${bed_file} --out-format SS > features.ss
    """
}

process aggregateSS {

    executor 'local'

    input:
        path "features.ss"
        val species_list

    output:
        path "all-features.ss"

    script:
        
    """
    msa_view --aggregate ${species_list} --unordered-ss --out-format SS features.ss* > all-features.ss
    """
}

process getTree {

    executor 'local'

    input:
        path hal_file

    output:
        path "tree.nwk"

    script:
    """
    halStats --tree ${hal_file} | sed -e 's/:[0-9.]*//g' > tree.nwk
    """
}

process runPhyloFit {

    publishDir './results', mode: 'copy'

    memory { 40.GB * task.attempt }
    time { 24.h * task.attempt }
    errorStrategy { task.exitStatus in [143,137,104,134,139,9,140] ? 'retry' : 'terminate' }
    maxRetries 3
    executor 'slurm'
    queue 'long'
    cpus 96

    input:
    path maf_file
	path tree

    output:
        path "ancestral-repeats.mod"

    script:
    """

    # Extract features from MAF
    msa_view ${maf_file} --out-format SS --unordered-ss > features.ss
    
    # Train phyloFit model
    /private/groups/colquittlab/opt/phast_new/bin/phyloFit --threads 96 --tree ${tree} --msa-format SS --out-root ancestral-repeats features.ss
    """
}

process rescaleByChromosome {

    publishDir './results/rescaled_models', mode: 'copy'

    memory { 20.GB * task.attempt }
    time { 2.h * task.attempt }
    errorStrategy { task.exitStatus in [143,137,104,134,139,9, 140] ? 'retry' : 'terminate' }
    maxRetries 3
    executor 'slurm'
    queue 'medium'
    cpus 1

    input:
        tuple val(set_name), path(maf_file)
        path model_file

    output:
        path "rescaled_${set_name}.mod"

    script:
        """
    # phyloFit rescale takes MAF directly and uses internal msa_view
    phyloFit --scale-only --msa-format MAF --init-model ${model_file} --out-root rescaled_${set_name} ${maf_file}
    """
}
    
workflow {

    // Extract ancestral genome and get repeats
    buildRepeatMaskerLib()
    extractGenomeFasta(params.hal_file, params.ancestral_species)
	faSplit(extractGenomeFasta.out)
	getAncestralRepeats(faSplit.out.flatten(), buildRepeatMaskerLib.out)
	generateAncestralRepeatsBed(getAncestralRepeats.out.flatten())
	combineAncestralRepeatsBed(generateAncestralRepeatsBed.out.collect())

    // Get ancestral chromosome sizes for filtering
    getAncestralChromSizes(params.hal_file, params.ancestral_species)

    // Filter combined repeats to ensure they are within ancestral genome ranges
    filterBedByChromSizes(combineAncestralRepeatsBed.out, getAncestralChromSizes.out)

    // Randomly select N positions from filtered ancestral repeats for base model training
    randomSelectBedPositions(filterBedByChromSizes.out, params.n_positions)

    // Extract tree from HAL only if params.tree is not provided
    if (params.tree == null || params.tree == '') {
	    getTree(params.hal_file)
	    tree_ch = getTree.out
    } else {
        tree_ch = channel.fromPath(params.tree)
    }

    // Lift over the randomly selected positions to TARGET genome for rescaling
    liftOverRepeatsBed(randomSelectBedPositions.out)

    // Split lifted bed by chromosome SET (not individual chromosomes)
    def chrom_sets_json = groovy.json.JsonOutput.toJson(params.chrom_sets)
	splitBedByChromosomeSet(liftOverRepeatsBed.out, chrom_sets_json)

    // Prepare channels for both ancestral and target MAF extraction
    // Ancestral: single bed file with 'ancestral' prefix using ancestral species
    ancestral_bed_ch = randomSelectBedPositions.out
        .map { bed -> tuple(params.hal_file, params.ancestral_species, 'ancestral', bed) }

    // Target: multiple bed files (one per chromosome set) with set name as prefix using target species
    target_beds_ch = splitBedByChromosomeSet.out.flatten()
	    .map { fname ->
            def set_name = fname.getBaseName().toString().split(".bed")[0]
            tuple(params.hal_file, params.target_species, set_name, fname)
        }
        .filter { _hal, _species, set_name, _bed -> set_name != 'unassigned' }

    // Combine both channels and extract MAFs in a single process call
    all_beds_ch = ancestral_bed_ch.mix(target_beds_ch)

    extractMafFromHal(all_beds_ch)

    // Split output: ancestral MAFs for phyloFit, target MAFs for rescaling
    extractMafFromHal.out
        .branch { row ->
            ancestral: row[0] == 'ancestral'
            target: row[0] != 'ancestral'
        }
        .set { maf_outputs }

    // Train base phyloFit model on ANCESTRAL genome
    runPhyloFit(
        maf_outputs.ancestral.map { _prefix, maf -> maf },
        tree_ch
    )

    // Prepare target MAFs for rescaling (already have set names as prefix)
    maf_outputs.target.view()

    // Rescale phyloFit model for each chromosome set on TARGET genome
    // Use .first() to convert queue channel to value channel that can be reused
    rescaleByChromosome(maf_outputs.target, runPhyloFit.out.first())
}
