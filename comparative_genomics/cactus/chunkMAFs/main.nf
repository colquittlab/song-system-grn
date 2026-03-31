#! /usr/bin/env nextflow

process extractGenomeFasta {

    executor 'local'

    input:
        path hal_file
        val species

    output:
        path "genome.fasta"

    script:
    """
    hal2fasta ${hal_file} ${species} > genome.fasta
    """
}

process splitFastaByChromosome {

    executor 'local'

    input:
        path fasta_file

    output:
        path "fasta_by_chrom/*.fa"

    script:
    """
    mkdir fasta_by_chrom
    faSplit byname ${fasta_file} fasta_by_chrom/
    """
}

process splitMAFByChromosome {

    executor 'slurm'
    queue 'medium'
    time '12:00:00'
    cpus 1
    memory '16G'

    input:
        path maf_file

    output:
        path "bychrom_*.maf"

    script:
    """
    python3 << 'PYEOF'
outfiles = {}
header_lines = []
current_block = []
ref_chrom = None

with open("${maf_file}") as f:
    for line in f:
        if line.startswith('#'):
            if not current_block:
                header_lines.append(line)
            continue
        if line.strip() == '':
            if current_block and ref_chrom:
                fname = f'bychrom_{ref_chrom}.maf'
                if ref_chrom not in outfiles:
                    outfiles[ref_chrom] = open(fname, 'w')
                    outfiles[ref_chrom].writelines(header_lines)
                outfiles[ref_chrom].write(''.join(current_block))
                outfiles[ref_chrom].write('\\n')
            current_block = []
            ref_chrom = None
            continue
        if line.startswith('a'):
            current_block = [line]
            ref_chrom = None
        elif line.startswith('s') and ref_chrom is None:
            src = line.split()[1]
            ref_chrom = src.split('.')[-1] if '.' in src else src
            current_block.append(line)
        else:
            current_block.append(line)

# Handle last block (no trailing newline)
if current_block and ref_chrom:
    fname = f'bychrom_{ref_chrom}.maf'
    if ref_chrom not in outfiles:
        outfiles[ref_chrom] = open(fname, 'w')
        outfiles[ref_chrom].writelines(header_lines)
    outfiles[ref_chrom].write(''.join(current_block))
    outfiles[ref_chrom].write('\\n')

for f in outfiles.values():
    f.close()
PYEOF
    """
}

process chunkMAFs {

    publishDir(
        path: "${params.output_dir}",
        mode: 'move',
        saveAs: { filename -> "${chrom}/${filename}" }
    )

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

    # Generate windows
    bedtools makewindows -g ${refseq_file}.fai -w ${params.chunk_size} > windows.bed

    # Split MAF by windows
    mafSplit windows.bed chunk_${root}_ ${maf_file}
    """
}

workflow {

    // Extract genome fasta
    extractGenomeFasta(params.hal_file, params.target_species)

    if (params.maf_file) {
        // Single MAF: pass directly to chunkMAFs with full genome fasta
        maf_fasta_ch = channel.fromPath(params.maf_file)
            .combine(extractGenomeFasta.out)
            .map { maf_file, fasta_file ->
                def chrom = maf_file.baseName.replaceAll(/\.maf$/, '')
                tuple(chrom, maf_file, fasta_file)
            }
    } else {
        // Multiple MAFs from directory: split genome fasta by chromosome and combine
        splitFastaByChromosome(extractGenomeFasta.out)
        fasta_ch = splitFastaByChromosome.out.flatten()
            .map { file ->
                def chrom = file.baseName.replaceAll(".fa", "")
                tuple(chrom, file)
            }

        mafs_ch = channel.fromPath("${params.maf_dir}/*.maf")
            .map { file ->
                def chrom = file.baseName.replaceAll(".maf", "")
                tuple(chrom, file)
            }

        if (params.chromosomes) {
            def chrom_list = params.chromosomes.split(',').collect { it.trim() }
            mafs_ch = mafs_ch.filter { chrom, _file -> chrom_list.contains(chrom) }
        }

        maf_fasta_ch = mafs_ch.combine(fasta_ch, by: 0)
    }

    // Chunk MAFs
    chunkMAFs(maf_fasta_ch)
}
