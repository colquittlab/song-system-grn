#!/bin/bash
#SBATCH --job-name=extract_maf
#SBATCH --output=logs/extract_maf_%j.log
#SBATCH --error=logs/extract_maf_%j.log
#SBATCH --time=24:00:00
#SBATCH --partition=long
#SBATCH --cpus-per-task=48
#SBATCH --mem=100G
#
# Extract and sort MAF files from a HAL alignment
#
# Usage:
#   Direct: ./extract_maf.sh <hal_file> <species> <output_dir> [cpus]
#   SLURM:  sbatch extract_maf.sh <hal_file> <species> <output_dir> [cpus]
#
# This script:
# 1. Extracts MAF for the entire genome using cactus-hal2maf
# 2. Sorts MAF by reference sequence order using taffy
# 3. Splits into per-chromosome MAF files
#
# Prerequisites:
# - Sufficient disk space for intermediate files

set -euo pipefail

# Activate cactus environment
source /private/home/bcolquit/.bashrc
mamba activate cactus

if [[ $# -lt 3 ]]; then
    echo "Usage: $0 <hal_file> <species> <output_dir> [cpus]"
    echo ""
    echo "Arguments:"
    echo "  hal_file    Path to HAL alignment file"
    echo "  species     Reference species name"
    echo "  output_dir  Directory for output MAF files"
    echo "  cpus        Number of CPUs for cactus-hal2maf (default: 48)"
    echo ""
    echo "Examples:"
    echo "  Direct run:"
    echo "    ./extract_maf.sh birds.hal Lonchura_striata ./maf"
    echo ""
    echo "  Submit to SLURM:"
    echo "    sbatch extract_maf.sh birds.hal Lonchura_striata ./maf"
    exit 1
fi

HAL_FILE="$1"
SPECIES="$2"
OUTPUT_DIR="$3"
CPUS="${4:-48}"

# Validate inputs
if [[ ! -f "$HAL_FILE" ]]; then
    echo "Error: HAL file not found: $HAL_FILE"
    exit 1
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Create temporary working directory
WORK_DIR=$(mktemp -d)
echo "Working directory: $WORK_DIR"

cleanup() {
    echo "Cleaning up temporary files..."
    rm -rf "$WORK_DIR"
}
trap cleanup EXIT

echo "=== Job Info ==="
echo "Job ID: ${SLURM_JOB_ID:-N/A}"
echo "Node: ${SLURMD_NODENAME:-$(hostname)}"
echo "Start time: $(date)"
echo ""

echo "=== Step 1: Extracting MAF from HAL ==="
echo "HAL file: $HAL_FILE"
echo "Species: $SPECIES"
echo "CPUs: $CPUS"

cactus-hal2maf \
    --chunkSize 100000 \\
        --batchCores 96 \\
        --batchCount 10 \\
        --noAncestors \\
        --filterGapCausingDupes \\
        --outType single \\
        --batchSystem slurm \\
        --maxLocalJobs 800 \\
        --slurmTime 200:00:00 \\
        --slurmPartition long \\
	--logFile avian.maf.gz.log \\
    --refGenome "$SPECIES" \
    "/js" \
    "$HAL_FILE" \
    "$WORK_DIR/genome.maf"

echo "=== Step 2: Getting chromosome order ==="
halStats --chromSizes "$SPECIES" "$HAL_FILE" | cut -f1 > "$WORK_DIR/genome.list"
echo "Found $(wc -l < "$WORK_DIR/genome.list") sequences"

echo "=== Step 3: Sorting MAF by reference sequence order ==="
taffy view "$WORK_DIR/genome.maf" \
    | taffy sort -n "$WORK_DIR/genome.list" \
    | taffy view -m > "$WORK_DIR/genome.sorted.maf"
rm "$WORK_DIR/genome.maf"

# echo "=== Step 4: Splitting MAF by chromosome ==="
# msa_split "$WORK_DIR/genome.sorted.maf" \
#     --in-format MAF \
#     --by-index \
#     --out-root "$OUTPUT_DIR/chrom" \
#     --out-format MAF

# Rename files to use chromosome names instead of indices
# echo "=== Step 5: Renaming output files ==="
# cd "$OUTPUT_DIR"
# i=0
# while read -r chrom; do
#     if [[ -f "chrom.${i}.maf" ]]; then
#         mv "chrom.${i}.maf" "${chrom}.maf"
#         echo "  chrom.${i}.maf -> ${chrom}.maf"
#     fi
#     ((i++))
# done < "$WORK_DIR/genome.list"

echo ""
echo "=== Done ==="
echo "End time: $(date)"
echo "MAF files written to: $OUTPUT_DIR"
echo "Total files: $(ls -1 "$OUTPUT_DIR"/*.maf 2>/dev/null | wc -l)"
