#!/bin/bash
#SBATCH --job-name=create_test_hal
#SBATCH --output=logs/create_test_hal_%j.log
#SBATCH --error=logs/create_test_hal_%j.log
#SBATCH --time=2:00:00
#SBATCH --partition=medium
#SBATCH --cpus-per-task=4
#SBATCH --mem=100G
#
# Create a subsetted test HAL file from a larger HAL alignment
#
# Usage:
#   Direct: ./create_test_hal.sh <hal_file> <output_hal> [scale_factor]
#   SLURM:  sbatch create_test_hal.sh <hal_file> <output_hal> [scale_factor]
#
# This script uses halLodExtract to create a downsampled HAL file
# suitable for pipeline testing.

set -euo pipefail

# Activate cactus environment
source /private/home/bcolquit/.bashrc
mamba activate cactus

if [[ $# -lt 2 ]]; then
    echo "Usage: $0 <hal_file> <output_hal> [scale_factor]"
    echo ""
    echo "Arguments:"
    echo "  hal_file      Path to input HAL alignment file"
    echo "  output_hal    Path for output test HAL file"
    echo "  scale_factor  Downsampling scale factor (default: 100)"
    echo "                Higher values = smaller/faster test file"
    echo ""
    echo "Examples:"
    echo "  ./create_test_hal.sh birds.hal test.hal"
    echo "  ./create_test_hal.sh birds.hal test.hal 50   # less downsampling"
    echo "  ./create_test_hal.sh birds.hal test.hal 200  # more downsampling"
    exit 1
fi

HAL_FILE="$1"
OUTPUT_HAL="$2"
SCALE_FACTOR="${3:-100}"

# Validate inputs
if [[ ! -f "$HAL_FILE" ]]; then
    echo "Error: HAL file not found: $HAL_FILE"
    exit 1
fi

# Create output directory if needed
OUTPUT_DIR=$(dirname "$OUTPUT_HAL")
mkdir -p "$OUTPUT_DIR"

echo "=== Job Info ==="
echo "Job ID: ${SLURM_JOB_ID:-N/A}"
echo "Node: ${SLURMD_NODENAME:-$(hostname)}"
echo "Start time: $(date)"
echo ""

echo "=== Parameters ==="
echo "HAL file: $HAL_FILE"
echo "Output HAL: $OUTPUT_HAL"
echo "Scale factor: $SCALE_FACTOR"
echo ""

#echo "=== Input HAL stats ==="
#halStats "$HAL_FILE" | head -10
#echo ""

echo "=== Creating downsampled HAL ==="
halLodExtract \
    --allSequences \
    "$HAL_FILE" \
    "$OUTPUT_HAL" \
    "$SCALE_FACTOR"

echo ""
echo "=== Output HAL stats ==="
halStats "$OUTPUT_HAL" | head -30

echo ""
echo "=== Done ==="
echo "End time: $(date)"
echo "Test HAL written to: $OUTPUT_HAL"
echo "Input size: $(du -h "$HAL_FILE" | cut -f1)"
echo "Output size: $(du -h "$OUTPUT_HAL" | cut -f1)"
