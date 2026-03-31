#! /bin/bash

HAL_FILE="$1"
SPECIES="$2"
RESTART="${3:-}"
OUTDIR="maf"

RESTART_FLAG=""
if [[ "$RESTART" == "restart" ]]; then
    RESTART_FLAG="--restart"
fi

mkdir -p $OUTDIR

cactus-hal2maf \
    --chunkSize 100000 \
    --batchCores 8 \
    --batchCount 50 \
    --batchMemory 256G \
    --noAncestors \
    --filterGapCausingDupes \
    --outType single \
    --batchSystem slurm \
    --maxLocalJobs 800 \
    --slurmTime 200:00:00 \
    --slurmPartition long \
    --logFile "$OUTDIR/cactus-hal2maf.log" \
    --refGenome "$SPECIES" \
    $RESTART_FLAG \
    "$OUTDIR/js" \
    "$HAL_FILE" \
    "$OUTDIR/genome.maf"
