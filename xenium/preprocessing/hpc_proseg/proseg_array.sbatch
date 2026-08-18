#!/bin/bash
#SBATCH --job-name=proseg
#SBATCH --array=1-10
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --output=logs/proseg_%A_%a.out
#SBATCH --error=logs/proseg_%A_%a.err
# NOTE on the two lines above: SLURM opens these files before this script's own
# body runs, so `logs/` (relative to wherever `sbatch` is invoked from) must
# already exist at submit time — nothing this script does can create it in
# time for its own job. It ships pre-created (with a `.keep` placeholder, since
# rclone/rsync can drop empty directories) by 00_make_manifest.R. If submission
# still fails with "couldn't open output file", run `mkdir -p logs` yourself
# before `sbatch` — most likely the sync dropped it despite the placeholder.
# ---------------------------------------------------------------------------
# EDIT THESE for your cluster:
#SBATCH --account=CHANGE_ME
#SBATCH --partition=CHANGE_ME
# ---------------------------------------------------------------------------
#
# One array task per SECTION — this is the coarse, embarrassingly-parallel cut
# (ten independent tissue pieces, no spatial dependency between them).
#
# Do NOT also try to spatially tile within a section. proseg already does that
# internally: it partitions space into quads (--quad-size, default 150 um) and
# updates non-adjacent quads concurrently (graph coloring), targeting
# --cells-per-chunk (default 100) cells per chunk. On the 165,094-cell section
# measured locally that is ~1,650 independent chunks — far more parallelism
# than -t exploited at 22 cores. That internal scheme is what -t parallelizes;
# hand-tiling a section would cut across the sampler's spatial dependencies and
# is not how the tool is meant to be split.
#
# So: give each task as many cores as one node provides (32 here — CHANGE to
# match your cluster's node size; the local run used only 22 and 1,650 chunks
# were nowhere near saturating that many workers, so a full node should help).
#
# Local reference point: OR52YW26_2_2, 96.7M transcripts, 165,094 cells,
# 22 threads -> ~45 min wall. --time=02:00:00 leaves headroom for slower nodes
# or sections with more transcripts; check `sacct` after the first run and
# tighten it.
#
# Submit (everything below lives in one self-contained directory — built by
# 00_make_manifest.R, moved as a single tree by rclone, no separate paths to
# keep in sync):
#   Rscript 00_make_manifest.R          # locally, once — stages manifest.tsv,
#                                        # src/<section>/transcripts.parquet,
#                                        # and scripts/ under one directory
#   rclone sync <that dir> remote:hpc_proseg
#   # on the cluster:
#   rclone copy remote:hpc_proseg /scratch/$USER/hpc_proseg
#   cd /scratch/$USER/hpc_proseg
#   sbatch scripts/proseg_array.sbatch .
#
# Resume after partial failure: this script is idempotent (skips a section
# whose output already has cell-metadata.parquet), so a plain re-submit of the
# same --array range only redoes what's missing.

set -euo pipefail

BASE="${1:?usage: sbatch proseg_array.sbatch <hpc_proseg_dir>}"
BASE="$(cd "$BASE" && pwd)"       # absolute, so relative manifest paths resolve
MANIFEST="$BASE/manifest.tsv"

# Belt-and-suspenders for everything this script's own body can still control
# (does NOT cover logs/, which SLURM needs before this line ever runs — see the
# --output/--error note above): results/ for the mkdir below, and re-create
# logs/ here too so a resubmitted array recovers if it was ever removed after
# the first successful submission.
mkdir -p "$BASE/results" "$BASE/logs"

# --- environment ------------------------------------------------------------
# proseg itself is a standalone binary once built and needs nothing at runtime
# beyond it, but it was built against this env's toolchain (conda-forge rust +
# compilers), so activate it here too in case that toolchain provides shared
# libs (libstdc++ etc.) not otherwise on the node's library path.
if command -v micromamba >/dev/null 2>&1; then
    eval "$(micromamba shell hook --shell bash)"
    micromamba activate xenium-proseg
else
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate xenium-proseg
fi

export PATH="$HOME/.cargo/bin:$PATH"   # `cargo install` puts proseg here regardless
if ! command -v proseg >/dev/null 2>&1; then
    echo "proseg not found. See README.md step 0 (micromamba create -f scripts/environment.yml; cargo install proseg)." >&2
    exit 1
fi

# tail avoids the header row; row N of the manifest is array index N
LINE=$(tail -n +2 "$MANIFEST" | awk -v i="$SLURM_ARRAY_TASK_ID" 'NR==i')
SECTION_ID=$(echo "$LINE" | cut -f2)
# transcripts_path in the manifest is relative to BASE (e.g.
# src/OR52YW26_1_4/transcripts.parquet) so the whole tree stays portable —
# resolve it against BASE here rather than baking an absolute path in earlier.
TRANSCRIPTS="$BASE/$(echo "$LINE" | cut -f3)"

OUT_DIR="$BASE/results/$SECTION_ID"
mkdir -p "$OUT_DIR"

if [ -f "$OUT_DIR/cell-metadata.parquet" ]; then
    echo "section $SECTION_ID already done — exiting"
    exit 0
fi

echo "section $SECTION_ID : $TRANSCRIPTS"
cd "$OUT_DIR"

# --output-counts-fmt pinned explicitly: proseg ignores the file extension you
# give --output-expected-counts and always writes gzipped MatrixMarket unless
# told otherwise. Naming the file .csv.gz here matches what it actually is —
# earlier local runs discovered a mismatch when the file was named .parquet.
proseg --xenium "$TRANSCRIPTS" xenium \
    -t "$SLURM_CPUS_PER_TASK" \
    --output-expected-counts expected-counts.csv.gz \
    --output-counts-fmt csv-gz \
    --output-cell-metadata cell-metadata.parquet \
    --output-cell-metadata-fmt parquet \
    --output-transcript-metadata transcript-metadata.parquet \
    --output-transcript-metadata-fmt parquet \
    --output-cell-polygons cell-polygons.geojson.gz
