#!/usr/bin/env bash
# Score one part of the motif list for the ra-arco-hvc-nc_hybrid cisTarget database (create_cistarget_motif_databases.py -p).
#
#   CT_ROOT=<database dir> ./score_part.sh PART NPARTS [THREADS]
#
# CT_ROOT holds lonStrDom2_1kb_bg_padding.fa, motifs.txt, singletons/ (the .cb motif files) and tools/ (the
# create_cisTarget_databases scripts and a static cbust binary), as copied by the transfer in README.md. Partial score databases
# go to $CT_ROOT/partial/. A part that already finished (its .done marker exists) is skipped, so a requeued array task resumes.
# CT_FASTA / CT_MOTIFS override the inputs (used for the local smoke test).
set -euo pipefail

ROOT="${CT_ROOT:?set CT_ROOT to the database directory}"
PART="${1:?part}"; NPARTS="${2:?number of parts}"; THREADS="${3:-40}"
PREFIX="${CT_PREFIX:-ra-arco-hvc-nc_hybrid}"
FASTA="${CT_FASTA:-$ROOT/lonStrDom2_1kb_bg_padding.fa}"
MOTIFS="${CT_MOTIFS:-$ROOT/motifs.txt}"

# Fail in seconds, not after queueing, if the active python lacks something the scripts import. tools/ carries its own copy of
# flatbuffers (pure python; prism's scenicplus env does not have it); numpy, pandas, pyarrow and numba must come from the env.
(cd "$ROOT/tools" && python -c "import numpy, pandas, pyarrow, numba, cistarget_db") \
    || { echo "python cannot import what the cisTarget scripts need (numpy, pandas, pyarrow, numba, tools/); is the scenicplus env active?" >&2; exit 1; }

mkdir -p "$ROOT/partial"
marker="$ROOT/partial/part_${PART}_of_${NPARTS}.done"
if [ -e "$marker" ]; then echo "part $PART/$NPARTS already done"; exit 0; fi

echo "$(date +%T) scoring part $PART of $NPARTS with $THREADS threads"
python "$ROOT/tools/create_cistarget_motif_databases.py" \
    -f "$FASTA" \
    -M "$ROOT/singletons" \
    -m "$MOTIFS" \
    -p "$PART" "$NPARTS" \
    -o "$ROOT/partial/$PREFIX" \
    --bgpadding 1000 \
    --cbust "$ROOT/tools/cbust" \
    -t "$THREADS"
touch "$marker"
echo "$(date +%T) part $PART/$NPARTS done"
