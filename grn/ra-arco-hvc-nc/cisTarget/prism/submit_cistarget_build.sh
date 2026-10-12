#!/usr/bin/env bash
# Submit the cisTarget database build for ra-arco-hvc-nc_hybrid on prism: an array of scoring jobs, then one job that combines
# the parts and makes the rankings once all of them have succeeded.
#
#   cd /private/groups/colquittlab/scenicplus/cistarget/ra-arco-hvc-nc_hybrid
#   ./prism/submit_cistarget_build.sh [NPARTS]        # default 6
#   DRY_RUN=1 ./prism/submit_cistarget_build.sh       # print the sbatch commands only
#
# The database directory (this script's parent) must hold lonStrDom2_1kb_bg_padding.fa, motifs.txt, singletons/, tools/ and
# prism/ (see prism/README.md). The finished files are written into that same directory, which is where the SCENIC+
# configs 36-39 read them from.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
NPARTS="${1:-6}"

for f in lonStrDom2_1kb_bg_padding.fa motifs.txt tools/create_cistarget_motif_databases.py tools/cbust tools/feather_v1_fbs singletons; do
    [ -e "$ROOT/$f" ] || { echo "missing $ROOT/$f -- transfer the database directory first (prism/README.md)" >&2; exit 1; }
done
[ "$(ls "$ROOT/singletons" | wc -l)" -ge "$(wc -l < "$ROOT/motifs.txt")" ] || { echo "singletons/ has fewer files than motifs.txt lists" >&2; exit 1; }
chmod +x "$ROOT/tools/cbust"

cd "$ROOT"
mkdir -p logs
if [ -n "${DRY_RUN:-}" ]; then
    echo "[dry run] sbatch --array=1-$NPARTS --export=ALL,CT_ROOT=$ROOT,NPARTS=$NPARTS prism/score_parts.sbatch"
    echo "[dry run] sbatch --dependency=afterok:<array job> --export=ALL,CT_ROOT=$ROOT,NPARTS=$NPARTS prism/combine_and_rank.sbatch"
    exit 0
fi
jid="$(sbatch --parsable --array=1-"$NPARTS" --export=ALL,CT_ROOT="$ROOT",NPARTS="$NPARTS" prism/score_parts.sbatch)"
echo "scoring array: $jid (tasks 1-$NPARTS)"
cid="$(sbatch --parsable --dependency=afterok:"$jid" --export=ALL,CT_ROOT="$ROOT",NPARTS="$NPARTS" prism/combine_and_rank.sbatch)"
echo "combine + rankings: $cid (after $jid)"
