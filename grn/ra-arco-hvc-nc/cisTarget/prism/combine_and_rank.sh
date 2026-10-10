#!/usr/bin/env bash
# Combine the partial score databases and create the rankings database (steps 2 and 3 of the partial workflow in the
# create_cisTarget_databases README).
#
#   CT_ROOT=<database dir> ./combine_and_rank.sh NPARTS
#
# Needs every part's .done marker. Writes, into $CT_ROOT, the files SCENIC+ reads for configs 36-39:
#   <prefix>.regions_vs_motifs.rankings.feather   (ctx_db_fname)
#   <prefix>.regions_vs_motifs.scores.feather     (dem_db_fname)
# plus <prefix>.motifs_vs_regions.scores.feather, which the old database also kept. The partial files and the
# motifs_vs_regions rankings (not read by SCENIC+) are removed at the end, only after the files above are written.
set -euo pipefail

ROOT="${CT_ROOT:?set CT_ROOT to the database directory}"
NPARTS="${1:?number of parts}"
PREFIX="${CT_PREFIX:-ra-arco-hvc-nc_hybrid}"

for p in $(seq 1 "$NPARTS"); do
    [ -e "$ROOT/partial/part_${p}_of_${NPARTS}.done" ] || { echo "part $p of $NPARTS is not finished" >&2; exit 1; }
done

# Run inside the database directory with relative paths: the tools split database file names with a regex that reads any
# ".something" in a parent directory name (e.g. a ".claude" or a "x.0.05" component) as a species field and mangles the output.
cd "$ROOT"

echo "$(date +%T) combining $NPARTS partial score databases"
python tools/combine_partial_motifs_or_tracks_vs_regions_or_genes_scores_cistarget_dbs.py -i partial -o .

echo "$(date +%T) creating rankings"
python tools/convert_motifs_or_tracks_vs_regions_or_genes_scores_to_rankings_cistarget_dbs.py -i "$PREFIX.motifs_vs_regions.scores.feather"

for f in "$PREFIX.regions_vs_motifs.rankings.feather" "$PREFIX.regions_vs_motifs.scores.feather" "$PREFIX.motifs_vs_regions.scores.feather"; do
    [ -s "$ROOT/$f" ] || { echo "missing $f; partial files kept" >&2; exit 1; }
done
rm -f "$ROOT/$PREFIX.motifs_vs_regions.rankings.feather"
rm -r "$ROOT/partial"
echo "$(date +%T) done"
ls -lh "$ROOT"/*.feather
