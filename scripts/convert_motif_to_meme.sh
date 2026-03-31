#!/usr/bin/env bash
# Usage: convert_motif_to_meme.sh [DIR]
#
# Recursively convert HOMER .motif and .cb files to MEME format using chen2meme.
# Skips files with "similar", "RV", or "known" in the path (reverse-complement
# logos, redundant motifs, and known-motif files).
#
# Arguments:
#   DIR   Directory to search (default: current directory)
#
# Output:
#   For each foo.motif or foo.cb a sibling file foo.meme is written in place.
#
# Requirements:
#   chen2meme (part of the MEME Suite) must be on PATH.
#
# Example:
#   cd peaks_homer/homer_dars
#   bash ~/ssd/rstudio/multiome/motor-pathway/scripts/convert_motif_to_meme.sh
#   bash ~/ssd/rstudio/multiome/motor-pathway/scripts/convert_motif_to_meme.sh ra_pos_vs_neg/

DIR="${1:-.}"

for f in $(find "$DIR" \( -name "*.motif" -o -name "*.cb" \) | grep -v "similar" | grep -v "RV" | grep -v "known"); do
  ext="${f##*.}"
  chen2meme < "$f" > "${f%.$ext}.meme"
done
