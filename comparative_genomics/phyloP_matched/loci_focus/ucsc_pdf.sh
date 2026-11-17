#!/bin/bash
# Render a UCSC Genome Browser view of the Colquitt lab hub (lonStrDom2, GCF_005870125.1) to PDF.
#
#   UCSC_API_KEY_FILE=~/.ucsc_api_key ucsc_pdf.sh <region> <out.pdf> [width_px] [track=vis ...]
#
# region   e.g. chr20:12945900-12947200
# track=vis  trackDb names with full|pack|dense|squish|hide; with any of these, all other tracks are
#            hidden first (hideTracks=1), otherwise the hub's default visibility is used.
# The API key (UCSC account > Hub Development > API key) bypasses the browser's bot check; keep it in a
# private file, never in this repository. A fresh session is used per call, so settings do not carry over.
set -euo pipefail
REGION=$1; OUT=$2; PIX=${3:-1200}; shift 2; [ $# -gt 0 ] && [[ $1 =~ ^[0-9]+$ ]] && shift
KEY=$(cat "${UCSC_API_KEY_FILE:?set UCSC_API_KEY_FILE to a file holding the UCSC API key}")
HUB="https://colquitt-track-hub.pbsci.ucsc.edu/lonStrDom2/hub.txt"
Q="hubUrl=$HUB&genome=GCF_005870125.1&position=$REGION&pix=$PIX&hgt.psOutput=on&apiKey=$KEY"
if [ $# -gt 0 ]; then Q="$Q&hideTracks=1"; for tv in "$@"; do Q="$Q&${tv}"; done; fi
TMP=$(mktemp -d); trap 'rm -rf "$TMP"' EXIT
curl -s -L --max-time 300 -c "$TMP/cj" -b "$TMP/cj" "https://genome.ucsc.edu/cgi-bin/hgTracks?$Q" -o "$TMP/page.html"
P=$(grep -oE "\.\./trash/hgt/hgt_genome_[^'\"]*\.pdf" "$TMP/page.html" | head -1 || true)
[ -n "$P" ] || { echo "no PDF link returned (check the key, region and track names)" >&2; exit 1; }
curl -s --max-time 120 -b "$TMP/cj" "https://genome.ucsc.edu/${P#../}" -o "$OUT"
echo "wrote $OUT ($(stat -c %s "$OUT") bytes)"
