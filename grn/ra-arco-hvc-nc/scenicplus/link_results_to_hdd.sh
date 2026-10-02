#!/usr/bin/env bash
# Keep SCENIC+ results off /ssd: make configN/outs and configN/Snakemake/logs symlinks into /hdd storage.
# Local-machine utility (the /hdd path below is this machine's); prism does not use it.
#
#   ./link_results_to_hdd.sh 3 4 5      # specific configs
#   ./link_results_to_hdd.sh            # every configN/ next to this script
#
# Run it BEFORE results are transferred in, so they land on /hdd. Each config needs ~80 GB. Transfers
# that write into this directory otherwise create real directories on whichever disk the repo is on.
# If outs/ is already a real directory it is copied to /hdd, verified (rsync dry run must list nothing),
# then removed and replaced by the symlink. A directory with a *.partial file is skipped, because a
# transfer is still writing into it -- rerun once that finishes.
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
STORE="${STORE:-/hdd/jupyter/brad/scenicplus/motor-pathway_multiome/motor-pathway_multiome_seurat_cellbender.0.05_preprocess_cr/ra-arco-hvc-nc_hybrid/results}"

if [ "$#" -gt 0 ]; then
    nums=("$@")
else
    mapfile -t nums < <(find "$here" -maxdepth 1 -type d -name 'config[0-9]*' | sed 's#.*/config##' | sort -n)
fi

for n in "${nums[@]}"; do
    for sub in outs Snakemake/logs; do
        src="$here/config$n/$sub"
        dst="$STORE/config$n/$(basename "$sub")"
        mkdir -p "$dst"
        if [ -L "$src" ]; then echo "config$n/$sub: already a symlink"; continue; fi
        if [ -d "$src" ]; then
            if find "$src" -name '*.partial' | grep -q .; then
                echo "config$n/$sub: has a .partial file (transfer in flight), skipped" >&2; continue
            fi
            echo "$(date +%T) config$n/$sub: copying $(du -sh --one-file-system "$src" | cut -f1)"
            rsync -a "$src/" "$dst/"
            left="$(rsync -an --itemize-changes "$src/" "$dst/" | wc -l)"
            if [ "$left" -ne 0 ]; then echo "VERIFY FAILED for $src ($left differences); source kept" >&2; exit 1; fi
            rm -r "$src"
        fi
        ln -s "$dst" "$src"
        echo "config$n/$sub -> $dst"
    done
done
