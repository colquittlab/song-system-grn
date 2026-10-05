#!/usr/bin/env bash
# Run ONCE on prism, before submitting the sweep configs (16-28): make shared/outs/ from config1's motif-independent
# outputs, so every sweep config reuses them instead of redoing the 3.5 h region-to-gene fit and the 35 GB
# accessibility+expression object.
#
#   ./prepare_shared_upstream.sh
#
# Shared: ACC_GEX.h5mu (depends only on the cisTopic object + expression), search_space.tsv (config1's search space,
# which the sweep keeps), region_to_gene_adj.tsv (depends only on those two). Not shared: motif enrichment, TF-to-gene
# (its TF list comes from the cistromes), eRegulons and activity scores, which differ per config.
#
# Hardlinks when config1 is on the same filesystem (no extra 35 GB), a copy otherwise. Files are made read-only, so a
# Snakemake rerun that decided to regenerate one fails loudly instead of overwriting config1's results.
# If config1/outs no longer has them on prism, restore them (the same three files are on /hdd under
# .../ra-arco-hvc-nc_hybrid/results/config1/outs/) or rerun config1's first rules; don't let a sweep job recompute them.
set -euo pipefail

ROOT=/private/groups/colquittlab/scenicplus/motor-pathway_multiome/motor-pathway_multiome_seurat_cellbender.0.05_preprocess_cr/ra-arco-hvc-nc_hybrid/scenicplus
mkdir -p "$ROOT/shared/outs"

for f in ACC_GEX.h5mu search_space.tsv region_to_gene_adj.tsv; do
    src="$ROOT/config1/outs/$f"
    dst="$ROOT/shared/outs/$f"
    if [ -e "$dst" ]; then echo "$f: already in shared/outs"; continue; fi
    [ -s "$src" ] || { echo "missing $src -- restore it before submitting the sweep (see header)" >&2; exit 1; }
    if ln "$src" "$dst" 2>/dev/null; then how=hardlinked; else cp -p "$src" "$dst"; how=copied; fi
    chmod a-w "$dst"
    echo "$f: $how ($(du -h "$dst" | cut -f1))"
done
ls -l "$ROOT/shared/outs"
