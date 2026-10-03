#!/usr/bin/env bash
# Assemble scplusmdata_slim.h5mu for every config in results/ that has its small outputs but no merged file yet.
#   ./assemble_slim_all.sh            # all such configs
#   ./assemble_slim_all.sh 16 17      # only these
# Needs AUCell_direct.h5mu, AUCell_extended.h5mu, eRegulon_direct.tsv and eRegulons_extended.tsv in configN/outs.
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
RES=/hdd/jupyter/brad/scenicplus/motor-pathway_multiome/motor-pathway_multiome_seurat_cellbender.0.05_preprocess_cr/ra-arco-hvc-nc_hybrid/results
PY=/home/brad/micromamba/envs/scenicplus/bin/python

if [ "$#" -gt 0 ]; then nums=("$@"); else mapfile -t nums < <(ls "$RES" | sed -n 's/^config//p' | sort -n); fi

for n in "${nums[@]}"; do
    outs="$RES/config$n/outs"
    if [ -e "$outs/scplusmdata.h5mu" ] || [ -e "$outs/scplusmdata_slim.h5mu" ]; then echo "config$n: already has a merged file"; continue; fi
    if [ ! -e "$outs/AUCell_extended.h5mu" ] || [ ! -e "$outs/eRegulons_extended.tsv" ]; then echo "config$n: small outputs not here yet, skipped"; continue; fi
    echo "config$n: assembling"
    "$PY" "$here/assemble_slim_scplusmdata.py" "$outs" 2>&1 | grep -v "Warning\|warnings.warn" | tail -1
done
