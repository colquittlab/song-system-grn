#!/usr/bin/env bash
# Figure 6 -- the analysis tables the panels draw from, rebuilt in order (results/mafb_domain/, maf_switch/data/). Run once before the panels,
# or after any change to the curation files (loci_focus/*_pvalb1.bed, excluded_*.bed, manual_lobes.bed, lobe_count.bed).
#   1. element_profiles.py           DAR elements of the MAFB domain (PVALB-1 comparator), lobes, accessibility-weighted phyloP, alignment depth
#   2. element_rule_classes.py       per-element summaries the class step reads
#   3. element_threshold_classes.py  phyloP classes (avian 0.52 / oscine 0.27 / strong 3.0)
#   4. element_ar_sites.py           FIMO scan for the AR half-site (AREh, panel F)
#   5. hal2maf                       the avian363 alignment (with ancestors) over every element +/- 100 bp (~1 h; skipped if present)
#   6. lobe_site_gain.py             regulator motif sites, presence across species, oscine gain and placement (panel F bold)
#   7. lobe_gain_null.py (optional)  shuffled / background-motif expectations (not drawn in the figure; FIG6_NULLS=1)
#   bash figures/figure6/make_inputs.sh
source "$(dirname "$0")/common.sh"
cd "$LOCI"
"$PY" element_profiles.py domain
"$PY" element_rule_classes.py
"$PY" element_threshold_classes.py
"$PY" element_ar_sites.py
MAF=$MAF_DATA/data/lobes_gap_pvalb1.maf
if [ ! -s "$MAF" ]; then
  HAL=${FIG6_HAL:-/private/warm-archive/colquittlab/assembly/cactus/avian363/363-avian-2020_lsd.hal}
  export PATH=${FIG6_CACTUS_BIN:-/private/groups/colquittlab/opt/cactus-bin-v3.1.2/bin}:$PATH
  tail -n +2 "$RESULTS/mafb_domain/element_threshold_classes.tsv" | awk -v OFS='\t' '{print $3,$4-100,$5+100}' | sort -k2,2n \
    | bedtools merge > "$MAF_DATA/data/lobes_windows_pvalb1.bed"
  hal2maf --refGenome Lonchura_striata --refTargets "$MAF_DATA/data/lobes_windows_pvalb1.bed" --maxRefGap 100 "$HAL" "$MAF"
fi
# the maf_switch scripts read data/ relative to the working directory; put the repository's modules first on the path
cd "$MAF_DATA"
export PYTHONPATH=$LOCI/maf_switch:$LOCI${PYTHONPATH:+:$PYTHONPATH}
LOBE_MAF=data/lobes_gap_pvalb1.maf "$PY" "$LOCI/maf_switch/lobe_site_gain.py"
if [ "${FIG6_NULLS:-0}" = 1 ]; then
  LOBE_MAF=data/lobes_gap_pvalb1.maf "$PY" "$LOCI/maf_switch/lobe_gain_null.py" scan
  "$PY" "$LOCI/maf_switch/lobe_gain_null_report.py"
fi
