# figures/figure6/common.sh -- paths and interpreters for the Figure 6 panel scripts (sourced by each).
#
# The Figure 6 analyses live in comparative_genomics/phyloP_matched/loci_focus (LOCI); their tables and plots are written under RESULTS, and the
# alignment extracts under MAF_DATA. Each panel script runs the analysis code for its panel and copies what it draws to FIG6_OUT
# (default figures/figure6/output/, not tracked). Every variable below can be set in the environment instead.
set -euo pipefail
REPO=$(git -C "$(dirname "${BASH_SOURCE[0]}")" rev-parse --show-toplevel)
LOCI=${FIG6_LOCI:-$REPO/comparative_genomics/phyloP_matched/loci_focus}
RESULTS=${FIG6_RESULTS:-/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus/results}
MAF_DATA=${FIG6_MAF_DATA:-/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus/maf_switch}
FIG6_OUT=${FIG6_OUT:-$REPO/figures/figure6/output}
PY=${FIG6_PY:-$HOME/miniforge3/envs/hic/bin/python}              # numpy, pandas, scipy, pysam, pyBigWig, matplotlib
PY_ALN=${FIG6_PY_ALN:-$HOME/miniforge3/envs/crested/bin/python}  # as PY, plus lxml (SVG clean-up) and FIMO in the same env
RSCRIPT=${FIG6_RSCRIPT:-$HOME/miniforge3/envs/r_env/bin/Rscript}  # data.table, ggplot2, ComplexHeatmap, circlize, cowplot
mkdir -p "$FIG6_OUT"
# copy <src> to FIG6_OUT/<name>, failing loudly if the analysis did not write it
put() { [ -s "$1" ] || { echo "missing output: $1" >&2; exit 1; }; cp "$1" "$FIG6_OUT/$2"; echo "  -> $FIG6_OUT/$2"; }
