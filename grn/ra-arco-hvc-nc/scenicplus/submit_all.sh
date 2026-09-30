#!/usr/bin/env bash
# Submit the SCENIC+ Snakemake workflow for every configN/ next to this script (one SLURM job each).
#
#   ./submit_all.sh            # all configN
#   ./submit_all.sh 1 3        # only config1 and config3
#   DRY_RUN=1 ./submit_all.sh  # print what would be submitted
#
# Each job runs the whole workflow inside its own allocation (see configN/run_snakemake.sbatch for the
# partition/memory/time). Two configs at once means two full allocations -- 2 x 40 cpus / 600G on the
# medium partition with the current settings -- running concurrently.
#
# Run from anywhere on prism; paths resolve from this script's own location. SLURM opens the
# #SBATCH --output file before the job body runs, so logs/ is created here, at submit time.
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"  # -P: resolve a symlinked scenicplus/ dir, or find() sees nothing

if [ "$#" -gt 0 ]; then
    cfgs=()
    for n in "$@"; do cfgs+=("$here/config$n"); done
else
    mapfile -t cfgs < <(find "$here" -maxdepth 1 -type d -name 'config[0-9]*' | sort -V)
fi
[ "${#cfgs[@]}" -gt 0 ] || { echo "no configN/ directories found in $here" >&2; exit 1; }

for cfg in "${cfgs[@]}"; do
    name="$(basename "$cfg")"
    if [ ! -f "$cfg/Snakemake/config/config.yaml" ] || [ ! -f "$cfg/run_snakemake.sbatch" ]; then
        echo "skip $name: missing Snakemake/config/config.yaml or run_snakemake.sbatch" >&2
        continue
    fi
    if [ -n "${DRY_RUN:-}" ]; then
        echo "[dry run] (cd $cfg/Snakemake && mkdir -p logs && sbatch ../run_snakemake.sbatch)"
        continue
    fi
    echo -n "$name: "
    (cd "$cfg/Snakemake" && mkdir -p logs && sbatch ../run_snakemake.sbatch)
done
