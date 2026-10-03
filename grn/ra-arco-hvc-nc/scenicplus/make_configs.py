#!/usr/bin/env python3
"""Generate the per-config SCENIC+ Snakemake directories for ra-arco-hvc-nc_hybrid.

Each configN/ mirrors the layout used on prism (and in the earlier runs):

    configN/Snakemake/config/config.yaml    all paths absolute, under PRISM_ROOT
    configN/Snakemake/workflow/Snakefile    copied from workflow/Snakefile (shared)
    configN/run_snakemake.sbatch            submit from configN/Snakemake

Only the parameters that differ between configs live in CONFIGS below; everything else comes
from BASE_PARAMS, so a sweep is a diff of two small dicts rather than of two 100-line YAMLs.
Re-run this script after editing; it overwrites configN/ but never touches outs/.

    python make_configs.py
"""
import math
import random
import shutil
from pathlib import Path

HERE = Path(__file__).resolve().parent
PRISM_ROOT = ("/private/groups/colquittlab/scenicplus/motor-pathway_multiome/"
              "motor-pathway_multiome_seurat_cellbender.0.05_preprocess_cr/ra-arco-hvc-nc_hybrid")
CISTARGET_DIR = "/private/groups/colquittlab/scenicplus/cistarget/ra-arco-hvc-nc_seurat-clustering"  # existing DB, unchanged
CISTARGET_PREFIX = "ra-arco-hvc-nc_seurat-clustering"
MOTIF_ANNOT = ("/private/groups/colquittlab/scenicplus/v10nr_clust_public/snapshots/"
               "motifs-v10-nr.hgnc-m0.00001-o0.0.tbl")
N_CPU = 40

# Parameter values of the earlier full-dataset run (old config3), relaxed relative to SCENIC+ defaults.
BASE_PARAMS = {
    "params_data_preparation": {
        "direct_annotation": "Direct_annot",
        "extended_annotation": "Orthology_annot",
        "search_space_upstream": "1000 150000",
        "search_space_downstream": "1000 150000",
        "search_space_extend_tss": "10 10",
    },
    "params_motif_enrichment": {
        "motif_similarity_fdr": 0.01,
        "annotations_to_use": "Direct_annot Orthology_annot Motif_similarity_annot",
        "dem_adj_pval_thr": 0.05,
        "dem_log2fc_thr": 1.0,
        "ctx_auc_threshold": 0.005,
        "ctx_nes_threshold": 3.0,
        "ctx_rank_threshold": 0.2,
    },
}

CONFIGS = {
    # config1: same thresholds as the earlier full-dataset run, so old vs new is a labels-only comparison
    1: {},
    # config2: SCENIC+ defaults for the motif-similarity FDR and cisTarget rank threshold
    2: {"params_motif_enrichment": {"motif_similarity_fdr": 0.001, "ctx_rank_threshold": 0.05}},

    # config3-11 sweep around config1. The two sets run so far fix the motif-similarity FDR and the cisTarget
    # rank threshold together, so 3-5 separate them; 6-9 vary one other parameter each; 10 changes the search
    # space; 11 stacks the relaxed values. Every value below appeared in an earlier config or as a commented
    # alternative in one (the previous sweep's own record was not kept), none is new.
    3: {"params_motif_enrichment": {"ctx_rank_threshold": 0.1}},
    4: {"params_motif_enrichment": {"ctx_rank_threshold": 0.05}},
    5: {"params_motif_enrichment": {"motif_similarity_fdr": 0.001}},
    6: {"params_motif_enrichment": {"ctx_nes_threshold": 2.0}},
    7: {"params_motif_enrichment": {"ctx_auc_threshold": 0.0025}},
    8: {"params_motif_enrichment": {"dem_adj_pval_thr": 0.1, "dem_log2fc_thr": 0.5}},
    9: {"params_data_preparation": {"extended_annotation": "Orthology_annot Motif_similarity_annot"}},
    # region_to_gene is rerun on a different search space, so this one is slower than the rest
    10: {"params_data_preparation": {"search_space_upstream": "1000 100000",
                                     "search_space_downstream": "1000 100000"}},
    11: {"params_motif_enrichment": {"motif_similarity_fdr": 0.001, "ctx_rank_threshold": 0.1,
                                     "ctx_nes_threshold": 2.0, "ctx_auc_threshold": 0.0025,
                                     "dem_adj_pval_thr": 0.1, "dem_log2fc_thr": 0.5}},

    # config12/13: CONTROLS, not part of the sweep. Same parameters as config1 / config2 but the expression input
    # is the legacy double-transformed matrix the earlier run used (see GEX_ANNDATA), with the new labels, cells and
    # region sets unchanged. Comparing 12 vs 1 and 13 vs 2 isolates the effect of the old .raw bug.
    12: {},
    13: {"params_motif_enrichment": {"motif_similarity_fdr": 0.001, "ctx_rank_threshold": 0.05}},
}

# Expression input per config (file under anndata_rna/ on prism). Default is raw counts in .raw, which is what SCENIC+
# expects. adata_doublenorm.h5ad has .raw = log1p(normalize_total(SCT data)), made by anndata_rna/make_adata_doublenorm.py.
GEX_DEFAULT = "adata.h5ad"
GEX_ANNDATA = {12: "adata_doublenorm.h5ad", 13: "adata_doublenorm.h5ad"}

# config14/15: number of LDA topics (the hand-picked step in pycisTopic.ipynb; 20 everywhere else). Same parameters and
# raw-count input as config1. Each topic count has its own cisTopic object (the selected model drives the imputed
# accessibility SCENIC+ builds from it) and its own region sets (topic sets and the DARs, which are computed from that
# model's imputed accessibility). Both made from the already-trained models by pycisTopic/pycistopic_topic_k.py.
CONFIGS[14] = {}
CONFIGS[15] = {}
N_TOPICS = {14: 15, 15: 30}
N_TOPICS_DEFAULT = 20


def cistopic_obj_name(n):
    k = N_TOPICS.get(n)
    return f"cistopic_obj_glut_k{k}.pkl" if k else "cistopic_obj_glut.pkl"


def region_sets_name(n):
    if n in REGION_SET_VARIANT:
        return REGION_SET_VARIANT[n]
    k = N_TOPICS.get(n)
    return f"region_sets_k{k}" if k else "region_sets"


# config16-28: sweep of the motif-stage thresholds in the looser regime, designed to need little recomputation.
# AR, for example, only gets an eRegulon when these are looser than config1 (configs 6, 7, 8, 10, 11), so this maps that
# regime instead of stepping one parameter at a time. A seeded Latin hypercube over six parameters (ranges extend from
# config1's strict values toward looser ones), plus the loosest corner as an anchor. Each parameter's range:
#   ctx_nes_threshold      1.5 - 3.0   lower is looser (config1 3.0, config6 2.0)
#   ctx_auc_threshold      0.001 - 0.01, log  NOT directional (a smaller AUC window is not simply looser); brackets
#                          config7's 0.0025, config1's 0.005 and the 0.01 other configs used as "default"
#   ctx_rank_threshold     0.2 - 0.4   higher is looser (config1 0.2)
#   dem_adj_pval_thr       0.05 - 0.2  higher is looser (config8 0.1)
#   dem_log2fc_thr         0.25 - 1.0  lower is looser (config8 0.5)
#   motif_similarity_fdr   0.01 - 0.05, log  higher is looser (config1 0.01; config5's 0.001 was stricter)
# Everything else (annotations, eRegulon-assembly filters, search space, input matrix, topic count) stays at config1's.
SWEEP_SPECS = [  # (parameter, low, high, log scale, round to a multiple of / None = 2 significant digits)
    ("ctx_nes_threshold", 1.5, 3.0, False, 0.05),
    ("ctx_auc_threshold", 0.001, 0.01, True, None),
    ("ctx_rank_threshold", 0.2, 0.4, False, 0.01),
    ("dem_adj_pval_thr", 0.05, 0.2, False, 0.01),
    ("dem_log2fc_thr", 0.25, 1.0, False, 0.05),
    ("motif_similarity_fdr", 0.01, 0.05, True, None),
]
SWEEP_N, SWEEP_SEED, SWEEP_FIRST = 12, 2026, 16
SWEEP_CORNER = {"ctx_nes_threshold": 1.5, "ctx_auc_threshold": 0.0025, "ctx_rank_threshold": 0.4,
                "dem_adj_pval_thr": 0.2, "dem_log2fc_thr": 0.25, "motif_similarity_fdr": 0.05}


def latin_hypercube(n, d, seed):
    """n points in d dimensions, exactly one per stratum in every dimension (stdlib only, deterministic)."""
    rng = random.Random(seed)
    cols = []
    for _ in range(d):
        perm = list(range(n))
        rng.shuffle(perm)
        cols.append([(perm[i] + rng.random()) / n for i in range(n)])
    return [[cols[j][i] for j in range(d)] for i in range(n)]


def sweep_value(u, lo, hi, log, step):
    v = math.exp(math.log(lo) + u * (math.log(hi) - math.log(lo))) if log else lo + u * (hi - lo)
    return float(f"{v:.2g}") if step is None else round(round(v / step) * step, 4)


for _i, _u in enumerate(latin_hypercube(SWEEP_N, len(SWEEP_SPECS), SWEEP_SEED)):
    CONFIGS[SWEEP_FIRST + _i] = {"params_motif_enrichment": {
        name: sweep_value(u, lo, hi, log, step) for u, (name, lo, hi, log, step) in zip(_u, SWEEP_SPECS)}}
CONFIGS[SWEEP_FIRST + SWEEP_N] = {"params_motif_enrichment": dict(SWEEP_CORNER)}

# The sweep configs share config1's accessibility+expression object, search space and region-to-gene fit (the only slow,
# motif-independent steps; see shared/ on prism, made by prepare_shared_upstream.sh) and skip the 35 GB merged file.
# config29-33: does removing interneuron region sets from motif enrichment change the homeodomain eRegulons?
# For homeodomain TFs (ALX4, EMX2, LHX2, ...) the cistrome regions peak in GABA-LGE although the target genes peak where the
# TF is expressed, because the shared homeodomain motif is enriched in the GABA region sets. These remove those sets (folders
# from pycisTopic/make_region_set_variants.py). Each arm runs at config1's strict thresholds and at config11's loose ones
# (where AR, ALX4 and others appear); everything else, including the cisTopic object and expression, is config1's.
#   noGABA  all 9 GABA DAR sets + topics 1, 15, 20          noLGE  GABA-LGE-1 and -2 DAR sets + Topic20
#   ctrl    size-matched random removal of non-GABA, non-focal sets (at config1's thresholds only)
REGION_SET_VARIANT = {29: "region_sets_noGABA", 30: "region_sets_noGABA", 31: "region_sets_noLGE", 32: "region_sets_noLGE",
                      33: "region_sets_ctrl"}
CONFIGS[29] = {}
CONFIGS[30] = {k: dict(v) for k, v in CONFIGS[11].items()}
CONFIGS[31] = {}
CONFIGS[32] = {k: dict(v) for k, v in CONFIGS[11].items()}
CONFIGS[33] = {}

SHARED_UPSTREAM = set(range(SWEEP_FIRST, SWEEP_FIRST + SWEEP_N + 1)) | set(REGION_SET_VARIANT)


def merged(overrides):
    out = {k: dict(v) for k, v in BASE_PARAMS.items()}
    for section, vals in overrides.items():
        out[section].update(vals)
    return out


def q(v):
    return f'"{v}"' if isinstance(v, str) else str(v)


def render(n, overrides):
    p = merged(overrides)
    outs = f"{PRISM_ROOT}/scenicplus/config{n}/outs"
    up = f"{PRISM_ROOT}/scenicplus/shared/outs" if n in SHARED_UPSTREAM else outs  # shared, motif-independent outputs
    tmp = f"{PRISM_ROOT}/scenicplus/config{n}/tmp"
    py = f"{PRISM_ROOT}/pycisTopic"
    me, dp = p["params_motif_enrichment"], p["params_data_preparation"]
    lines = f"""input_data:
  cisTopic_obj_fname: "{py}/{cistopic_obj_name(n)}"
  GEX_anndata_fname: "{PRISM_ROOT}/anndata_rna/{GEX_ANNDATA.get(n, GEX_DEFAULT)}"
  region_set_folder: "{py}/{region_sets_name(n)}"
  ctx_db_fname: "{CISTARGET_DIR}/{CISTARGET_PREFIX}.regions_vs_motifs.rankings.feather"
  dem_db_fname: "{CISTARGET_DIR}/{CISTARGET_PREFIX}.regions_vs_motifs.scores.feather"
  path_to_motif_annotations: "{MOTIF_ANNOT}"

output_data:
  combined_GEX_ACC_mudata: "{up}/ACC_GEX.h5mu"
  dem_result_fname: "{outs}/dem_results.hdf5"
  ctx_result_fname: "{outs}/ctx_results.hdf5"
  output_fname_dem_html: "{outs}/dem_results.html"
  output_fname_ctx_html: "{outs}/ctx_results.html"
  cistromes_direct: "{outs}/cistromes_direct.h5ad"
  cistromes_extended: "{outs}/cistromes_extended.h5ad"
  tf_names: "{outs}/tfs.txt"
  genome_annotation: "{PRISM_ROOT}/scenicplus/genome_annotation.tsv"
  chromsizes: "{PRISM_ROOT}/scenicplus/chromsizes.tsv"
  search_space: "{up}/search_space.tsv"
  tf_to_gene_adjacencies: "{outs}/tf_to_gene_adj.tsv"
  region_to_gene_adjacencies: "{up}/region_to_gene_adj.tsv"
  eRegulons_direct: "{outs}/eRegulon_direct.tsv"
  eRegulons_extended: "{outs}/eRegulons_extended.tsv"
  AUCell_direct: "{outs}/AUCell_direct.h5mu"
  AUCell_extended: "{outs}/AUCell_extended.h5mu"
  scplus_mdata: "{outs}/scplusmdata.h5mu"

params_general:
  temp_dir: "{tmp}"
  n_cpu: {N_CPU}
  seed: 666{"" if n not in SHARED_UPSTREAM else chr(10) + "  build_scplus_mudata: False"}

params_data_preparation:
  bc_transform_func: "\\"lambda x: f'{{x}}'\\""
  is_multiome: True
  key_to_group_by: ""
  nr_cells_per_metacells: 10
  direct_annotation: {q(dp["direct_annotation"])}
  extended_annotation: {q(dp["extended_annotation"])}
  species: "hsapiens"
  biomart_host: "http://www.ensembl.org"
  search_space_upstream: {q(dp["search_space_upstream"])}
  search_space_downstream: {q(dp["search_space_downstream"])}
  search_space_extend_tss: {q(dp["search_space_extend_tss"])}

params_motif_enrichment:
  species: "homo_sapiens"
  annotation_version: "v10nr_clust"
  motif_similarity_fdr: {me["motif_similarity_fdr"]}
  orthologous_identity_threshold: 0.0
  annotations_to_use: {q(me["annotations_to_use"])}
  fraction_overlap_w_dem_database: 0.4
  dem_max_bg_regions: 500
  dem_balance_number_of_promoters: True
  dem_promoter_space: 1_000
  dem_adj_pval_thr: {me["dem_adj_pval_thr"]}
  dem_log2fc_thr: {me["dem_log2fc_thr"]}
  dem_mean_fg_thr: 0.0
  dem_motif_hit_thr: 3.0
  fraction_overlap_w_ctx_database: 0.4
  ctx_auc_threshold: {me["ctx_auc_threshold"]}
  ctx_nes_threshold: {me["ctx_nes_threshold"]}
  ctx_rank_threshold: {me["ctx_rank_threshold"]}

params_inference:
  tf_to_gene_importance_method: "GBM"
  region_to_gene_importance_method: "GBM"
  region_to_gene_correlation_method: "SR"
  order_regions_to_genes_by: "importance"
  order_TFs_to_genes_by: "importance"
  gsea_n_perm: 1000
  quantile_thresholds_region_to_gene: "0.85 0.90 0.95"
  top_n_regionTogenes_per_gene: "5 10 15"
  top_n_regionTogenes_per_region: ""
  min_regions_per_gene: 0
  rho_threshold: 0.05
  min_target_genes: 10
"""
    return lines


SBATCH = """#!/bin/bash
#SBATCH --job-name=scenicplus_hybrid_config{n}
#SBATCH --partition=medium
#SBATCH --cpus-per-task={cpus}
#SBATCH --mem=600G
#SBATCH --time={time}
#SBATCH --output=logs/scenicplus_config{n}_%j.out
#SBATCH --error=logs/scenicplus_config{n}_%j.err
# No --account: prism doesn't need one here. --mem/--time are set by hand, not measured; tf_to_gene and
# region_to_gene (the GBM steps) dominate wall time, so check `sacct` after the first run.
# The whole workflow runs inside this one allocation (--cores {cpus}); it is not a per-rule submit.
#
# Submit from configN/Snakemake (so relative config/ and workflow/ resolve), after `mkdir -p logs`:
#   cd {prism_cfg}/Snakemake && mkdir -p logs && sbatch ../run_snakemake.sbatch
set -euo pipefail
set +u   # conda's activate scripts read unset variables
source /private/home/${{USER}}/miniforge3/etc/profile.d/conda.sh
conda activate scenicplus
set -u
{precheck}mkdir -p {prism_cfg}/outs {prism_cfg}/tmp
snakemake --cores {cpus} --rerun-incomplete --printshellcmds
"""


SHARED_PRECHECK = """# These come from config1 via prepare_shared_upstream.sh; without them Snakemake would silently redo the 3.5 h
# region-to-gene fit and the 35 GB object, per config.
for f in ACC_GEX.h5mu search_space.tsv region_to_gene_adj.tsv; do
    test -s {shared}/$f || {{ echo "missing shared input {shared}/$f -- run prepare_shared_upstream.sh first" >&2; exit 1; }}
done
"""


def write_parameter_table():
    """config_parameters.tsv: one row per config, every swept parameter, so analyses can label configs."""
    rows = []
    for n, overrides in CONFIGS.items():
        p = merged(overrides)
        row = {"config": f"config{n}"}
        for section in BASE_PARAMS:
            row.update(p[section])
        row["gex_anndata"] = GEX_ANNDATA.get(n, GEX_DEFAULT)
        row["n_topics"] = N_TOPICS.get(n, N_TOPICS_DEFAULT)
        row["shared_upstream"] = "yes" if n in SHARED_UPSTREAM else "no"
        row["region_sets"] = region_sets_name(n)
        changed = [k for section in BASE_PARAMS for k, v in p[section].items()
                   if v != merged(CONFIGS[1])[section][k]]
        if row["gex_anndata"] != GEX_DEFAULT:
            changed.append("gex_anndata")
        if row["n_topics"] != N_TOPICS_DEFAULT:
            changed.append("n_topics")
        if n in REGION_SET_VARIANT:
            changed.append("region_sets")
        row["changed_from_config1"] = ",".join(changed) or "-"
        rows.append(row)
    cols = list(rows[0])
    with open(HERE / "config_parameters.tsv", "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")


def main():
    write_parameter_table()
    src_snakefile = HERE / "workflow" / "Snakefile"
    for n, overrides in CONFIGS.items():
        d = HERE / f"config{n}"
        (d / "Snakemake" / "config").mkdir(parents=True, exist_ok=True)
        (d / "Snakemake" / "workflow").mkdir(parents=True, exist_ok=True)
        (d / "Snakemake" / "config" / "config.yaml").write_text(render(n, overrides))
        shutil.copy(src_snakefile, d / "Snakemake" / "workflow" / "Snakefile")
        (d / "run_snakemake.sbatch").write_text(
            SBATCH.format(n=n, cpus=N_CPU, prism_cfg=f"{PRISM_ROOT}/scenicplus/config{n}",
                          time="08:00:00" if n in SHARED_UPSTREAM else "12:00:00",
                          precheck=SHARED_PRECHECK.format(shared=f"{PRISM_ROOT}/scenicplus/shared/outs") if n in SHARED_UPSTREAM else ""))
        print("wrote", d)


if __name__ == "__main__":
    main()
