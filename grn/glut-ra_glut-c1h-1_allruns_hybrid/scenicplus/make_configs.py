#!/usr/bin/env python3
"""Per-config SCENIC+ Snakemake directories for glut-ra_glut-c1h-1_allruns_hybrid (Glut-CACNA1H-RA vs Glut-CACNA1H-1, all runs (the hybrid labels of Glut-RA vs Glut-Arco-1, run3 included), hybrid labels).

Reuses the YAML/sbatch rendering and the shared Snakefile of ../../ra-arco-hvc-nc/scenicplus/make_configs.py (importing it only
defines its tables; nothing is written until main()), so the runs cannot drift apart. Only the root and the configs differ.
Derived from ../../ra-arco-hvc-nc_glut_hybrid/scenicplus/make_configs.py by ../../make_hybrid_subset_dirs.py.

    config1  the earlier run's parameters (glut-ra_glut-arco-1 (config18)): motif_similarity_fdr 0.001, dem_adj_pval_thr 0.05, dem_log2fc_thr 1.0, ctx_auc_threshold 0.005, ctx_nes_threshold 3.0, ctx_rank_threshold 0.05; extended_annotation 'Orthology_annot Motif_similarity_annot'.
             10 topics, raw counts in .raw.
    config2  the all-cell hybrid run's strict config1 thresholds (base parameters; = the all-cell config37's), reusing config1's accessibility +
             expression object, search space and region-to-gene fit (none depends on the motif database): see prepare_shared_upstream.sh.
    config3  config2 with a larger region-to-gene search space: up to 250 kb upstream and downstream of each gene (search_space_upstream /
             search_space_downstream "1000 250000"; config1/2 use "1000 150000", the SCENIC+ default).
    config4  the same at 500 kb ("1000 500000").
             The search space changes the search-space table, the region-to-gene fit and everything after it, so 3 and 4 compute their own; they
             reuse only config1's accessibility + expression object (ACC_GEX.h5mu in shared/outs, which does not depend on the search space).
             Motif enrichment and TF-to-gene are recomputed as in config2 (they do not depend on the search space either, but the workflow has no
             switch to share them). Compare 2 -> 3 -> 4 for the effect of the search space alone.
    All use the hybrid cisTarget database (cistarget/ra-arco-hvc-nc_hybrid), built on the consensus regions; its region list is identical to the
    cisTopic object's 499,348 regions, so no region is dropped or shifted by the 40 % overlap match the older shared database needs.

    python make_configs.py
"""
import importlib.util
import shutil
from pathlib import Path

HERE = Path(__file__).resolve().parent
ALL = HERE.parent.parent / "ra-arco-hvc-nc" / "scenicplus" / "make_configs.py"
spec = importlib.util.spec_from_file_location("hybrid_all_make_configs", ALL)
base = importlib.util.module_from_spec(spec)
spec.loader.exec_module(base)

base.PRISM_ROOT = ("/private/groups/colquittlab/scenicplus/motor-pathway_multiome/"
                   "motor-pathway_multiome_seurat_cellbender.0.05_preprocess_cr/glut-ra_glut-c1h-1_allruns_hybrid")
base.N_TOPICS = {}              # one topic count (10); region_sets
base.REGION_SET_VARIANT = {}
base.SHARED_UPSTREAM = {2}
base.N_CPU_MOTIF = {}
base.cistopic_obj_name = lambda n: "cistopic_obj_ra_c1h1_all.pkl"   # the all-cell script names it cistopic_obj_glut.pkl

SPEC = {'params_motif_enrichment': {'motif_similarity_fdr': 0.001, 'dem_adj_pval_thr': 0.05, 'dem_log2fc_thr': 1.0, 'ctx_auc_threshold': 0.005, 'ctx_nes_threshold': 3.0, 'ctx_rank_threshold': 0.05}, 'params_data_preparation': {'extended_annotation': 'Orthology_annot Motif_similarity_annot'}}
SEARCH_SPACE = {3: "1000 250000", 4: "1000 500000"}
CONFIGS = {1: SPEC, 2: {}}
CONFIGS.update({n: {"params_data_preparation": {"search_space_upstream": ss, "search_space_downstream": ss}} for n, ss in SEARCH_SPACE.items()})
# configs that take only ACC_GEX.h5mu from shared/outs (not the search space or the region-to-gene fit)
SHARED_ACC_GEX_ONLY = set(SEARCH_SPACE)
# wall time: the region-to-gene GBM grows with the number of candidate regions per gene, roughly with the search-space width (not measured here)
TIME = {1: "12:00:00", 2: "08:00:00", 3: "12:00:00", 4: "12:00:00"}   # the medium partition allows at most 12 h
SHARED_ACC_GEX_PRECHECK = """# ACC_GEX.h5mu comes from config1 via prepare_shared_upstream.sh (the search space and region-to-gene fit are this config's own);
# without it Snakemake would silently rebuild the accessibility + expression object.
test -s {shared}/ACC_GEX.h5mu || {{ echo "missing shared input {shared}/ACC_GEX.h5mu -- run prepare_shared_upstream.sh first" >&2; exit 1; }}
"""

HYBRID_DB = ("/private/groups/colquittlab/scenicplus/cistarget/ra-arco-hvc-nc_hybrid", "ra-arco-hvc-nc_hybrid")
CONFIG_DB = {n: HYBRID_DB for n in CONFIGS}   # every config here uses the hybrid cisTarget database


def write_parameter_table():
    rows = []
    ref = base.merged({})
    for n, overrides in CONFIGS.items():
        p = base.merged(overrides)
        row = {"config": f"config{n}"}
        for section in base.BASE_PARAMS:
            row.update(p[section])
        row.update(gex_anndata=base.GEX_DEFAULT, n_topics=10, region_sets="region_sets", cistarget_db=CONFIG_DB[n][1],
                   shared_upstream="yes" if n in base.SHARED_UPSTREAM else ("ACC_GEX only" if n in SHARED_ACC_GEX_ONLY else "no"))
        row["changed_from_hybrid_all_config1"] = ",".join(
            k for s in base.BASE_PARAMS for k, v in p[s].items() if v != ref[s][k]) or "-"
        rows.append(row)
    cols = list(rows[0])
    (HERE / "config_parameters.tsv").write_text(
        "\t".join(cols) + "\n" + "".join("\t".join(str(r[c]) for c in cols) + "\n" for r in rows))


def render(n, overrides):
    """base.render(), with ACC_GEX.h5mu taken from shared/outs (and no full MuData, as in config2) for the configs in SHARED_ACC_GEX_ONLY."""
    y = base.render(n, overrides)
    if n in SHARED_ACC_GEX_ONLY:
        own = f'combined_GEX_ACC_mudata: "{base.PRISM_ROOT}/scenicplus/config{n}/outs/ACC_GEX.h5mu"'
        assert own in y, "ACC_GEX path not found in the rendered config"
        y = y.replace(own, f'combined_GEX_ACC_mudata: "{base.shared_outs(n)}/ACC_GEX.h5mu"')
        assert "build_scplus_mudata" not in y
        y = y.replace("  seed: 666", "  seed: 666\n  build_scplus_mudata: False", 1)
    return y


def precheck(n):
    shared = base.shared_outs(n)
    if n in base.SHARED_UPSTREAM:
        return base.SHARED_PRECHECK.format(shared=shared, args="")
    if n in SHARED_ACC_GEX_ONLY:
        return SHARED_ACC_GEX_PRECHECK.format(shared=shared)
    return ""


def main():
    write_parameter_table()
    src_snakefile = ALL.parent / "workflow" / "Snakefile"
    for n, overrides in CONFIGS.items():
        d = HERE / f"config{n}"
        (d / "Snakemake" / "config").mkdir(parents=True, exist_ok=True)
        (d / "Snakemake" / "workflow").mkdir(parents=True, exist_ok=True)
        base.CISTARGET_DIR, base.CISTARGET_PREFIX = CONFIG_DB[n]   # render() reads these module globals
        (d / "Snakemake" / "config" / "config.yaml").write_text(render(n, overrides))
        shutil.copy(src_snakefile, d / "Snakemake" / "workflow" / "Snakefile")
        (d / "run_snakemake.sbatch").write_text(
            base.SBATCH.format(n=n, cpus=base.N_CPU, prism_cfg=f"{base.PRISM_ROOT}/scenicplus/config{n}", time=TIME[n],
                               precheck=precheck(n)).replace("scenicplus_hybrid_config", "scenicplus_ra_c1h1_all_hybrid_config"))
        print("wrote", d)


if __name__ == "__main__":
    main()
