#!/usr/bin/env python3
"""Per-config SCENIC+ Snakemake directories for ra-arco-hvc-nc_gaba_hybrid (GABAergic neurons, hybrid labels).

Reuses the YAML/sbatch rendering and the shared Snakefile of ../../ra-arco-hvc-nc/scenicplus/make_configs.py (importing it only
defines its tables; nothing is written until main()), so the runs cannot drift apart. Only the root and the configs differ.
Derived from ../../ra-arco-hvc-nc_glut_hybrid/scenicplus/make_configs.py by ../../make_hybrid_subset_dirs.py.

    config1  the earlier run's parameters (ra-arco-hvc-nc_gaba (config13)): motif_similarity_fdr 0.01, dem_adj_pval_thr 0.05, dem_log2fc_thr 1.0, ctx_auc_threshold 0.01, ctx_nes_threshold 2.0, ctx_rank_threshold 0.1.
             20 topics, raw counts in .raw.
    config2  the all-cell hybrid run's strict config1 thresholds (base parameters; = the all-cell config37's), reusing config1's accessibility +
             expression object, search space and region-to-gene fit (none depends on the motif database): see prepare_shared_upstream.sh.
    Both use the hybrid cisTarget database (cistarget/ra-arco-hvc-nc_hybrid), built on the consensus regions; its region list is identical to the
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
                   "motor-pathway_multiome_seurat_cellbender.0.05_preprocess_cr/ra-arco-hvc-nc_gaba_hybrid")
base.N_TOPICS = {}              # one topic count (20); region_sets
base.REGION_SET_VARIANT = {}
base.SHARED_UPSTREAM = {2}
base.N_CPU_MOTIF = {}
base.cistopic_obj_name = lambda n: "cistopic_obj_gaba.pkl"   # the all-cell script names it cistopic_obj_glut.pkl

SPEC = {'params_motif_enrichment': {'motif_similarity_fdr': 0.01, 'dem_adj_pval_thr': 0.05, 'dem_log2fc_thr': 1.0, 'ctx_auc_threshold': 0.01, 'ctx_nes_threshold': 2.0, 'ctx_rank_threshold': 0.1}}
CONFIGS = {1: SPEC, 2: {}}

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
        row.update(gex_anndata=base.GEX_DEFAULT, n_topics=20, region_sets="region_sets", cistarget_db=CONFIG_DB[n][1],
                   shared_upstream="yes" if n in base.SHARED_UPSTREAM else "no")
        row["changed_from_hybrid_all_config1"] = ",".join(
            k for s in base.BASE_PARAMS for k, v in p[s].items() if v != ref[s][k]) or "-"
        rows.append(row)
    cols = list(rows[0])
    (HERE / "config_parameters.tsv").write_text(
        "\t".join(cols) + "\n" + "".join("\t".join(str(r[c]) for c in cols) + "\n" for r in rows))


def main():
    write_parameter_table()
    src_snakefile = ALL.parent / "workflow" / "Snakefile"
    for n, overrides in CONFIGS.items():
        d = HERE / f"config{n}"
        (d / "Snakemake" / "config").mkdir(parents=True, exist_ok=True)
        (d / "Snakemake" / "workflow").mkdir(parents=True, exist_ok=True)
        base.CISTARGET_DIR, base.CISTARGET_PREFIX = CONFIG_DB[n]   # render() reads these module globals
        (d / "Snakemake" / "config" / "config.yaml").write_text(base.render(n, overrides))
        shutil.copy(src_snakefile, d / "Snakemake" / "workflow" / "Snakefile")
        (d / "run_snakemake.sbatch").write_text(
            base.SBATCH.format(n=n, cpus=base.N_CPU, prism_cfg=f"{base.PRISM_ROOT}/scenicplus/config{n}",
                               time="08:00:00" if n in base.SHARED_UPSTREAM else "12:00:00",
                               precheck=base.SHARED_PRECHECK.format(shared=f"{base.PRISM_ROOT}/scenicplus/shared/outs")
                               if n in base.SHARED_UPSTREAM else "").replace("scenicplus_hybrid_config", "scenicplus_gaba_hybrid_config"))
        print("wrote", d)


if __name__ == "__main__":
    main()
