#!/usr/bin/env python3
"""Per-config SCENIC+ Snakemake directories for ra-arco-hvc-nc_glut_hybrid (glutamatergic subset, hybrid labels).

Reuses the YAML/sbatch rendering and the shared Snakefile of ../../ra-arco-hvc-nc/scenicplus/make_configs.py (importing it only
defines its tables; nothing is written until main()), so the two runs cannot drift apart. Only the root and the configs differ.

    config1  the earlier glut run's parameters (ra-arco-hvc-nc_glut/scenicplus/scenicplus_config.yaml, "config15"): the relaxed
             motif stage -- motif_similarity_fdr 0.001, dem 0.1 / 0.5, ctx auc 0.0025, nes 2.0, rank 0.2 (hybrid-all config11 but
             with rank 0.2). 20 topics, raw counts in .raw, cisTarget database shared with the all-cell run.
    config2  control: the all-cell hybrid run's strict config1 thresholds (base parameters), so glut-only vs all-cell can also be
             compared at identical thresholds.
    config3  config2's parameters (= the all-cell config37's) with the hybrid cisTarget database, built on this object's own consensus
             regions (cistarget/ra-arco-hvc-nc_hybrid; its region list is identical to the cisTopic object's 499,348 regions, so no
             region is dropped or shifted by the 40 % overlap match that config1/2's older shared database needs).
    config4  config1's (glut spec) thresholds with the hybrid database.
    Everything else stays at 20 topics, raw counts, region_sets (ArchR song-pair DARs). Configs 3/4 reuse config1's accessibility +
    expression object, search space and region-to-gene fit (none depends on the motif database): see prepare_shared_upstream.sh.

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
                   "motor-pathway_multiome_seurat_cellbender.0.05_preprocess_cr/ra-arco-hvc-nc_glut_hybrid")
base.N_TOPICS = {}              # one topic count (20); cistopic_obj_glut.pkl / region_sets
base.REGION_SET_VARIANT = {}
base.SHARED_UPSTREAM = {3, 4}
base.N_CPU_MOTIF = {}

GLUT_SPEC = {"params_motif_enrichment": {"motif_similarity_fdr": 0.001, "dem_adj_pval_thr": 0.1, "dem_log2fc_thr": 0.5,
                                         "ctx_auc_threshold": 0.0025, "ctx_nes_threshold": 2.0, "ctx_rank_threshold": 0.2}}
CONFIGS = {1: GLUT_SPEC, 2: {}, 3: {}, 4: GLUT_SPEC}

OLD_DB = (base.CISTARGET_DIR, base.CISTARGET_PREFIX)
HYBRID_DB = ("/private/groups/colquittlab/scenicplus/cistarget/ra-arco-hvc-nc_hybrid", "ra-arco-hvc-nc_hybrid")
CONFIG_DB = {1: OLD_DB, 2: OLD_DB, 3: HYBRID_DB, 4: HYBRID_DB}


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
                               if n in base.SHARED_UPSTREAM else "").replace("scenicplus_hybrid_config", "scenicplus_glut_hybrid_config"))
        print("wrote", d)


if __name__ == "__main__":
    main()
