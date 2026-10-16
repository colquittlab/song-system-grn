#!/usr/bin/env python3
"""Derive the per-subset SCENIC+ run directories (hybrid labels) from ra-arco-hvc-nc_glut_hybrid/.

ra-arco-hvc-nc_glut_hybrid/ is the reference: glutamatergic neurons, hybrid labels, 20 topics, ArchR song-pair DARs, configs 1-4
(config1 = the older subset run's motif thresholds, config2 = the all-cell run's strict thresholds, config3/4 = the same two with the
hybrid cisTarget database). Each subset below changes only what the older run of that subset specified: the cell filter, the topic count,
the DAR thresholds, which ArchR song pairs apply, and config1's motif/data-preparation parameters. Everything else is copied from the
reference files by exact-string replacement (each replacement is asserted, so a change in the reference breaks this loudly instead of
silently producing a different script).

    python make_hybrid_subset_dirs.py [subset ...]   # (re)writes the named subset directories (default: all); nothing under /hdd is touched
"""
import json
import re
import shutil
from pathlib import Path

HERE = Path(__file__).resolve().parent
REF = HERE / "ra-arco-hvc-nc_glut_hybrid"
REF_ROOT = "ra-arco-hvc-nc_glut_hybrid"

SUBSETS = {
    "gaba": dict(
        dirname="ra-arco-hvc-nc_gaba_hybrid", label="GABAergic neurons", old="ra-arco-hvc-nc_gaba (config13)",
        filter_doc=("clusters containing \"GABA\", minus the immature stage. The earlier run excluded it with `~cluster.str.contains('Pre')`;\n"
                    "the hybrid labels renamed GABA-Pre to GABA-Im (snrna/naming/hybrid_division_naming.qmd), so a substring match would silently\n"
                    "keep it. The list is explicit and asserted instead."),
        filter_code='''PRECURSORS = ["GABA-Im"]   # GABA-Pre in the older labels

adata = sc.read_h5ad(SRC)
assert np.all(adata.X.data == np.round(adata.X.data)), "source X is not integer counts"

cl = adata.obs["cluster"].astype(str)
keep = cl.str.startswith("GABA") & ~cl.isin(PRECURSORS)
adata = adata[keep.values].copy()
assert set(PRECURSORS).isdisjoint(adata.obs["cluster"].astype(str)), "immature cluster survived the filter"
adata.obs["cluster"] = adata.obs["cluster"].astype(str)''',
        K=20, k_note="the earlier GABA run's hand-picked choice (it trained 2-40)",
        dar_adjp=0.05, dar_log2fc="np.log2(1.5)", pairs=[],
        spec_motif=dict(motif_similarity_fdr=0.01, dem_adj_pval_thr=0.05, dem_log2fc_thr=1.0,
                        ctx_auc_threshold=0.01, ctx_nes_threshold=2.0, ctx_rank_threshold=0.1),
        spec_dataprep={}, spec_name="config13 of the earlier GABA run",
    ),
    "ra_c1h1": dict(
        dirname="glut-ra_glut-c1h-1_hybrid", label="Glut-CACNA1H-RA vs Glut-CACNA1H-1 (the hybrid labels of Glut-RA vs Glut-Arco-1)",
        old="glut-ra_glut-arco-1 (config18)",
        filter_doc=("the two clusters Glut-CACNA1H-RA and Glut-CACNA1H-1 (older labels Glut-RA, Glut-Arco-1), with the run3 samples excluded as in\n"
                    "the earlier run (`~sample_short.str.contains('run3')`). That halves the cells (see the README); set EXCLUDE_RUN3=0 to keep them."),
        filter_code='''CLUSTERS = ["Glut-CACNA1H-RA", "Glut-CACNA1H-1"]
EXCLUDE_RUN3 = os.environ.get("EXCLUDE_RUN3", "1") == "1"   # the earlier run: ~sample_short.str.contains('run3')

adata = sc.read_h5ad(SRC)
assert np.all(adata.X.data == np.round(adata.X.data)), "source X is not integer counts"

cl = adata.obs["cluster"].astype(str)
keep = cl.isin(CLUSTERS)
if EXCLUDE_RUN3:
    keep &= ~adata.obs["sample_short"].astype(str).str.contains("run3")
adata = adata[keep.values].copy()
assert set(adata.obs["cluster"].astype(str)) == set(CLUSTERS), "a cluster is missing after the filter"
adata.obs["cluster"] = adata.obs["cluster"].astype(str)''',
        K=10, k_note="the earlier RA vs Arco-1 run's hand-picked choice (it trained 2-40)",
        dar_adjp=0.1, dar_log2fc="np.log2(1.5)", pairs=["Glut-CACNA1H-RA_VS_Glut-CACNA1H-1"],
        spec_motif=dict(motif_similarity_fdr=0.001, dem_adj_pval_thr=0.05, dem_log2fc_thr=1.0,
                        ctx_auc_threshold=0.005, ctx_nes_threshold=3.0, ctx_rank_threshold=0.05),
        spec_dataprep={"extended_annotation": "Orthology_annot Motif_similarity_annot"}, spec_name="config18 of the earlier RA vs Arco-1 run",
    ),
    "ra_c1h1_all": dict(
        dirname="glut-ra_glut-c1h-1_allruns_hybrid",
        label="Glut-CACNA1H-RA vs Glut-CACNA1H-1, all runs (the hybrid labels of Glut-RA vs Glut-Arco-1, run3 included)",
        old="glut-ra_glut-arco-1 (config18)",
        filter_doc=("the two clusters Glut-CACNA1H-RA and Glut-CACNA1H-1 (older labels Glut-RA, Glut-Arco-1) from ALL runs. The earlier run excluded the run3\n"
                    "samples; glut-ra_glut-c1h-1_hybrid keeps that exclusion (745 cells), this variant does not (1,966 cells)."),
        filter_code='''CLUSTERS = ["Glut-CACNA1H-RA", "Glut-CACNA1H-1"]

adata = sc.read_h5ad(SRC)
assert np.all(adata.X.data == np.round(adata.X.data)), "source X is not integer counts"

cl = adata.obs["cluster"].astype(str)
keep = cl.isin(CLUSTERS)   # all runs, run3 included
adata = adata[keep.values].copy()
assert set(adata.obs["cluster"].astype(str)) == set(CLUSTERS), "a cluster is missing after the filter"
adata.obs["cluster"] = adata.obs["cluster"].astype(str)''',
        K=10, k_note="the earlier RA vs Arco-1 run's hand-picked choice on its 745 cells (it trained 2-40); not re-selected for the larger cell set",
        dar_adjp=0.1, dar_log2fc="np.log2(1.5)", pairs=["Glut-CACNA1H-RA_VS_Glut-CACNA1H-1"],
        spec_motif=dict(motif_similarity_fdr=0.001, dem_adj_pval_thr=0.05, dem_log2fc_thr=1.0,
                        ctx_auc_threshold=0.005, ctx_nes_threshold=3.0, ctx_rank_threshold=0.05),
        spec_dataprep={"extended_annotation": "Orthology_annot Motif_similarity_annot"}, spec_name="config18 of the earlier RA vs Arco-1 run",
    ),
}


def rep(text, old, new, count=1):
    assert text.count(old) == count, f"expected {count} x {old[:70]!r}, found {text.count(old)}"
    return text.replace(old, new)


def make_adata_script(s, tag, root):
    t = (REF / "anndata_rna" / "make_glut_adata.py").read_text()
    head_end = t.index('import os\n')
    doc = f'''#!/usr/bin/env python3
"""{s["label"]}: SCENIC+ expression input, hybrid labels (the hybrid counterpart of {s["old"].split(" ")[0]}/anndata_rna/anndata.ipynb).

Same steps as ../../ra-arco-hvc-nc/anndata_rna/anndata.ipynb (which reads data/adata.h5ad from export_h5ad_hybrid.R), with the cell filter of the
earlier run: {s["filter_doc"]}

    adata.h5ad   .raw = raw integer counts (what SCENIC+ reads); X = log-normalized, HVG, scaled (display only)

    python make_{tag}_adata.py
"""
'''
    body = t[head_end:]
    a = body.index('PRECURSORS = [')
    b = body.index('adata.obsm["archr_umap"]')
    body = body[:a] + s["filter_code"] + '\nprint(adata.obs["cluster"].value_counts().to_string())\n\n' + body[b:]
    body = rep(body, 'OUT_DIR = B + "ra-arco-hvc-nc_glut_hybrid/anndata_rna/"', f'OUT_DIR = B + "{root}/anndata_rna/"')
    return doc + body


def make_pycistopic_script(s, tag, root):
    t = (REF / "pycisTopic" / "pycistopic_glut.py").read_text()
    body = t[t.index("import os\n"):]
    pairs_doc = ("  4. song-pair DARs: the ArchR ones from the all-cell run, copied (" + ", ".join(p.replace("_VS_", " vs ") for p in s["pairs"]) + ")"
                 if s["pairs"] else "  4. no song-pair DARs for this subset (the earlier run had none)")
    doc = f'''#!/usr/bin/env python3
"""{s["label"]}: pycisTopic stage, hybrid labels (counterpart of {s["old"].split(" ")[0]}/pycisTopic/pycisTopic.ipynb).

Starts from the all-cell hybrid cisTopic object (../ra-arco-hvc-nc_hybrid/pycisTopic/cistopic_obj.pkl: same QC, same consensus regions
and so the same cisTarget database) and does what the earlier notebook did after QC, on this subset's cells only:

  1. subset to the cells of anndata_rna/adata.h5ad that have an ATAC profile (make_{tag}_adata.py)
  2. train the mallet LDA model with {s["K"]} topics ({s["k_note"]}; the grid is not retrained)
  3. topic region sets (otsu, top 3k), imputed accessibility, highly variable regions, DARs for every cluster (adj. p {s["dar_adjp"]}, |log2FC| >= log2(1.5))
{pairs_doc}

Parameters are copied from ../../ra-arco-hvc-nc/pycisTopic/pycistopic_topic_k.py (itself copied from the notebook). Everything
expensive is cached; PYCISTOPIC_REDO=1 recomputes the LDA models. Derived from ../../ra-arco-hvc-nc_glut_hybrid/pycisTopic/pycistopic_glut.py by
../../make_hybrid_subset_dirs.py.

    python pycistopic_{tag}.py [n_cpu]
"""
'''
    body = rep(body, "K = 20\nN_TOPICS = [K]  # the earlier glut run trained 2-40 and picked 20 by hand; only the chosen model is needed (~30 min)",
               f"K = {s['K']}\nN_TOPICS = [K]  # {s['k_note']}; only the chosen model is needed")
    body = rep(body, 'MAIN = B + "ra-arco-hvc-nc_glut_hybrid/"', f'MAIN = B + "{root}/"')
    body = rep(body, 'TMP = "/hdd/brad/ray_glut_hybrid"', f'TMP = "/hdd/brad/ray_{tag}"')
    body = rep(body, "cistopic_obj_glut_nomodel.pkl", f"cistopic_obj_{tag}_nomodel.pkl")
    body = rep(body, 'WORK + "cistopic_obj_glut.pkl"', f'WORK + "cistopic_obj_{tag}.pkl"', count=2)
    body = rep(body, "glut RNA cells have no ATAC profile", "RNA cells have no ATAC profile")
    body = rep(body, "adjpval_thr=0.05, log2fc_thr=np.log2(1.5), n_cpu=N_CPU, _temp_dir=TMP, split_pattern=\"-\")\nprint(\"DARs:\"",
               f"adjpval_thr={s['dar_adjp']}, log2fc_thr={s['dar_log2fc']}, n_cpu=N_CPU, _temp_dir=TMP, split_pattern=\"-\")\nprint(\"DARs:\"")
    body = rep(body, "# the object SCENIC+ reads: glut cells only, K-topic model selected", "# the object SCENIC+ reads: this subset's cells only, K-topic model selected")
    a = body.index("# ArchR song-pair DARs from the all-cell hybrid run")
    b = body.index("from pycisTopic.diff_features")
    if s["pairs"]:
        pairs = ", ".join(json.dumps(p) for p in s["pairs"])
        body = body[:a] + (f"# ArchR song-pair DARs from the all-cell hybrid run (read-only source). Only these may be in DARs_song-pairs: SCENIC+ runs motif\n"
                           f"# enrichment on every .bed in the folder.\nARCHR_PAIRS = [{pairs}]\n"
                           'ARCHR_PAIRS_DIR = B + "ra-arco-hvc-nc_hybrid/pycisTopic/region_sets_k40_archr/DARs_song-pairs/"\n\n') + body[b:]
    else:
        body = body[:a] + body[b:]
        a = body.index("# 4. ArchR song-pair DARs, copied")
        b = body.index('for name, d in [("region_bin_topics_top3k"')
        body = body[:a] + "# 4. no song-pair DARs for this subset\n" + body[b:]
    return doc + body


def make_notebook(script_text, tag, s):
    import nbformat
    from nbformat.v4 import new_code_cell, new_markdown_cell, new_notebook
    doc = re.search(r'"""(.*?)"""', script_text, re.S).group(1).strip()
    body = script_text[script_text.index('"""', script_text.index('"""') + 3) + 3:].lstrip("\n")
    body = rep(body, "N_CPU = int(sys.argv[1]) if len(sys.argv) > 1 else 40", 'N_CPU = int(os.environ.get("N_CPU", 40))  # kernel-safe: sys.argv is not ours in Jupyter')
    body = rep(body, "import sys\n", "")
    doc = doc.replace(f"    python pycistopic_{tag}.py [n_cpu]", f"Run top to bottom. Set N_CPU in the environment to change the worker count (default 40); `PYCISTOPIC_REDO=1` retrains the LDA model.\n`pycistopic_{tag}.py` in this folder is the same code as a script.")
    marks = [("from pycisTopic.diff_features", "## Imports"), ("adata = ad.read_h5ad", "## Cells (from `anndata_rna/make_" + tag + "_adata.py`)"),
             ("# 1. subset", "## 1. Subset the all-cell hybrid cisTopic object to this subset's cells with ATAC"),
             ("# 2. LDA", f"## 2. LDA model ({s['K']} topics)"), ("# 3. topic sets", "## 3. Topic region sets, imputed accessibility, per-cluster DARs"),
             ("markers_dict_all = find_diff_features", "### DARs for every cluster (pycisTopic)"), ("def write_sets", "## Save region sets")]
    if s["pairs"]:
        marks.append(("# 4. ArchR song-pair", "## 4. Song-pair DARs: ArchR, copied from the all-cell run"))
    marks.append(("# the object SCENIC+ reads", "## Save the cisTopic object SCENIC+ reads"))
    lines = body.splitlines(keepends=True)
    idx = [0] + [next(k for k, l in enumerate(lines) if l.startswith(m)) for m, _ in marks] + [len(lines)]
    nb = new_notebook()
    nb.cells.append(new_markdown_cell(f"# {s['dirname']}: pycisTopic ({s['label']}, hybrid labels)\n\n" + doc))
    nb.cells.append(new_markdown_cell("## Setup"))
    nb.cells.append(new_code_cell("".join(lines[idx[0]:idx[1]]).rstrip() + "\n"))
    for j, (m, head) in enumerate(marks):
        nb.cells.append(new_markdown_cell(head))
        nb.cells.append(new_code_cell("".join(lines[idx[j + 1]:idx[j + 2]]).rstrip() + "\n"))
    nb.metadata["kernelspec"] = {"display_name": "scenicplus", "language": "python", "name": "python3"}
    return nbformat.writes(nb)


def make_configs_script(s, tag, root):
    mot = s["spec_motif"]
    spec = {"params_motif_enrichment": dict(mot)}
    if s["spec_dataprep"]:
        spec["params_data_preparation"] = dict(s["spec_dataprep"])
    return f'''#!/usr/bin/env python3
"""Per-config SCENIC+ Snakemake directories for {root} ({s["label"]}, hybrid labels).

Reuses the YAML/sbatch rendering and the shared Snakefile of ../../ra-arco-hvc-nc/scenicplus/make_configs.py (importing it only
defines its tables; nothing is written until main()), so the runs cannot drift apart. Only the root and the configs differ.
Derived from ../../ra-arco-hvc-nc_glut_hybrid/scenicplus/make_configs.py by ../../make_hybrid_subset_dirs.py.

    config1  the earlier run's parameters ({s["old"]}): {", ".join(f"{k} {v}" for k, v in mot.items())}{
        "; " + ", ".join(f"{k} {v!r}" for k, v in s["spec_dataprep"].items()) if s["spec_dataprep"] else ""}.
             {s["K"]} topics, raw counts in .raw, cisTarget database shared with the all-cell run.
    config2  control: the all-cell hybrid run's strict config1 thresholds (base parameters).
    config3  config2's parameters (= the all-cell config37's) with the hybrid cisTarget database, built on the consensus regions
             (cistarget/ra-arco-hvc-nc_hybrid; its region list is identical to the cisTopic object's 499,348 regions, so no region is dropped
             or shifted by the 40 % overlap match that config1/2's older shared database needs).
    config4  config1's thresholds with the hybrid database.
    Configs 3/4 reuse config1's accessibility + expression object, search space and region-to-gene fit (none depends on the motif database):
    see prepare_shared_upstream.sh.

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
                   "motor-pathway_multiome_seurat_cellbender.0.05_preprocess_cr/{root}")
base.N_TOPICS = {{}}              # one topic count ({s["K"]}); region_sets
base.REGION_SET_VARIANT = {{}}
base.SHARED_UPSTREAM = {{3, 4}}
base.N_CPU_MOTIF = {{}}
base.cistopic_obj_name = lambda n: "cistopic_obj_{tag}.pkl"   # the all-cell script names it cistopic_obj_glut.pkl

SPEC = {spec!r}
CONFIGS = {{1: SPEC, 2: {{}}, 3: {{}}, 4: SPEC}}

OLD_DB = (base.CISTARGET_DIR, base.CISTARGET_PREFIX)
HYBRID_DB = ("/private/groups/colquittlab/scenicplus/cistarget/ra-arco-hvc-nc_hybrid", "ra-arco-hvc-nc_hybrid")
CONFIG_DB = {{1: OLD_DB, 2: OLD_DB, 3: HYBRID_DB, 4: HYBRID_DB}}


def write_parameter_table():
    rows = []
    ref = base.merged({{}})
    for n, overrides in CONFIGS.items():
        p = base.merged(overrides)
        row = {{"config": f"config{{n}}"}}
        for section in base.BASE_PARAMS:
            row.update(p[section])
        row.update(gex_anndata=base.GEX_DEFAULT, n_topics={s["K"]}, region_sets="region_sets", cistarget_db=CONFIG_DB[n][1],
                   shared_upstream="yes" if n in base.SHARED_UPSTREAM else "no")
        row["changed_from_hybrid_all_config1"] = ",".join(
            k for s in base.BASE_PARAMS for k, v in p[s].items() if v != ref[s][k]) or "-"
        rows.append(row)
    cols = list(rows[0])
    (HERE / "config_parameters.tsv").write_text(
        "\\t".join(cols) + "\\n" + "".join("\\t".join(str(r[c]) for c in cols) + "\\n" for r in rows))


def main():
    write_parameter_table()
    src_snakefile = ALL.parent / "workflow" / "Snakefile"
    for n, overrides in CONFIGS.items():
        d = HERE / f"config{{n}}"
        (d / "Snakemake" / "config").mkdir(parents=True, exist_ok=True)
        (d / "Snakemake" / "workflow").mkdir(parents=True, exist_ok=True)
        base.CISTARGET_DIR, base.CISTARGET_PREFIX = CONFIG_DB[n]   # render() reads these module globals
        (d / "Snakemake" / "config" / "config.yaml").write_text(base.render(n, overrides))
        shutil.copy(src_snakefile, d / "Snakemake" / "workflow" / "Snakefile")
        (d / "run_snakemake.sbatch").write_text(
            base.SBATCH.format(n=n, cpus=base.N_CPU, prism_cfg=f"{{base.PRISM_ROOT}}/scenicplus/config{{n}}",
                               time="08:00:00" if n in base.SHARED_UPSTREAM else "12:00:00",
                               precheck=base.SHARED_PRECHECK.format(shared=f"{{base.PRISM_ROOT}}/scenicplus/shared/outs")
                               if n in base.SHARED_UPSTREAM else "").replace("scenicplus_hybrid_config", "scenicplus_{tag}_hybrid_config"))
        print("wrote", d)


if __name__ == "__main__":
    main()
'''


def make_readme(s, tag, root):
    mot = s["spec_motif"]
    cells_note = {
        "gaba": "all `GABA-*` clusters minus `GABA-Im` (the renamed `GABA-Pre`): 8 clusters, 3,275 RNA cells",
        "ra_c1h1": "`Glut-CACNA1H-RA` and `Glut-CACNA1H-1` with the run3 samples excluded as in the earlier run: 745 RNA cells (258 RA, 487 CACNA1H-1). "
                   "Without the exclusion (`EXCLUDE_RUN3=0`) it is 1,966 (560 and 1,406): that is `glut-ra_glut-c1h-1_allruns_hybrid`.",
        "ra_c1h1_all": "`Glut-CACNA1H-RA` and `Glut-CACNA1H-1` from all runs, run3 included: 1,966 RNA cells (560 RA, 1,406 CACNA1H-1). The earlier run excluded run3 "
                       "(`glut-ra_glut-c1h-1_hybrid`, 745 cells). The topic count (10) was picked for those 745 cells and is not re-selected here.",
    }[tag]
    pairs = (f"DARs for the ArchR song pair copied from the all-cell run (`{s['pairs'][0].replace('_VS_', '` vs `')}`)" if s["pairs"]
             else "no song-pair DARs (the earlier run had none)")
    return f'''# {root}

SCENIC+ on **{s["label"]}**, hybrid labels (`cluster_hybrid`). Derived from `ra-arco-hvc-nc_glut_hybrid` (see that README) by
`../make_hybrid_subset_dirs.py`, changing only what the earlier run ({s["old"]}) specified.

| | earlier run | this directory |
|---|---|---|
| cells | older labels | {cells_note} |
| `.raw` | older export (not re-checked here) | raw counts |
| LDA topics | {s["K"]} (hand-picked from 2-40) | {s["K"]} (only this model trained) |
| DARs (every cluster) | adj. p {s["dar_adjp"]}, abs(log2FC) >= log2(1.5) | same; {pairs} |
| motif stage | {", ".join(f"{k} {v}" for k, v in mot.items())}{"; " + ", ".join(f"{k} {v!r}" for k, v in s["spec_dataprep"].items()) if s["spec_dataprep"] else ""} | **config1 = this**, config2 = the all-cell strict thresholds; config3/4 = the same two with the hybrid cisTarget database |

## Steps (same as the glut run)

1. `anndata_rna/make_{tag}_adata.py` -> `anndata_rna/adata.h5ad` (`.raw` = raw counts).
2. `pycisTopic/pycistopic_{tag}.py` (script) or `pycisTopic/pycisTopic.ipynb` (the same code, interactive): subset the all-cell hybrid cisTopic object to
   the cells with ATAC, train the {s["K"]}-topic model, write topic sets and DARs. Output `cistopic_obj_{tag}.pkl` and `region_sets/`.
3. `scenicplus/make_configs.py` -> config1-4. Configs 3/4 reuse config1's accessibility+expression object, search space and region-to-gene fit:
   run `scenicplus/submit_all.sh 1 2` on prism, then `scenicplus/prepare_shared_upstream.sh`, then `scenicplus/submit_all.sh 3 4`.
   Helpers (`submit_all.sh`, `link_results_to_hdd.sh`, `prepare_shared_upstream.sh`) are the glut run's with this root.

Data live under `/hdd/jupyter/brad/scenicplus/motor-pathway_multiome/.../{root}/`, outside git. Analysis templates are not made yet (they need results).
'''


def build(tag, s):
    root = s["dirname"]
    D = HERE / root
    for sub in ("anndata_rna", "pycisTopic", "scenicplus"):
        (D / sub).mkdir(parents=True, exist_ok=True)
    (D / "anndata_rna" / f"make_{tag}_adata.py").write_text(make_adata_script(s, tag, root))
    pc = make_pycistopic_script(s, tag, root)
    (D / "pycisTopic" / f"pycistopic_{tag}.py").write_text(pc)
    (D / "pycisTopic" / "pycisTopic.ipynb").write_text(make_notebook(pc, tag, s))
    (D / "scenicplus" / "make_configs.py").write_text(make_configs_script(s, tag, root))
    for name in ("submit_all.sh", ".gitignore", "chromsizes.tsv", "genome_annotation.tsv"):
        shutil.copy2(REF / "scenicplus" / name, D / "scenicplus" / name)
    for name in ("link_results_to_hdd.sh", "prepare_shared_upstream.sh"):
        text = (REF / "scenicplus" / name).read_text().replace(REF_ROOT, root)
        assert root in text
        (D / "scenicplus" / name).write_text(text)
        (D / "scenicplus" / name).chmod((REF / "scenicplus" / name).stat().st_mode & 0o777)
    (D / "README.md").write_text(make_readme(s, tag, root))
    print("wrote", D)


if __name__ == "__main__":
    import sys
    wanted = sys.argv[1:] or list(SUBSETS)
    for tag in wanted:
        build(tag, SUBSETS[tag])
