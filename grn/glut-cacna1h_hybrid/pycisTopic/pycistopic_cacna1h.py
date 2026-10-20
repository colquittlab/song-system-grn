#!/usr/bin/env python3
"""all Glut-CACNA1H subtypes (Glut-CACNA1H-RA, -1 and -2): pycisTopic stage, hybrid labels (a new subset, no earlier run; modelled on glut-ra_glut-arco-1).

Starts from the all-cell hybrid cisTopic object (../ra-arco-hvc-nc_hybrid/pycisTopic/cistopic_obj.pkl: same QC, same consensus regions
and so the same cisTarget database) and does what the earlier notebook did after QC, on this subset's cells only:

  1. subset to the cells of anndata_rna/adata.h5ad that have an ATAC profile (make_cacna1h_adata.py)
  2. train the mallet LDA model with 10 topics (10 as for the RA vs C1H-1 subsets: there is no earlier run of this subset and no topic grid was trained; the grid is not retrained)
  3. topic region sets (otsu, top 3k), imputed accessibility, highly variable regions, DARs for every cluster (adj. p 0.1, |log2FC| >= log2(1.5))
  4. song-pair DARs: the ArchR ones from the all-cell run, copied (Glut-CACNA1H-RA vs Glut-CACNA1H-1) to pycisTopic/archr_song_pair/, NOT into region_sets/ (a single set has no DEM background; DARs_all already has RA-up and C1H-1-up)

Parameters are copied from ../../ra-arco-hvc-nc/pycisTopic/pycistopic_topic_k.py (itself copied from the notebook). Everything
expensive is cached; PYCISTOPIC_REDO=1 recomputes the LDA models. Derived from ../../ra-arco-hvc-nc_glut_hybrid/pycisTopic/pycistopic_glut.py by
../../make_hybrid_subset_dirs.py.

    python pycistopic_cacna1h.py [n_cpu]
"""
import os
import pickle
import sys

import matplotlib
matplotlib.use("Agg")  # evaluate_models plots; with a stale DISPLAY it hung in futex_wait with no CPU use
import anndata as ad
import numpy as np

N_CPU = int(sys.argv[1]) if len(sys.argv) > 1 else 40
K = 10
N_TOPICS = [K]  # 10 as for the RA vs C1H-1 subsets: there is no earlier run of this subset and no topic grid was trained; only the chosen model is needed
REDO = os.environ.get("PYCISTOPIC_REDO") == "1"

B = ("/hdd/jupyter/brad/scenicplus/motor-pathway_multiome/"
     "motor-pathway_multiome_seurat_cellbender.0.05_preprocess_cr/")
SRC_OBJ = B + "ra-arco-hvc-nc_hybrid/pycisTopic/cistopic_obj.pkl"
MAIN = B + "glut-cacna1h_hybrid/"
WORK = MAIN + "pycisTopic/"
SETS = WORK + "region_sets/"
# Ray's Unix socket path must stay under 107 bytes: short temp dir on /hdd (see pycistopic_topic_k.py).
TMP = "/hdd/brad/ray_cacna1h"
os.makedirs(WORK, exist_ok=True)
os.makedirs(TMP, exist_ok=True)

# ArchR song-pair DARs from the all-cell hybrid run (read-only source). Only these may be in DARs_song-pairs: SCENIC+ runs motif
# enrichment on every .bed in the folder.
ARCHR_PAIRS = ["Glut-CACNA1H-RA_VS_Glut-CACNA1H-1"]
ARCHR_PAIRS_DIR = B + "ra-arco-hvc-nc_hybrid/pycisTopic/region_sets_k40_archr/DARs_song-pairs/"

from pycisTopic.diff_features import (find_diff_features, find_highly_variable_features,  # noqa: E402
                                      impute_accessibility, normalize_scores)
from pycisTopic.lda_models import evaluate_models, run_cgs_models_mallet  # noqa: E402
from pycisTopic.topic_binarization import binarize_topics  # noqa: E402
from pycisTopic.utils import region_names_to_coordinates  # noqa: E402

adata = ad.read_h5ad(MAIN + "anndata_rna/adata.h5ad", backed="r")
cells = adata.obs_names.tolist()
labels = adata.obs["cluster"].astype(str)

# 1. subset
obj_fname = WORK + "cistopic_obj_cacna1h_nomodel.pkl"
if os.path.exists(obj_fname) and not REDO:
    obj = pickle.load(open(obj_fname, "rb"))
else:
    full = pickle.load(open(SRC_OBJ, "rb"))
    # RNA cells that failed ATAC QC have no cisTopic profile (as in the earlier glut run, which subset the cisTopic object to
    # the labelled cells); SCENIC+ intersects the two modalities itself, so the RNA input keeps them.
    have = set(full.cell_names)
    both = [c for c in cells if c in have]
    print(f"{len(cells) - len(both)} of {len(cells)} RNA cells have no ATAC profile; cisTopic subset = {len(both)}", flush=True)
    assert len(both) > 0.9 * len(cells), "unexpectedly many RNA cells missing from the cisTopic object"
    obj = full.subset(cells=both, copy=True)
    del full
    obj.cell_data["cluster"] = obj.cell_data["cluster"].astype(str)
    assert (obj.cell_data.loc[both, "cluster"].values == labels.loc[both].values).all(), "cisTopic and adata labels disagree"
    pickle.dump(obj, open(obj_fname, "wb"))
print(f"cisTopic: {len(obj.cell_names)} cells, {len(obj.region_names)} regions", flush=True)

# 2. LDA
models_fname = WORK + "models.pkl"
if os.path.exists(models_fname) and not REDO:
    models = pickle.load(open(models_fname, "rb"))
else:
    os.environ["MALLET_MEMORY"] = "200G"
    tmp_path = os.path.join(TMP, "pycistopic_mallet")
    os.makedirs(tmp_path, exist_ok=True)
    models = run_cgs_models_mallet(
        obj, n_topics=N_TOPICS, n_cpu=N_CPU, n_iter=500, random_state=555, alpha=50, alpha_by_topic=True,
        eta=0.1, eta_by_topic=False, tmp_path=tmp_path, save_path=tmp_path,
        mallet_path="/home/brad/ssd/repos/Mallet-202108/bin/mallet")
    pickle.dump(models, open(models_fname, "wb"))
print("models:", [m.n_topic for m in models], flush=True)

model = evaluate_models(models, select_model=K, return_model=True)
obj.add_LDA_model(model)
assert obj.selected_model.n_topic == K

# 3. topic sets, imputed accessibility, DARs
region_bin_topics_top_3k = binarize_topics(obj, method="ntop", ntop=3_000, plot=False)
region_bin_topics_otsu = binarize_topics(obj, method="otsu", plot=False)
imputed_acc_obj = impute_accessibility(obj, selected_cells=None, selected_regions=None, scale_factor=10**6)
normalized_imputed_acc_obj = normalize_scores(imputed_acc_obj, scale_factor=10**4)
variable_regions = find_highly_variable_features(
    normalized_imputed_acc_obj, min_disp=0.05, min_mean=0.0125, max_mean=3, max_disp=np.inf,
    n_bins=20, n_top_features=None, plot=False)
print("variable regions:", len(variable_regions), flush=True)

markers_dict_all = find_diff_features(
    obj, imputed_acc_obj, variable="cluster", var_features=variable_regions, contrasts=None,
    adjpval_thr=0.1, log2fc_thr=np.log2(1.5), n_cpu=N_CPU, _temp_dir=TMP, split_pattern="-")
print("DARs:", {k: len(v) for k, v in markers_dict_all.items()}, flush=True)


def write_sets(subdir, sets):
    d = os.path.join(SETS, subdir)
    os.makedirs(d, exist_ok=True)
    for name, df in sets.items():
        region_names_to_coordinates(df.index).sort_values(["Chromosome", "Start", "End"]).to_csv(
            os.path.join(d, f"{name}.bed"), sep="\t", header=False, index=False)


write_sets("Topics_otsu", region_bin_topics_otsu)
write_sets("Topics_top_3k", region_bin_topics_top_3k)
write_sets("DARs_all", markers_dict_all)
# 4. ArchR song-pair DARs, copied
import glob
import shutil
pairs_out = os.path.join(WORK, "archr_song_pair")  # NOT under region_sets/: a one-set folder has an empty DEM background and fails
os.makedirs(pairs_out, exist_ok=True)
for f in glob.glob(os.path.join(pairs_out, "*.bed")):
    os.remove(f)
for nm in ARCHR_PAIRS:
    shutil.copy(ARCHR_PAIRS_DIR + nm + ".bed", pairs_out)
assert sorted(os.listdir(pairs_out)) == sorted(n + ".bed" for n in ARCHR_PAIRS)
have = set(obj.region_names)
for nm in ARCHR_PAIRS:
    bed = [l.split("	")[:3] for l in open(os.path.join(pairs_out, nm + ".bed")).read().splitlines()]
    n_in = sum(f"{c}:{s0}-{e}" in have for c, s0, e in bed)
    print(f"song-pair DARs {nm}: {len(bed)} regions, {n_in} in this cisTopic object", flush=True)
for name, d in [("region_bin_topics_top3k", region_bin_topics_top_3k), ("region_bin_topics_otsu", region_bin_topics_otsu),
                ("markers_dict_all", markers_dict_all)]:
    pickle.dump(d, open(SETS + name + ".pkl", "wb"))

# the object SCENIC+ reads: this subset's cells only, K-topic model selected
pickle.dump(obj, open(WORK + "cistopic_obj_cacna1h.pkl", "wb"))
print("wrote", WORK + "cistopic_obj_cacna1h.pkl", "and", SETS, flush=True)
