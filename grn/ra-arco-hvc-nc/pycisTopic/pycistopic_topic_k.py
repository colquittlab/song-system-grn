#!/usr/bin/env python3
"""Redo everything downstream of LDA model selection for a different number of topics, without retraining.

pycisTopic.ipynb trains models for several topic counts (models.pkl) and then picks one by hand (20). That choice
feeds three things SCENIC+ sees: the topic region sets, the imputed accessibility (so the DARs and the ACC matrix
SCENIC+ builds from the cisTopic object), and the cisTopic object itself. This script repeats exactly the notebook's
steps after selection for another topic count, from the same trained models, and writes

    cistopic_obj_glut_k<K>.pkl                       cisTopic object with the K-topic model selected
    region_sets_k<K>/{Topics_otsu,Topics_top_3k,DARs_all,DARs_song-pairs}/*.bed

Parameters are copied from the notebook (binarization, imputation scale factors, HVF filter, DAR thresholds,
contrasts); only K differs. The production cistopic_obj_glut.pkl / region_sets are not touched.

    python pycistopic_topic_k.py K [n_cpu]
"""
import os
import pickle
import sys

import numpy as np

K = int(sys.argv[1])
N_CPU = int(sys.argv[2]) if len(sys.argv) > 2 else 20

MAIN = ("/hdd/jupyter/brad/scenicplus/motor-pathway_multiome/"
        "motor-pathway_multiome_seurat_cellbender.0.05_preprocess_cr/ra-arco-hvc-nc_hybrid/")
WORK = MAIN + "pycisTopic/"
OUT_OBJ = WORK + f"cistopic_obj_glut_k{K}.pkl"
OUT_SETS = WORK + f"region_sets_k{K}/"
# Ray puts a Unix socket under its temp dir and those paths cannot exceed 107 bytes, so this must be SHORT (the long
# project path fails with "AF_UNIX path length cannot exceed 107 bytes"). On /hdd rather than /tmp (small root disk),
# and separate per K so parallel runs cannot collide.
TMP = f"/hdd/brad/ray_k{K}"
os.makedirs(TMP, exist_ok=True)

from pycisTopic.diff_features import (find_diff_features, find_highly_variable_features,  # noqa: E402
                                      impute_accessibility, normalize_scores)
from pycisTopic.topic_binarization import binarize_topics  # noqa: E402
from pycisTopic.utils import region_names_to_coordinates  # noqa: E402

models = pickle.load(open(WORK + "models.pkl", "rb"))
assert any(m.n_topic == K for m in models), f"no {K}-topic model in models.pkl: {[m.n_topic for m in models]}"
model = next(m for m in models if m.n_topic == K)

obj = pickle.load(open(WORK + "cistopic_obj_glut.pkl", "rb"))
obj.add_LDA_model(model)
assert obj.selected_model.n_topic == K, obj.selected_model.n_topic
obj.cell_data["cluster"] = obj.cell_data["cluster"].astype(str)
print(f"K={K}: {obj.selected_model.n_topic} topics, {len(obj.cell_names)} cells, {len(obj.region_names)} regions", flush=True)

# topic region sets (notebook cells 47-48)
region_bin_topics_top_3k = binarize_topics(obj, method="ntop", ntop=3_000, plot=False)
region_bin_topics_otsu = binarize_topics(obj, method="otsu", plot=False)
print("binarized topics:", len(region_bin_topics_otsu), len(region_bin_topics_top_3k), flush=True)

# DARs from the imputed accessibility of THIS model (notebook cells 57-63)
imputed_acc_obj = impute_accessibility(obj, selected_cells=None, selected_regions=None, scale_factor=10**6)
normalized_imputed_acc_obj = normalize_scores(imputed_acc_obj, scale_factor=10**4)
variable_regions = find_highly_variable_features(
    normalized_imputed_acc_obj, min_disp=0.05, min_mean=0.0125, max_mean=3, max_disp=np.inf,
    n_bins=20, n_top_features=None, plot=False)
print("variable regions:", len(variable_regions), flush=True)

markers_dict_all = find_diff_features(
    obj, imputed_acc_obj, variable="cluster", var_features=variable_regions, contrasts=None,
    adjpval_thr=0.05, log2fc_thr=np.log2(1.5), n_cpu=N_CPU, _temp_dir=TMP, split_pattern="-")
contrasts = [[["Glut-CACNA1H-RA"], ["Glut-CACNA1H-1"]],
             [["Glut-CACNA1H-RA"], ["Glut-CACNA1H-2"]],
             [["Glut-DACH2-HVCra"], ["Glut-DACH2-HVCx"]]]
markers_dict_pairs = find_diff_features(
    obj, imputed_acc_obj, variable="cluster", var_features=variable_regions, contrasts=contrasts,
    adjpval_thr=0.05, log2fc_thr=np.log2(1.5), n_cpu=N_CPU, _temp_dir=TMP, split_pattern="-")
print("DARs:", {k: len(v) for k, v in list(markers_dict_all.items())[:3]}, "...", len(markers_dict_all), "sets;",
      {k: len(v) for k, v in markers_dict_pairs.items()}, flush=True)


def write_sets(subdir, sets):
    d = os.path.join(OUT_SETS, subdir)
    os.makedirs(d, exist_ok=True)
    for name, df in sets.items():
        region_names_to_coordinates(df.index).sort_values(["Chromosome", "Start", "End"]).to_csv(
            os.path.join(d, f"{name}.bed"), sep="\t", header=False, index=False)


write_sets("Topics_otsu", region_bin_topics_otsu)
write_sets("Topics_top_3k", region_bin_topics_top_3k)
write_sets("DARs_all", markers_dict_all)
write_sets("DARs_song-pairs", markers_dict_pairs)

pickle.dump(obj, open(OUT_OBJ, "wb"))
print("wrote", OUT_OBJ, "and", OUT_SETS, flush=True)
