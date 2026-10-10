#!/usr/bin/env python3
"""Topic-number grid for the glut cells (hybrid labels): train mallet LDA models for 2-40 topics and collect their metrics.

pycistopic_glut.py trains only the 20-topic model (the earlier glut run's choice). This trains the grid that run used for model
selection, on the same 9,063-cell cisTopic object and with the same mallet settings, one process per topic count so the models run
side by side instead of one after another. It writes ONLY new files; models.pkl, region_sets/ and cistopic_obj_glut.pkl are not touched.

    python pycistopic_glut_grid.py K [threads]     # train one model -> models_grid/model_K.pkl
    python pycistopic_glut_grid.py merge           # collect -> models_grid.pkl and models_grid_metrics.csv

Start from cistopic_obj_glut_nomodel.pkl (written by pycistopic_glut.py / pycisTopic.ipynb, step 1).
"""
import os
import pickle
import sys

import matplotlib
matplotlib.use("Agg")
import pandas as pd

B = ("/hdd/jupyter/brad/scenicplus/motor-pathway_multiome/"
     "motor-pathway_multiome_seurat_cellbender.0.05_preprocess_cr/ra-arco-hvc-nc_glut_hybrid/pycisTopic/")
GRID = [2, 4, 6, 8, 10, 15, 20, 30, 40]
OUT = B + "models_grid/"
os.makedirs(OUT, exist_ok=True)

if sys.argv[1] == "merge":
    models = []
    for k in GRID:
        f = OUT + f"model_{k}.pkl"
        assert os.path.exists(f), f"model_{k}.pkl missing"
        models.extend(pickle.load(open(f, "rb")))
    pickle.dump(models, open(B + "models_grid.pkl", "wb"))
    rows = [{"n_topic": m.n_topic, **{k: float(v) for k, v in m.metrics.items()
                                      if not hasattr(v, "__len__") or isinstance(v, (int, float))}} for m in models]
    pd.DataFrame(rows).to_csv(B + "models_grid_metrics.csv", index=False)
    print(pd.DataFrame(rows).to_string(index=False))
    sys.exit()

K = int(sys.argv[1])
threads = int(sys.argv[2]) if len(sys.argv) > 2 else 5
# one tmp dir per K (mallet writes a corpus and state files there); Ray is not used by mallet training
tmp = f"/hdd/brad/ray_glut_grid/k{K}"
os.makedirs(tmp, exist_ok=True)
os.environ["MALLET_MEMORY"] = "48G"
from pycisTopic.lda_models import run_cgs_models_mallet  # noqa: E402

obj = pickle.load(open(B + "cistopic_obj_glut_nomodel.pkl", "rb"))
print(f"K={K}: {len(obj.cell_names)} cells, {len(obj.region_names)} regions, {threads} threads", flush=True)
models = run_cgs_models_mallet(
    obj, n_topics=[K], n_cpu=threads, n_iter=500, random_state=555, alpha=50, alpha_by_topic=True,
    eta=0.1, eta_by_topic=False, tmp_path=tmp, save_path=tmp,
    mallet_path="/home/brad/ssd/repos/Mallet-202108/bin/mallet")
pickle.dump(models, open(OUT + f"model_{K}.pkl", "wb"))
print(f"K={K}: done", flush=True)
