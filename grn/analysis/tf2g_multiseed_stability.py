#!/usr/bin/env python3
"""Seed stability of the SCENIC+ TF-to-gene step for AR and NR3C2 (RA vs C1H-1, all runs, config2).

SCENIC+ fits one early-stopped gradient-boosting model per target gene (arboreto GRNBoost2, `infer_partial_network`, SGBM_KWARGS, seed 666) with
the TFs as predictors. The importances are noisy: reseeding alone moves a TF's share by ~30% and can change its rank for a gene substantially. This repeats the
step with several seeds on the raw-count expression matrix used by the run, and asks which AR / NR3C2 links are stable.

  baseline   every gene, N_SEEDS seeds, all TFs as predictors
  variants   AR's and NR3C2's eRegulon targets plus MAFB, N_SEEDS seeds, once without NR3C2 and once without AR (drop-one)

The TF list is approximated by the TFs of the run's cistromes that are in the expression matrix (tfs.txt was not transferred). Per (gene, seed, model) we keep AR's and
NR3C2's importance and rank; for MAFB the whole predictor list.  Outputs: OUT_DIR/*.csv (large, outside git) and tf2g_stability_*.csv next to this script (small).

    python tf2g_multiseed_stability.py [n_seeds=10] [n_cpu=40]
"""
import multiprocessing as mp
import os
import sys
import time
import warnings

import anndata
import mudata
import numpy as np
import pandas as pd
from arboreto.core import EARLY_STOP_WINDOW_LENGTH, SGBM_KWARGS, infer_partial_network
from scipy.stats import rankdata

warnings.filterwarnings("ignore")
HERE = os.path.dirname(os.path.abspath(__file__))
N_SEEDS = int(sys.argv[1]) if len(sys.argv) > 1 else 10
N_CPU = int(sys.argv[2]) if len(sys.argv) > 2 else 40
SEEDS = [666 + i for i in range(N_SEEDS)]
ROOT = ("/hdd/jupyter/brad/scenicplus/motor-pathway_multiome/motor-pathway_multiome_seurat_cellbender.0.05_preprocess_cr/"
        "glut-ra_glut-c1h-1_allruns_hybrid/results/config2/")
OUT_DIR = ROOT + "tf2g_stability/"
os.makedirs(OUT_DIR, exist_ok=True)

rna = mudata.read(ROOT + "outs/scplusmdata_perturb.h5mu/scRNA_counts")
df = rna.to_df()
GENES = list(df.columns)
GIDX = {g: i for i, g in enumerate(GENES)}
X = df.to_numpy(dtype=np.float32)
cis = [anndata.read_h5ad(ROOT + f"outs/cistromes_{k}.h5ad").var_names for k in ("direct", "extended")]
TFS = sorted((set(cis[0]) | set(cis[1])) & set(GENES))
TF_COLS = {t: GIDX[t] for t in TFS}
er = pd.concat([pd.read_table(ROOT + "outs/" + f, usecols=["TF", "Gene"]) for f in ("eRegulon_direct.tsv", "eRegulons_extended.tsv")])
clean = lambda s: [g for g in sorted(s) if g in GIDX and "LOC" not in g]
AR_T, NR_T = clean(set(er[er.TF == "AR"].Gene) - {"AR"}), clean(set(er[er.TF == "NR3C2"].Gene) - {"NR3C2"})
SUBSET = sorted(set(AR_T) | set(NR_T) | {"MAFB"})
print(f"{len(GENES)} genes, {len(TFS)} TF predictors, {N_SEEDS} seeds; variant subset {len(SUBSET)} genes", flush=True)


def fit(task):
    gene, seed, drop, tag = task
    use = [t for t in TFS if t != gene and t not in drop]
    a = infer_partial_network(regressor_type="GBM", regressor_kwargs=SGBM_KWARGS, tf_matrix=X[:, [TF_COLS[t] for t in use]], tf_matrix_gene_names=use,
                              target_gene_name=gene, target_gene_expression=X[:, GIDX[gene]], include_meta=False,
                              early_stop_window_length=EARLY_STOP_WINDOW_LENGTH, seed=seed)
    s = a.set_index("TF").importance
    rk = s.rank(ascending=False)
    row = dict(gene=gene, seed=seed, tag=tag, total=float(s.sum()), n_pred=len(s))
    for t in ("AR", "NR3C2"):
        row[f"{t}_imp"] = float(s.get(t, 0.0)) if t in use else np.nan      # 0 = predictor present but not selected; nan = removed
        row[f"{t}_rank"] = float(rk.get(t, np.nan))
    full = s.reset_index().assign(gene=gene, seed=seed, tag=tag, rank=rk.values) if gene == "MAFB" else None
    return row, full


if __name__ == "__main__":
    tasks = [(g, s, (), "full") for g in GENES for s in SEEDS]
    tasks += [(g, s, ("NR3C2",), "noNR3C2") for g in SUBSET for s in SEEDS] + [(g, s, ("AR",), "noAR") for g in SUBSET for s in SEEDS]
    print(len(tasks), "fits", flush=True)
    rows, mafb = [], []
    t0 = time.time()
    with mp.get_context("fork").Pool(N_CPU) as pool:
        for k, (row, full) in enumerate(pool.imap_unordered(fit, tasks, chunksize=64)):
            rows.append(row)
            if full is not None:
                mafb.append(full)
            if (k + 1) % 20000 == 0:
                print(f"  {k + 1}/{len(tasks)} fits, {time.time() - t0:.0f}s", flush=True)
    R = pd.DataFrame(rows)
    R.to_csv(OUT_DIR + "tf2g_ar_nr3c2_by_gene_seed.csv", index=False)
    M = pd.concat(mafb, ignore_index=True)
    M.to_csv(OUT_DIR + "tf2g_mafb_predictors_by_seed.csv", index=False)

    # seed-independent correlation (Spearman across cells) of AR and NR3C2 with every gene, as SCENIC+ adds to the adjacencies
    rx = np.apply_along_axis(rankdata, 0, X).astype(np.float32)
    rx = (rx - rx.mean(0)) / (rx.std(0) + 1e-9)
    rho = pd.DataFrame({t: (rx[:, GIDX[t]] @ rx) / X.shape[0] for t in ("AR", "NR3C2")}, index=GENES)
    rho.to_csv(OUT_DIR + "rho_ar_nr3c2_by_gene.csv")
    print(f"done: {len(R)} fits in {time.time() - t0:.0f}s -> {OUT_DIR}", flush=True)
