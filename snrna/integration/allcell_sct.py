"""SCTransform-normalized inputs for the all-cell cross-species integrations
(finch x chicken = gg, finch x mouse = yao), mirroring the glut-only SCT pipeline so all
integrations share one normalization.

  gsi-export  --pair gg|yao      OMA 1:1 orthologs; robust genes on the complete data, cluster
                                  subsampling and UMI thinning exactly as run_gsi.py, raw counts
                                  written for SCTransform (gg_glutonly_hybrid/gsi/gsi_sct_export_inputs.py)
  gsi-run     --pair gg|yao      GSI Spearman correlation on SCT-corrected cluster means
                                  (gg_glutonly_hybrid/gsi/run_gsi_sct.py) ->
                                  composite_scoring/results/gsi_corr_<pair>_adult_hybrid_SCT.csv
  samap-export   --species bf|gg|mm   full raw counts (all cells, all genes) for SCTransform
  samap-assemble --species bf|gg|mm   SCT-corrected counts + original obs -> h5ad for SAMap

SCTransform itself is run between export and the next step with the toolkit's
R/sct_transform_export.R (<prefix> -> <prefix>_sct); see run_allcell_sct.sh.
"""
import argparse, importlib.util, sys
from pathlib import Path
import anndata as ad
import numpy as np
import pandas as pd
import scipy.io as sio
import scipy.sparse as sp
from scipy.io import mmwrite

INT = Path(__file__).resolve().parent
LIVE = Path("/private/groups/colquittlab/song-system-grn/snrna/integration")   # data live here
CS = LIVE / "composite_scoring"
SCT = CS / "data" / "sct"
H5 = {"bf": LIVE / "datasets/snrna-bf-adult_snrna-gg-adult/data/bf_adult_hybrid.h5ad",
      "gg": LIVE / "datasets/snrna-bf-adult_snrna-gg-adult/data/gg_adult.h5ad",
      "mm": Path("/private/groups/colquittlab/saturn/snrna-bf-dev_snrna-yao2023/data/mm_yao.h5ad")}
LABEL = {"bf": "celltype_hybrid", "gg": "cluster", "mm": "cluster"}
PAIR = {"gg": ("gg", CS / "data/orthologs_oma_bfgg_adult_aliased.csv"),
        "yao": ("mm", CS / "data/orthologs_oma_bfmm_adult_aliased.csv")}
RENAME = {"Glut-DACH2-HVCra-Pre": "Glut-DACH2-HVCra-Int"}


def _load(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m); return m


def gsi_export(pair: str):
    ex = _load(INT / "gg_glutonly_hybrid/gsi/gsi_sct_export_inputs.py", "gsi_sct_export_inputs")
    ref, orth_csv = PAIR[pair]
    out = SCT / f"gsi_{pair}"; out.mkdir(parents=True, exist_ok=True)
    orth = pd.read_csv(orth_csv)
    q = ad.read_h5ad(H5["bf"], backed="r"); r = ad.read_h5ad(H5[ref], backed="r")
    orth = orth[orth.finch.isin(set(q.var_names)) & orth.mouse.isin(set(r.var_names))]
    print(f"{pair}: {len(orth):,} OMA 1:1 orthologs present in both datasets")
    orth.to_csv(out / "orthologs_used.csv", index=False)
    okq = ex.export_species(H5["bf"], LABEL["bf"], list(orth.finch), "finch (query)", out / "finch")
    okr = ex.export_species(H5[ref], LABEL[ref], list(orth.mouse), f"{ref} (reference)", out / "ref")
    both = okq & okr
    print(f"genes robust in BOTH species: {both.sum():,}/{len(both):,}")
    np.save(out / "both_robust_mask.npy", both)


def gsi_run(pair: str):
    rg = _load(INT / "gg_glutonly_hybrid/gsi/run_gsi_sct.py", "run_gsi_sct")
    d = SCT / f"gsi_{pair}"
    orth = pd.read_csv(d / "orthologs_used.csv"); both = np.load(d / "both_robust_mask.npy")
    Mq = rg.load_sct_means(d / "finch_sct", list(orth.finch), d / "finch_labels.csv")
    Mr = rg.load_sct_means(d / "ref_sct", list(orth.mouse), d / "ref_labels.csv")
    Mq, Mr = Mq.loc[np.asarray(orth.finch)[both]], Mr.loc[np.asarray(orth.mouse)[both]]
    Gq, Gr = rg.gsi(Mq), rg.gsi(Mr)
    shared = min(len(Gq), len(Gr)); Gq, Gr = Gq.iloc[:shared], Gr.iloc[:shared]
    Rq, Rr = Gq.rank(axis=0).values, Gr.rank(axis=0).values
    Zq = (Rq - Rq.mean(0)) / np.where(Rq.std(0) == 0, 1, Rq.std(0))
    Zr = (Rr - Rr.mean(0)) / np.where(Rr.std(0) == 0, 1, Rr.std(0))
    C = pd.DataFrame(Zq.T @ Zr / Zq.shape[0], index=Gq.columns, columns=Gr.columns).rename(index=RENAME)
    out = CS / "results" / f"gsi_corr_{pair}_adult_hybrid_SCT.csv"
    C.round(5).to_csv(out)
    print(f"GSI (SCT) {C.shape}, {shared:,} genes -> {out}; range {C.values.min():+.3f}..{C.values.max():+.3f}")


def samap_export(species: str):
    out = SCT / "samap"; out.mkdir(parents=True, exist_ok=True)
    a = ad.read_h5ad(H5[species])
    X = a.X.tocsr() if sp.issparse(a.X) else sp.csr_matrix(a.X)
    prefix = out / species
    mmwrite(str(prefix) + "_matrix.mtx", X.T.tocoo())     # genes x cells
    pd.Series(a.var_names).to_csv(str(prefix) + "_genes.tsv", index=False, header=False)
    pd.Series(a.obs_names).to_csv(str(prefix) + "_cells.tsv", index=False, header=False)
    print(f"{species}: wrote {X.shape[0]:,} cells x {X.shape[1]:,} genes -> {prefix}_*")


def samap_assemble(species: str):
    prefix = SCT / "samap" / f"{species}_sct"
    X = sio.mmread(str(prefix) + "_matrix.mtx").tocsr().T.tocsr()     # cells x genes
    genes = pd.read_csv(str(prefix) + "_genes.tsv", header=None)[0].values
    cells = pd.read_csv(str(prefix) + "_cells.tsv", header=None)[0].astype(str).values
    obs = ad.read_h5ad(H5[species], backed="r").obs.copy()
    obs.index = obs.index.astype(str); obs = obs.loc[cells]
    obs["cluster"] = obs[LABEL[species]].astype(str).replace(RENAME)   # SAMap key column
    out = ad.AnnData(X=X, obs=obs, var=pd.DataFrame(index=genes)); out.var_names_make_unique()
    out.write_h5ad(SCT / "samap" / f"{species}_sct.h5ad", compression="gzip")
    print(f"{species}: {out.shape}, {out.obs.cluster.nunique()} labels -> {SCT / 'samap' / f'{species}_sct.h5ad'}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("step", choices=["gsi-export", "gsi-run", "samap-export", "samap-assemble"])
    ap.add_argument("--pair", choices=list(PAIR))
    ap.add_argument("--species", choices=list(H5))
    a = ap.parse_args()
    {"gsi-export": lambda: gsi_export(a.pair), "gsi-run": lambda: gsi_run(a.pair),
     "samap-export": lambda: samap_export(a.species), "samap-assemble": lambda: samap_assemble(a.species)}[a.step]()
