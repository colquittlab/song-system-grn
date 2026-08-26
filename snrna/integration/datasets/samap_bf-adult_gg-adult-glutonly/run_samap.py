"""SAMap: developing finch (v3 reclustered) x developing mouse (La Manno et al. 2021).

Run as an independent check on the SATURN integration. The two methods differ in ways
that matter for this project:

  * homology comes from reciprocal blastp bit scores here, vs ESM protein-language-model
    embeddings in SATURN;
  * SAMap iteratively reweights the gene-gene homology graph using expression
    correlation in mapped cells, a feedback loop SATURN has no direct analogue for;
  * critically, SAMap never consumes the cell-type labels -- they are used only to
    report cluster alignment scores afterwards. SATURN's triplet loss does use them, so
    for judging the identity of an UNKNOWN finch cluster SAMap is the less circular test.

Alignment score (get_mapping_scores) is the cross-species mutual-nearest-neighbour edge
fraction per cluster pair, which is what the SATURN pipeline reconstructs via KNN label
transfer + mixing_index.
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from samalg import SAM
from samap.mapping import SAMAP
from samap.analysis import get_mapping_scores

BASE = Path("/private/groups/colquittlab/song-system-grn/snrna/integration/datasets/samap_bf-adult_gg-adult-glutonly")
LM = Path("/private/groups/colquittlab/song-system-grn/snrna/integration/datasets/snrna-bf-adult_snrna-gg-adult-glutonly/data")
GG = Path("/private/groups/colquittlab/song-system-grn/snrna/integration/datasets/snrna-bf-adult_snrna-gg-adult-glutonly/data")


def load_sam(h5ad: Path, label_col: str, name: str, preprocess: bool) -> SAM:
    a = sc.read_h5ad(h5ad)
    a.obs[label_col] = a.obs[label_col].astype(str)
    # SAMap expects raw-ish counts and runs its own SAM preprocessing/feature weighting.
    sam = SAM(counts=a)
    if preprocess:
        sam.preprocess_data(sum_norm="cell_median", norm="log",
                            thresh_low=0.0, thresh_high=0.96, min_expression=1)
        sam.run(preprocessing="StandardScaler", npcs=100, weight_PCs=False,
                k=20, n_genes=3000, seed=0)
    print(f"{name}: {sam.adata.shape}, {a.obs[label_col].nunique()} labels in '{label_col}'")
    return sam


def main(out_dir: Path, numiters: int, crossk: int, ncpus: int, seed: int, finch_h5ad: str, chicken_h5ad: str):
    out_dir.mkdir(parents=True, exist_ok=True)
    sams = {
        "bf": load_sam(Path(finch_h5ad), "cluster", "finch", True),
        "gg": load_sam(Path(chicken_h5ad), "cluster", "chicken", True),
    }
    keys = {"bf": "cluster", "gg": "cluster"}

    sm = SAMAP(sams, f_maps="/private/groups/colquittlab/song-system-grn/snrna/integration/datasets/samap_bf-adult_gg-adult/maps/", keys=keys)
    print(f"\nrunning SAMap... NUMITERS={numiters} crossK={crossk}")
    sm.run(NUMITERS=numiters, crossK=crossk, ncpus=ncpus, umap=True)

    # cluster-level alignment scores
    D, MappingTable = get_mapping_scores(sm, keys, n_top=0)
    D.to_csv(out_dir / "samap_alignment_pairs.csv")
    MappingTable.to_csv(out_dir / "samap_mapping_table.csv")
    print(f"\nmapping table: {MappingTable.shape}")

    # per-finch-cluster best mouse match
    MT = MappingTable.copy()
    bf_rows = [i for i in MT.index if str(i).startswith("bf_")]
    mm_cols = [c for c in MT.columns if str(c).startswith("gg_")]
    sub = MT.loc[bf_rows, mm_cols]
    best = pd.DataFrame({
        "finch_cluster": [str(i)[3:] for i in sub.index],
        "samap_top": [str(sub.columns[j])[3:] for j in np.argmax(sub.values, axis=1)],
        "samap_score": sub.values.max(axis=1).round(4),
    })
    order = np.argsort(-sub.values, axis=1)
    best["samap_top2"] = [str(sub.columns[o[1]])[3:] if sub.shape[1] > 1 else "" for o in order]
    best["samap_score2"] = [round(sub.values[i, o[1]], 4) if sub.shape[1] > 1 else np.nan
                            for i, o in enumerate(order)]
    best = best.sort_values("samap_score", ascending=False)
    best.to_csv(out_dir / "samap_finch_calls.csv", index=False)

    pd.set_option("display.width", 200)
    print("\n=== SAMap best mouse match per finch cluster ===")
    print(best.to_string(index=False))

    sm.samap.adata.write_h5ad(out_dir / "samap_joint.h5ad")
    try:
        import dill
        with open(out_dir / "samap_object.pkl", "wb") as f:
            dill.dump(sm, f)
    except Exception as e:
        print(f"(could not pickle SAMAP object: {e})")
    print(f"\nwrote -> {out_dir}")


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--out_dir", type=Path, default=BASE / "results")
    p.add_argument("--numiters", type=int, default=3)
    p.add_argument("--crossk", type=int, default=20)
    p.add_argument("--ncpus", type=int, default=32)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--finch_h5ad", type=str, default=str(LM / "bf_adult_glut.h5ad"))
    p.add_argument("--chicken_h5ad", type=str, default=str(GG / "gg_adult_ex.h5ad"))
    main(**vars(p.parse_args()))
