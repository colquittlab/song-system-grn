"""SAMap on SCTransform-corrected input, glut-only chicken hybrid-label comparison --
variant of run_samap.py for the manuscript's own normalization method.

SAMap's preprocess_data() normally does its own sum-normalization + log-transform
(sum_norm="cell_median", norm="log") as an integral part of the algorithm -- the
downstream iterative graph refinement is tuned to operate on that representation. Here
that step is replaced upstream (samap_sct_export_inputs.py + R/sct_transform_export.R
independently SCTransform each species' FULL count matrix), so preprocess_data is called
with sum_norm=None, norm=None -- confirmed via inspecting SAM.preprocess_data's source
that this is the genuine no-op (sum_norm=1 would be WRONG: preprocess_data treats a float
sum_norm as a target total-count-per-cell to rescale to, which would undo SCTransform's
own depth correction; only sum_norm=None skips that step entirely). The remaining
preprocess_data arguments (min_expression, thresh_low/thresh_high) still apply the same
gene-presence filtering SAM normally does, now on the SCT-corrected values.

Both species now come from SCT-corrected h5ads (built by samap_sct_assemble_h5ad.py),
unlike run_samap.py where only --finch_h5ad was ever swapped and chicken stayed fixed --
here BOTH sides changed normalization, so both are CLI args.
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from samalg import SAM
from samap.mapping import SAMAP
from samap.analysis import get_mapping_scores

BASE = Path("/private/groups/colquittlab/song-system-grn/snrna/integration/gg_glutonly_hybrid/samap")
SCT_DIR = Path("/private/groups/colquittlab/song-system-grn/snrna/integration/gg_glutonly_hybrid/samap/data")


def load_sam_sct(h5ad: Path, label_col: str, name: str) -> SAM:
    a = sc.read_h5ad(h5ad)
    a.obs[label_col] = a.obs[label_col].astype(str)
    sam = SAM(counts=a)
    # sum_norm=None, norm=None: skip SAM's own normalization entirely (data is already
    # SCT-corrected) -- see module docstring for why sum_norm=1 would be wrong here.
    sam.preprocess_data(sum_norm=None, norm=None, thresh_low=0.0, thresh_high=0.96, min_expression=1)
    sam.run(preprocessing="StandardScaler", npcs=100, weight_PCs=False, k=20, n_genes=3000, seed=0)
    print(f"{name}: {sam.adata.shape}, {a.obs[label_col].nunique()} labels in '{label_col}'")
    return sam


def main(out_dir: Path, numiters: int, crossk: int, ncpus: int, seed: int,
         finch_h5ad: str, chicken_h5ad: str):
    out_dir.mkdir(parents=True, exist_ok=True)
    sams = {
        "bf": load_sam_sct(Path(finch_h5ad), "cluster", "finch"),
        "gg": load_sam_sct(Path(chicken_h5ad), "cluster", "chicken"),
    }
    keys = {"bf": "cluster", "gg": "cluster"}

    sm = SAMAP(sams, f_maps="/private/groups/colquittlab/song-system-grn/snrna/integration/datasets/samap_bf-adult_gg-adult/maps/", keys=keys)
    print(f"\nrunning SAMap (SCT input)... NUMITERS={numiters} crossK={crossk}")
    sm.run(NUMITERS=numiters, crossK=crossk, ncpus=ncpus, umap=True)

    D, MappingTable = get_mapping_scores(sm, keys, n_top=0)
    D.to_csv(out_dir / "samap_alignment_pairs.csv")
    MappingTable.to_csv(out_dir / "samap_mapping_table.csv")
    print(f"\nmapping table: {MappingTable.shape}")

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
    print("\n=== SAMap (SCT input) best chicken match per finch cluster ===")
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
    p.add_argument("--out_dir", type=Path, default=BASE / "results" / "SCT")
    p.add_argument("--numiters", type=int, default=3)
    p.add_argument("--crossk", type=int, default=20)
    p.add_argument("--ncpus", type=int, default=32)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--finch_h5ad", type=str, default=str(SCT_DIR / "finch_sct.h5ad"))
    p.add_argument("--chicken_h5ad", type=str, default=str(SCT_DIR / "chicken_sct.h5ad"))
    main(**vars(p.parse_args()))
