"""SATURN input CSV for the Glut(finch)-only x Excitatory(chicken)-only subset."""
from pathlib import Path
import anndata as ad, pandas as pd, torch

BASE = Path("/private/groups/colquittlab/song-system-grn/snrna/integration/datasets/snrna-bf-adult_snrna-gg-adult-glutonly")
DATA = BASE / "data"
E = Path("/private/groups/colquittlab/saturn/embeddings")
MIN_MARGIN = 0.005

ROWS = [
    ("finch", DATA / "bf_adult_glut.h5ad", E / "GCF_005870125.1_lonStrDom2_translated_cds.gene_symbol_to_embedding_ESM1b.pt"),
    ("chicken", DATA / "gg_adult_ex.h5ad", E / "Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.pep.all.gene_symbol_to_embedding_ESM1b.pt"),
]

if __name__ == "__main__":
    embs = {s: set(torch.load(p)) for s, _, p in ROWS}
    for sp, h5, _ in ROWS:
        genes = set(ad.read_h5ad(h5, backed="r").var_names)
        own = len(genes & embs[sp]) / len(genes)
        other = max(len(genes & e) / len(genes) for s, e in embs.items() if s != sp)
        print(f"{sp:8s} {len(genes):6,} genes | own {own:6.1%} | other {other:6.1%} | margin {own-other:+.1%}")
        assert own > other + MIN_MARGIN, f"{sp}: embedding columns may be transposed"
    df = pd.DataFrame({"path": [str(h) for _, h, _ in ROWS], "species": [s for s, _, _ in ROWS],
                       "embedding_path": [str(p) for _, _, p in ROWS]})
    out = DATA / "bf_gg_glutonly_run.csv"
    df.to_csv(out, index=False)
    print(f"\nwrote {out}")
    clusters = pd.concat([pd.DataFrame({"species": s, "cluster": sorted(ad.read_h5ad(h, backed="r").obs["cluster"].unique())})
                          for s, h, _ in ROWS])
    clusters.to_csv(DATA / "bf_gg_glutonly_clusters.csv", index=False)
    print(f"wrote {DATA / 'bf_gg_glutonly_clusters.csv'} ({len(clusters)} labels)")
