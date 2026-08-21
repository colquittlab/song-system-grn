"""SATURN input CSV for the adult song-system finch (hybrid celltype_hybrid labels,
lonStrDom2) x Yao et al. 2023 mouse telencephalon (ABC Atlas, GRCm39).

Finch side reuses bf_adult_hybrid.h5ad directly (no copy) -- same object used for the
bf-adult x gg-adult-hybrid pairing. Mouse side reuses mm_yao.h5ad directly, unchanged
from the existing bf-mature(v3 dev) x Yao pipeline.

Embedding choice: this finch object's var_names are lonStrDom2 (GCF_005870125.1) gene
IDs, NOT lonStrDom3 -- use the lonStrDom2 ESM1b embedding (matches bf-adult x gg-adult's
make_run_csv.py), not the lonStrDom3_plus_v2fallback embedding used by the unrelated
v3-reclustered developmental finch pipeline (different genome build, would silently
mismatch almost every gene).
"""
from pathlib import Path
import anndata as ad, pandas as pd, torch

BASE = Path("/private/groups/colquittlab/song-system-grn/snrna/integration/datasets/snrna-bf-adult_snrna-yao")
FINCH_H5AD = Path("/private/groups/colquittlab/song-system-grn/snrna/integration/datasets/snrna-bf-adult_snrna-gg-adult/data/bf_adult_hybrid.h5ad")
MOUSE_H5AD = Path("/private/groups/colquittlab/saturn/snrna-bf-dev_snrna-yao2023/data/mm_yao.h5ad")
E = Path("/private/groups/colquittlab/saturn/embeddings")
MIN_MARGIN = 0.005

ROWS = [
    ("finch", FINCH_H5AD, E / "GCF_005870125.1_lonStrDom2_translated_cds.gene_symbol_to_embedding_ESM1b.pt"),
    ("mouse", MOUSE_H5AD, E / "Mus_musculus.GRCm39.pep.all.gene_symbol_to_embedding_ESM1b.pt"),
]

if __name__ == "__main__":
    embs = {s: set(torch.load(p)) for s, _, p in ROWS}
    for sp, h5, _ in ROWS:
        genes = set(ad.read_h5ad(h5, backed="r").var_names)
        own = len(genes & embs[sp]) / len(genes)
        other = max(len(genes & e) / len(genes) for s, e in embs.items() if s != sp)
        print(f"{sp:8s} {len(genes):6,} genes | own {own:6.1%} | other {other:6.1%} | margin {own-other:+.1%}")
        assert own > other + MIN_MARGIN, f"{sp}: embedding columns may be transposed"

    df = pd.DataFrame({"path": [str(h) for _, h, _ in ROWS],
                       "species": [s for s, _, _ in ROWS],
                       "embedding_path": [str(p) for _, _, p in ROWS]})
    out = BASE / "data" / "bf_yao_adult_hybrid_run.csv"
    df.to_csv(out, index=False)
    print(f"\nwrote {out}")

    clusters = pd.concat([
        pd.DataFrame({"species": s, "cluster": sorted(ad.read_h5ad(h, backed="r").obs["cluster"].unique())})
        for s, h, _ in ROWS])
    clusters.to_csv(BASE / "data" / "bf_yao_adult_hybrid_clusters.csv", index=False)
    print(f"wrote {BASE / 'data' / 'bf_yao_adult_hybrid_clusters.csv'} ({len(clusters)} labels)")
