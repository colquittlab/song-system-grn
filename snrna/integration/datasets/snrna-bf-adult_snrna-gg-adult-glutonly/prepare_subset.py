"""Glutamatergic/excitatory-only subset of the adult finch x adult chicken pair, with
the manuscript's cluster refinements applied.

Reviewer-requested test: does an RPCA finding (song nucleus Glut types -- HVC, RA, LMAN,
Area X -- fail to integrate onto any single well-defined chicken excitatory type, while
non-song Glut types integrate well) replicate under other integration methods restricted
to the same Glut-vs-excitatory subset. Both sides use the SAME h5ads/embeddings/BLAST
maps as the full adult comparison -- only the cell subset/labels change.

finch:   cluster starts with "Glut-", EXCLUDING Glut-Nido-2, Glut-Nido-3, and the three
         Glut-Pre-* clusters (dropped in the manuscript's later cluster refinement);
         Glut-Nido-1 is relabelled to Glut-NC-1 (collapsed together in that refinement,
         not merely excluded).
chicken: class == "Excitatory neurons", EXCLUDING Ex_TCF7L2.
"""
import anndata as ad
import pandas as pd

FULL = "/private/groups/colquittlab/song-system-grn/snrna/integration/datasets/snrna-bf-adult_snrna-gg-adult/data"
OUT = "/private/groups/colquittlab/song-system-grn/snrna/integration/datasets/snrna-bf-adult_snrna-gg-adult-glutonly/data"

DROP_FINCH = ["Glut-Nido-2", "Glut-Nido-3", "Glut-Pre-1", "Glut-Pre-2", "Glut-Pre-3"]
COLLAPSE_FINCH = {"Glut-Nido-1": "Glut-NC-1"}
DROP_CHICKEN = ["Ex_TCF7L2"]

b = ad.read_h5ad(f"{FULL}/bf_adult.h5ad")
cl = b.obs["cluster"].astype(str)
glut = b[cl.str.startswith("Glut-") & ~cl.isin(DROP_FINCH)].copy()
glut.obs["cluster"] = glut.obs["cluster"].astype(str).replace(COLLAPSE_FINCH)
print(f"finch: {b.n_obs} -> {glut.n_obs} cells, {glut.obs['cluster'].nunique()} clusters "
      f"(dropped {DROP_FINCH}, collapsed {COLLAPSE_FINCH})")
print(glut.obs["cluster"].value_counts().to_string())
glut.write_h5ad(f"{OUT}/bf_adult_glut.h5ad", compression="gzip")

g = ad.read_h5ad(f"{FULL}/gg_adult.h5ad")
gcl = g.obs["cluster"].astype(str)
ex = g[(g.obs["class"].astype(str) == "Excitatory neurons") & ~gcl.isin(DROP_CHICKEN)].copy()
ex.obs["cluster"] = ex.obs["cluster"].astype(str)
print(f"\nchicken: {g.n_obs} -> {ex.n_obs} cells, {ex.obs['cluster'].nunique()} clusters "
      f"(dropped {DROP_CHICKEN})")
ex.write_h5ad(f"{OUT}/gg_adult_ex.h5ad", compression="gzip")
