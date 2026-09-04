"""Glutamatergic/excitatory-only subset of the adult finch x adult chicken pair, using the
UPDATED hybrid-labeled finch data (celltype_hybrid, 47 types) in place of the original flat
'cluster' scheme (49 types) -- repeat of prepare_subset.py's finch-side filtering logic
against bf_adult_hybrid.h5ad instead of bf_adult.h5ad. Chicken side is untouched (same
gg_adult_ex.h5ad reused from the non-hybrid glutonly run).

finch: 'cluster' (== celltype_hybrid) starts with "Glut-", EXCLUDING Glut-Im, Glut-NB,
       Glut-NSC, Glut-GABA, Glut-SATB2-1 -- the five hybrid types the hybrid dataset's own
       docstring flags as "new developmental/ambiguous populations ... with no clear
       analogue in the original scheme" (analogous to how prepare_subset.py dropped the
       developmental Glut-Pre-1/2/3 and Glut-Nido-2/3 populations from the old scheme).
       No collapsing is needed here -- unlike the old scheme's Glut-Nido-1 -> Glut-NC-1
       merge, celltype_hybrid's song-nucleus populations are already the manuscript's
       refined resolution.
chicken: unchanged -- reuses gg_adult_ex.h5ad (class == "Excitatory neurons", Ex_TCF7L2
         dropped), written by the original prepare_subset.py.
"""
import anndata as ad

FULL = "/private/groups/colquittlab/song-system-grn/snrna/integration/datasets/snrna-bf-adult_snrna-gg-adult/data"
OUT = "/private/groups/colquittlab/song-system-grn/snrna/integration/datasets/snrna-bf-adult_snrna-gg-adult-glutonly/data"

DROP_FINCH_HYBRID = ["Glut-Im", "Glut-NB", "Glut-NSC", "Glut-GABA", "Glut-SATB2-1"]

b = ad.read_h5ad(f"{FULL}/bf_adult_hybrid.h5ad")
cl = b.obs["cluster"].astype(str)
glut = b[cl.str.startswith("Glut-") & ~cl.isin(DROP_FINCH_HYBRID)].copy()
print(f"finch (hybrid): {b.n_obs} -> {glut.n_obs} cells, {glut.obs['cluster'].nunique()} clusters "
      f"(dropped {DROP_FINCH_HYBRID})")
print(glut.obs["cluster"].value_counts().to_string())
glut.write_h5ad(f"{OUT}/bf_adult_glut_hybrid.h5ad", compression="gzip")

print(f"\nchicken side unchanged -- reusing {OUT}/gg_adult_ex.h5ad")
g = ad.read_h5ad(f"{OUT}/gg_adult_ex.h5ad", backed="r")
print(f"chicken: {g.n_obs} cells, {g.obs['cluster'].nunique()} clusters")
