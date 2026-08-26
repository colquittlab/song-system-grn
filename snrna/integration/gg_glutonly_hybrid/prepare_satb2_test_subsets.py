"""Two follow-up subsets testing whether Glut-DACH2-HVCra's strong SCT-corrected match to
Ex_SATB2_ZNF385B/Ex_SATB2_SOX6 reflects HVCra's own biology, or is an artifact of
Glut-SATB2-1 (100% identical to the OLD scheme's Glut-NR-5, confirmed by cell-barcode
crosswalk) being excluded from the glut-only hybrid finch subset -- RPCA has no "reject"
option, so removing a population's true match doesn't remove the pull on the chicken
side, it just reassigns it to whatever finch population is next-best.

1. bf_adult_glut_hybrid_plusSATB2.h5ad: the standard glut-only hybrid finch subset
   (Glut-* excluding Im/NB/NSC/GABA) but WITH Glut-SATB2-1 added back in -- if HVCra's
   match is an artifact, restoring its "competitor" for the SATB2-family chicken types
   should visibly dilute HVCra's score. Chicken side unchanged (gg_adult_ex.h5ad).

2. gg_adult_ex_noMeso.h5ad: the standard chicken excitatory reference with the
   mesopallial/SATB2-family types removed entirely (any cluster matching Ex_SATB2* or
   Ex_KIAA1217*) -- if HVCra's SATB2-family match is real (not just "closest available
   option"), removing those targets should force HVCra onto a clearly worse/lower-scoring
   match rather than finding an equally strong new one. Finch side unchanged (the
   standard bf_adult_glut_hybrid.h5ad, WITHOUT Glut-SATB2-1, matching the main analysis).
"""
import anndata as ad

FULL = "/private/groups/colquittlab/song-system-grn/snrna/integration/datasets/snrna-bf-adult_snrna-gg-adult/data"
GLUTONLY = "/private/groups/colquittlab/song-system-grn/snrna/integration/datasets/snrna-bf-adult_snrna-gg-adult-glutonly/data"

# --- 1. finch: add Glut-SATB2-1 back in ---
DROP_FINCH_MINUS_SATB2 = ["Glut-Im", "Glut-NB", "Glut-NSC", "Glut-GABA"]  # SATB2-1 NOT dropped this time

b = ad.read_h5ad(f"{FULL}/bf_adult_hybrid.h5ad")
cl = b.obs["cluster"].astype(str)
glut = b[cl.str.startswith("Glut-") & ~cl.isin(DROP_FINCH_MINUS_SATB2)].copy()
print(f"finch +SATB2-1: {b.n_obs} -> {glut.n_obs} cells, {glut.obs['cluster'].nunique()} clusters "
      f"(dropped {DROP_FINCH_MINUS_SATB2}, KEPT Glut-SATB2-1)")
print(glut.obs["cluster"].value_counts().to_string())
glut.write_h5ad(f"{GLUTONLY}/bf_adult_glut_hybrid_plusSATB2.h5ad", compression="gzip")

# --- 2. chicken: remove mesopallial/SATB2-family types ---
g = ad.read_h5ad(f"{GLUTONLY}/gg_adult_ex.h5ad")
gcl = g.obs["cluster"].astype(str)
is_meso = gcl.str.startswith("Ex_SATB2") | gcl.str.startswith("Ex_KIAA1217")
print(f"\nchicken clusters being removed as mesopallial: {sorted(gcl[is_meso].unique())}")
ex_nomeso = g[~is_meso].copy()
print(f"chicken no-meso: {g.n_obs} -> {ex_nomeso.n_obs} cells, {ex_nomeso.obs['cluster'].nunique()} clusters "
      f"(was {gcl.nunique()})")
ex_nomeso.write_h5ad(f"{GLUTONLY}/gg_adult_ex_noMeso.h5ad", compression="gzip")
