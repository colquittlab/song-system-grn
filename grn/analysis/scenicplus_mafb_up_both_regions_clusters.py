"""Summaries behind the region x cell type accessibility heatmap: RA vs PVALB-1 accessibility correlation across the regions, regions open in both (z > 1) against
independence, and the correlation tree cut into 8 clusters. Run: python scenicplus_mafb_up_both_regions_clusters.py"""
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform

d = pd.read_csv("/ssd/brad/rstudio/multiome/song-system-grn/.claude/worktrees/grn-k40-nosongpairs/grn/analysis/scenicplus_mafb_up_both_regions_accessibility_config37.csv").set_index("region")
Z = d[[c for c in d.columns if c.startswith("z_")]]
Z.columns = [c[2:] for c in Z.columns]
ra, c1, pv1, pv2, sst, lamp = "Glut-CACNA1H-RA", "Glut-CACNA1H-1", "GABA-MGE-PVALB-1", "GABA-MGE-PVALB-2", "GABA-MGE-SST-1", "GABA-MGE-LAMP5"
print("regions:", len(Z))
print("Spearman z(RA) vs z(PVALB-1) across these regions: %.2f ; vs PVALB-2: %.2f ; RA vs C1H-1: %.2f" % (
    Z[ra].corr(Z[pv1], method="spearman"), Z[ra].corr(Z[pv2], method="spearman"), Z[ra].corr(Z[c1], method="spearman")))
for t in (1.0, 0.5):
    A, B = Z[ra] > t, Z[pv1] > t
    print(f"z > {t}: RA {int(A.sum())}, PVALB-1 {int(B.sum())}, both {int((A & B).sum())} (expected if independent {A.sum() * B.sum() / len(Z):.1f})")
A, B = Z[ra] > 1, Z[pv1] > 1
both = d[A & B]
print("regions open (z>1) in RA and PVALB-1: DAR flags -> RA DAR", int(both.DAR_RA.sum()), "| PV1 DAR", int(both.DAR_PV1.sum()), "| both DAR", int((both.DAR_RA & both.DAR_PV1).sum()), "of", len(both))
print("genes among those regions:", len(set(";".join(both.genes).split(";"))), "->", ", ".join(sorted(set(";".join(both.genes).split(";"))))[:300])
print("PVALB-1 open (z>1) but NOT a PV1 DAR:", int((B & (d.DAR_PV1 == 0)).sum()), "of", int(B.sum()), "| RA open but not an RA DAR:", int((A & (d.DAR_RA == 0)).sum()), "of", int(A.sum()))

# the correlation tree of the heatmap, cut into clusters
D = 1 - np.corrcoef(Z.values)
np.fill_diagonal(D, 0)
lk = linkage(squareform(np.clip(D, 0, None), checks=False), "average")
k = 8
lab = fcluster(lk, k, "maxclust")
rows = []
for c in sorted(set(lab)):
    m = lab == c
    z = Z[m]
    rows.append(dict(cluster=c, n_regions=int(m.sum()), n_genes=len(set(";".join(d.genes[m]).split(";"))), RA=z[ra].mean().round(2), C1H1=z[c1].mean().round(2),
                     PVALB1=z[pv1].mean().round(2), PVALB2=z[pv2].mean().round(2), SST1=z[sst].mean().round(2), LAMP5=z[lamp].mean().round(2),
                     DAR_RA=int(d.DAR_RA[m].sum()), DAR_PV1=int(d.DAR_PV1[m].sum()), top=z.mean().sort_values(ascending=False).index[0]))
pd.set_option("display.width", 220)
print("\ncorrelation tree cut into", k, "clusters (mean z by cell type):")
print(pd.DataFrame(rows).to_string(index=False))

for c in (6, 4, 5):
    m = lab == c
    gs = pd.Series(";".join(d.genes[m]).split(";")).value_counts()
    print(f"\ncluster {c}: genes (regions each):", ", ".join(f"{g}({n})" for g, n in gs.items()))
