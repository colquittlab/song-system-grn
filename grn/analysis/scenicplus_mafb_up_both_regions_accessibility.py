#!/usr/bin/env python3
"""Pseudobulk accessibility (regions x cell types) of the regions behind the genes that are up in BOTH RA and PVALB-1, whatever their DAR status.

The DAR overlap (scenicplus_mafb_up_both_dar_overlap.py) is limited by the power of the PVALB-1 vs LAMP5 test (309 vs 167 cells). Here no
significance filter is applied: every region linked to the gene set is kept and its accessibility is looked at across all cell types, so
patterns that fall short of a DAR threshold can still show up when the regions are clustered (scenicplus_mafb_region_accessibility_heatmap.R).

Gene set: MAFB +/+ targets of config37 with log2FC > 1 in both RA vs C1H-1 and PVALB-1 vs LAMP5. Regions: the region-gene links of MAFB's +/+
eRegulons (config37), on the 499,348-region consensus of the 40-topic cisTopic object (its fragment matrix is the source of the counts).
Accessibility: summed fragments per cell type, normalized to CPM, log2(CPM + 1), then z-scored per region across the cell types (>= 50 cells).

    python scenicplus_mafb_up_both_regions_accessibility.py
"""
import pickle
import re
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

HERE = Path(__file__).resolve().parent
STORE = Path("/hdd/jupyter/brad/scenicplus/motor-pathway_multiome/motor-pathway_multiome_seurat_cellbender.0.05_preprocess_cr/ra-arco-hvc-nc_hybrid")
ERE = Path("~/ssd/rstudio/multiome/motor-pathway/scenicplus/motor-pathway_scenicplus_v2_hybrid_all/config37/scenicplus_eRegulons.txt").expanduser()
THR = 1.0
ORDER = ["Glut-DACH2-HVCra", "Glut-DACH2-HVCx", "Glut-DACH2-1", "Glut-DACH2-2", "Glut-DACH2-3", "Glut-DACH2-4", "Glut-CACNA1H-RA", "Glut-CACNA1H-1",
         "Glut-CACNA1H-2", "Glut-Im", "Glut-NB", "Glut-NSC", "GABA-LGE-1", "GABA-LGE-2", "GABA-MGE-SST-1", "GABA-MGE-PVALB-1", "GABA-MGE-PVALB-2",
         "GABA-MGE-LAMP5", "GABA-MGE-LHX8", "GABA-CGE", "GABA-Im", "Astro", "Epen", "Oligo", "OPC", "Micro", "Endo"]   # project hybrid order


def bed_names(path):
    b = pd.read_csv(path, sep="\t", header=None, names=["c", "s", "e"])
    return set(b.c + ":" + b.s.astype(str) + "-" + b.e.astype(str))


DAR = {"DAR_RA": bed_names(STORE / "pycisTopic/region_sets_k40_archr/DARs_song-pairs/Glut-CACNA1H-RA_VS_Glut-CACNA1H-1.bed"),
       "DAR_PV1": bed_names(STORE / "pycisTopic/dars_extra/GABA-MGE-PVALB-1_VS_GABA-MGE-LAMP5.bed"),
       "DAR_PV2": bed_names(STORE / "pycisTopic/dars_extra/GABA-MGE-PVALB-2_VS_GABA-MGE-LAMP5.bed")}
W = pd.read_csv(HERE / "scenicplus_mafb_pvalb_vs_mafb_low_mge_config37.csv").set_index("gene")
G = set(W.index[(W.lfc_RA_vs_C1H1 > THR) & (W.lfc_PV1_vs_low > THR)])
e = pd.read_csv(ERE, sep="\t", usecols=["TF", "Gene", "Region", "eRegulon_name"])
e = e[(e.TF == "MAFB") & e.eRegulon_name.str.endswith("+/+") & e.Gene.isin(G)].drop_duplicates(["Gene", "Region"])
meta = e.groupby("Region").Gene.apply(lambda s: ";".join(sorted(s))).rename("genes").to_frame()
print(f"{len(G)} genes up in both (log2FC > {THR}); {len(meta)} linked regions")

o = pickle.load(open(STORE / "pycisTopic/cistopic_obj_glut_k40.pkl", "rb"))
cl = o.cell_data.loc[o.cell_names, "cluster"].astype(str).values
n = pd.Series(cl).value_counts()
clusters = [c for c in ORDER if c in n.index and n[c] >= 50]
fm = o.fragment_matrix.tocsr()
idx = {c: i for i, c in enumerate(clusters)}
onehot = sp.csr_matrix((np.ones(len(cl)), (np.arange(len(cl)), [idx.get(c, len(clusters)) for c in cl])), shape=(len(cl), len(clusters) + 1))[:, : len(clusters)]
tot = np.asarray((fm @ onehot).sum(0)).ravel()                      # fragments per cluster, over all regions, for CPM
rid = pd.Series(np.arange(len(o.region_names)), index=o.region_names)
ix = rid.reindex(meta.index)
assert ix.notna().all(), f"{int(ix.isna().sum())} regions are not in the cisTopic object"
sub = (fm[ix.astype(int).values] @ onehot).toarray().astype(float)    # regions x clusters
L = np.log2(sub / tot[None, :] * 1e6 + 1)
sd = L.std(1)
Z = np.where(sd[:, None] > 0, (L - L.mean(1, keepdims=True)) / np.where(sd[:, None] > 0, sd[:, None], 1), 0.0)
out = meta.copy()
for k, v in DAR.items():
    out[k] = [int(r in v) for r in out.index]
out["n_genes"] = out.genes.str.count(";") + 1
for j, c in enumerate(clusters):
    out[f"z_{c}"] = np.round(Z[:, j], 3)
out.index.name = "region"
out.to_csv(HERE / "scenicplus_mafb_up_both_regions_accessibility_config37.csv")
print("clusters:", len(clusters), "| regions with no fragments in any cluster:", int((sd == 0).sum()))
print(out[["DAR_RA", "DAR_PV1", "DAR_PV2"]].sum().to_dict())
