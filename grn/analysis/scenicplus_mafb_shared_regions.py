#!/usr/bin/env python3
"""Are the regions behind MAFB's targets accessible in both RA and PVALB-2? (companion to scenicplus_mafb_shared_program.R)

Reuse of a GRN can mean the same genes under the same enhancers or the same genes under different ones. The gene-level test cannot
tell these apart; this one looks at the regions. Each region gets an accessibility specificity z-score for RA (Glut-CACNA1H-RA) and
PVALB-2 (GABA-MGE-PVALB-2) across the neuron clusters, from depth-normalized pseudobulk fragment counts (CPM, log2). The MAFB +/+
linked regions are compared with random regions MATCHED on mean accessibility and on how variable it is across neuron clusters, so
"open in every neuron" (promoters, housekeeping) does not count as shared:
    S1  regions with z > 1 in both    S2  mean z_RA * z_PVALB-2    S3  Spearman(z_RA, z_PVALB-2)
Regions are also split by what their linked gene does (specific in both / RA only / PVALB-2 only / neither; from the R script's
table) to ask whether the genes that ARE shared sit under shared regions. Regions are the 517k consensus of configs 1-33, read
from the 30-topic cisTopic object (same regions and cells; only its fragment matrix is used).

    python scenicplus_mafb_shared_regions.py [partner cluster, default GABA-MGE-PVALB-2]
"""
import pickle
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

HERE = Path(__file__).resolve().parent
S = "/hdd/jupyter/brad/scenicplus/motor-pathway_multiome/motor-pathway_multiome_seurat_cellbender.0.05_preprocess_cr/ra-arco-hvc-nc_hybrid/"
OUT = Path("~/ssd/rstudio/multiome/motor-pathway/scenicplus/motor-pathway_scenicplus_v2_hybrid_all/").expanduser()
RA = "Glut-CACNA1H-RA"
PV = sys.argv[1] if len(sys.argv) > 1 else "GABA-MGE-PVALB-2"          # partner cell type; others write files with a suffix
TAG = PV.replace("GABA-MGE-", "").replace("-", "")
SUF = "" if PV == "GABA-MGE-PVALB-2" else f"_vs_{TAG}"
NPERM = 2000
rng = np.random.default_rng(2026)

o = pickle.load(open(S + "pycisTopic/cistopic_obj_glut_k30.pkl", "rb"))
cl = o.cell_data.loc[o.cell_names, "cluster"].astype(str).values
fm = o.fragment_matrix.tocsr()
regions = np.array(o.region_names)
n_cells = pd.Series(cl).value_counts()
clusters = [c for c in sorted(n_cells.index) if n_cells[c] >= 50 and re.match(r"^(Glut|GABA)-", c)]
print("partner:", PV, "| neuron clusters:", len(clusters), "| cells:", int(n_cells[clusters].sum()))
onehot = sp.csr_matrix((np.ones(len(cl)), (np.arange(len(cl)), [clusters.index(c) if c in clusters else len(clusters) for c in cl])),
                       shape=(len(cl), len(clusters) + 1))[:, : len(clusters)]
counts = (fm @ onehot).toarray().astype(np.float64)               # regions x clusters
cpm = counts / counts.sum(0, keepdims=True) * 1e6
L = np.log2(cpm + 1)
mu, sd = L.mean(1), L.std(1)
ok = (sd > 0) & (counts.sum(1) >= 20)
Z = np.zeros_like(L)
Z[ok] = (L[ok] - mu[ok, None]) / sd[ok, None]
zra, zpv = Z[:, clusters.index(RA)], Z[:, clusters.index(PV)]
qm = pd.qcut(mu, 10, labels=False, duplicates="drop")
qs = pd.qcut(sd, 4, labels=False, duplicates="drop")
stratum = np.where(ok, qm * 10 + qs, -1)
strata_idx = {k: np.where(stratum == k)[0] for k in np.unique(stratum) if k >= 0}
rid = pd.Series(np.arange(len(regions)), index=regions)
print("regions usable:", int(ok.sum()), "of", len(regions), "| strata:", len(strata_idx))


def stats(ix):
    a, b = zra[ix], zpv[ix]
    both = (a > 1) & (b > 1)
    rho = pd.Series(a).corr(pd.Series(b), method="spearman")
    return both.sum(), (a * b).mean(), rho


def enrich(ix, nperm=NPERM):
    ix = np.unique(ix[ok[ix]])
    if len(ix) < 20:
        return None
    obs = np.array(stats(ix), dtype=float)
    need = pd.Series(stratum[ix]).value_counts()
    null = np.empty((nperm, 3))
    for i in range(nperm):
        pick = np.concatenate([rng.choice(strata_idx[k], size=n, replace=len(strata_idx[k]) < n) for k, n in need.items()])
        null[i] = stats(pick)
    return dict(n_regions=len(ix), S1=int(obs[0]), S1_null=null[:, 0].mean(), S1_ratio=obs[0] / null[:, 0].mean(),
                S1_p=(1 + (null[:, 0] >= obs[0]).sum()) / (nperm + 1),
                S2=obs[1], S2_null=null[:, 1].mean(), S2_z=(obs[1] - null[:, 1].mean()) / null[:, 1].std(),
                S3=obs[2], S3_null=null[:, 2].mean(), S3_z=(obs[2] - null[:, 2].mean()) / null[:, 2].std())


def load(cfg):
    p = OUT / cfg / "scenicplus_eRegulons.txt"
    if not p.exists():
        return None
    d = pd.read_csv(p, sep="\t", usecols=["TF", "Gene", "Region", "eRegulon_name"])
    return d[d.eRegulon_name.str.endswith("+/+")]


rows = []
for cfg in [f"config{n}" for n in list(range(1, 12)) + list(range(16, 29))]:
    d = load(cfg)
    if d is None:
        continue
    ix = rid.reindex(d[d.TF == "MAFB"].Region.unique()).dropna().astype(int).values
    r = enrich(ix)
    if r:
        rows.append(dict(set=cfg, **r))
print(pd.DataFrame(rows).round(3).to_string(index=False))

# gene classes (from the R script): do shared genes sit under shared regions?
g = pd.read_csv(OUT / f"mafb_gene_specificity_RA_{'PVALB2' if SUF == '' else TAG}.csv").set_index("gene")
g = g.rename(columns={"z_PVALB2": "z_partner", "expr_PVALB2": "expr_partner"})
d1 = load("config1")
m = d1[d1.TF == "MAFB"].drop_duplicates(["Gene", "Region"]).copy()
m = m[m.Gene.isin(g.index)]
za, zb = g.loc[m.Gene, "z_RA"].values, g.loc[m.Gene, "z_partner"].values
m["cls"] = np.select([(za > 1) & (zb > 1), (za > 1) & (zb <= 0), (zb > 1) & (za <= 0)], ["shared", "RA only", "partner only"], "neither")
print("\nconfig1 MAFB +/+ links by what the linked gene does (z across neuron clusters):")
cls_rows = []
for c, sub in m.groupby("cls"):
    ix = rid.reindex(sub.Region.unique()).dropna().astype(int).values
    r = enrich(ix)
    if r:
        cls_rows.append(dict(gene_class=c, n_genes=sub.Gene.nunique(), **r))
        ixu = np.unique(ix[ok[ix]])
        cls_rows[-1]["frac_open_in_RA"] = (zra[ixu] > 1).mean()
        cls_rows[-1]["frac_open_in_partner"] = (zpv[ixu] > 1).mean()
cls = pd.DataFrame(cls_rows)
print(cls.round(3).to_string(index=False))

# the regions open in both: which genes do they link to, and are those genes themselves specific in either cell type?
mm = d1[d1.TF == "MAFB"].drop_duplicates(["Gene", "Region"]).copy()
mm["rix"] = rid.reindex(mm.Region).values
mm = mm.dropna(subset=["rix"])
mm["rix"] = mm.rix.astype(int)
mm = mm[ok[mm.rix.values]]
mm["z_region_RA"], mm["z_region_PVALB2"] = zra[mm.rix.values], zpv[mm.rix.values]
mm["open_both"] = (mm.z_region_RA > 1) & (mm.z_region_PVALB2 > 1)
mm["z_gene_RA"] = g.reindex(mm.Gene).z_RA.values
mm["z_gene_PVALB2"] = g.reindex(mm.Gene).z_partner.values
both = mm[mm.open_both]
print(f"\nconfig1: {both.Region.nunique()} MAFB regions open in both, linked to {both.Gene.nunique()} genes")
print("  linked genes' own specificity (neuron z): RA median %.2f, partner median %.2f; share with z>1 in RA %.2f, in partner %.2f, both %.2f" % (
    both.z_gene_RA.median(), both.z_gene_PVALB2.median(), (both.z_gene_RA > 1).mean(), (both.z_gene_PVALB2 > 1).mean(),
    ((both.z_gene_RA > 1) & (both.z_gene_PVALB2 > 1)).mean()))
gb = both.groupby("Gene").agg(n_open_both_regions=("Region", "nunique")).join(g[["z_RA", "z_partner", "expr_RA", "expr_partner"]])
print(gb.sort_values("n_open_both_regions", ascending=False).head(15).round(2).to_string())
gb.sort_values("n_open_both_regions", ascending=False).to_csv(HERE / f"scenicplus_mafb_genes_under_shared_regions{SUF}.csv")

# reference: every other TF's +/+ regions in config1
ref = []
for tf, sub in d1.groupby("TF"):
    ix = rid.reindex(sub.Region.unique()).dropna().astype(int).values
    if len(np.unique(ix[ok[ix]])) >= 100:
        r = enrich(ix, nperm=300)
        if r:
            ref.append(dict(TF=tf, **r))
ref = pd.DataFrame(ref)
mf = ref[ref.TF == "MAFB"].iloc[0]
print(f"\nconfig1: MAFB vs {len(ref) - 1} other TFs with >= 100 +/+ regions")
for k in ("S1_ratio", "S2_z", "S3_z"):
    other = ref.loc[ref.TF != "MAFB", k]
    print(f"  {k}: MAFB {mf[k]:.2f}, rank {int((ref[k] > mf[k]).sum()) + 1} of {len(ref)}; others median {other.median():.2f} (IQR {other.quantile(.25):.2f} to {other.quantile(.75):.2f})")

pd.DataFrame(rows).to_csv(HERE / f"scenicplus_mafb_shared_regions{SUF}.csv", index=False)
cls.to_csv(HERE / f"scenicplus_mafb_shared_regions_by_gene_class{SUF}.csv", index=False)
ref.to_csv(HERE / f"scenicplus_shared_regions_all_tfs_config1{SUF}.csv", index=False)
