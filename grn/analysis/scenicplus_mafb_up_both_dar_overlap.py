#!/usr/bin/env python3
"""Do the genes that are up in BOTH RA (vs C1H-1) and PVALB-1 (vs MAFB-low MGE) use the same differentially accessible regions?

Gene set G: MAFB +/+ targets of the main config (config37) with log2FC > THR in both contrasts (directional rule, no padj; THR = 1 as in the
scatter labels, with 0.5 and 0 as sensitivity). Regions: the region-gene links of MAFB's +/+ eRegulons (config37), i.e. the regions behind each
target. DARs: ArchR bias-matched Wilcoxon on the consensus peak matrix, FDR <= 0.05 and Log2FC >= 0.585, regions MORE accessible in the
foreground: RA vs C1H-1 (make_song_pair_dars_archr.R) and PVALB-1 / PVALB-2 vs LAMP5 (make_interneuron_dars_archr.R).

For G, per contrast pair:
  1. regions: how many of G's linked regions are DARs in RA, in PVALB-1, and in both, against the count independence would give;
  2. genes: each gene is "same region" (one region is a DAR in both), "different regions" (an RA-DAR region and a PVALB-1-DAR region, but none
     shared), "RA only", "PVALB only" or "neither";
  3. context: the same numbers for other MAFB targets (RA-up-only, PVALB-up-only, all targets), for 5,000 random gene sets of the same size
     drawn from the testable MAFB targets, and the genome-wide overlap of the two DAR sets.

    python scenicplus_mafb_up_both_dar_overlap.py
"""
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
STORE = Path("/hdd/jupyter/brad/scenicplus/motor-pathway_multiome/motor-pathway_multiome_seurat_cellbender.0.05_preprocess_cr/ra-arco-hvc-nc_hybrid")
OUT = Path("~/ssd/rstudio/multiome/motor-pathway/scenicplus/motor-pathway_scenicplus_v2_hybrid_all/config37/scenicplus_eRegulons.txt").expanduser()
MAIN_REGIONS = 499348   # consensus regions in the cisTopic object of configs 34-39
rng = np.random.default_rng(2026)


def bed_names(path):
    b = pd.read_csv(path, sep="\t", header=None, names=["c", "s", "e"])
    return set(b.c + ":" + b.s.astype(str) + "-" + b.e.astype(str))


DAR = {"RA": bed_names(STORE / "pycisTopic/region_sets_k40_archr/DARs_song-pairs/Glut-CACNA1H-RA_VS_Glut-CACNA1H-1.bed"),
       "PV1": bed_names(STORE / "pycisTopic/dars_extra/GABA-MGE-PVALB-1_VS_GABA-MGE-LAMP5.bed"),
       "PV2": bed_names(STORE / "pycisTopic/dars_extra/GABA-MGE-PVALB-2_VS_GABA-MGE-LAMP5.bed")}
print("DARs:", {k: len(v) for k, v in DAR.items()})

W = pd.read_csv(HERE / "scenicplus_mafb_pvalb_vs_mafb_low_mge_config37.csv").set_index("gene")
e = pd.read_csv(OUT, sep="\t", usecols=["TF", "Gene", "Region", "eRegulon_name"])
e = e[(e.TF == "MAFB") & e.eRegulon_name.str.endswith("+/+")].drop_duplicates(["Gene", "Region"])
e = e[e.Gene.isin(W.index)]
regs = e.groupby("Gene").Region.apply(set).to_dict()
U = set(e.Region)
print("testable MAFB +/+ targets with linked regions:", len(regs), "| linked regions:", len(U))


def region_stats(genes, a, b):
    R = set().union(*[regs[g] for g in genes if g in regs]) if genes else set()
    n = len(R)
    na, nb = len(R & DAR[a]), len(R & DAR[b])
    both = len(R & DAR[a] & DAR[b])
    exp = na * nb / n if n else np.nan
    uni = len((R & DAR[a]) | (R & DAR[b]))
    return dict(n_genes=len(genes), n_regions=n, n_DAR_RA=na, n_DAR_partner=nb, n_DAR_both=both, expected_if_independent=round(exp, 1),
                obs_over_expected=round(both / exp, 2) if exp else np.nan, jaccard=round(both / uni, 3) if uni else np.nan,
                frac_regions_both=round(both / n, 3) if n else np.nan)


def gene_class(g, a, b):
    r = regs.get(g, set())
    ra, pb = r & DAR[a], r & DAR[b]
    if ra & pb:
        return "same region (a region is a DAR in both)"
    if ra and pb:
        return "different regions (RA-DAR and partner-DAR regions, none shared)"
    if ra:
        return "RA DARs only"
    if pb:
        return "partner DARs only"
    return "no DARs among linked regions"


summary, per_gene = [], []
pairs = (("RA", "PV1", "PV1_vs_low"), ("RA", "PV2", "PV2_vs_low"))
for a, b, col in pairs:
    for thr in (1.0, 0.5, 0.0):
        G = [g for g in regs if W.loc[g, "lfc_RA_vs_C1H1"] > thr and W.loc[g, f"lfc_{col}"] > thr]
        row = region_stats(G, a, b)
        summary.append(dict(partner=b, rule=f"both log2FC > {thr}", set="up in both", **row))
        if thr == 1.0:
            cls = pd.Series({g: gene_class(g, a, b) for g in G})
            for g in G:
                r = regs[g]
                per_gene.append(dict(partner=b, gene=g, lfc_RA=round(W.loc[g, "lfc_RA_vs_C1H1"], 2), lfc_partner=round(W.loc[g, f"lfc_{col}"], 2), n_regions=len(r),
                                     n_RA_DARs=len(r & DAR[a]), n_partner_DARs=len(r & DAR[b]), n_both_DARs=len(r & DAR[a] & DAR[b]), gene_class=cls[g]))
            print(f"\nRA + {b}, genes up in both (log2FC > 1): {len(G)}")
            print(cls.value_counts().to_string())
    # context sets at THR = 1
    ra_only = [g for g in regs if W.loc[g, "lfc_RA_vs_C1H1"] > 1 and W.loc[g, f"lfc_{col}"] <= 0]
    pb_only = [g for g in regs if W.loc[g, f"lfc_{col}"] > 1 and W.loc[g, "lfc_RA_vs_C1H1"] <= 0]
    for lab, S in (("RA up only (log2FC > 1; partner <= 0)", ra_only), ("partner up only (log2FC > 1; RA <= 0)", pb_only), ("all testable MAFB targets", list(regs))):
        summary.append(dict(partner=b, rule="", set=lab, **region_stats(S, a, b)))
    # random gene sets of the size of G
    G1 = [g for g in regs if W.loc[g, "lfc_RA_vs_C1H1"] > 1 and W.loc[g, f"lfc_{col}"] > 1]
    allg = np.array(list(regs))
    fr, nb = [], []
    for _ in range(5000):
        S = list(rng.choice(allg, size=len(G1), replace=False))
        st = region_stats(S, a, b)
        fr.append(st["frac_regions_both"])
        nb.append(st["obs_over_expected"])
    obs = region_stats(G1, a, b)
    summary.append(dict(partner=b, rule="", set=f"random {len(G1)} MAFB targets (mean of 5,000)", n_genes=len(G1), n_regions=np.nan, n_DAR_RA=np.nan, n_DAR_partner=np.nan,
                        n_DAR_both=np.nan, expected_if_independent=np.nan, obs_over_expected=round(np.nanmean(nb), 2), jaccard=np.nan, frac_regions_both=round(np.nanmean(fr), 3)))
    p = (1 + sum(f >= obs["frac_regions_both"] for f in fr)) / 5001
    print(f"  fraction of linked regions that are DARs in both: G {obs['frac_regions_both']} vs random sets {np.mean(fr):.3f} (95% {np.percentile(fr, 95):.3f}); one-sided p = {p:.4f}")
    summary[-1]["p_frac_both_vs_random"] = round(p, 4)
    # genome-wide overlap of the two DAR sets
    both = len(DAR[a] & DAR[b])
    exp = len(DAR[a]) * len(DAR[b]) / MAIN_REGIONS
    summary.append(dict(partner=b, rule="", set="all consensus regions (genome-wide)", n_genes=np.nan, n_regions=MAIN_REGIONS, n_DAR_RA=len(DAR[a]), n_DAR_partner=len(DAR[b]),
                        n_DAR_both=both, expected_if_independent=round(exp, 1), obs_over_expected=round(both / exp, 2), jaccard=round(both / len(DAR[a] | DAR[b]), 3),
                        frac_regions_both=round(both / MAIN_REGIONS, 4)))

S = pd.DataFrame(summary)
pd.set_option("display.width", 250, "display.max_columns", 20)
print("\n" + S.to_string(index=False))
S.to_csv(HERE / "scenicplus_mafb_up_both_dar_overlap_config37.csv", index=False)
pd.DataFrame(per_gene).to_csv(HERE / "scenicplus_mafb_up_both_dar_overlap_genes_config37.csv", index=False)
g = pd.DataFrame(per_gene)
print("\nper gene (RA + PVALB-1):")
print(g[g.partner == "PV1"].sort_values("n_both_DARs", ascending=False).to_string(index=False))
