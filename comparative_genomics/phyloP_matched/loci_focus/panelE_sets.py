#!/usr/bin/env python
"""Rebuild the Fig. 6E region set with the hybrid-label inputs, and the filter ablation.

Fig. 6E (comparative_genomics/deeptools/deeptools_mafb.qmd) draws phyloP at
high-accessibility-difference bases in MAFB-regulating regions. Its regions are built by
filter_bed():
  1. SnapATAC2 cluster peaks: Glut-RA.bed (RA set) and GABA-4-1.bed (G4 set)
  2. overlapping a MAFB peak-to-gene link (ArchR, r >= 0.4; RA and GABA links pooled;
     linked peak resized to 500 bp about its centre)
  3. overlapping a DAR of the same contrast with FDR < 0.1 and Log2FC > 1
     (RA: Glut-RA vs Glut-Arco-1, `ra`; G4: GABA-MGE-PVALB-2 vs GABA-LGE-1, `pvalb2_lge1`)
  4. peak midpoint 20-200 kb from the MAFB TSS (chr20:13,369,534, minus strand)
  5. no overlap with an exon
Inputs here are the hybrid-label versions: peaks2geneslinks_hybrid, peaks_glut_hybrid and
peaks_gaba_hybrid DAR objects (exported to TSV by export_dar.R), and the stringtie GTF.

Ablation: the loci_focus domain DARs (RA_over_Arco / GABA4_over_GABA1, FDR < 0.1, no fold-change
cutoff, elements up in both dropped) with Fig. 6E's filters added one at a time -- p2g link,
Log2FC > 1, TSS window, exon exclusion -- and all together.

Writes results/panelE/<set>.bed (chrom, start, end, name, class) and filter_counts.tsv.
"""
import collections, os, re, sys
import numpy as np, pandas as pd

F = "/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus"
sys.path.insert(0, F)
from span import elem_coords, in_domain, result

GRN = "/private/home/bcolquit/group/song-system-grn/multiome/archr"
T = "/private/warm-archive/colquittlab/assembly/tracks_lark"
GTF = ("/private/home/bcolquit/group/assembly/lonStrDom2/GCF_005870125.1_lonStrDom2_genomic_stringtie_"
       "FPKMthresh0.5_FPKMmaxfrac0.1_minhits5_ucsc_sort_gene-name.gtf")
GENE, P2G_R, LFC, FDR = "MAFB", 0.4, 1.0, 0.1
TSS_EXCL, TSS_WIN = 20_000, 200_000
CLS = {"RA": ("RA_over_Arco", "Glut-RA", "ra"), "G4": ("GABA4_over_GABA1", "GABA-4-1", "pvalb2_lge1")}
OUT = os.path.dirname(result("x", "panelE"))


def overlaps(q, iv):
    """q: DataFrame chrom/start/end (0-based half-open); iv: DataFrame same. -> bool array"""
    out = np.zeros(len(q), bool)
    for c, g in q.groupby("chrom"):
        h = iv[iv.chrom == c]
        if h.empty: continue
        hs = h.sort_values("start"); S, E = hs.start.values, hs.end.values
        emax = np.maximum.accumulate(E)
        for i, (s, e) in zip(g.index, zip(g.start.values, g.end.values)):
            j = np.searchsorted(S, e)                     # intervals starting before e
            if j and emax[j - 1] > s: out[q.index.get_loc(i)] = True
    return out


def main():
    # MAFB TSS and exons
    tss = None; ex = []
    for l in open(GTF):
        x = l.split("\t")
        if len(x) < 9: continue
        if x[2] == "exon": ex.append((x[0], int(x[3]) - 1, int(x[4])))
        elif x[2] == "gene" and f'gene_name "{GENE}"' in x[8]:
            tss = (x[0], int(x[4]) if x[6] == "-" else int(x[3]) - 1)
    exons = pd.DataFrame(ex, columns=["chrom", "start", "end"])
    print(f"  {GENE} TSS {tss[0]}:{tss[1]:,}; exons {len(exons):,}")

    # MAFB p2g links (pooled), linked peak from peak_id, resized to 500 bp
    links = []
    for k in ("ra", "gaba"):
        d = pd.read_csv(f"{GRN}/peaks2geneslinks_hybrid/peak2geneLinks_{k}_corCutOff-{P2G_R}.bed", sep="\t", header=None)
        d = d[d[11] == GENE]
        for pid in d[16]:
            m = re.match(r"(.+)_(\d+)-(\d+)$", pid); c, s, e = m.group(1), int(m.group(2)), int(m.group(3))
            mid = (s + e) // 2; links.append((c, mid - 250, mid + 250))
    links = pd.DataFrame(sorted(set(links)), columns=["chrom", "start", "end"])
    print(f"  {GENE} p2g-linked peaks (r >= {P2G_R}, RA+GABA pooled): {len(links)}")

    # DARs from the hybrid ArchR objects (rowData start is 1-based)
    dar = {}
    for k, (_, _, con) in CLS.items():
        d = pd.read_csv(f"{T}/panelE_inputs/dar_{con}.tsv", sep="\t")
        d["start"] = d.start - 1
        dar[k] = d

    def tss_ok(df):
        mid = (df.start + df.end) // 2; dist = (mid - tss[1]).abs()
        return ((df.chrom == tss[0]) & (dist > TSS_EXCL) & (dist <= TSS_WIN)).values

    counts = []; sets = {}
    # ---- exact Fig. 6E set
    for k, (cl, snap, con) in CLS.items():
        p = pd.read_csv(f"{T}/snapatac_peaks/{snap}.bed", sep="\t", header=None, usecols=[0, 1, 2], names=["chrom", "start", "end"])
        p = p.reset_index(drop=True); n0 = len(p)
        p = p[overlaps(p, links)]; n1 = len(p)
        hi = dar[k][(dar[k].FDR < FDR) & (dar[k].Log2FC > LFC)]
        p = p[overlaps(p, hi)]; n2 = len(p)
        p = p[tss_ok(p)]; n3 = len(p)
        p = p[~overlaps(p, exons)]; n4 = len(p)
        counts.append(dict(set="panelE", cls=k, start=n0, p2g=n1, dar_fdr_lfc=n2, tss_20_200kb=n3, no_exon=n4))
        sets.setdefault("panelE", []).append(p.assign(name=[f"{snap}:{c}:{s}-{e}" for c, s, e in zip(p.chrom, p.start, p.end)], cls=cl))
    # ---- ablation on the loci_focus domain DARs
    EL = elem_coords(); cl_of = {}
    for l in open(f"{F}/allelem_cls.bed"):
        x = l.rstrip("\n").split("\t"); cl_of[x[3]] = x[4]
    dom = pd.DataFrame([(c, s, e, n, cl_of[n]) for n, (c, s, e) in EL.items() if in_domain((c, s, e))],
                       columns=["chrom", "start", "end", "name", "cls"]).reset_index(drop=True)
    f_p2g = overlaps(dom, links)
    f_lfc = np.zeros(len(dom), bool)
    for k, (cl, _, _) in CLS.items():
        hi = dar[k][(dar[k].FDR < FDR) & (dar[k].Log2FC > LFC)]
        m = (dom.cls == cl).values; f_lfc[m] = overlaps(dom[m], hi)
    f_tss = tss_ok(dom); f_ex = ~overlaps(dom, exons)
    ab = {"domain": np.ones(len(dom), bool), "domain+p2g": f_p2g, "domain+lfc1": f_lfc,
          "domain+tss": f_tss, "domain+noexon": f_ex, "domain+all": f_p2g & f_lfc & f_tss & f_ex}
    for nm, m in ab.items():
        sets[nm] = [dom[m]]
        for k, (cl, _, _) in CLS.items():
            counts.append(dict(set=nm, cls=k, n=int(((dom.cls == cl).values & m).sum())))
    for nm, parts in sets.items():
        d = pd.concat(parts)[["chrom", "start", "end", "name", "cls"]].sort_values(["chrom", "start"])
        d.to_csv(f"{OUT}/{nm}.bed", sep="\t", header=False, index=False)
    C = pd.DataFrame(counts); C.to_csv(f"{OUT}/filter_counts.tsv", sep="\t", index=False)
    print(C.to_string(index=False))
    print(f"  wrote {OUT}/<set>.bed and filter_counts.tsv")


if __name__ == "__main__":
    main()
