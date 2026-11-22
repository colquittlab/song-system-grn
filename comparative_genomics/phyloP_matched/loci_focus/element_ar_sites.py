#!/usr/bin/env python
"""AR motif sites in the curated MAFB-domain elements: enrichment by DAR class and phyloP at the sites.

Elements are the curated lobes (element_threshold_classes.tsv; RA, G4, shared = "both"), each scanned over its
scored extent. Motifs and settings as in the earlier AR scans (s_homer.sh, s_tf3.sh): FIMO p < 1e-4 with
  full    HOMER de-novo ARE-full (motif2, 12 bp; found in RA vs Arco DARs genome-wide -- RA enrichment is
          partly built in)
  half    HOMER de-novo ARE-half (motif4, 8 bp; same origin)
  jaspar  JASPAR MA0007.3 (16 bp; independent of these DARs)
Overlapping hits of one motif in a lobe are merged into one site (both strand reports of a palindromic ARE),
represented by its best-scoring hit (exactly the motif width).
Enrichment: sites per kb by class, RA vs G4 / G4 vs shared / RA vs shared by an exact Poisson rate test
(binomial on the pair's site counts given the bases scanned), and the share of lobes with >= 1 site (Fisher).
Constraint: per site, mean whole-tree ("avian") and oscine-lineage phyloP over its bases; compared between
classes (Mann-Whitney) and with the non-motif bases of the site's own lobe (paired Wilcoxon on site minus lobe).
Run on all lobes and on lobes with every 50-bp window aligned in >= 0.8 of oscine and non-oscine genomes.

    python element_ar_sites.py
Writes results/mafb_domain/element_ar_sites.tsv (one row per site).
"""
import os, subprocess, tempfile
import numpy as np, pandas as pd, pyBigWig, pysam
from scipy import stats

R = "/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus/results/mafb_domain"
FA = "/private/groups/colquittlab/assembly/lonStrDom2/GCF_005870125.1_lonStrDom2_genomic_ucsc_only.fna"
PB = "/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/perbase/tracks"
FIMO = "/private/home/bcolquit/miniforge3/envs/crested/bin/fimo"
H = "/private/home/bcolquit/group/motifs/homer/ra_pos_vs_neg/homerResults_manual"
MOT = {"full": f"{H}/motif2.meme", "half": f"{H}/motif4.meme",
       "jaspar": "/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus/ar_jaspar.meme"}
PAIRS = [("RA", "G4"), ("G4", "both"), ("RA", "both")]
NAME = {"RA": "RA", "G4": "G4", "both": "shared"}


def scan(E, fa):
    with tempfile.TemporaryDirectory(dir=os.environ.get("TMPDIR")) as d:
        with open(f"{d}/e.fa", "w") as fh:
            # plain index names: FIMO reads "name:chrom:start-end" headers as coordinates and rewrites the hits
            for i, r in enumerate(E.itertuples()): fh.write(f">e{i}\n{fa.fetch(r.chrom, r.ext_start, r.ext_end).upper()}\n")
        out = []
        for m, f in MOT.items():
            t = subprocess.run([FIMO, "--verbosity", "1", "--thresh", "1e-4", "--text", f, f"{d}/e.fa"], capture_output=True, text=True, check=True).stdout
            h = pd.read_csv(pd.io.common.StringIO(t), sep="\t")
            h = h.rename(columns={"p-value": "p"}).assign(motif=m)
            h["element"] = E.element.values[h.sequence_name.str[1:].astype(int).values] if len(h) else []
            out.append(h[["motif", "element", "start", "stop", "strand", "score", "p"]])
    H_ = pd.concat(out)
    sites = []
    for (m, e), g in H_.groupby(["motif", "element"]):        # merge overlapping hits; keep the best hit per site
        g = g.sort_values("start"); grp, end = [], -1
        for r in g.itertuples():
            if r.start > end: grp.append([])
            grp[-1].append(r); end = max(end, r.stop)
        for c in grp:
            b = max(c, key=lambda r: r.score); sites.append(dict(motif=m, element=e, s=b.start - 1, e=b.stop, strand=b.strand, score=b.score, p=b.p))
    return pd.DataFrame(sites)


def main():
    E = pd.read_csv(f"{R}/element_threshold_classes.tsv", sep="\t")
    fa = pysam.FastaFile(FA)
    S = scan(E, fa).merge(E[["element", "cls", "chrom", "ext_start", "ext_end", "thr_class", "osc_depth_min50", "nonosc_depth_min50"]], on="element")
    S["start"], S["end"] = S.ext_start + S.s, S.ext_start + S.e
    bw = {"avian": pyBigWig.open(f"{PB}/wholetree_perbase_lsd2.bw"), "oscine": pyBigWig.open(f"{PB}/subtree222_perbase_lsd2.bw")}
    v = lambda b, c, a, z: np.nan_to_num(np.array(b.values(c, int(a), int(z)), float))
    for t, b in bw.items():
        S[t] = [v(b, c, a, z).mean() for c, a, z in zip(S.chrom, S.start, S.end)]
        bg = []
        for r in S.itertuples():                              # lobe bases outside every site of this motif
            x = v(b, r.chrom, r.ext_start, r.ext_end); m = np.ones(len(x), bool)
            for a, z in zip(S[(S.element == r.element) & (S.motif == r.motif)].s, S[(S.element == r.element) & (S.motif == r.motif)].e): m[a:z] = False
            bg.append(x[m].mean() if m.any() else np.nan)
        S[f"{t}_lobe_bg"] = bg
    S.drop(columns=["s", "e"]).to_csv(f"{R}/element_ar_sites.tsv", sep="\t", index=False, float_format="%.3f")
    E["bp"] = E.ext_end - E.ext_start
    for lab, keep in (("all lobes", np.ones(len(E), bool)), ("depth >= 0.8", (E.osc_depth_min50 >= .8) & (E.nonosc_depth_min50 >= .8))):
        e = E[keep]; s = S[S.element.isin(e.element)]
        print(f"\n===== {lab}: lobes RA {sum(e.cls == 'RA')}, G4 {sum(e.cls == 'G4')}, shared {sum(e.cls == 'both')}")
        for m in MOT:
            sm = s[s.motif == m]
            n = {k: (sm.cls == k).sum() for k in NAME}; bp = {k: e[e.cls == k].bp.sum() for k in NAME}
            has = {k: e[e.cls == k].element.isin(sm.element).sum() for k in NAME}; nl = {k: (e.cls == k).sum() for k in NAME}
            print(f"  {m:6s} sites/kb: " + ", ".join(f"{NAME[k]} {1000 * n[k] / bp[k]:.2f} ({n[k]})" for k in NAME)
                  + " | lobes with a site: " + ", ".join(f"{NAME[k]} {has[k]}/{nl[k]}" for k in NAME))
            for a, b in PAIRS:
                pr = stats.binomtest(n[a], n[a] + n[b], bp[a] / (bp[a] + bp[b])).pvalue if n[a] + n[b] else np.nan
                pf = stats.fisher_exact([[has[a], nl[a] - has[a]], [has[b], nl[b] - has[b]]])[1]
                rr = (n[a] / bp[a]) / (n[b] / bp[b]) if n[b] else np.inf
                print(f"      {NAME[a]} vs {NAME[b]}: rate ratio {rr:.2f}, Poisson p = {pr:.3g}; lobes-with-site Fisher p = {pf:.3g}")
            for t in ("avian", "oscine"):
                d = sm[t] - sm[f"{t}_lobe_bg"]
                line = f"      {t:6s} phyloP at sites (median; site - own lobe, Wilcoxon p): "
                for k in NAME:
                    x = sm[sm.cls == k]; dd = (x[t] - x[f"{t}_lobe_bg"]).dropna()
                    line += f"{NAME[k]} {x[t].median():+.2f} ({dd.median():+.2f}, p={stats.wilcoxon(dd).pvalue:.2g}, n={len(dd)})  " if len(dd) >= 5 else f"{NAME[k]} n={len(dd)}  "
                print(line)
                cmp_ = []
                for a, b in PAIRS:
                    xa, xb = sm[sm.cls == a][t], sm[sm.cls == b][t]
                    if len(xa) >= 3 and len(xb) >= 3: cmp_.append(f"{NAME[a]} vs {NAME[b]} MW p = {stats.mannwhitneyu(xa, xb).pvalue:.2g}")
                print("        " + "; ".join(cmp_))


if __name__ == "__main__":
    main()
