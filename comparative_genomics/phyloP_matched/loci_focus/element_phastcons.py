#!/usr/bin/env python
"""Clade phastCons (oscine vs non-oscine, chr20; phastcons_clade/v2) at the MAFB-domain elements.

phastCons was run separately on the oscine (146 spp.) and non-oscine (217 spp.) sub-alignments
(TC 0.25, EL 20, per-clade rho). Elements from these alignments are short (median a few bp), so
oscine-specific elements -- oscine elements with no overlap with any non-oscine element -- are
merged across gaps of <= 10 bp and kept if >= 20 bp. Per lobe, over its full extent
(element_rule_classes.tsv), weighted by accessibility as for the rule scores:
  pc_osc / pc_nonosc   mean posterior in each clade
  pc_diff              mean of oscine - non-oscine posterior
  pc_oscspec           % of bases in merged oscine-specific elements
Each is compared with the visual oscine rating (Spearman; best lobe per rated region, as in
element_rule_classes.py) and summarised by rule class.

    python element_phastcons.py
Writes results/mafb_domain/element_phastcons.tsv and, in phastcons_clade/v2/results/hub/,
phastCons_oscspec_merged_chr20.bb (merged oscine-specific elements) and
phastCons_osc_minus_nonosc_chr20.bw (posterior difference).
"""
import subprocess
import numpy as np, pandas as pd, pyBigWig
from scipy import stats

F = "/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus"
R = f"{F}/results/mafb_domain"
V = "/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/phastcons_clade/v2/results"
FR = "/private/warm-archive/colquittlab/assembly/tracks_lark/fragment"
UCSC = "/private/groups/colquittlab/opt/ucsc_utilities"
SIZES = "/private/warm-archive/colquittlab/assembly/cactus/avian363/vgp577/fa/lonStrDom2.sizes"
UP = {"RA": ["Glut-RA"], "G4": ["GABA-4-1"], "both": ["Glut-RA", "GABA-4-1"]}
GAP, MINLEN = 10, 20
RANK = {"neutral": 0, "weak": 1, "moderate": 2, "strong": 3}


def merged_oscspec():
    e = pd.read_csv(f"{V}/phastCons_oscine_specific_chr20.bed", sep="\t", header=None, usecols=[0, 1, 2], names=["c", "s", "e"]).sort_values("s")
    out, cs, ce = [], None, None
    for s, t in zip(e.s, e.e):
        if cs is not None and s - ce <= GAP: ce = max(ce, t)
        else:
            if cs is not None: out.append((cs, ce))
            cs, ce = s, t
    out.append((cs, ce))
    M = pd.DataFrame(out, columns=["s", "e"]); M = M[M.e - M.s >= MINLEN]
    print(f"  oscine-specific elements: {len(e)} -> {len(M)} after merging <= {GAP} bp gaps and keeping >= {MINLEN} bp "
          f"({(M.e - M.s).sum() / 1e3:.0f} kb, median {int((M.e - M.s).median())} bp)")
    return M


def tracks(M):
    import os
    H = f"{V}/hub"; os.makedirs(H, exist_ok=True)
    bed = f"{H}/phastCons_oscspec_merged_chr20.bed"
    M.assign(c="chr20", n=[f"oscspec_{i + 1}" for i in range(len(M))])[["c", "s", "e", "n"]].to_csv(bed, sep="\t", header=False, index=False)
    subprocess.run([f"{UCSC}/bedToBigBed", "-type=bed4", bed, SIZES, bed[:-4] + ".bb"], check=True)
    o, n = pyBigWig.open(f"{V}/phastCons_osc_chr20.bw"), pyBigWig.open(f"{V}/phastCons_nonosc_chr20.bw")
    L = o.chroms("chr20"); d = np.nan_to_num(np.array(o.values("chr20", 0, L))) - np.nan_to_num(np.array(n.values("chr20", 0, L)))
    w = pyBigWig.open(f"{H}/phastCons_osc_minus_nonosc_chr20.bw", "w"); w.addHeader([("chr20", L)])
    w.addEntries("chr20", 0, values=d.astype(float).tolist(), span=1, step=1); w.close()
    print(f"  wrote {H}/phastCons_oscspec_merged_chr20.bb and phastCons_osc_minus_nonosc_chr20.bw")


def main():
    M = merged_oscspec(); tracks(M)
    S = pd.read_csv(f"{R}/element_rule_classes.tsv", sep="\t")
    o, n = pyBigWig.open(f"{V}/phastCons_osc_chr20.bw"), pyBigWig.open(f"{V}/phastCons_nonosc_chr20.bw")
    fb = {x: pyBigWig.open(f"{FR}/{x}-fragment-bin1-smoothNone.bw") for x in ("Glut-RA", "GABA-4-1")}
    v = lambda bw, c, s, e: np.nan_to_num(np.array(bw.values(c, s, e), float))
    rows = []
    for r in S.itertuples(index=False):
        s, e = int(r.ext_start), int(r.ext_end)
        w = sum(v(fb[x], r.chrom, s, e) / max(v(fb[x], r.chrom, s, e).max(), 1e-9) for x in UP[r.cls])
        w = np.clip(w, 0, None); w = w if w.sum() > 0 else np.ones(e - s)
        po, pn = v(o, r.chrom, s, e), v(n, r.chrom, s, e)
        inside = np.zeros(e - s, bool)
        for a, b in zip(M.s, M.e):
            if a < e and b > s: inside[max(a, s) - s:min(b, e) - s] = True
        avg = lambda x: (w * x).sum() / w.sum()
        rows.append(dict(element=r.element, pc_osc=avg(po), pc_nonosc=avg(pn), pc_diff=avg(po - pn), pc_oscspec=100 * avg(inside)))
    P = S.merge(pd.DataFrame(rows), on="element")
    P[["element", "cls", "rule", "chrom", "ext_start", "ext_end", "osc_ext_pct1", "whole_ext_block", "pc_osc", "pc_nonosc", "pc_diff", "pc_oscspec"]].to_csv(
        f"{R}/element_phastcons.tsv", sep="\t", index=False, float_format="%.3f")
    print("\n  by rule class (mean):")
    print(P.groupby("rule")[["pc_osc", "pc_nonosc", "pc_diff", "pc_oscspec", "osc_ext_pct1"]].mean().round(3)
          .reindex(["deep+oscine", "deep", "oscine", "weak"]).to_string())
    m = pd.read_csv(f"{F}/results/manual_phylop/manual_phylop.tsv", sep="\t", dtype=str).fillna("")
    m = m[m.oscine.isin(RANK)]; Q = []
    for x in m.itertuples(index=False):
        a, b = x.region.strip().split(":")[-1].replace(",", "").split("-"); a, b = int(a) - 1, int(b)
        h = P[(P.ext_start < b) & (P.ext_end > a) & P.cls.isin([x[3], "both"])]
        if len(h): Q.append(dict(rating=RANK[x.oscine], **{k: h[k].max() for k in ("osc_ext_pct1", "pc_diff", "pc_oscspec", "pc_osc")}))
    Q = pd.DataFrame(Q)
    print(f"\n  Spearman with the visual oscine rating ({len(Q)} regions): " +
          "   ".join(f"{k} {stats.spearmanr(Q.rating, Q[k])[0]:+.2f}" for k in ("osc_ext_pct1", "pc_diff", "pc_oscspec", "pc_osc")))
    print(f"  rated moderate/strong vs neutral/weak, pc_oscspec median: "
          f"{Q[Q.rating >= 2].pc_oscspec.median():.1f} vs {Q[Q.rating < 2].pc_oscspec.median():.1f}")


if __name__ == "__main__":
    main()
