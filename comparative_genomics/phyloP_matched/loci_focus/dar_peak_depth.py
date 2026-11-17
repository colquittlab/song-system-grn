#!/usr/bin/env python
"""Are weakly constrained DAR peaks songbird-specific sequence? Per-clade alignment depth vs phyloP.

For each genome-wide DAR peak (results/genome_wide/dar_peak_repeats.tsv.gz; element_repeats.py), alignment
depth is taken from the 363-way per-element depth scan (depth_scan/out; element_depth.py): the mean
fraction of oscine (Anc222 clade, 145) and non-oscine (217) genomes with a non-gap base, over each
500-bp ArchR peak; a SnapATAC peak takes the overlap-weighted mean over the ArchR peaks it overlaps.
  oscine-specific sequence   oscine depth >= 0.5 and non-oscine depth < 0.1
  shared                     non-oscine depth >= 0.5
The oscine-lineage phyloP track is phyloP --subtree birdAnc222 (LRT, CONACC): oscine rate given the
rest of the tree, so sequence absent outside oscines has nothing to be contrasted with. The clade
phastCons posterior on the oscine-only sub-alignment (phastcons_clade/v2; chr20 only) measures
conservation within oscines without that contrast.

    python dar_peak_depth.py
Writes results/genome_wide/dar_peak_depth.tsv.gz.
"""
import glob
from collections import defaultdict
import numpy as np, pandas as pd, pyBigWig
from scipy import stats

R = "/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus/results/genome_wide"
DS = "/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/depth_scan/out"
PC = "/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/phastcons_clade/v2/results"


def main():
    Q = pd.read_csv(f"{R}/dar_peak_repeats.tsv.gz", sep="\t")
    cols = ["eid", "chrom", "start", "end", "bp", "in_depth", "out_depth"]
    A = pd.concat(pd.read_csv(f, sep="\t", header=None, usecols=range(7), names=cols, comment="e") for f in glob.glob(f"{DS}/*.tsv"))
    A = A.drop_duplicates("eid")
    idx = {c: g.sort_values("start")[["start", "end", "in_depth", "out_depth"]].values for c, g in A.groupby("chrom")}
    ind, outd = [], []
    for c, s, e in zip(Q.chrom, Q.start, Q.end):
        a = idx.get(c)
        if a is None: ind.append(np.nan); outd.append(np.nan); continue
        i = max(0, np.searchsorted(a[:, 0], s - 600)); j = np.searchsorted(a[:, 0], e)
        w = np.clip(np.minimum(a[i:j, 1], e) - np.maximum(a[i:j, 0], s), 0, None)
        if w.sum() == 0: ind.append(np.nan); outd.append(np.nan); continue
        ind.append((w * a[i:j, 2]).sum() / w.sum()); outd.append((w * a[i:j, 3]).sum() / w.sum())
    Q["osc_depth"], Q["nonosc_depth"] = ind, outd
    Q = Q.dropna(subset=["osc_depth"])
    Q["seq"] = np.select([(Q.osc_depth >= 0.5) & (Q.nonosc_depth < 0.1), Q.nonosc_depth >= 0.5], ["oscine-specific", "shared"], "partial")
    o, n = pyBigWig.open(f"{PC}/phastCons_osc_chr20.bw"), pyBigWig.open(f"{PC}/phastCons_nonosc_chr20.bw")
    c20 = Q.chrom == "chr20"
    Q.loc[c20, "pc_osc"] = [o.stats(c, s, e, type="mean")[0] for c, s, e in zip(Q[c20].chrom, Q[c20].start, Q[c20].end)]
    Q.loc[c20, "pc_nonosc"] = [n.stats(c, s, e, type="mean")[0] for c, s, e in zip(Q[c20].chrom, Q[c20].start, Q[c20].end)]
    Q.to_csv(f"{R}/dar_peak_depth.tsv.gz", sep="\t", index=False, float_format="%.3f")

    for cl in ("RA", "G4"):
        q = Q[Q.cls == cl]; w3 = q.whole.quantile([1 / 3, 2 / 3]).values; o3 = q.oscine.quantile([1 / 3, 2 / 3]).values
        q = q.assign(grp=np.select([(q.whole <= w3[0]) & (q.oscine <= o3[0]), (q.whole >= w3[1]) | (q.oscine >= o3[1])], ["weak", "constrained"], "mid"))
        print(f"\n== {cl} DAR peaks (n {len(q)})")
        t = pd.crosstab(q.grp, q.seq, normalize="index").mul(100).round(1)[["shared", "partial", "oscine-specific"]]
        print(t.assign(n=q.grp.value_counts()).reindex(["weak", "mid", "constrained"]).to_string())
        x = q[q.grp != "mid"]; a = (x.grp == "weak"); b = (x.seq == "oscine-specific")
        orr, p = stats.fisher_exact([[(a & b).sum(), (a & ~b).sum()], [(~a & b).sum(), (~a & ~b).sum()]])
        print(f"   oscine-specific sequence, weak vs constrained: OR {orr:.2f}, p {p:.2g}")
        te = x.te_frac > 0
        print(f"   of weak peaks touching a TE, oscine-specific: {100 * ((x.seq == 'oscine-specific') & te & a).sum() / max(1, (te & a).sum()):.1f}%  "
              f"(weak without TE: {100 * ((x.seq == 'oscine-specific') & ~te & a).sum() / max(1, (~te & a).sum()):.1f}%)")
        print("   median phyloP by sequence class:")
        print(q.groupby("seq")[["whole", "oscine", "nonosc_depth", "osc_depth"]].median().round(2).reindex(["shared", "partial", "oscine-specific"]).to_string())
    s = Q[Q.chrom == "chr20"].dropna(subset=["pc_osc"])
    print("\n== chr20: oscine-only phastCons posterior vs oscine subtree phyloP, by sequence class")
    print(s.groupby("seq").agg(n=("pc_osc", "size"), pc_osc=("pc_osc", "median"), pc_nonosc=("pc_nonosc", "median"),
                              osc_phyloP=("oscine", "median"), whole_phyloP=("whole", "median")).round(2).to_string())
    os_ = s[s.seq == "oscine-specific"]
    if len(os_) > 5:
        print(f"   oscine-specific peaks: Spearman(pc_osc, oscine phyloP) {stats.spearmanr(os_.pc_osc, os_.oscine)[0]:+.2f}; "
              f"{(os_.pc_osc > 0.5).sum()} of {len(os_)} with oscine phastCons posterior > 0.5")


if __name__ == "__main__":
    main()
