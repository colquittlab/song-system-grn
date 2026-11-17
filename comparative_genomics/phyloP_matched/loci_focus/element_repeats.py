#!/usr/bin/env python
"""Sequence composition of accessible elements by constraint: are weakly constrained ones repeat-derived?

Elements:
  domain   MAFB-domain lobes over their scored extent (element_rule_classes.tsv), with the
           element_cluster_k<k>.tsv group
  genome   all SnapATAC2 peaks of the up cluster overlapping an own-contrast ArchR DAR (FDR < 0.1,
           Log2FC > 0), genome-wide: Glut-RA peaks for RA, GABA-4-1 peaks for G4
Per element:
  te_frac / te_class   fraction of bases in a RepeatMasker TE (lonStrDom2_TE.bed; TEs only, run with
                       -nolow) and the class covering most bases
  gc, cpg_oe           GC fraction; CpG observed / expected
  dust_frac            fraction masked by dustmasker (low complexity)
  trf_frac             fraction in a TRF tandem repeat (2 7 7 80 10 50 500)
  whole / oscine       mean matched per-base phyloP (whole tree / oscine lineage Anc222); for
                       domain lobes, the accessibility-weighted full-extent means
Genome-wide, elements are split into thirds of whole-tree and of oscine mean phyloP; "weak" = lower
third on both, "constrained" = upper third on either. Composition is compared between them
(Mann-Whitney; Fisher for any-TE overlap and for each TE class), within each class.

    python element_repeats.py [k=4]
Writes results/mafb_domain/element_repeats.tsv and results/genome_wide/dar_peak_repeats.tsv.gz.
"""
import os, subprocess, sys, tempfile
from collections import defaultdict
import numpy as np, pandas as pd, pyBigWig, pysam
from scipy import stats

F = "/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus"
R = f"{F}/results"
T = "/private/warm-archive/colquittlab/assembly/tracks_lark"
PB = "/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/perbase/tracks"
FA = "/private/groups/colquittlab/assembly/lonStrDom2/GCF_005870125.1_lonStrDom2_genomic_ucsc_only.fna"
TE = "/private/warm-archive/colquittlab/assembly/cactus/lonStrDom2_repeats/lonStrDom2_TE.bed"
TRF = "/private/groups/colquittlab/opt/trf409.linux64"
CLS = {"RA": ("Glut-RA", "ra"), "G4": ("GABA-4-1", "pvalb2_lge1")}


def load_te():
    te = defaultdict(list)
    for l in open(TE):
        c, s, e, k = l.split("\t")[:4]; te[c].append((int(s), int(e), k.split("/")[0] if "/" in k else k))
    out = {}
    for c, v in te.items():
        v.sort(); out[c] = (np.array([x[0] for x in v]), np.array([x[1] for x in v]), [x[2] for x in v])
    return out


def te_overlap(te, c, s, e):
    if c not in te: return 0.0, ""
    S, E, K = te[c]; i = np.searchsorted(E, s, "right"); cov = np.zeros(e - s, bool); by = defaultdict(int)
    # TE ends are not sorted by start order in general; scan a window back for long elements
    j0 = max(0, np.searchsorted(S, s - 20000)); j1 = np.searchsorted(S, e)
    for j in range(j0, j1):
        if E[j] > s and S[j] < e:
            a, b = max(S[j], s) - s, min(E[j], e) - s; cov[a:b] = True; by[K[j]] += b - a
    return cov.mean(), (max(by, key=by.get) if by else "")


def masked_frac(seqs, tool):
    """fraction of each sequence masked by dustmasker or TRF"""
    with tempfile.TemporaryDirectory(dir=os.environ.get("TMPDIR", "/data/tmp")) as d:
        fa = f"{d}/s.fa"
        with open(fa, "w") as fh:
            for i, s in enumerate(seqs): fh.write(f">s{i}\n{s}\n")
        mask = [np.zeros(len(s), bool) for s in seqs]
        if tool == "dust":
            out = subprocess.run(["dustmasker", "-in", fa, "-outfmt", "interval"], capture_output=True, text=True, check=True).stdout
            cur = None
            for l in out.splitlines():
                if l.startswith(">"): cur = int(l[2:]); continue
                a, b = map(int, l.replace(" ", "").split("-")); mask[cur][a:b + 1] = True
        else:
            out = subprocess.run([TRF, fa, "2", "7", "7", "80", "10", "50", "500", "-h", "-ngs"], capture_output=True, text=True, cwd=d).stdout
            cur = None
            for l in out.splitlines():
                if l.startswith("@"): cur = int(l[2:].split()[0]); continue
                x = l.split()
                if cur is not None and len(x) > 2: mask[cur][int(x[0]) - 1:int(x[1])] = True
    return np.array([m.mean() if len(m) else np.nan for m in mask])


def composition(df, fa, te):
    seqs = [fa.fetch(c, s, e).upper() for c, s, e in zip(df.chrom, df.start, df.end)]
    gc = np.array([(x.count("G") + x.count("C")) / max(1, len(x) - x.count("N")) for x in seqs])
    def oe(x):
        c, g, n = x.count("C"), x.count("G"), len(x)
        return x.count("CG") * n / (c * g) if c and g else np.nan
    t = [te_overlap(te, c, s, e) for c, s, e in zip(df.chrom, df.start, df.end)]
    return df.assign(te_frac=[a for a, _ in t], te_class=[b for _, b in t], gc=gc, cpg_oe=[oe(x) for x in seqs],
                     dust_frac=masked_frac(seqs, "dust"), trf_frac=masked_frac(seqs, "trf"))


def compare(a, b, label):
    rows = []
    for f in ("te_frac", "gc", "cpg_oe", "dust_frac", "trf_frac"):
        x, y = a[f].dropna(), b[f].dropna()
        rows.append(dict(feature=f, weak=x.mean(), constrained=y.mean(), p=stats.mannwhitneyu(x, y).pvalue))
    x1, y1 = (a.te_frac > 0).sum(), (b.te_frac > 0).sum()
    orr, p = stats.fisher_exact([[x1, len(a) - x1], [y1, len(b) - y1]])
    rows.append(dict(feature="any TE (fraction of elements)", weak=x1 / len(a), constrained=y1 / len(b), p=p))
    print(f"\n  {label}: weak n={len(a)} vs constrained n={len(b)}")
    print(pd.DataFrame(rows).to_string(index=False, float_format=lambda v: f"{v:.3g}"))
    for k in sorted(set(a.te_class) | set(b.te_class) - {""}):
        if not k: continue
        x1, y1 = (a.te_class == k).sum(), (b.te_class == k).sum()
        if x1 + y1 >= 5:
            orr, p = stats.fisher_exact([[x1, len(a) - x1], [y1, len(b) - y1]])
            print(f"    {k:8s} weak {x1:5d} ({100 * x1 / len(a):4.1f}%)  constrained {y1:5d} ({100 * y1 / len(b):4.1f}%)  OR {orr:5.2f}  p {p:.2g}")


def main():
    k = int(sys.argv[1]) if len(sys.argv) > 1 else 4
    fa, te = pysam.FastaFile(FA), load_te()
    # --- domain lobes
    S = pd.read_csv(f"{R}/mafb_domain/element_rule_classes.tsv", sep="\t")
    G = pd.read_csv(f"{R}/mafb_domain/element_cluster_k{k}.tsv", sep="\t")[["element", "group"]]
    S = S.merge(G, on="element").rename(columns={"ext_start": "start", "ext_end": "end"})
    S = composition(S, fa, te)
    S[["element", "cls", "group", "chrom", "start", "end", "whole_ext_mean", "osc_ext_mean", "acc_ra_ext", "acc_g4_ext",
       "te_frac", "te_class", "gc", "cpg_oe", "dust_frac", "trf_frac"]].to_csv(f"{R}/mafb_domain/element_repeats.tsv", sep="\t", index=False, float_format="%.3f")
    print("  MAFB domain lobes by cluster group (mean):")
    print(S.groupby("group").agg(n=("element", "size"), any_te=("te_frac", lambda v: (v > 0).mean()), te_frac=("te_frac", "mean"),
                                 gc=("gc", "mean"), cpg_oe=("cpg_oe", "mean"), dust=("dust_frac", "mean"), trf=("trf_frac", "mean")).round(3).to_string())
    print("  TE-overlapping lobes:"); print(S[S.te_frac > 0][["element", "group", "te_class", "te_frac", "whole_ext_mean", "osc_ext_mean"]].round(2).to_string(index=False))
    for f in ("te_frac", "gc", "cpg_oe", "dust_frac", "trf_frac"):
        print(f"    Spearman {f:9s} vs whole {stats.spearmanr(S[f], S.whole_ext_mean, nan_policy='omit')[0]:+.2f}   vs oscine {stats.spearmanr(S[f], S.osc_ext_mean, nan_policy='omit')[0]:+.2f}")
    # --- genome-wide DAR peaks
    bw = {"whole": pyBigWig.open(f"{PB}/wholetree_perbase_lsd2.bw"), "oscine": pyBigWig.open(f"{PB}/subtree222_perbase_lsd2.bw")}
    out = []
    for cl, (up, con) in CLS.items():
        d = pd.read_csv(f"{T}/panelE_inputs/dar_{con}.tsv", sep="\t"); d = d[(d.FDR < 0.1) & (d.Log2FC > 0)]; d["start"] = d.start - 1
        p = pd.read_csv(f"{T}/snapatac_peaks/{up}.bed", sep="\t", header=None, usecols=[0, 1, 2], names=["chrom", "start", "end"])
        keep = []
        for c, pp in p.groupby("chrom"):
            dd = d[d.chrom == c].sort_values("start")
            if dd.empty: continue
            ds, de = dd.start.values, dd.end.values; mx = np.maximum.accumulate(de)
            for s, e in zip(pp.start, pp.end):
                i = np.searchsorted(ds, e)            # DARs starting before peak end
                if i and mx[i - 1] > s: keep.append((c, s, e))
        q = pd.DataFrame(keep, columns=["chrom", "start", "end"]).assign(cls=cl)
        for t, b in bw.items():
            ch = b.chroms()
            q[t] = [b.stats(c, s, e, type="mean")[0] if c in ch else None for c, s, e in zip(q.chrom, q.start, q.end)]
        out.append(q)
    Q = pd.concat(out).dropna(subset=["whole", "oscine"]).reset_index(drop=True)
    print(f"\n  genome-wide DAR peaks with phyloP: " + ", ".join(f"{c} {n}" for c, n in Q.cls.value_counts().items()))
    Q = composition(Q, fa, te)
    Q.to_csv(f"{R}/genome_wide/dar_peak_repeats.tsv.gz", sep="\t", index=False, float_format="%.3f")
    for cl in ("RA", "G4"):
        q = Q[Q.cls == cl]; w3 = q.whole.quantile([1 / 3, 2 / 3]).values; o3 = q.oscine.quantile([1 / 3, 2 / 3]).values
        weak = q[(q.whole <= w3[0]) & (q.oscine <= o3[0])]; con = q[(q.whole >= w3[1]) | (q.oscine >= o3[1])]
        compare(weak, con, f"{cl} DAR peaks genome-wide")
        print("    any-TE fraction by whole-tree phyloP decile: " + " ".join(f"{v:.3f}" for v in q.groupby(pd.qcut(q.whole, 10, labels=False)).te_frac.apply(lambda v: (v > 0).mean())))


if __name__ == "__main__":
    main()
