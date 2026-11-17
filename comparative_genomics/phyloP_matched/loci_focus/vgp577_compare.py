#!/usr/bin/env python
"""Does the VGP 577-way change element-level phyloP? Compare with the 363-way per-base tracks.

Elements (results/vgp577/elements_lsd2.bed): the MAFB-domain lobes (scored extents) and a random
3,000 RA + 3,000 G4 genome-wide DAR peaks. They are lifted lonStrDom2 -> bLonStr1 (liftOver chain,
minMatch 0.5) -> zebra finch GCA_048771995.1 (halLiftover through vgp-577way-v1.hal; run_lift.sh).
Fragments are kept on the zebra finch sequence carrying most of an element's lifted bases, within
3x the element's length of each other. Per element:
  vgp_birds   mean VGP zebra-finch-referenced phyloP, bird subtree model (rBirdsAnc0; 137 birds)
  vgp_all     mean VGP phyloP, whole 577-vertebrate model
  ours_whole  mean 363-way matched per-base whole-tree phyloP over the lonStrDom2 element
  ours_osc    mean 363-way oscine-lineage (subtree birdAnc222) phyloP
and the lifted fraction of the element.

    python vgp577_compare.py
Writes results/vgp577/vgp577_compare.tsv.
"""
import numpy as np, pandas as pd, pyBigWig
from scipy import stats

W = "/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus/results/vgp577"
G = "/private/home/ghickey/dev/work/vgp-cactus/577way/conservation/GCA_048771995.1/vgp-577way-v1-zebra_finch.single.phyloP"
PB = "/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/perbase/tracks"


def main():
    E = pd.read_csv(f"{W}/elements_lsd2.bed", sep="\t", header=None, names=["chrom", "start", "end", "id"])
    Z = pd.read_csv(f"{W}/elements_zf_frag.bed", sep="\t", header=None, usecols=[0, 1, 2, 3], names=["c", "s", "e", "id"])
    vb, va = pyBigWig.open(f"{G}.rBirdsAnc0.bw"), pyBigWig.open(f"{G}.bw")
    ow, oo = pyBigWig.open(f"{PB}/wholetree_perbase_lsd2.bw"), pyBigWig.open(f"{PB}/subtree222_perbase_lsd2.bw")
    L = dict(zip(E.id, E.end - E.start)); rows = []
    for i, z in Z.groupby("id"):
        z = z.assign(bp=z.e - z.s); c = z.groupby("c").bp.sum().idxmax(); z = z[z.c == c].sort_values("s")
        mid = z.s.iloc[np.argmax(z.bp.values)]; z = z[(z.s - mid).abs() <= 3 * L[i]]
        v = {k: [] for k in ("b", "a")}
        for s, e in zip(z.s, z.e):
            for k, bw in (("b", vb), ("a", va)):
                if c in bw.chroms(): v[k].append(np.array(bw.values(c, int(s), int(e)), float))
        get = lambda k: np.nanmean(np.concatenate(v[k])) if v[k] and np.isfinite(np.concatenate(v[k])).any() else np.nan
        rows.append(dict(id=i, zf=f"{c}:{z.s.min() + 1}-{z.e.max()}", lifted_frac=z.bp.sum() / L[i], vgp_birds=get("b"), vgp_all=get("a")))
    V = E.merge(pd.DataFrame(rows), on="id", how="left")
    V["ours_whole"] = [ow.stats(c, s, e, type="mean")[0] for c, s, e in zip(V.chrom, V.start, V.end)]
    V["ours_osc"] = [oo.stats(c, s, e, type="mean")[0] for c, s, e in zip(V.chrom, V.start, V.end)]
    V["set"] = V.id.str.split("|").str[0]; V["cls"] = V.id.str.split("|").str[1].str.split(":").str[0]
    V.to_csv(f"{W}/vgp577_compare.tsv", sep="\t", index=False, float_format="%.3f")
    for st, v in V.groupby("set"):
        ok = v.dropna(subset=["vgp_birds", "ours_whole"]); ok = ok[ok.lifted_frac >= 0.5]
        print(f"\n== {st}: {len(v)} elements, {len(ok)} lifted to zebra finch (>= 50% of bases) with VGP phyloP")
        print(f"   Spearman ours whole vs VGP birds {stats.spearmanr(ok.ours_whole, ok.vgp_birds)[0]:+.2f}, "
              f"vs VGP 577 {stats.spearmanr(ok.ours_whole, ok.vgp_all)[0]:+.2f}; ours oscine vs VGP birds {stats.spearmanr(ok.ours_osc, ok.vgp_birds)[0]:+.2f}")
        for t, ours in (("ours_whole", "whole"),):
            q = ok[t].quantile([1 / 3, 2 / 3]).values; r = ok.vgp_birds.quantile([1 / 3, 2 / 3]).values
            a = pd.cut(ok[t], [-np.inf, *q, np.inf], labels=["low", "mid", "high"]); b = pd.cut(ok.vgp_birds, [-np.inf, *r, np.inf], labels=["low", "mid", "high"])
            print("   thirds, ours (rows) vs VGP birds (cols):"); print(pd.crosstab(a, b).to_string())
        print("   by class, median:"); print(ok.groupby("cls")[["ours_whole", "vgp_birds", "vgp_all"]].median().round(2).to_string())
        if st == "gw":
            for cl, g in ok.groupby("cls"): pass
            x, y = ok[ok.cls == "RA"], ok[ok.cls == "G4"]
            for k in ("ours_whole", "vgp_birds", "vgp_all"):
                print(f"   G4 vs RA {k:10s}: median {y[k].median():+.2f} vs {x[k].median():+.2f}, Mann-Whitney p {stats.mannwhitneyu(y[k], x[k]).pvalue:.2g}")
        unl = v[~v.id.isin(ok.id)]
        print(f"   not compared ({len(unl)}): median ours whole {unl.ours_whole.median():+.2f} vs compared {ok.ours_whole.median():+.2f}")


if __name__ == "__main__":
    main()
