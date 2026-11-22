#!/usr/bin/env python
"""Per-regulator matrix panel for the motif-definition sensitivity analysis.

For each MAFB regulator, every defensible matrix, tagged by where it comes from:

  cluster    every member of the regulator's top SCENIC+ motif (the cluster
             build_scenic_motifs.py takes the medoid of); `medoid` marks the one
             the main analyses use
  homer      HOMER known motifs, taken from the cistarget collection's homer__
             entries: named for the TF or a direct alias (own), or for close
             paralogs when the regulator has none of its own or is reported at
             family level (family), or a composite with the TF (composite)
  denovo     HOMER de novo motifs from the Glut-RA vs negative-set run
             (ra_pos_vs_neg). Discovered on Glut-RA DARs, so a Glut-RA vs GABA-4
             contrast is not neutral for them -- tagged so they can be set apart.

NR3C2 shares AR's cluster (metacluster_185.1) and is not scanned separately.
Matrices of 4-30 columns, written from counts at full width. Duplicates (the same
matrix reached by two routes) are kept once, with every route recorded.

Writes motif_panel.meme and results/regulator_motifs/motif_panel.json.
"""
import collections, glob, json, os, re, sys
import numpy as np

F = "/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus"
V = "/private/groups/colquittlab/scenicplus/v10nr_clust_public/singletons"
HM = "/private/home/bcolquit/group/motifs/homer/ra_pos_vs_neg/homerResults_manual"
sys.path.insert(0, F)
from span import result
from build_scenic_motifs import cb_split, norm, consensus

CHOICE = json.load(open(result("tf_scenic_choice.json", "regulator_motifs")))
TFS = [t for t in CHOICE if t != "NR3C2"]

HOMER = {   # tf -> [(homer motif id, relation)]
    "AR": [("homer__AGAACAGNCTGTTCTT_ARE", "own"), ("homer__CCAGGAACAG_AR-halfsite", "own"),
           ("homer__AGTAAACAAAAAAGAACANA_FOXA1_AR", "composite")],
    "KLF12": [("homer__VDGGGYGGGGCY_KLF1", "family"), ("homer__NRGCCCCRCCCHBNN_KLF3", "family"),
              ("homer__GCCACACCCA_Klf4", "family"), ("homer__DGGGYGKGGC_KLF5", "family"),
              ("homer__MKGGGYGTGGCC_KLF6", "family"), ("homer__GCCACRCCCACY_Klf9", "family"),
              ("homer__GGGGGTGTGTCC_KLF10", "family"), ("homer__RGKGGGCGKGGC_KLF14", "family")],
    "NFIA": [("homer__TTGCCAAG_NF1-halfsite", "family"), ("homer__NNTGTTTATTTTGGCA_NF1_FOXA1", "composite")],
    "STAT1": [("homer__NATTTCCNGGAAAT_STAT1", "own")],
    "FOSL2": [("homer__NATGASTCABNN_Fosl2", "own"), ("homer__GGATGACTCATC_Fra2", "own"),
              ("homer__NNATGASTCATH_Fra1", "family")],
    "MAF": [("homer__TGCTGACTCA_MafA", "family"), ("homer__WNTGCTGASTCAGCANWTTY_MafB", "family"),
            ("homer__HWWGTCAGCAWWTTT_MafF", "family"), ("homer__GCTGASTCAGCA_MafK", "family")],
    "SOX6": [("homer__CCATTGTTNY_Sox6", "own")],
    "TCF12": [("homer__ACAGCTGCTG_Tcf12", "own"), ("homer__VCAGCTGBNN_HEB", "own")],
    "MEIS2": [("homer__VGCTGWCAVB_Meis1", "family")],
}
DENOVO = {"AR": ["motif2", "motif4"], "FOSL2": ["motif1"], "NFIA": ["motif3"],
          "TCF12": ["motif7"], "MAF": ["motif8"], "KLF12": ["motif19"]}


def homer_counts(p):
    rows = [[float(x) for x in l.split()] for l in open(p) if not l.startswith(">") and l.strip()]
    return np.array(rows) * 100


def main():
    lib = {}
    for p in glob.glob(f"{V}/*.cb"):
        for n, a in cb_split(p): lib.setdefault(n, a)
    panel = collections.OrderedDict()

    def add(tf, mid, cnt, source, relation=""):
        if not (4 <= len(cnt) <= 30): return
        key = (tf, mid)
        if key in panel:
            panel[key]["routes"].append(f"{source}:{relation}" if relation else source); return
        p = norm(cnt)
        panel[key] = dict(tf=tf, matrix=mid, source=source, relation=relation, width=len(p),
                          consensus=consensus(p), routes=[f"{source}:{relation}" if relation else source],
                          medoid=(mid == CHOICE[tf]["matrix"]), p=p)
    for tf in TFS:
        c = CHOICE[tf]
        for n, a in cb_split(f"{V}/{c['cluster']}.cb"):
            add(tf, n, a, "cluster")
        for mid, rel in HOMER.get(tf, []):
            if mid in lib: add(tf, mid, lib[mid], "homer", rel)
            else: print(f"  missing {mid}")
        for m in DENOVO.get(tf, []):
            p = f"{HM}/{m}.motif"
            lab = open(p).readline().split("\t")[1].strip()
            add(tf, f"homer_denovo__{m}_{lab}", homer_counts(p), "denovo")
    meme = ["MEME version 4", "", "ALPHABET= ACGT", "", "strands: + -", "",
            "Background letter frequencies", "A 0.25 C 0.25 G 0.25 T 0.25", ""]
    out = []
    for i, ((tf, mid), d) in enumerate(panel.items()):
        uid = f"P{i:04d}"                                    # short, FIMO-safe id
        d["id"] = uid
        meme += [f"MOTIF {uid} {tf}|{d['source']}|{mid}", "",
                 f"letter-probability matrix: alength= 4 w= {d['width']} nsites= 20 E= 0", ""]
        meme += [" " + "  ".join(f"{x:.6f}" for x in r) for r in d["p"]] + [""]
        out.append({k: v for k, v in d.items() if k != "p"})
    open(f"{F}/motif_panel.meme", "w").write("\n".join(meme) + "\n")
    json.dump(out, open(result("motif_panel.json", "regulator_motifs"), "w"), indent=1)
    by = collections.Counter((d["tf"], d["source"]) for d in out)
    print(f"  {len(out)} matrices")
    for tf in TFS:
        print(f"    {tf:7s} cluster {by[(tf, 'cluster')]:3d}   homer {by[(tf, 'homer')]:2d}   denovo {by[(tf, 'denovo')]}")


if __name__ == "__main__":
    main()
