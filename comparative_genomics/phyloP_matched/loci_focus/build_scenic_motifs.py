#!/usr/bin/env python
"""One matrix per MAFB regulator, chosen from the motif evidence SCENIC+ itself used.

Source runs (identical parameters): glut-ra_glut-c1h-1_allruns_hybrid config4 for the
Glut-RA regulators, AR included; ra-arco-hvc-nc_gaba_hybrid config4 for GABA-4.

Rule, the same for every TF:
  1. take the enriched motifs (cisTarget + DEM) annotated to the TF at the level of
     its MAFB eRegulon -- direct for all but SATB1, whose MAFB link is extended-only;
  2. take the highest-NES one;
  3. if it is a metacluster, use its medoid: the member matrix with the highest mean
     similarity to the other members (column Pearson correlation, best offset and
     strand, overlap >= 60% of the shorter matrix).

One matrix per TF, rather than the union of a cluster's members, so every TF is
scanned at the same per-base false-positive rate; a union over a 35-matrix
cluster would carry ~35x the background hits of a single-matrix TF.

Matrices are written from the .cb counts at full width. (The earlier tf13.meme and
ar_jaspar.meme conversions dropped each matrix's first column.)

Writes tf_scenic.meme and results/regulator_motifs/tf_scenic_choice.json.
"""
import json, os, re, sys
import numpy as np

F = "/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus"
V = "/private/groups/colquittlab/scenicplus/v10nr_clust_public/singletons"
sys.path.insert(0, F)
from span import result

EV = json.load(open(result("scenic_motif_evidence.json", "regulator_motifs")))
REG = {"AR": "Glut-RA", "KLF12": "Glut-RA", "NR3C2": "Glut-RA", "NFIA": "Glut-RA",
       "LIN28B": "Glut-RA", "STAT1": "Glut-RA", "FOSL2": "Glut-RA",
       "MAF": "GABA-4", "SOX6": "GABA-4", "TCF12": "GABA-4", "ARX": "GABA-4",
       "SATB1": "GABA-4", "MEIS2": "GABA-4"}


def cb_split(p):
    out = []; nm = None; rows = []
    for l in open(p):
        if l.startswith(">"):
            if nm is not None and rows: out.append((nm, np.array(rows)))
            nm = l[1:].split()[0]; rows = []
        else:
            q = l.split()
            if len(q) == 4:
                try: rows.append([float(x) for x in q])
                except ValueError: pass
    if nm is not None and rows: out.append((nm, np.array(rows)))
    return [(n, a) for n, a in out]


def norm(c):
    return (c + 1e-9) / (c + 1e-9).sum(1, keepdims=True)


def sim(a, b):
    """Mean column Pearson r at the best offset/strand; overlap >= 60% of the shorter."""
    best = -1.0
    for bb in (b, b[::-1, ::-1]):
        la, lb = len(a), len(bb); mn = max(4, int(np.ceil(.6 * min(la, lb))))
        for off in range(-(lb - mn), la - mn + 1):
            i0, j0 = max(0, off), max(0, -off)
            n = min(la - i0, lb - j0)
            if n < mn: continue
            x = a[i0:i0 + n]; y = bb[j0:j0 + n]
            xc = x - x.mean(1, keepdims=True); yc = y - y.mean(1, keepdims=True)
            r = (xc * yc).sum(1) / np.sqrt((xc ** 2).sum(1) * (yc ** 2).sum(1) + 1e-12)
            best = max(best, float(r.mean()))
    return best


IUPAC = {frozenset("AG"): "R", frozenset("CT"): "Y", frozenset("GC"): "S", frozenset("AT"): "W",
         frozenset("GT"): "K", frozenset("AC"): "M"}


def consensus(p):
    s = ""
    for r in p:
        ic = 2 + (r * np.log2(np.clip(r, 1e-9, 1))).sum()
        o = np.argsort(r)[::-1]
        c = "ACGT"[o[0]] if r[o[0]] >= .6 else (
            IUPAC.get(frozenset("ACGT"[o[0]] + "ACGT"[o[1]]), "n") if r[o[0]] + r[o[1]] >= .75 else "n")
        s += c if ic >= 1 else c.lower()
    return s


def main():
    choice = {}; meme = ["MEME version 4", "", "ALPHABET= ACGT", "", "strands: + -", "",
                         "Background letter frequencies", "A 0.25 C 0.25 G 0.25 T 0.25", ""]
    for tf in REG:
        ev = EV[tf]["motifs"]
        top = ev[0]["motif"]
        mats = [(n, a) for n, a in cb_split(f"{V}/{top}.cb") if 4 <= len(a) <= 30]
        if len(mats) == 1:
            name, cnt = mats[0]; how = "single matrix"; msim = None
        else:
            P = [norm(a) for _, a in mats]
            S = np.array([[sim(P[i], P[j]) if i != j else np.nan for j in range(len(P))]
                          for i in range(len(P))])
            k = int(np.nanargmax(np.nanmean(S, 1)))
            name, cnt = mats[k]; how = f"medoid of {len(mats)}"; msim = float(np.nanmean(S[k]))
        p = norm(cnt)
        choice[tf] = dict(regulon=REG[tf], run=EV[tf]["run"], level=EV[tf]["level"], cluster=top,
                          nes=ev[0]["nes"], matrix=name, how=how, mean_sim=msim, width=len(p),
                          consensus=consensus(p), runner_up=[m["motif"] for m in ev[1:3]])
        meme += [f"MOTIF {tf} {name}", "",
                 f"letter-probability matrix: alength= 4 w= {len(p)} nsites= 20 E= 0", ""]
        meme += [" " + "  ".join(f"{x:.6f}" for x in r) for r in p] + [""]
        print(f"  {tf:7s} {top:34s} NES {ev[0]['nes']:5.1f}  {how:14s} -> {name:42s} {consensus(p)}")
    open(f"{F}/tf_scenic.meme", "w").write("\n".join(meme) + "\n")
    json.dump(choice, open(result("tf_scenic_choice.json", "regulator_motifs"), "w"), indent=1)
    print(f"\n  wrote tf_scenic.meme ({len(choice)} matrices) and results/regulator_motifs/tf_scenic_choice.json")


if __name__ == "__main__":
    main()
