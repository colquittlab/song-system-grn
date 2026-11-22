#!/usr/bin/env python
"""Rebuild every site-gain table on the corrected presence call.

Presence comes from presence.project(): one reference-ordered sequence per
genome, each reference position taken once, so duplicated spans and paralogous
rows can no longer splice into chimeric sequence. A genome or ancestor is
COVERED when >=80% of the motif's reference positions carry a base; uncovered
genomes are missing, not absent -- dropped from the frequencies, and on the
tree no event is called on any branch touching one.

Writes, in the schemas the downstream scripts already read (originals kept as
*_concat.json):
  data/site_gain.json, data/tf_site_gain.json      avian363 frequencies
  trees_gain/counts.json, trees_tf/counts.json     74-way events, every site
  trees_gain/manifest.json, trees_tf/manifest.json trees drawn
and redraws the tree PDFs.

New fields on the 74-way records:
  gain_nodes   the branches (by child node) where the site was gained
  where        each gain's location: oscine stem / within oscines /
               above the oscine stem / non-oscine lineage
  tier         oscine-stem   exactly one gain, on the branch into the oscine
                             MRCA (Anc54) -- the oscine-specific gains
               single-other  exactly one gain, somewhere else
               two-gain      two independent gains
               turnover      more than two
               no-gain       none observed
"""
import collections, glob, json, os, shutil, sys
import numpy as np

sys.path.insert(0, "."); sys.path.insert(0, "..")
import presence as PR
import sys as _s; _s.path.insert(0, "/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus"); from span import per_site  # results/per_site/
from treelib import parse_newick, leaves
from gain_trees import draw
from recount_cov import LO, TH, best, sites, PAD, OSC363, OUT363, OSC74

STEM = "Anc54"
COLL = {"AR": "gain_trees_ar", "TF": "gain_trees_tf"}
FLOOR_OSC, FLOOR_OUT = 40, 60          # avian363 coverage floor, as before
FLOOR_OSC74, FLOOR_OUT74 = 8, 20       # 74-way floor, as before

tree_str = open("data74/tree74.nw").read().strip()
tree = parse_newick(tree_str)
full = {}; par = {}
def reg(n):
    full[n["name"]] = n
    for c in n["ch"]: par[c["name"]] = n["name"]; reg(c)
reg(tree)
edges = [(p, c) for c, p in par.items()]


def where(n):
    L = {x["name"] for x in leaves(full[n])}
    if n == STEM: return "oscine stem"
    if L <= OSC74: return "within oscines"
    if not L & OSC74: return "non-oscine lineage"
    return "above the oscine stem"


def tier(gn):
    if not gn: return "no-gain"
    if len(gn) == 1: return "oscine-stem" if gn[0] == STEM else "single-other"
    return "two-gain" if len(gn) == 2 else "turnover"


def call(B, st):
    """-> {genome: True/False/None}; None = no orthologous sequence at the motif."""
    s, e, k = st["s"], st["e"], st["motif"]
    w0, w1 = s - PAD, e + PAD
    out = {}
    for g, (sq, nmot, _) in PR.project(PR.overlapping(*B, w0, w1), w0, w1, s, e).items():
        if not PR.covered(nmot, s, e):
            out[g] = None; continue
        v = best(sq, k) if len(sq) >= len(LO[k]) else None
        out[g] = v is not None and v >= TH[k]
    return out


def backup(path):
    b = path.replace(".json", "_concat.json")
    if os.path.exists(path) and not os.path.exists(b):
        shutil.copy(path, b)


def main():
    S = sites()
    for fam, maf363, maf74, gfile, tdir in (
            ("AR", "data/ar_windows_gap.maf", "data74/ar_windows_gap_74.maf",
             "data/site_gain.json", "trees_gain"),
            ("TF", "data/tfsc_windows_gap.maf", "data74/tfsc_windows_gap_74.maf",
             "data/tf_site_gain.json", "trees_tf")):
        print(f"  {fam}: loading alignments", flush=True)
        B3 = PR.load(maf363); B7 = PR.load(maf74)
        for f in (gfile, f"{tdir}/counts.json", f"{tdir}/manifest.json"):
            backup(f)
        # trees are written to results/per_site/; the manifest says which are current
        gain, counts, man = [], [], []
        for st in [x for x in S if x["fam"] == fam]:
            keep = {j: st[j] for j in ("motif", "name", "cls", "elem")}
            # ---- avian363
            d = call(B3, st)
            o = [g for g, v in d.items() if g in OSC363 and v is not None]
            u = [g for g, v in d.items() if g in OUT363 and v is not None]
            if len(o) >= FLOOR_OSC and len(u) >= FLOOR_OUT:
                po, pu = sum(d[g] for g in o), sum(d[g] for g in u)
                gain.append(dict(**keep, s=st["s"], e=st["e"], n_osc=len(o), n_out=len(u),
                                 p_osc=po, p_out=pu, f_osc=po / len(o), f_out=pu / len(u),
                                 n_osc_nocov=sum(1 for g, v in d.items()
                                                 if g in OSC363 and v is None),
                                 n_out_nocov=sum(1 for g, v in d.items()
                                                 if g in OUT363 and v is None)))
            # ---- 74-way
            p7 = call(B7, st)
            ext = [g for g in p7 if not g.startswith("Anc") and p7[g] is not None]
            oo = [g for g in ext if g in OSC74]; uu = [g for g in ext if g not in OSC74]
            if len(oo) < FLOOR_OSC74 or len(uu) < FLOOR_OUT74:
                continue
            fo = sum(p7[g] for g in oo) / len(oo); fu = sum(p7[g] for g in uu) / len(uu)
            gn = []; nl = 0
            for a, b in edges:
                pa, pb = p7.get(a), p7.get(b)
                if pa is None or pb is None or pa == pb: continue
                if pb: gn.append(b)
                else: nl += 1
            pf = bool(fo >= .6 and fu <= .35)
            counts.append(dict(**keep, f_osc=round(fo, 3), f_out=round(fu, 3),
                               gains=len(gn), losses=nl, n_osc=len(oo), n_out=len(uu),
                               pass_freq=pf, gain_nodes=gn, where=[where(n) for n in gn],
                               tier=tier(gn)))
            if pf and len(gn) <= 2:
                stem = f'{st["motif"]}_{st["name"]}'
                lab = {"full": "ARE-full", "jaspar": "JASPAR AR"}.get(st["motif"], st["motif"])
                N, g2, l2 = draw(p7, tree_str, f"{per_site(COLL[fam], stem)}/{stem}", OSC74,
                                 f'{st["name"]}  ({lab})  chr20:{st["s"] + 1:,}-{st["e"]:,}'
                                 f'   oscine {sum(p7[g] for g in oo)}/{len(oo)}'
                                 f'  outgroup {sum(p7[g] for g in uu)}/{len(uu)}')
                assert (g2, l2) == (len(gn), nl), (stem, g2, l2, len(gn), nl)
                man.append(dict(file=stem + ".pdf", **keep, f_osc=round(fo, 3),
                                f_out=round(fu, 3), tips=N, gains=len(gn), losses=nl,
                                tier=tier(gn)))
        json.dump(gain, open(gfile, "w"), indent=1)
        json.dump(counts, open(f"{tdir}/counts.json", "w"), indent=1)
        json.dump(man, open(f"{tdir}/manifest.json", "w"), indent=1)
        t = collections.Counter(c["tier"] for c in counts)
        print(f"    {fam}: {len(gain)} sites scored (avian363), {len(counts)} on the 74-way,"
              f" {len(man)} trees   tiers: " + ", ".join(f"{k} {v}" for k, v in sorted(t.items())),
              flush=True)


if __name__ == "__main__":
    main()
