#!/usr/bin/env python
"""Site-gain placement and the gain-rate test on the avian363 tree.

The 74-way tree has no basal oscines -- all 17 of its oscines are Passerides or
Corvus -- so its "oscine stem" is the core-oscine stem. avian363 samples the
basal lineages (lyrebird, scrubbird, treecreeper, bowerbird, fairy-wren,
honeyeaters, logrunner, Australasian babbler), so a gain can be placed on the
true oscine stem or a few nodes inside it.

Song and the song system are unlikely to have arisen in one genetic event, so a
gain a few nodes inside the stem is as much of interest as one exactly on it.
Placement is therefore reported, not just pass/fail:

  stem          gain on the branch into the oscine MRCA
  near-stem     gain on one of the deep oscine backbone nodes between the oscine
                MRCA and the core oscines -- the nodes still holding >= 90% of
                oscines (1-4 nodes inside the stem; coverage falls from 92% to 71%
                at the next node, where the Corvides split, so the boundary is a
                natural one)
  within        a shallower gain inside oscines
  inherited     present in oscines because the site predates them
  none          no gain in or above oscines

Per site the PRINCIPAL oscine gain -- the one with most oscines below it -- sets
the tier, unless the site is present at the oscine MRCA: then it is stem (gained on
the stem branch) or inherited, and gains inside oscines are regains after a loss.
(Before 2026-11-19 a regain inside oscines outranked a present MRCA, labelling such
sites near-stem or within.) Independent gains outside oscines are counted as homoplasy, not used to
discard the site; extra gains inside oscines (re-gains) are counted too.

Presence is the corrected call (presence.py: one reference-projected sequence per
node, missing below 80% motif coverage); ancestral states are cactus's
reconstructions. Coverage floor: >= 40 oscines and >= 60 non-oscines with
sequence.

Rate test: gains per unit branch length, summed over sites, for (a) the stem
branch and (b) the stem plus the near-stem backbone, each against the rest of the
root-to-Lonchura path (the ascertainment-matched comparison) and against branches
off that path. AR and the other regulators separately, in the domain and +/-1.5 Mb.

Writes trees363/counts_ar.json, trees363/counts_tf.json, data363/rate_test363.json.
"""
import collections, json, os, re, sys
import numpy as np
from scipy import stats

sys.path.insert(0, "."); sys.path.insert(0, "..")
import presence as PR
from treelib import parse_newick, leaves
from span import elem_coords, in_domain
from recount_cov import LO, TH, best

F = "/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus"
PAD = 30
OSC = {l.strip() for l in open("data/osc363.txt") if l.strip()}
OUT = {l.strip() for l in open("data/out363.txt") if l.strip()}
EL = elem_coords()

tree = parse_newick(open("data/tree363.nw").read())
full = {}; par = {}
def _reg(n):
    full[n["name"]] = n
    for c in n["ch"]: par[c["name"]] = n["name"]; _reg(c)
_reg(tree)
TIPS = {}
def _tips(n):
    s = {n["name"]} if not n["ch"] else set().union(*[_tips(c) for c in n["ch"]])
    TIPS[n["name"]] = frozenset(s); return s
_tips(tree)
INTERNAL = {n for n, v in full.items() if v["ch"]}
STEM = min((n for n in full if OSC <= TIPS[n]), key=lambda n: len(TIPS[n]))
path = []; cur = "Lonchura_striata"
while cur in par: path.append(cur); cur = par[cur]
PATH = set(path)
BACKBONE = [n for n in path if n != STEM and TIPS[n] <= OSC and len(TIPS[n] & OSC) >= .9 * len(OSC)]
DEPTH = {n: i + 1 for i, n in enumerate(sorted(BACKBONE, key=lambda n: -len(TIPS[n])))}
EDGES = [(p, c) for c, p in par.items()]


def where(g):
    if g == STEM: return "stem"
    if g in DEPTH: return "near-stem"
    T = TIPS[g]
    if T <= OSC: return "within"
    if not T & OSC: return "non-oscine"
    return "above"


def classify(pres):
    gains = [c for p, c in EDGES if pres.get(p) is not None and pres.get(c) is not None
             and pres[p] != pres[c] and pres[c]]
    losses = sum(1 for p, c in EDGES if pres.get(p) is not None and pres.get(c) is not None
                 and pres[p] != pres[c] and not pres[c])
    loc = {g: where(g) for g in gains}
    osc_g = [g for g in gains if loc[g] in ("stem", "near-stem", "within")]
    if STEM in gains:
        pg = STEM; tier = "stem"
    elif pres.get(STEM):
        # present at the oscine MRCA without a stem gain: the site predates oscines, and any gain inside oscines
        # is a regain after a loss (counted in other_osc_gains), not the origin of the site
        pg = None; tier = "inherited"
    elif osc_g:
        pg = max(osc_g, key=lambda g: len(TIPS[g] & OSC)); tier = loc[pg]
    elif any(loc[g] == "above" for g in gains):
        pg = None; tier = "inherited"
    else:
        pg = None; tier = "none"
    return dict(tier=tier, principal=pg, depth=(0 if pg == STEM else DEPTH.get(pg)),
                osc_cov=(len(TIPS[pg] & OSC) / len(OSC) if pg else None),
                gains=len(gains), losses=losses,
                nonosc_gains=sum(1 for g in gains if loc[g] == "non-oscine"),
                other_osc_gains=len(osc_g) - (pg is not None),
                gain_nodes=[[g, loc[g]] for g in gains])


def call(B, st):
    s, e, k = st["s"], st["e"], st["motif"]
    w0, w1 = s - PAD, e + PAD
    out = {}
    for g, (sq, nmot, _) in PR.project(PR.overlapping(*B, w0, w1), w0, w1, s, e).items():
        if not PR.covered(nmot, s, e): out[g] = None; continue
        v = best(sq, k) if len(sq) >= len(LO[k]) else None
        out[g] = v is not None and v >= TH[k]
    return out


def sites(fam):
    rows = []
    if fam == "AR":
        for bed, mo in ((f"{F}/tracks/mafb_1p5mb/arefull_sites_MAFB_lsd2.bed", "full"),
                        (f"{F}/tracks/mafb_1p5mb/jasparAR_sites_MAFB_lsd2.bed", "jaspar")):
            for l in open(bed):
                x = l.rstrip("\n").split("\t")
                rows.append(dict(motif=mo, name=x[3], s=int(x[1]), e=int(x[2]), cls=x[9], elem=x[10]))
    else:
        for l in open(f"{F}/tracks/mafb_1p5mb/tf_sites_MAFB_lsd2.bed"):
            x = l.rstrip("\n").split("\t")
            rows.append(dict(motif=x[9], name=x[3], s=int(x[1]), e=int(x[2]), cls=x[10], elem=x[11]))
    return rows


def rate(recs, keep):
    G = collections.Counter(); L = collections.Counter(); BL = collections.Counter(); n = 0
    for r in recs:
        if not keep(r): continue
        n += 1
        pres = r["_pres"]
        for p, c in EDGES:
            a, b = pres.get(p), pres.get(c)
            if a is None or b is None: continue
            BL[c] += full[c]["bl"]
            if a != b: (G if b else L)[c] += 1
    def tot(keys):
        return (sum(G[k] for k in keys), sum(L[k] for k in keys), sum(BL[k] for k in keys))
    stem = [STEM]; back = [STEM] + list(DEPTH)
    rest = [k for k in BL if k in PATH and k not in back]
    off = [k for k in BL if k not in PATH]
    out = dict(n_sites=n)
    for lab, keys in (("stem", stem), ("stem+backbone", back)):
        g0, l0, b0 = tot(keys); g1, l1, b1 = tot(rest); g2, l2, b2 = tot(off)
        e1 = b0 * g1 / b1 if b1 else np.nan; e2 = b0 * g2 / b2 if b2 else np.nan
        out[lab] = dict(branches=len(keys), gains=g0, losses=l0, length=b0,
                        rest_gains=g1, rest_losses=l1, rest_length=b1, rest_branches=len(rest),
                        off_gains=g2, off_losses=l2, off_length=b2,
                        exp_matched=e1, fold_matched=(g0 / e1 if e1 else np.nan),
                        p_matched=float(stats.poisson.sf(g0 - 1, e1)) if e1 == e1 else None,
                        exp_off=e2, fold_off=(g0 / e2 if e2 else np.nan),
                        p_off=float(stats.poisson.sf(g0 - 1, e2)) if e2 == e2 else None)
    return out


def main():
    os.makedirs("trees363", exist_ok=True); os.makedirs("data363", exist_ok=True)
    print(f"  oscine MRCA {STEM}; near-stem backbone {[(n, len(TIPS[n] & OSC)) for n in sorted(DEPTH, key=DEPTH.get)]}")
    RATE = {}
    fams = sys.argv[1:] or ["AR", "TF"]
    for fam, maf in (("AR", "data/ar_windows_gap.maf"), ("TF", "data/tfsc_windows_gap_anc.maf")):
        if fam not in fams: continue
        B = PR.load(maf)
        recs = []
        for st in sites(fam):
            pres = call(B, st)
            o = [g for g, v in pres.items() if g in OSC and v is not None]
            u = [g for g, v in pres.items() if g in OUT and v is not None]
            if len(o) < 40 or len(u) < 60: continue
            r = dict(motif=st["motif"], name=st["name"], cls=st["cls"], elem=st["elem"],
                     dom=in_domain(EL[st["elem"]]), n_osc=len(o), n_out=len(u),
                     f_osc=sum(pres[g] for g in o) / len(o), f_out=sum(pres[g] for g in u) / len(u),
                     **classify(pres))
            r["_pres"] = pres; recs.append(r)
        json.dump([{k: v for k, v in r.items() if k != "_pres"} for r in recs],
                  open(f"trees363/counts_{fam.lower()}.json", "w"), indent=1)
        t = collections.Counter(r["tier"] for r in recs)
        print(f"  {fam}: {len(recs)} sites   " + ", ".join(f"{k} {v}" for k, v in sorted(t.items())), flush=True)
        RATE[fam] = {sp: rate(recs, keep) for sp, keep in (("domain", lambda r: r["dom"]), ("1p5mb", lambda r: True))}
    if os.path.exists("data363/rate_test363.json"):
        old = json.load(open("data363/rate_test363.json"))["rate"]
        RATE = {**{k: v for k, v in old.items() if k not in RATE}, **RATE}
    json.dump(dict(stem=STEM, backbone={n: DEPTH[n] for n in DEPTH},
                   path=sorted(PATH), rate=RATE), open("data363/rate_test363.json", "w"), indent=1, default=float)
    for fam, R in RATE.items():
        if fam not in fams: continue
        for sp, d in R.items():
            for lab in ("stem", "stem+backbone"):
                x = d[lab]
                print(f"  {fam:2s} {sp:6s} {lab:14s} sites {d['n_sites']:3d}  gains {x['gains']:3d} over {x['length']:.2f}"
                      f"   matched: exp {x['exp_matched']:.2f}, {x['fold_matched']:.2f}x, p = {x['p_matched']:.3g}"
                      f"   off-path: {x['fold_off']:.1f}x")
    print("  wrote trees363/counts_{ar,tf}.json and data363/rate_test363.json")


if __name__ == "__main__":
    main()
