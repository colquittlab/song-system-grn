#!/usr/bin/env python
"""Site gain recounted with missing data treated as missing, not absent.

The original presence call (site_gain.py, gain_trees.py, tf_*) scans each
genome's degapped sequence across a +/-30 bp window and calls the site absent
if the best hit misses threshold. A genome with no orthologous sequence across
the motif itself -- only a few aligned bases elsewhere in the window -- is then
scored "absent", which inflates outgroup absence, and an ancestor with no
sequence produces a "gain" on the branch below it that is really a transition
from missing to present.

Presence and coverage both come from presence.project(): one reference-ordered
sequence per genome, each reference position taken once, so duplicated
reference spans can no longer splice into chimeric sequence. A genome or
ancestor is COVERED at a site when >=80% of the motif's reference positions
carry a base; uncovered genomes are dropped from the frequencies and uncovered
nodes are missing on the tree, so no event is called on a branch touching one.

Writes data/site_cov.json (avian363 frequencies) and
trees_cov.json (74-way events), keyed by (matrix, site).
"""
import collections, json, re, sys
import numpy as np
from scipy import stats

sys.path.insert(0, "."); sys.path.insert(0, "..")
from treelib import parse_newick
from span import F, elem_coords, in_domain
import presence as PR

HM = "/private/home/bcolquit/group/motifs/homer/ra_pos_vs_neg/homerResults_manual"
PC = 0.001; PAD = 30; COV = 0.8
OSC363 = {l.strip() for l in open("data/osc363.txt") if l.strip()}
OUT363 = {l.strip() for l in open("data/out363.txt") if l.strip()}
OSC74 = {l.strip() for l in open("data74/osc74.txt") if l.strip()}


def _p_meme(f):
    out = {}; nm = None; rows = []; want = 0
    for l in open(f):
        if l.startswith("MOTIF"):
            if nm and rows: out[nm] = np.array(rows)
            nm = l.split()[1]; rows = []; want = 0
        elif l.startswith("letter-probability"):
            want = int(re.search(r"w=\s*(\d+)", l).group(1))
        elif want and len(l.split()) == 4 and len(rows) < want:
            try: rows.append([float(x) for x in l.split()])
            except ValueError: pass
    if nm and rows: out[nm] = np.array(rows)
    return out


P = {"full": np.array([[float(x) for x in l.split()] for l in open(f"{HM}/motif2.motif")
                       if not l.startswith(">") and l.strip()]),
     "jaspar": list(_p_meme(f"{F}/ar_jaspar.meme").values())[0]}
P.update(_p_meme(f"{F}/tf_scenic.meme"))   # SCENIC+ regulator matrices
LO = {k: np.log2((v + PC) / (1 + 4 * PC) / .25) for k, v in P.items()}
RC = {k: v[::-1, ::-1] for k, v in LO.items()}


def thresh(lo, n=200000, q=99.99):
    rng = np.random.default_rng(0); w = len(lo)
    s = rng.integers(0, 4, (n, w))
    return float(np.percentile(np.array([lo[k][s[:, k]] for k in range(w)]).sum(0), q))


TH = {k: thresh(v) for k, v in LO.items()}
CODE = np.full(256, 4, np.int8)
for i, c in enumerate("ACGT"):
    CODE[ord(c)] = i; CODE[ord(c.lower())] = i


def best(seq, k):
    w = len(LO[k])
    a = CODE[np.frombuffer(seq.encode(), np.uint8)]
    if a.size < w: return None
    sl = np.lib.stride_tricks.sliding_window_view(a, w)
    ok = (sl < 4).all(1)
    if not ok.any(): return None
    sl = sl[ok]; idx = np.arange(w)
    return float(max(LO[k][idx, sl].sum(1).max(), RC[k][idx, sl].sum(1).max()))


def per_site(B, s, e, k):
    """-> {genome: (present, covered)} from the reference projection."""
    w0, w1 = s - PAD, e + PAD
    bl = PR.overlapping(*B, w0, w1)
    out = {}
    for g, (sq, nmot, _) in PR.project(bl, w0, w1, s, e).items():
        v = best(sq, k) if len(sq) >= len(LO[k]) else None
        out[g] = (v is not None and v >= TH[k], PR.covered(nmot, s, e))
    return out


def sites():
    rows = []
    for bed, mo in ((f"{F}/tracks/mafb_1p5mb/arefull_sites_MAFB_lsd2.bed", "full"),
                    (f"{F}/tracks/mafb_1p5mb/jasparAR_sites_MAFB_lsd2.bed", "jaspar")):
        for l in open(bed):
            x = l.rstrip("\n").split("\t")
            rows.append(dict(motif=mo, name=x[3], s=int(x[1]), e=int(x[2]),
                             cls=x[9], elem=x[10], fam="AR"))
    for l in open(f"{F}/tracks/mafb_1p5mb/tf_sites_MAFB_lsd2.bed"):
        x = l.rstrip("\n").split("\t")
        rows.append(dict(motif=x[9], name=x[3], s=int(x[1]), e=int(x[2]),
                         cls=x[10], elem=x[11], fam="TF"))
    return rows


def main():
    S = sites(); el = elem_coords()
    tree = parse_newick(open("data74/tree74.nw").read())
    edges = []
    def walk(n):
        for c in n["ch"]: edges.append((n["name"], c["name"])); walk(c)
    walk(tree)
    stem = "Anc54"

    res363 = {}; res74 = {}
    for fam, maf363, maf74 in (("AR", "data/ar_windows_gap.maf", "data74/ar_windows_gap_74.maf"),
                               ("TF", "data/tf_windows_gap.maf", "data74/tf_windows_gap_74.maf")):
        print(f"  loading {fam} alignments", flush=True)
        B3 = PR.load(maf363); B7 = PR.load(maf74)
        for st in [x for x in S if x["fam"] == fam]:
            key = f'{st["motif"]}|{st["name"]}'
            # --- avian363 frequencies, covered genomes only
            d = per_site(B3, st["s"], st["e"], st["motif"])
            o = [g for g in d if g in OSC363]; u = [g for g in d if g in OUT363]
            oc = [g for g in o if d[g][1]]; uc = [g for g in u if d[g][1]]
            res363[key] = dict(
                motif=st["motif"], name=st["name"], cls=st["cls"], elem=st["elem"],
                dom=in_domain(el[st["elem"]]),
                n_osc=len(o), n_out=len(u),
                p_osc=sum(d[g][0] for g in o), p_out=sum(d[g][0] for g in u),
                n_osc_cov=len(oc), n_out_cov=len(uc),
                p_osc_cov=sum(d[g][0] for g in oc), p_out_cov=sum(d[g][0] for g in uc),
                out_absent_nocov=sum(1 for g in u if not d[g][1] and not d[g][0]))
            # --- 74-way events, uncovered nodes as missing
            d7 = per_site(B7, st["s"], st["e"], st["motif"])
            pres_old = {g: v[0] for g, v in d7.items()}
            pres_new = {g: (v[0] if v[1] else None) for g, v in d7.items()}
            def events(pres):
                g_ = l_ = 0; stem_g = 0; where = []
                for a, b in edges:
                    pa, pb = pres.get(a), pres.get(b)
                    if pa is None or pb is None or pa == pb: continue
                    if pb: g_ += 1; where.append(b); stem_g += (b == stem)
                    else: l_ += 1
                return g_, l_, stem_g, where
            ext = [g for g in d7 if not g.startswith("Anc")]
            oo = [g for g in ext if g in OSC74]; uu = [g for g in ext if g not in OSC74]
            fo = (sum(d7[g][0] for g in oo) / len(oo)) if oo else 0
            fu = (sum(d7[g][0] for g in uu) / len(uu)) if uu else 0
            uuc = [g for g in uu if d7[g][1]]
            fuc = (sum(d7[g][0] for g in uuc) / len(uuc)) if uuc else float("nan")
            go, lo_, so, _ = events(pres_old)
            gn, ln, sn, wn = events(pres_new)
            res74[key] = dict(f_osc=fo, f_out=fu, f_out_cov=fuc, n_out=len(uu),
                              n_out_cov=len(uuc), gains_old=go, losses_old=lo_,
                              gains=gn, losses=ln, stem_gain=sn, gain_nodes=wn,
                              pass_freq=bool(len(oo) >= 8 and len(uu) >= 20
                                             and fo >= .6 and fu <= .35))
        print(f"    {fam}: {sum(1 for x in S if x['fam'] == fam)} sites", flush=True)

    json.dump(res363, open("data/site_cov.json", "w"), indent=1)
    json.dump(res74, open("trees_cov.json", "w"), indent=1)
    print("  wrote data/site_cov.json and trees_cov.json")


if __name__ == "__main__":
    main()
