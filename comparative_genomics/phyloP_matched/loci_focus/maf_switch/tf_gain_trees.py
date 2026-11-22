#!/usr/bin/env python
"""Branch gains/losses on the 74-way tree for the 13 other MAFB regulators.

gain_trees.py run on the SCENIC+-nominated inputs instead of AR. Cactus stores
a reconstructed sequence for every ancestral node, so presence is read from
each node's own sequence and a gain or loss is OBSERVED on a branch rather
than inferred from the tips.

Counts are recorded for every site with coverage (counts.json). Trees are
drawn only for sites that clear the candidate filter, to keep the figure set
readable.
"""
import sys, collections, json, re
import numpy as np

sys.path.insert(0, ".")
from treelib import parse_newick, leaves
from gain_trees import count_events, draw
import sys as _s; _s.path.insert(0, "/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus"); from span import per_site  # results/per_site/

F = "/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus"
PC = 0.001; PAD = 30
OSC = {l.strip() for l in open("data74/osc74.txt") if l.strip()}
tree = open("data74/tree74.nw").read().strip()

LO, W = {}, {}
name = None; rows = []; want = 0
for l in open(f"{F}/tf13.meme"):
    if l.startswith("MOTIF"):
        if name and rows:
            LO[name] = np.log2((np.array(rows) + PC) / (1 + 4 * PC) / .25)
        name = l.split()[1]; rows = []; want = 0
    elif l.startswith("letter-probability"):
        want = int(re.search(r"w=\s*(\d+)", l).group(1)); W[name] = want
    elif want and len(l.split()) == 4 and len(rows) < want:
        try: rows.append([float(x) for x in l.split()])
        except ValueError: pass
if name and rows:
    LO[name] = np.log2((np.array(rows) + PC) / (1 + 4 * PC) / .25)
RC = {k: m[::-1, ::-1] for k, m in LO.items()}


def thresh(lo, n=200000, q=99.99):
    rng = np.random.default_rng(0); w = len(lo)
    s = rng.integers(0, 4, (n, w))
    return float(np.percentile(np.array([lo[k][s[:, k]] for k in range(w)]).sum(0), q))


TH = {k: thresh(v) for k, v in LO.items()}
CODE = np.full(256, 4, np.int8)
for i, c in enumerate("ACGT"):
    CODE[ord(c)] = i; CODE[ord(c.lower())] = i


def best(seq, k):
    w = W[k]
    a = CODE[np.frombuffer(seq.encode(), np.uint8)]
    if a.size < w:
        return -99.0
    sl = np.lib.stride_tricks.sliding_window_view(a, w)
    ok = (sl < 4).all(1)
    if not ok.any():
        return -99.0
    sl = sl[ok]; idx = np.arange(w)
    return float(max(LO[k][idx, sl].sum(1).max(), RC[k][idx, sl].sum(1).max()))


blocks = []; cur = []
for l in open("data74/tf_windows_gap_74.maf"):
    if l.startswith("a"):
        if cur: blocks.append(cur)
        cur = []
    elif l.startswith("s\t"):
        f = l.rstrip("\n").split("\t"); cur.append((f[1], int(f[2]), int(f[3]), f[6]))
if cur: blocks.append(cur)
print(f"  blocks {len(blocks):,}")

sites = []
for l in open(f"{F}/tracks/mafb_1p5mb/tf_sites_MAFB_lsd2.bed"):
    x = l.rstrip("\n").split("\t")
    sites.append(dict(motif=x[9], name=x[3], s=int(x[1]), e=int(x[2]),
                      cls=x[10], elem=x[11]))

man = []; counts = []
for st in sites:
    k = st["motif"]; w0, w1 = st["s"] - PAD, st["e"] + PAD
    seq = collections.defaultdict(str)
    for b in blocks:
        ref = [r for r in b if r[0].startswith("Lonchura_striata.")]
        if not ref: continue
        coord = []; p = ref[0][1]
        for c in ref[0][3]:
            if c == "-": coord.append(None)
            else: coord.append(p); p += 1
        real = [c for c in coord if c is not None]
        if not real or real[-1] < w0 or real[0] >= w1: continue
        keep = [i for i, c in enumerate(coord) if c is None or w0 <= c < w1]
        if not keep: continue
        for sp, _, _, t in b:
            seq[sp.split(".")[0]] += "".join(t[i] for i in keep if i < len(t)).replace("-", "")
    pres = {g: (best(s, k) >= TH[k]) for g, s in seq.items() if len(s) >= W[k]}
    ext = [g for g in pres if not g.startswith("Anc")]
    o = [g for g in ext if g in OSC]; u = [g for g in ext if g not in OSC]
    if len(o) < 8 or len(u) < 20: continue
    fo = sum(pres[g] for g in o) / len(o); fu = sum(pres[g] for g in u) / len(u)
    ng, nl = count_events(pres, tree)
    counts.append(dict(**{j: st[j] for j in ("motif", "name", "cls", "elem")},
                       f_osc=round(fo, 3), f_out=round(fu, 3), gains=ng, losses=nl,
                       n_osc=len(o), n_out=len(u),
                       pass_freq=bool(fo >= 0.6 and fu <= 0.35)))
    if not (fo >= 0.6 and fu <= 0.35) or ng > 2:
        continue
    nm = f"{k}_{st['name']}"
    N, g2, l2 = draw(pres, tree, f"{per_site('gain_trees_tf', nm)}/{nm}", OSC,
                     f"{st['name']}  ({k})"
                     f"  chr20:{st['s'] + 1:,}-{st['e']:,}   oscine {sum(pres[g] for g in o)}/{len(o)}"
                     f"  outgroup {sum(pres[g] for g in u)}/{len(u)}")
    man.append(dict(file=nm + ".pdf", **{j: st[j] for j in ("motif", "name", "cls", "elem")},
                    f_osc=round(fo, 3), f_out=round(fu, 3), tips=N, gains=g2, losses=l2))
    print(f"    {nm:<46} osc {fo:.2f} out {fu:.2f}   {g2} gains, {l2} losses")

json.dump(man, open("trees_tf/manifest.json", "w"), indent=1)
json.dump(counts, open("trees_tf/counts.json", "w"), indent=1)
print(f"  sites with coverage: {len(counts)}   trees drawn: {len(man)}")
c = collections.Counter(x["motif"] for x in man)
for m in sorted(c):
    print(f"    {m:<8} {c[m]}")
