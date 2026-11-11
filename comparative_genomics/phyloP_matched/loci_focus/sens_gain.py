#!/usr/bin/env python
"""Motif-definition sensitivity, site gain at the MAFB locus: one summary per matrix.

For every matrix in the panel, the procedure rebuild_gain.py applies to the main
set:
  sites      the matrix's own FIMO hits in the +/-1.5 Mb DARs; overlapping hits in
             one element collapse to the best-scoring one, exactly W bp
  presence   presence.project() over a +/-30 bp window (one sequence per genome,
             each reference position once); a genome with bases at <80% of the
             motif's positions is missing, not absent; present if the best hit
             clears the matrix's 99.99th-percentile uniform-null score
  gain       fraction of covered oscines minus fraction of covered non-oscines
             (avian363; >=40 and >=60 covered)
  tier       74-way branch reconstruction; oscine-stem = exactly one gain, on the
             branch into the oscine MRCA (Anc54)

Per matrix and span (domain, +/-1.5 Mb): sites, mean gain, candidates (gain >=
0.50 and BH q < 0.05 within the matrix's sites in that span), oscine-stem gains.

Alignments: maf_switch/data/dar1p5mb_gap.maf and data74/dar1p5mb_gap_74.maf, every
+/-1.5 Mb DAR +/-30 bp, so every matrix's sites are covered by one extraction.

Writes results/regulator_motifs/sens_gain.json (per matrix) and
results/regulator_motifs/sens_gain_sites.json (per site).
"""
import collections, glob, json, re, sys
import numpy as np
from scipy import stats

F = "/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus"
M = f"{F}/maf_switch"
sys.path.insert(0, F); sys.path.insert(0, M)
from span import result, elem_coords, in_domain, wide_elements
import presence as PR
from treelib import parse_newick, leaves

PC = 0.001; PAD = 30; STEM = "Anc54"
PANEL = {d["id"]: d for d in json.load(open(result("motif_panel.json", "regulator_motifs")))}
OSC363 = {l.strip() for l in open(f"{M}/data/osc363.txt") if l.strip()}
OUT363 = {l.strip() for l in open(f"{M}/data/out363.txt") if l.strip()}
OSC74 = {l.strip() for l in open(f"{M}/data74/osc74.txt") if l.strip()}
EL = elem_coords(); WIDE = wide_elements()


def load_panel():
    P = {}; nm = None; rows = []; w = 0
    for l in open(f"{F}/motif_panel.meme"):
        if l.startswith("MOTIF"):
            if nm and rows: P[nm] = np.array(rows)
            nm = l.split()[1]; rows = []; w = 0
        elif l.startswith("letter-probability"):
            w = int(re.search(r"w=\s*(\d+)", l).group(1))
        elif w and len(l.split()) == 4 and len(rows) < w:
            try: rows.append([float(x) for x in l.split()])
            except ValueError: pass
    if nm and rows: P[nm] = np.array(rows)
    return P


def thresh(lo, n=200000, q=99.99):
    rng = np.random.default_rng(0); w = len(lo)
    s = rng.integers(0, 4, (n, w))
    return float(np.percentile(np.array([lo[k][s[:, k]] for k in range(w)]).sum(0), q))


CODE = np.full(256, 4, np.int8)
for i, c in enumerate("ACGT"):
    CODE[ord(c)] = i; CODE[ord(c.lower())] = i


def best(seq, lo, rc):
    w = len(lo)
    a = CODE[np.frombuffer(seq.encode(), np.uint8)]
    if a.size < w: return None
    sl = np.lib.stride_tricks.sliding_window_view(a, w)
    ok = (sl < 4).all(1)
    if not ok.any(): return None
    sl = sl[ok]; idx = np.arange(w)
    return float(max(lo[idx, sl].sum(1).max(), rc[idx, sl].sum(1).max()))


def bh(p):
    p = np.asarray(p, float); o = np.argsort(p); q = np.empty(len(p))
    q[o] = np.minimum(np.minimum.accumulate((p[o] * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1], 1)
    return q


def sites():
    start = {}
    for l in open(f"{F}/am_dars.bed"):
        x = l.rstrip("\n").split("\t"); start[x[3]] = int(x[1])
    per = collections.defaultdict(list)
    for f in sorted(glob.glob(f"{F}/sens/chunk*_hits.txt")):
        for i, l in enumerate(open(f)):
            if not i or not l.strip(): continue
            p = l.rstrip("\n").split("\t")
            if len(p) < 7 or p[2] not in WIDE or p[2] not in start: continue
            s0 = start[p[2]]
            per[(p[0], p[2])].append((s0 + int(p[3]) - 1, s0 + int(p[4]), float(p[6])))
    out = []
    for (m, e), v in per.items():
        v.sort(); groups = []
        for h in v:
            if groups and h[0] < groups[-1][-1][1]: groups[-1].append(h)
            else: groups.append([h])
        for j, g in enumerate(groups):
            a, b, sc = max(g, key=lambda x: x[2])
            out.append(dict(id=m, elem=e, s=a, e=b, k=j + 1))
    return out


def main():
    P = load_panel()
    LO = {k: np.log2((v + PC) / (1 + 4 * PC) / .25) for k, v in P.items()}
    RC = {k: v[::-1, ::-1] for k, v in LO.items()}
    TH = {k: thresh(v) for k, v in LO.items()}
    S = sites()
    print(f"  {len(S):,} sites across {len({s['id'] for s in S})} matrices", flush=True)
    tree = parse_newick(open(f"{M}/data74/tree74.nw").read())
    par = {}
    def idx(n):
        for c in n["ch"]: par[c["name"]] = n["name"]; idx(c)
    idx(tree)
    edges = [(p, c) for c, p in par.items()]
    B3 = PR.load(f"{M}/data/dar1p5mb_gap.maf"); print("  avian363 loaded", flush=True)
    B7 = PR.load(f"{M}/data74/dar1p5mb_gap_74.maf"); print("  74-way loaded", flush=True)

    def call(B, st):
        w0, w1 = st["s"] - PAD, st["e"] + PAD; k = st["id"]
        out = {}
        for g, (sq, nmot, _) in PR.project(PR.overlapping(*B, w0, w1), w0, w1, st["s"], st["e"]).items():
            if not PR.covered(nmot, st["s"], st["e"]): out[g] = None; continue
            v = best(sq, LO[k], RC[k]) if len(sq) >= len(LO[k]) else None
            out[g] = v is not None and v >= TH[k]
        return out

    rows = []
    for n, st in enumerate(S):
        if n and n % 1000 == 0: print(f"    {n:,}/{len(S):,} sites", flush=True)
        d = call(B3, st)
        o = [g for g, v in d.items() if g in OSC363 and v is not None]
        u = [g for g, v in d.items() if g in OUT363 and v is not None]
        if len(o) < 40 or len(u) < 60: continue
        po, pu = sum(d[g] for g in o), sum(d[g] for g in u)
        r = dict(id=st["id"], elem=st["elem"], s=st["s"], e=st["e"], dom=in_domain(EL[st["elem"]]),
                 cls=None, po=po, no=len(o), pu=pu, nu=len(u), gain=po / len(o) - pu / len(u),
                 p=float(stats.fisher_exact([[po, len(o) - po], [pu, len(u) - pu]])[1]), tier=None)
        p7 = call(B7, st)
        ext = [g for g in p7 if not g.startswith("Anc") and p7[g] is not None]
        if sum(1 for g in ext if g in OSC74) >= 8 and sum(1 for g in ext if g not in OSC74) >= 20:
            gn = [c for a, c in edges if p7.get(a) is not None and p7.get(c) is not None
                  and p7[a] != p7[c] and p7[c]]
            r["tier"] = ("no-gain" if not gn else ("oscine-stem" if gn == [STEM] else
                         ("single-other" if len(gn) == 1 else ("two-gain" if len(gn) == 2 else "turnover"))))
        rows.append(r)
    json.dump(rows, open(result("sens_gain_sites.json", "regulator_motifs"), "w"))

    summ = []
    for mid, meta in PANEL.items():
        rec = dict(id=mid, **{k: meta[k] for k in ("tf", "matrix", "source", "relation", "medoid", "width")})
        for sp, keep in (("domain", lambda r: r["dom"]), ("1p5mb", lambda r: True)):
            R = [r for r in rows if r["id"] == mid and keep(r)]
            if not R:
                rec[sp] = None; continue
            q = bh([r["p"] for r in R])
            rec[sp] = dict(n=len(R), mean_gain=float(np.mean([r["gain"] for r in R])),
                           cand=int(sum(1 for r, qq in zip(R, q) if r["gain"] >= .5 and qq < .05)),
                           stem=int(sum(1 for r in R if r["tier"] == "oscine-stem")),
                           n_tree=int(sum(1 for r in R if r["tier"] is not None)))
        summ.append(rec)
    json.dump(summ, open(result("sens_gain.json", "regulator_motifs"), "w"), indent=1)
    print(f"  wrote sens_gain.json ({len(summ)} matrices) and sens_gain_sites.json ({len(rows):,} sites)")


if __name__ == "__main__":
    main()
