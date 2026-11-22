#!/usr/bin/env python
"""Site gain/loss for the 13 SCENIC+-nominated MAFB regulators, avian363.

Identical procedure to site_gain.py, run on the other predicted inputs instead
of AR. For each motif site called on lonStrDom2, every species' aligned
sequence across a +/-30 bp window is degapped back to that species' own
sequence and rescanned with the PWM on both strands; the species has the site
if its best hit clears the 99.99th percentile of that matrix's uniform-sequence
null, the stringency FIMO's p < 1e-4 implies. Indels are therefore a mechanism
of gain, not a confound -- which is the whole point of the gapped extraction.

The scan is vectorised here because there are 578 sites rather than 81.
"""
import sys, collections, json, re
import numpy as np

sys.path.insert(0, "."); sys.path.insert(0, "..")
F = "/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus"
PC = 0.001; PAD = 30
OSC = {l.strip() for l in open("data/osc363.txt") if l.strip()}
OUT = {l.strip() for l in open("data/out363.txt") if l.strip()}

# --- matrices -----------------------------------------------------------
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
assert all(len(LO[k]) == W[k] for k in LO), {k: (len(LO[k]), W[k]) for k in LO}

RC = {}            # reverse-complement matrix, so one pass scores both strands
for k, m in LO.items():
    RC[k] = m[::-1, ::-1]


def thresh(lo, n=200000, q=99.99):
    rng = np.random.default_rng(0); w = len(lo)
    s = rng.integers(0, 4, (n, w))
    v = np.array([lo[k][s[:, k]] for k in range(w)]).sum(0)
    return float(np.percentile(v, q))


TH = {k: thresh(v) for k, v in LO.items()}
print("  thresholds (99.99th pct of the uniform null):")
for k in sorted(TH):
    print(f"    {k:<8} w {W[k]:>2}   {TH[k]:+.2f}")

CODE = np.full(256, 4, np.int8)
for i, c in enumerate("ACGT"):
    CODE[ord(c)] = i; CODE[ord(c.lower())] = i


def best(seq, k):
    """Best log-odds over both strands, or -99 if nothing scoreable."""
    w = W[k]
    a = CODE[np.frombuffer(seq.encode(), np.uint8)]
    if a.size < w:
        return -99.0
    sl = np.lib.stride_tricks.sliding_window_view(a, w)       # (n, w)
    ok = (sl < 4).all(1)
    if not ok.any():
        return -99.0
    sl = sl[ok]
    idx = np.arange(w)
    f = LO[k][idx, sl].sum(1)
    r = RC[k][idx, sl].sum(1)
    return float(max(f.max(), r.max()))


# --- alignment ----------------------------------------------------------
blocks = []; cur = []
for l in open("data/tf_windows_gap.maf"):
    if l.startswith("a"):
        if cur: blocks.append(cur)
        cur = []
    elif l.startswith("s\t"):
        f = l.rstrip("\n").split("\t"); cur.append((f[1], int(f[2]), int(f[3]), f[6]))
if cur: blocks.append(cur)
byref = collections.defaultdict(list)
for b in blocks:
    ref = [r for r in b if r[0].startswith("Lonchura_striata.")]
    if ref: byref[(ref[0][1], ref[0][1] + ref[0][2])].append(b)
print(f"  alignment blocks: {len(blocks):,} over {len(byref)} reference spans")

sites = []
for l in open(f"{F}/tracks/mafb_1p5mb/tf_sites_MAFB_lsd2.bed"):
    x = l.rstrip("\n").split("\t")
    sites.append(dict(motif=x[9], name=x[3], s=int(x[1]), e=int(x[2]),
                      cls=x[10], elem=x[11]))

spans = sorted(byref)
res = []
for n, st in enumerate(sites):
    k = st["motif"]; w0, w1 = st["s"] - PAD, st["e"] + PAD
    seq = collections.defaultdict(str)
    for bs, be in spans:
        if be <= w0 or bs >= w1: continue
        for b in byref[(bs, be)]:
            ref = [r for r in b if r[0].startswith("Lonchura_striata.")][0]
            coord = []; p = ref[1]
            for c in ref[3]:
                if c == "-": coord.append(None)
                else: coord.append(p); p += 1
            keep = [i for i, c in enumerate(coord) if c is None or w0 <= c < w1]
            if not keep: continue
            for sp, _, _, t in b:
                g = sp.split(".")[0]
                if g.startswith("birdAnc"): continue
                seq[g] += "".join(t[i] for i in keep if i < len(t)).replace("-", "")
    o = [g for g in seq if g in OSC]; u = [g for g in seq if g in OUT]
    if len(o) < 40 or len(u) < 60: continue
    po = sum(1 for g in o if best(seq[g], k) >= TH[k])
    pu = sum(1 for g in u if best(seq[g], k) >= TH[k])
    res.append(dict(**{j: st[j] for j in ("motif", "name", "cls", "elem", "s", "e")},
                    n_osc=len(o), n_out=len(u), p_osc=po, p_out=pu,
                    f_osc=po / len(o), f_out=pu / len(u)))
    if (n + 1) % 50 == 0:
        print(f"    {n + 1}/{len(sites)} sites", flush=True)

json.dump(res, open("data/tf_site_gain.json", "w"), indent=1)
print(f"  sites scored: {len(res)} of {len(sites)}")
c = collections.Counter(r["motif"] for r in res)
for m in sorted(c):
    g = [r["f_osc"] - r["f_out"] for r in res if r["motif"] == m]
    print(f"    {m:<8} {c[m]:>4} sites   mean gain {np.mean(g):+.3f}"
          f"   max {max(g):+.3f}")
