#!/usr/bin/env python
"""External expectations for oscine motif-site gain in the curated MAFB-domain lobes.

lobe_site_gain.py and plot_tf_gain_rank.py compare the regulators with one another, so a gain common to every
motif in these sequences (reference ascertainment, or oscine-young DNA) sets their expectation and every matrix
looks ordinary. Here the expectation comes from motifs that should carry no regulator-specific signal, scanned and
called exactly as the regulators are (same lobes, threshold rule, presence call, coverage floor, tree placement):
  shuffled     NSHUF column permutations of each regulator matrix (same width, information content and base
               composition; the binding preference is destroyed)
  background   NBG matrices drawn at random from the SCENIC+ motif collection (v10nr_clust singletons, width
               MINW-MAXW), "motifs in general" in the same DNA
  regulator    the 12 SCENIC+ regulators (NR3C2 = AR dropped) and the AR references (full, jaspar)
Per site: the corrected presence call (presence.py projection, missing below 80% motif coverage), oscine and
outgroup presence fractions (floors 40 / 60), gain, one-sided Fisher p, the avian363 tier (tiers363.classify),
and the branch-rate components: gains and covered branch length on the oscine stem, the stem + near-stem
backbone, the rest of the root-to-Lonchura path, and off the path (as tiers363.rate).

    cd maf_switch && python lobe_gain_null.py scan      # sites + calls -> results/mafb_domain/lobe_gain_null_sites.tsv.gz
    python lobe_gain_null.py report                     # tables, figure
"""
import collections, glob, os, re, sys
import numpy as np, pandas as pd, pysam
from scipy import stats

sys.path.insert(0, "."); sys.path.insert(0, "..")
import presence as PR
from recount_cov import LO as RLO, PC, PAD, CODE, thresh
from tiers363 import classify, OSC, OUT, STEM, EDGES, DEPTH, PATH, full

R = "/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus/results/mafb_domain"
FA = "/private/groups/colquittlab/assembly/lonStrDom2/GCF_005870125.1_lonStrDom2_genomic_ucsc_only.fna"
V = "/private/groups/colquittlab/scenicplus/v10nr_clust_public/singletons"
MAF = "data/lobes_gap.maf"
OUTF = os.environ.get("OUTF", f"{R}/lobe_gain_null_sites.tsv.gz")
NSHUF, NBG = int(os.environ.get("NSHUF", "20")), int(os.environ.get("NBG", "1000"))
MINW, MAXW = 6, 30
FLOOR_OSC, FLOOR_OUT, DEPTH_MIN = 40, 60, 0.8
NPROC = int(os.environ.get("NPROC", "48"))
DROP = {"NR3C2"}
BACK = {STEM} | set(DEPTH)


def lo_of(p):
    p = p / p.sum(1, keepdims=True)
    return np.log2((p + PC) / (1 + 4 * PC) / .25)


def read_cb(f):
    rows = [[float(x) for x in l.split()] for l in open(f) if not l.startswith(">") and len(l.split()) == 4]
    return np.array(rows) if rows else None


def matrices():
    """-> {name: (set, parent, lo)}"""
    M = {}
    for k, lo in RLO.items():                         # the regulators' own log-odds, as lobe_site_gain.py scans them
        if k in DROP: continue
        M[k] = ("regulator", k, lo)
        rng = np.random.default_rng(sum(map(ord, k)))          # fixed per matrix (str hash is salted)
        for i in range(NSHUF):
            M[f"{k}~shuf{i:02d}"] = ("shuffled", k, lo[rng.permutation(len(lo))])
    fs = sorted(glob.glob(f"{V}/*.cb"))
    rng = np.random.default_rng(0); n = 0
    for f in rng.permutation(fs):
        p = read_cb(f)
        if p is None or not MINW <= len(p) <= MAXW or (p.sum(1) <= 0).any(): continue
        M["bg:" + os.path.basename(f)[:-3]] = ("background", "", lo_of(p)); n += 1
        if n == NBG: break
    return M


def _th(item):
    k, lo = item
    return k, thresh(lo)


def scan(seq, lo, th):
    w = len(lo); a = CODE[np.frombuffer(seq.encode(), np.uint8)]
    if a.size < w: return []
    sl = np.lib.stride_tricks.sliding_window_view(a, w); ok = (sl < 4).all(1)
    idx = np.arange(w); s2 = np.where(sl < 4, sl, 0); rc = lo[::-1, ::-1]
    sc = np.where(ok, np.maximum(lo[idx, s2].sum(1), rc[idx, s2].sum(1)), -np.inf)
    out = []
    for i in np.flatnonzero(sc >= th):
        if out and i < out[-1][0] + w:
            if sc[i] > out[-1][1]: out[-1] = (i, sc[i])
        else: out.append((i, sc[i]))
    return out


def best(seq, lo):
    w = len(lo); a = CODE[np.frombuffer(seq.encode(), np.uint8)]
    if a.size < w: return None
    sl = np.lib.stride_tricks.sliding_window_view(a, w); ok = (sl < 4).all(1)
    if not ok.any(): return None
    sl = sl[ok]; idx = np.arange(w); rc = lo[::-1, ::-1]
    return float(max(lo[idx, sl].sum(1).max(), rc[idx, sl].sum(1).max()))


B = None; MT = None; TH = None


def call_site(job):
    k, s, e = job[0], job[1], job[2]
    lo, th = MT[k][2], TH[k]
    w0, w1 = s - PAD, e + PAD
    d = {}
    for g, (sq, nmot, _) in PR.project(PR.overlapping(*B, w0, w1), w0, w1, s, e).items():
        if not PR.covered(nmot, s, e): d[g] = None; continue
        v = best(sq, lo) if len(sq) >= len(lo) else None
        d[g] = v is not None and v >= th
    o = [g for g in OSC if d.get(g) is not None]; u = [g for g in OUT if d.get(g) is not None]
    c = classify(d)
    gs = collections.Counter(); bl = collections.Counter()
    for p_, c_ in EDGES:
        a, b = d.get(p_), d.get(c_)
        if a is None or b is None: continue
        cat = "back" if c_ in BACK else ("rest" if c_ in PATH else "off")
        bl[cat] += full[c_]["bl"]
        if c_ == STEM: bl["stem"] += full[c_]["bl"]
        if b and not a:
            gs[cat] += 1
            if c_ == STEM: gs["stem"] += 1
    return job + (len(o), len(u), sum(d[g] for g in o), sum(d[g] for g in u), c["tier"], d.get(STEM),
                  gs["stem"], gs["back"], gs["rest"], gs["off"], bl["stem"], bl["back"], bl["rest"], bl["off"])


def run_scan():
    global B, MT, TH
    import multiprocessing as mp
    MT = matrices()
    print(f"  matrices: " + ", ".join(f"{k} {v}" for k, v in collections.Counter(x[0] for x in MT.values()).items()), flush=True)
    with mp.Pool(NPROC) as pool:
        TH = dict(pool.map(_th, [(k, v[2]) for k, v in MT.items()], chunksize=8))
    el = pd.read_csv(f"{R}/element_threshold_classes.tsv", sep="\t")
    fa = pysam.FastaFile(FA); jobs = []
    for r in el.itertuples(index=False):
        seq = fa.fetch(r.chrom, r.ext_start, r.ext_end).upper()
        dok = bool(r.osc_depth_min50 >= DEPTH_MIN and r.nonosc_depth_min50 >= DEPTH_MIN)
        for k, (st, par, lo) in MT.items():
            for i, sc in scan(seq, lo, TH[k]):
                jobs.append((k, r.ext_start + i, r.ext_start + i + len(lo), st, par, r.element, r.cls, dok, round(float(sc), 3)))
    print(f"  {len(jobs)} sites: " + ", ".join(f"{k} {v}" for k, v in collections.Counter(j[3] for j in jobs).items()), flush=True)
    print("  loading alignment", flush=True)
    B = PR.load(MAF)
    rows = []
    with mp.Pool(NPROC) as pool:                      # forked after B is loaded: workers share it
        for i, r in enumerate(pool.imap_unordered(call_site, jobs, chunksize=20), 1):
            rows.append(r)
            if i % 2000 == 0: print(f"    {i}/{len(jobs)}", flush=True)
    cols = ["matrix", "start", "end", "set", "parent", "element", "cls", "depth_ok", "score", "n_osc", "n_out",
            "p_osc", "p_out", "tier", "stem_state", "g_stem", "g_back", "g_rest", "g_off",
            "bl_stem", "bl_back", "bl_rest", "bl_off"]
    d = pd.DataFrame(rows, columns=cols).sort_values(["set", "matrix", "start"])
    d["floor_ok"] = (d.n_osc >= FLOOR_OSC) & (d.n_out >= FLOOR_OUT)
    d["gain"] = d.p_osc / d.n_osc.where(d.n_osc > 0) - d.p_out / d.n_out.where(d.n_out > 0)
    d["fisher_p"] = [stats.fisher_exact([[a, n - a], [b, m - b]], alternative="greater")[1] if ok else np.nan
                     for a, n, b, m, ok in zip(d.p_osc, d.n_osc, d.p_out, d.n_out, d.floor_ok)]
    d.to_csv(OUTF, sep="\t", index=False, float_format="%.5g")
    print(f"  wrote {OUTF}")


if __name__ == "__main__":
    if sys.argv[1:2] == ["scan"]: run_scan()
    elif sys.argv[1:2] == ["report"]:
        import lobe_gain_null_report; lobe_gain_null_report.main()
