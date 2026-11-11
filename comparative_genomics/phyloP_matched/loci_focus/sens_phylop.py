#!/usr/bin/env python
"""Motif-definition sensitivity, paired phyloP: one row per matrix in the panel.

The statistic is tf_core.diffs, unchanged: per element, motif phyloP minus the
phyloP of the rest of the same DAR (motif covering >=50% of its bases, >=50 bases
of remainder); Glut-RA minus GABA-4 is the DoD; p from a permutation of the class
labels (5,000; Welch's t above 2,000 elements, where the two agree). Minimum elements per class: 40 genome-wide, 8 at the locus -- the
thresholds tf_tables.py uses.

Outlier z: each matrix placed against the fixed genome-wide archetype background
(182 archetypes, pairs_all.json), with the same weighted fit and dispersion
estimate allmotif_outlier.py uses.

Check: each regulator's medoid matrix is in the panel, so its rows must reproduce
results/*/tf_tab_*.json exactly; the script reports any difference.

Writes results/regulator_motifs/sens_phylop.json.
"""
import collections, glob, json, sys
import numpy as np

F = "/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus"
sys.path.insert(0, F)
from span import result, elem_coords, in_domain, wide_elements

NP = 5000
rng = np.random.default_rng(0)
PANEL = {d["id"]: d for d in json.load(open(result("motif_panel.json", "regulator_motifs")))}
_cls = {l.split("\t")[3]: l.rstrip("\n").split("\t")[4] for l in open(f"{F}/allelem_cls.bed")}
EL = elem_coords(); WIDE = wide_elements()
DOM = {e for e in WIDE if in_domain(EL[e])}
SPANS = {"gw": (None, 40), "MAFB": (WIDE, 8), "domain": (DOM, 8)}
CL = ("RA_over_Arco", "GABA4_over_GABA1")


def dars(tree):
    d = {}
    for l in open(f"{F}/{tree}_am_dars.tab"):
        p = l.split("\t"); d[p[0]] = (int(p[1]), int(p[2]), float(p[3]))
    return d


def values(tree):
    """-> {(matrix, element): motif - rest}, exactly as tf_core.diffs."""
    dar = dars(tree)
    agg = collections.defaultdict(lambda: [0, 0, 0.0])
    for f in sorted(glob.glob(f"{F}/sens/chunk*_{tree}.tab")):
        for l in open(f):
            p = l.split("\t"); m, e, _ = p[0].split("@")
            a = agg[(m, e)]; a[0] += int(p[1]); a[1] += int(p[2]); a[2] += float(p[3])
    out = {}
    for (m, e), (msz, mcov, msum) in agg.items():
        if e not in dar or e not in _cls: continue
        dsz, dcov, dsum = dar[e]
        rcov = dcov - mcov; rsum = dsum - msum
        if mcov < 0.5 * msz or rcov < 50: continue
        out[(m, e)] = msum / mcov - rsum / rcov
    return out


def perm(a, b):
    """Permutation p on the class labels; above 2,000 elements, where the two
    agree, Welch's t (a 5,000 x n permutation matrix would dominate the run)."""
    o = a.mean() - b.mean(); c = np.concatenate([a, b]); k = len(a)
    if len(c) > 2000:
        from scipy import stats
        return float(o), float(stats.ttest_ind(a, b, equal_var=False).pvalue)
    idx = np.argsort(rng.random((NP, len(c))), axis=1)
    pc = c[idx]; d = pc[:, :k].mean(1) - pc[:, k:].mean(1)
    return float(o), float((np.sum(np.abs(d) >= abs(o)) + 1) / (NP + 1))


def wls(x, y, w):
    W = w.sum(); mx = (w * x).sum() / W; my = (w * y).sum() / W
    b = (w * (x - mx) * (y - my)).sum() / (w * (x - mx) ** 2).sum(); return my - b * mx, b


def paule_mandel(x, y, s2):
    n = len(x); lo, hi = 0.0, max(1e-12, np.var(y, ddof=1) * 50)
    def Q(t2):
        w = 1 / (s2 + t2); a, b = wls(x, y, w); r = y - (a + b * x)
        return (r ** 2 / (s2 + t2)).sum() - (n - 2)
    if Q(lo) <= 0: return 0.0
    for _ in range(200):
        m = (lo + hi) / 2
        if Q(m) > 0: lo = m
        else: hi = m
    return (lo + hi) / 2


def main():
    V = {t: values(t) for t in ("subtree222", "wholetree")}
    by = {t: collections.defaultdict(dict) for t in V}
    for t, d in V.items():
        for (m, e), v in d.items(): by[t][m][e] = v
    # archetype background for the outlier z (genome-wide)
    arch = [r["gw"] for r in json.load(open(result("pairs_all.json", "compare_mafb_spans")))
            if r["set"] == "arch" and r.get("gw") and r["gw"]["na"] >= 40 and r["gw"]["nb"] >= 40]
    X = np.array([a["x"] for a in arch]); Y = np.array([a["y"] for a in arch])
    b0, a0 = np.polyfit(X, Y, 1)
    S2 = np.array([a["vy"] + b0 * b0 * a["vx"] - 2 * b0 * a["cxy"] for a in arch])
    t2 = paule_mandel(X, Y, S2); w = 1 / (S2 + t2); aw, bw = wls(X, Y, w)
    Wt = w.sum(); mxw = (w * X).sum() / Wt; sxx = (w * (X - mxw) ** 2).sum()

    rows = []
    for mid, meta in PANEL.items():
        rec = dict(id=mid, **{k: meta[k] for k in ("tf", "matrix", "source", "relation", "medoid", "width", "consensus")})
        for sp, (keep, mn) in SPANS.items():
            for t in ("subtree222", "wholetree"):
                d = by[t].get(mid, {})
                arr = {c: np.array([v for e, v in d.items() if _cls.get(e) == c and (keep is None or e in keep)])
                       for c in CL}
                a, b = arr[CL[0]], arr[CL[1]]
                key = f"{t}|{sp}"
                if len(a) < mn or len(b) < mn:
                    rec[key] = None; continue
                dd, p = perm(a, b)
                rec[key] = dict(ra=float(a.mean()), g4=float(b.mean()), dod=dd, p=p, na=len(a), nb=len(b))
        # outlier z, genome-wide, elements aligned across the two trees
        dO, dW = by["subtree222"].get(mid, {}), by["wholetree"].get(mid, {})
        pr = {c: [(dO[e], dW[e]) for e in dO if e in dW and _cls.get(e) == c] for c in CL}
        if len(pr[CL[0]]) >= 40 and len(pr[CL[1]]) >= 40:
            A = np.array(pr[CL[0]]); Bm = np.array(pr[CL[1]]); na, nb = len(A), len(Bm)
            y = A[:, 0].mean() - Bm[:, 0].mean(); x = A[:, 1].mean() - Bm[:, 1].mean()
            vy = A[:, 0].var(ddof=1) / na + Bm[:, 0].var(ddof=1) / nb
            vx = A[:, 1].var(ddof=1) / na + Bm[:, 1].var(ddof=1) / nb
            cxy = np.cov(A[:, 0], A[:, 1])[0, 1] / na + np.cov(Bm[:, 0], Bm[:, 1])[0, 1] / nb
            vline = 1 / Wt + (x - mxw) ** 2 / sxx
            z = (y - (aw + bw * x)) / np.sqrt(vy + bw * bw * vx - 2 * bw * cxy + t2 + vline)
            from scipy import stats
            rec["outlier_gw"] = dict(x=float(x), y=float(y), res=float(y - (aw + bw * x)), z=float(z),
                                     p=float(2 * stats.norm.sf(abs(z))))
        else:
            rec["outlier_gw"] = None
        rows.append(rec)
    json.dump(rows, open(result("sens_phylop.json", "regulator_motifs"), "w"), indent=1)

    # check: the medoid rows reproduce the main tables
    print("  medoid check against tf_tables (oscine DoD; gw / MAFB / domain):")
    for sp, d in (("gw", "genome_wide"), ("MAFB", "mafb_1p5mb"), ("domain", "mafb_domain")):
        tab = {r["tf"]: r for r in json.load(open(f"{F}/results/{d}/tf_tab_subtree222_{sp}.json"))}
        bad = []
        for r in rows:
            if not r["medoid"]: continue
            name = "KLF" if r["tf"] == "KLF12" else r["tf"]
            v = r.get(f"subtree222|{sp}")
            if name in tab and v and abs(v["dod"] - tab[name]["dod"]) > 1e-9: bad.append((name, v["dod"], tab[name]["dod"]))
            if (name in tab) != bool(v): bad.append((name, "present" if v else "absent", "in tab" if name in tab else "not in tab"))
        print(f"    {sp:7s} {'all match' if not bad else bad}")
    print(f"  wrote results/regulator_motifs/sens_phylop.json ({len(rows)} matrices)")


if __name__ == "__main__":
    main()
