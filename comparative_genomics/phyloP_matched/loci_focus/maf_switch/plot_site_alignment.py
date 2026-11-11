#!/usr/bin/env python
"""Compact alignment view of one candidate motif site, with a tree.

Zooms to a candidate site and shows the aligned sequence of a subsample of
species chosen from informative clades, with a pruned cladogram to the left
(--phylogram: scaled by branch length).

Layout, left to right:
  tree       pruned topology of the chosen alignment (--aln 74 | 363) for the
             species shown (--phylogram: scaled by branch length). Internal edges
             show the reconstructed ancestor's state (red present, grey absent,
             dotted no sequence) -- the call the gain placement uses. Each tip
             edge is coloured by the share of the genomes its row stands for
             that carry the site (see "genomes"): grey none, pale to full red
             1-100%, dotted when none has sequence.
  rows       each row shows the most typical of the genomes it stands for (see
             "genomes"): the member with the fewest mismatches to the others over
             the displayed window, so a species whose own sequence is atypical
             (e.g. a genus-level loss) is replaced by one that represents its
             group. The tree is the full tree pruned to the species actually
             shown. --as-listed shows the subsample's listed species instead.
  sequence   reference (lonStrDom2) in full; every other row shows '.' for a
             match to the reference and the base where it differs. Gaps are
             '-', no alignment is blank. Columns where the reference has a gap
             (insertions in other species) are kept and tinted, because an
             indel can create a site. The motif columns are banded, with the
             matrix logo above in reference orientation.
  score      the species' best PWM score in the +/-30 bp analysis window, and
             a filled dot when it clears the 99.99th-percentile threshold. The
             sequence scored is presence.project()'s reference projection --
             the call the gain statistics use -- so figure and statistic agree.
             A species with bases at fewer than 80% of the motif's reference
             positions is marked "no seq": missing, not absent, and the tree
             draws no event on any branch that touches it.
  genomes    how many genomes each shown species stands for, and how many of
             them carry the site (bar: red present, open absent, pale no
             sequence; label k/n over genomes with sequence), so a basal
             lineage of two genomes is not read with the weight of a family of
             thirty. A genome goes to its unambiguous nearest shown relative;
             one branching off an internal edge of the pruned tree goes to the
             nearest member of the smallest clade bracket containing it.
             Genomes in no bracket are counted in the caption.

Below the tips, two reconstructed ancestors bracket the gain branch, so the
substitution that created the site can be read directly.

usage:
  plot_site_alignment.py <site name> <matrix> [--aln 74|363] [--phylogram] [--paper WIDTH_MM [--row-mm R]] [--flank N] [--out STEM]
    matrix: full | jaspar | any tf_scenic.meme name (AR, FOSL2, MAF, ...)
  e.g. plot_site_alignment.py RA_Glut-HVC-2.394967_1 jaspar
"""
import argparse, collections, re, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams.update({"font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Liberation Sans", "DejaVu Sans"],
    "pdf.fonttype": 42, "ps.fonttype": 42})
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, PathPatch
from matplotlib.textpath import TextPath
from matplotlib.font_manager import FontProperties
from matplotlib.transforms import Affine2D

sys.path.insert(0, ".")
from treelib import parse_newick, leaves
import presence as PR
import sys as _s; _s.path.insert(0, "/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus"); from span import per_site  # results/per_site/

F = "/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus"
HM = "/private/home/bcolquit/group/motifs/homer/ra_pos_vs_neg/homerResults_manual"
PC = 0.001; PAD = 30; COV = 0.8
INK, INK2, INK3, RULE = "#171B22", "#3A4350", "#6E7887", "#DDE1E7"
PRES, ABS, GAINC, LOSSC = "#C8161A", "#B7BEC8", "#C7521B", "#2C6E7F"
PRES_LO = "#F0C4C2"              # tip-edge ramp: PRES_LO (few genomes with the site) -> PRES (all)
def ramp(f):
    """tip-edge colour for the fraction f of represented genomes carrying the site"""
    if f <= 0: return ABS
    lo, hi = np.array(matplotlib.colors.to_rgb(PRES_LO)), np.array(matplotlib.colors.to_rgb(PRES))
    return tuple(lo + (hi - lo) * f)
BASE = {"A": "#2E8B57", "C": "#2459A8", "G": "#C08A00", "T": "#C8161A"}
OSCC = "#8A1C1F"
MINUS = "\u2212"
MONO = FontProperties(family="DejaVu Sans Mono")

# informative clades, top to bottom; the reference is always first
SUBSAMPLE = [
    ("oscines", ["Lonchura_striata", "Taeniopygia_guttata", "Vidua_chalybeata",
                 "Melospiza_georgiana", "Catharus_ustulatus", "Sylvia_atricapilla",
                 "Parus_atricapillus", "Corvus_moneduloides"]),
    ("suboscine", ["Chiroxiphia_lanceolata"]),
    ("NZ wren", ["Acanthisitta_chloris"]),
    ("parrots", ["Ara_ararauna", "Strigops_habroptila"]),
    ("falcons", ["Falco_peregrinus"]),
    ("seriema", ["Cariama_cristata"]),
    ("other Neoaves", ["Aquila_chrysaetos", "Calypte_anna", "Streptopelia_turtur"]),
    ("Galloanserae", ["Gallus_gallus", "Anas_platyrhynchos"]),
]
# avian363: richer at the nodes that matter -- the lyrebird (Menura), sister to all
# other oscines, decides whether a gain is on the oscine stem; both suboscine
# radiations (New World Manacus/Tyrannus/Furnarius, Old World Smithornis)
SUBSAMPLE_363 = [
    ("oscines", ["Lonchura_striata", "Taeniopygia_guttata", "Vidua_chalybeata",
                 "Melospiza_melodia", "Catharus_fuscescens", "Sylvia_atricapilla",
                 "Parus_major", "Corvus_moneduloides", "Ptilonorhynchus_violaceus",
                 "Malurus_elegans", "Menura_novaehollandiae"]),
    ("suboscines", ["Manacus_manacus", "Tyrannus_savana", "Furnarius_figulus", "Smithornis_capensis"]),
    ("NZ wren", ["Acanthisitta_chloris"]),
    ("parrots", ["Melopsittacus_undulatus", "Nestor_notabilis"]),
    ("falcons", ["Falco_peregrinus"]),
    ("seriema", ["Cariama_cristata"]),
    ("other Neoaves", ["Aquila_chrysaetos", "Calypte_anna", "Columba_livia"]),
    ("Galloanserae", ["Gallus_gallus", "Anas_platyrhynchos"]),
]
ALN = {
    "74": dict(label="74-way", unit="relative time", tree="data74/tree74.nw", osc="data74/osc74.txt",
               maf_ar="data74/ar_windows_gap_74.maf", maf_tf="data74/tfsc_windows_gap_74.maf",
               subsample=SUBSAMPLE, collection="site_alignments", prefix="aln"),
    "363": dict(label="avian363", unit="substitutions/site", tree="data/tree363.nw", osc="data/osc363.txt",
                maf_ar="data/ar_windows_gap.maf", maf_tf="data/tfsc_windows_gap.maf",
                subsample=SUBSAMPLE_363, collection="site_alignments_363", prefix="aln363"),
}
# successive ancestors of the oscine MRCA, named by the clade they found
CLADE_UP = ["oscines", "Eupasseres", "Passeriformes", "Psittacopasserae",
            "Eufalconimorphae", "Australaves"]


# ---------------------------------------------------------------- matrices --
def _p_homer(f):
    return np.array([[float(x) for x in l.split()] for l in open(f)
                     if not l.startswith(">") and l.strip()])


def _p_meme(f, name=None):
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
    return out[name] if name else out


def load_matrix(key):
    if key == "full":
        return _p_homer(f"{HM}/motif2.motif"), "ARE-full"
    if key == "jaspar":
        return list(_p_meme(f"{F}/ar_jaspar.meme").values())[0], "JASPAR AR"
    return _p_meme(f"{F}/tf_scenic.meme", key), key


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


# ------------------------------------------------------------------- site --
def find_site(name, key, cfg):
    if key in ("full", "jaspar"):
        bed = f"{F}/tracks/mafb_1p5mb/{'arefull' if key == 'full' else 'jasparAR'}_sites_MAFB_lsd2.bed"
        maf = cfg["maf_ar"]
        for l in open(bed):
            x = l.rstrip("\n").split("\t")
            if x[3] == name:
                return dict(chrom=x[0], s=int(x[1]), e=int(x[2]), strand=x[5],
                            cls=x[9], elem=x[10]), maf
    else:
        maf = cfg["maf_tf"]
        for l in open(f"{F}/tracks/mafb_1p5mb/tf_sites_MAFB_lsd2.bed"):
            x = l.rstrip("\n").split("\t")
            if x[3] == name and x[9] == key:
                return dict(chrom=x[0], s=int(x[1]), e=int(x[2]), strand=x[5],
                            cls=x[10], elem=x[11]), maf
    sys.exit(f"site {name!r} with matrix {key!r} not found")


# ------------------------------------------------------------------- tree --
def prune(node, keep):
    """Pruned copy holding only `keep` tips; unary nodes collapsed. Each
    pruned node records the full-tree node it stands for (its MRCA)."""
    if not node["ch"]:
        return dict(name=node["name"], full=node["name"], ch=[],
                    tips=frozenset([node["name"]])) if node["name"] in keep else None
    kids = [k for k in (prune(c, keep) for c in node["ch"]) if k]
    if not kids: return None
    if len(kids) == 1: return kids[0]
    return dict(name=node["name"], full=node["name"], ch=kids,
                tips=frozenset().union(*[k["tips"] for k in kids]))


def order_tips(n, rank):
    """Rotate children so the reference's lineage comes first."""
    if not n["ch"]: return
    for c in n["ch"]: order_tips(c, rank)
    n["ch"].sort(key=lambda c: min(rank.get(t, 99) for t in c["tips"]))


def tip_order(n):
    return [n["name"]] if not n["ch"] else sum((tip_order(c) for c in n["ch"]), [])


# --------------------------------------------------------------- drawing --
def logo(ax, p, xs, y0, h, strand):
    """Information-content logo, one column per entry of xs, in reference
    orientation. xs need not be contiguous: insertion columns can fall
    between motif bases."""
    if strand == "-":
        p = p[::-1, ::-1]               # reverse complement: A<->T, C<->G
    ic = 2 + (p * np.log2(np.clip(p, 1e-9, 1))).sum(1)
    for i, (col, tot) in enumerate(zip(p, ic)):
        y = y0
        for j in np.argsort(col):
            hh = col[j] * tot / 2 * h
            if hh < .004 * h: continue
            tp = TextPath((0, 0), "ACGT"[j], size=1, prop=FontProperties(family="DejaVu Sans",
                                                                          weight="bold"))
            bb = tp.get_extents()
            tr = (Affine2D().translate(-bb.x0, -bb.y0)
                  .scale(.92 / bb.width, hh / bb.height)
                  .translate(xs[i] - .46, y))
            ax.add_patch(PathPatch(tr.transform_path(tp), color=BASE["ACGT"[j]], lw=0))
            y += hh
    return ic


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("site"); ap.add_argument("matrix")
    ap.add_argument("--flank", type=int, default=12)
    ap.add_argument("--out", default=None)
    ap.add_argument("--aln", choices=sorted(ALN), default="74",
                    help="alignment: 74-way (default) or avian363")
    ap.add_argument("--phylogram", action="store_true",
                    help="scale the tree by branch length (default: cladogram). The 74-way tree is "
                         "dated (relative time, ultrametric); the avian363 tree is substitution-scaled, "
                         "so its tips do not align and dotted leaders run to the labels")
    ap.add_argument("--paper", type=float, default=None, metavar="WIDTH_MM",
                    help="figure-ready version: exact width in mm, all text Arial 5-6 pt, "
                         "svg with live text (written next to the screen version, suffix _paper)")
    ap.add_argument("--tree-scale", type=float, default=None,
                    help="tree panel width relative to the default (default 1; 0.5 with --paper)")
    ap.add_argument("--row-mm", type=float, default=2.3, help="row height in mm for --paper")
    ap.add_argument("--as-listed", action="store_true",
                    help="show the subsample's listed species instead of each clade's most typical member")
    a = ap.parse_args()

    cfg = ALN[a.aln]
    st, maf = find_site(a.site, a.matrix, cfg)
    p, mlab = load_matrix(a.matrix)
    lo = np.log2((p + PC) / (1 + 4 * PC) / .25)
    rc = lo[::-1, ::-1]
    TH = thresh(lo)

    tree = parse_newick(open(cfg["tree"]).read())
    OSC = {l.strip() for l in open(cfg["osc"]) if l.strip()}
    INTERNAL = set()
    def _int(n):
        if n["ch"]: INTERNAL.add(n["name"])
        for c in n["ch"]: _int(c)
    _int(tree)
    # the oscine MRCA: the smallest clade holding every oscine in the alignment
    def _mrca(n):
        L = {x["name"] for x in leaves(n)}
        if not OSC <= L: return None
        for c in n["ch"]:
            r = _mrca(c)
            if r: return r
        return n["name"]
    STEM = _mrca(tree)

    # presence at every node: the same reference projection the statistics use
    BLK = PR.load(maf)
    pw0, pw1 = st["s"] - PAD, st["e"] + PAD
    proj = PR.project(PR.overlapping(*BLK, pw0, pw1), pw0, pw1, st["s"], st["e"])
    score = {g: (best(sq, lo, rc) if len(sq) >= len(lo) else None)
             for g, (sq, _, _) in proj.items()}
    covered = {g: PR.covered(n, st["s"], st["e"]) for g, (_, n, _) in proj.items()}
    # missing data is missing, not absent: uncovered nodes carry no presence state
    pres = {g: ((v is not None and v >= TH) if covered[g] else None)
            for g, v in score.items()}
    pseq = proj

    # gain/loss branches on the FULL tree
    par = {}
    def idx(n):
        for c in n["ch"]: par[c["name"]] = n; idx(c)
    idx(tree)
    full = {}
    def reg(n):
        full[n["name"]] = n
        for c in n["ch"]: reg(c)
    reg(tree)
    events = []
    for c, pn in par.items():
        a_, b_ = pres.get(pn["name"]), pres.get(c)
        if a_ is not None and b_ is not None and a_ != b_:
            events.append((c, "+" if b_ else MINUS))

    # subsample, pruned tree, row order
    SUB = cfg["subsample"]
    want = [g for _, gs in SUB for g in gs]
    keep = {g for g in want if g in pseq}
    rank = {g: i for i, g in enumerate(want)}
    pt = prune(tree, keep); order_tips(pt, rank)
    tips = tip_order(pt)
    clade_of = {g: c for c, gs in SUB for g in gs}

    # what each shown species stands for. A genome goes to a shown species when that
    # species is its unambiguous nearest relative (the first ancestor holding any shown
    # species holds only one). A genome branching off an internal edge of the pruned tree
    # goes to the smallest clade bracket whose shown members' MRCA contains it, and there
    # to the bracket member fewest nodes away (ties: shortest branch length). Genomes in
    # no bracket (e.g. palaeognaths) are left unrepresented and counted in the caption.
    TIPS_ = {}
    def _tp(n):
        TIPS_[n["name"]] = ({n["name"]} if not n["ch"] else set().union(*[_tp(c) for c in n["ch"]]))
        return TIPS_[n["name"]]
    _tp(tree)
    def _down(n, t, d=0, b=0.0):
        if n["name"] == t: return (d, b)
        for c in n["ch"]:
            if t in TIPS_[c["name"]]: return _down(c, t, d + 1, b + c["bl"])
    def _mrca_of(S):
        n = tree
        while True:
            c = next((c for c in n["ch"] if S <= TIPS_[c["name"]]), None)
            if c is None: return n
            n = c
    BR = {}
    for g in keep: BR.setdefault(clade_of.get(g, ""), set()).add(g)
    BR_T = {c: TIPS_[_mrca_of(S)["name"]] for c, S in BR.items()}
    def _nearest(g, S):
        A = par.get(g)
        while A is not None and not (TIPS_[A["name"]] & S): A = par.get(A["name"])
        return None if A is None else (A, TIPS_[A["name"]] & S)
    rep_of = collections.defaultdict(list); unrep = []
    for g in (x for x in full if x not in INTERNAL and x in pres):
        if g in keep: rep_of[g].append(g); continue
        A, S = _nearest(g, keep)
        if len(S) == 1:
            rep_of[next(iter(S))].append(g); continue
        inb = [c for c in BR if g in BR_T[c]]
        if not inb: unrep.append(g); continue
        c = min(inb, key=lambda c: len(BR_T[c]))
        A, S = _nearest(g, BR[c])
        rep_of[min(S, key=lambda t: _down(A, t))].append(g)

    # each row shows the most typical of the genomes it stands for (rep_of), not
    # necessarily the species listed in the subsample: the member whose sequence over the
    # displayed window has the fewest mismatches to the others (gaps and missing bases
    # count; ties keep the listed species). The pruned tree is rebuilt from the full tree
    # for the species actually shown; a swap that would split a clade bracket is refused.
    fw0, fw1 = st["s"] - a.flank, st["e"] + a.flank
    def wseq(g):
        return {r: ch.upper() for r, ch in proj[g][2] if r is not None and fw0 <= r < fw1}
    def contiguous(tp):
        seen_, last = set(), None
        for t in tp:
            c = clade_of.get(t, "")
            if c != last and c in seen_: return False
            seen_.add(c); last = c
        return True
    shown_for = {}
    if not a.as_listed:
        for t in sorted(keep):
            if t == "Lonchura_striata": continue
            C = [g for g in rep_of[t] if g in proj and covered.get(g)]
            if len(C) < 3: continue
            W = {g: wseq(g) for g in C}
            D = {g: sum(sum(1 for r in range(fw0, fw1) if W[g].get(r, " ") != W[h].get(r, " "))
                        for h in C if h != g) for g in C}
            m = min(C, key=lambda g: (D[g], g != t))
            if m == t: continue
            k2 = (keep - {t}) | {m}; rank[m] = rank[t]; clade_of[m] = clade_of[t]
            p2 = prune(tree, k2); order_tips(p2, rank); t2 = tip_order(p2)
            if not contiguous(t2): continue
            keep, pt, tips = k2, p2, t2
            rep_of[m] = rep_of.pop(t); shown_for[m] = t

    # map events onto pruned edges by the selected-tip set under the branch
    def tipset(name):
        return frozenset(x["name"] for x in leaves(full[name])) & keep
    pnodes = {}
    def collect(n):
        pnodes[n["tips"]] = n
        for c in n["ch"]: collect(c)
    collect(pt)
    edge_ev = collections.defaultdict(list); hidden = collections.Counter()
    for c, kind in events:
        T = tipset(c)
        if T and T in pnodes: edge_ev[T].append((kind, c))
        else: hidden[kind] += 1

    tips_pre = tips
    def tipset_pre(name):
        return {x["name"] for x in leaves(full[name])}
    # the two ancestors bracketing the first gain on the reference's lineage
    path = []; cur = "Lonchura_striata"
    while cur in par: path.append(cur); cur = par[cur]["name"]
    gpath = [c for c, k in events if k == "+" and c in path]
    up = []; cur = STEM
    for lab in CLADE_UP:
        up.append((cur, lab)); cur = par[cur]["name"] if cur in par else None
        if cur is None: break
    cname = dict(up)
    anc_rows = []
    def genus(t): return t.split("_")[0]
    def describe(node):
        if node in cname:
            return f"{cname[node]} ancestor"
        T = [t for t in tips_pre if t in tipset_pre(node)]
        if len(T) >= 2:
            return f"MRCA {genus(T[0])}–{genus(T[-1])}"
        return node
    if gpath:
        g = gpath[-1]                                    # deepest gain on the path
        for node in (par[g]["name"], g):
            anc_rows.append((node, describe(node)))

    # display columns -- each reference position once, from one block
    w0, w1 = st["s"] - a.flank, st["e"] + a.flank
    rows_g = tips + [n for n, _ in anc_rows]
    dbl = PR.overlapping(*BLK, w0, w1)
    bpos = [{c: i for i, c in enumerate(b_[2]) if c is not None and w0 <= c < w1} for b_ in dbl]
    def bestrow(k, g):
        rows = dbl[k][3].get(g)
        if not rows: return None
        return max(rows, key=lambda t: sum(1 for i in bpos[k].values()
                                           if i < len(t) and t[i] != "-"))
    brow = {(k, g): bestrow(k, g) for k in range(len(dbl)) for g in rows_g}
    assign = {}
    for r in range(w0, w1):
        cand = [k for k in range(len(dbl)) if r in bpos[k]]
        if not cand: continue
        assign[r] = max(cand, key=lambda k: (
            sum(1 for g in rows_g if brow[(k, g)] is not None
                and brow[(k, g)][bpos[k][r]] != "-"), -k))
    mat = []
    rs = sorted(assign)
    def colof(k, i):
        c = {}
        for g in rows_g:
            t = brow[(k, g)]
            c[g] = " " if t is None else (t[i].upper() if i < len(t) else " ")
        return c
    for n_, r in enumerate(rs):
        k = assign[r]; i = bpos[k][r]
        mat.append((r, colof(k, i)))
        if n_ + 1 < len(rs) and assign[rs[n_ + 1]] == k:
            for ii in range(i + 1, bpos[k][rs[n_ + 1]]):
                if dbl[k][2][ii] is None:
                    mat.append((None, colof(k, ii)))
    # drop insertion columns that are gaps in every row shown
    mat = [(r, c) for r, c in mat if r is not None or any(v not in "- " for v in c.values())]
    ncol = len(mat)
    mcols = [k for k, (r, _) in enumerate(mat) if r is not None and st["s"] <= r < st["e"]]
    nins = (sum(1 for k in range(mcols[0], mcols[-1] + 1) if mat[k][0] is None)
            if mcols else 0)

    # ---------------------------------------------------------- figure --
    nrow = len(tips) + (len(anc_rows) + .6 if anc_rows else 0)
    cw, rh = .118, .158                                   # inches per column / row
    tree_w, lab_w, sc_w, rp_w, cl_w = 1.5, 1.35, .62, .9, .95
    tree_w *= a.tree_scale if a.tree_scale is not None else (.5 if a.paper else 1.0)
    if a.paper:                       # paper text is large relative to the layout: room for labels
        lab_w, rp_w, cl_w = 1.6, .8, .82
    seq_w = ncol * cw
    W = tree_w + lab_w + seq_w + sc_w + rp_w + cl_w + .25
    logo_h, top_h, bot_h = .62, .55, .9
    H = logo_h + top_h + nrow * rh + bot_h
    fig = plt.figure(figsize=(W, H))
    L = lambda x: x / W; B = lambda y: y / H
    ax_t = fig.add_axes([L(.1), B(bot_h), L(tree_w), B(nrow * rh)])
    ax_s = fig.add_axes([L(.1 + tree_w + lab_w), B(bot_h), L(seq_w), B(nrow * rh)])
    ax_l = fig.add_axes([L(.1 + tree_w + lab_w), B(bot_h + nrow * rh + .05), L(seq_w), B(logo_h - .05)])
    ax_r = fig.add_axes([L(.1 + tree_w + lab_w + seq_w + .08), B(bot_h), L(sc_w + rp_w + cl_w), B(nrow * rh)])
    for ax in (ax_t, ax_s, ax_r):
        ax.set_ylim(nrow - .5, -.5)
    for ax in (ax_t, ax_s, ax_l, ax_r):
        ax.set_axis_off()

    ypos = {g: i for i, g in enumerate(tips)}
    for j, (n, _) in enumerate(anc_rows):
        ypos[n] = len(tips) + .6 + j

    # tree -- cladogram (depth from the root) or phylogram (summed branch length; a
    # pruned edge carries the length of every full-tree branch it stands for)
    depth = {}
    def dep(n, d):
        depth[n["tips"]] = d
        for c in n["ch"]: dep(c, d + 1)
    dep(pt, 0)
    maxd = max(depth.values())
    droot = {}
    def _dr(n, x):
        droot[n["name"]] = x
        for c in n["ch"]: _dr(c, x + c["bl"])
    _dr(tree, 0.0)
    d0 = droot[pt["full"]]; dmax = max(droot[t] for t in tips) - d0
    def xy(n):
        if a.phylogram:
            x = (droot[n["full"]] - d0) / dmax
        else:
            x = 1.0 if not n["ch"] else depth[n["tips"]] / maxd
        y = ypos[n["name"]] if not n["ch"] else np.mean([xy(c)[1] for c in n["ch"]])
        return x, y
    node_xy = {}
    gk = 1.5 if a.paper else 1.0      # key spacing in rows: paper text is larger relative to a row
    # internal edges: the reconstructed ancestor's state (red present, grey absent,
    # dotted no sequence). Tip edges: the share of the genomes that row stands for
    # carrying the site; dotted when none of them has sequence
    def frac(n):
        G = [g for t in n["tips"] for g in rep_of.get(t, [t]) if pres.get(g) is not None]
        return (sum(pres[g] for g in G) / len(G)) if G else None
    def sty(f):
        if f is None: return dict(color="#D3D8DF", lw=1.0, ls=(0, (1.2, 1.4)))
        if isinstance(f, bool): return dict(color=PRES if f else ABS, lw=1.3)
        return dict(color=ramp(f), lw=1.3)
    def state(n):
        return frac(n) if not n["ch"] else pres.get(n["full"])
    def draw(n):
        x, y = xy(n); node_xy[n["full"]] = (x, y)
        for c in n["ch"]:
            cx, cy = xy(c)
            ax_t.plot([x, x], [y, cy], solid_capstyle="butt", zorder=2, **sty(state(n)))
            ax_t.plot([x, cx], [cy, cy], solid_capstyle="butt", zorder=2, **sty(state(c)))
            draw(c)
    draw(pt)
    ax_t.set_xlim(-.04, 1.02)
    if a.phylogram:
        for t in tips:                       # leaders from short tips to their labels
            tx, ty = xy(pnodes[frozenset([t])])
            if tx < .995: ax_t.plot([tx, 1.0], [ty, ty], color="#C9CED6", lw=.6, ls=(0, (1, 1.6)), zorder=1)
        # scale bar: a round length near a fifth of the tree's depth
        raw = dmax / 5; mag = 10 ** np.floor(np.log10(raw)); val = min((1, 2, 5, 10), key=lambda m: abs(m * mag - raw)) * mag
        yb = len(tips) - .5 + (len(anc_rows) + .6 if anc_rows else 0) + .2 * gk
        ax_t.plot([0, val / dmax], [yb, yb], color=INK2, lw=1.1, solid_capstyle="butt", clip_on=False)
        ax_t.text(0, yb + .35 * gk, f"{val:g} {cfg['unit']}", ha="left", va="top",
                  fontsize=6.2, color=INK3, clip_on=False)
    # key: edge colour ramp
    yf = len(tips) - .5 + (len(anc_rows) + .6 if anc_rows else 0) + (1.35 * gk if a.phylogram else .55 * gk)
    ax_t.text(0, yf, "internal edges: reconstructed ancestor\ntip edges: share of genomes with the site",
              ha="left", va="top", fontsize=5.6, color=INK3, clip_on=False)
    yk = yf + 1.45 * gk; ts = tree_w / 1.5; k0, kw = .07 / ts, .42 / ts      # key keeps its physical size
    ax_t.add_patch(Rectangle((0, yk - .16), .045 / ts, .32, fc=ABS, ec="none", clip_on=False))
    for i in range(60):
        ax_t.add_patch(Rectangle((k0 + kw * i / 60, yk - .16), kw / 60 * 1.05, .32, lw=0,
                                 fc=ramp((i + 1) / 60), ec="none", clip_on=False))
    for xv, tl in ((.0225 / ts, "0"), (k0 + kw / 2, "50%"), (k0 + kw, "100%")):
        ax_t.text(xv, yk + .4 * gk, tl, ha="center", va="top", fontsize=5.4, color=INK3, clip_on=False)
    # mark the reconstructed ancestors shown as rows
    for j, (n, lab) in enumerate(anc_rows):
        tset = tipset(n)
        pn = pnodes.get(tset)
        if pn is None: continue
        x, y = xy(pn)
        # beside the node, below its incoming edge, so it never covers an event glyph
        ax_t.plot([x], [y], "o", ms=3.4, color=INK2, zorder=5)
        ax_t.text(x + .022, y + .42, "ab"[j], ha="left", va="center", fontsize=6.4,
                  color=INK2, fontweight="bold", zorder=5)

    # species labels
    lx = L(.1 + tree_w + .06)
    for g in tips:
        ax_s.text(-.012 * ncol - .6, ypos[g], g.replace("_", " "), ha="right", va="center",
                  fontsize=7, style="italic", color=OSCC if g in OSC else INK,
                  fontweight="bold" if g == "Lonchura_striata" else "normal",
                  clip_on=False)
    for j, (n, lab) in enumerate(anc_rows):
        ax_s.text(-.012 * ncol - .6, ypos[n], f"{'ab'[j]}  {lab}", ha="right", va="center",
                  fontsize=6.6, color=INK2, clip_on=False)

    # sequence
    ax_s.set_xlim(-.5, ncol - .5)
    span = (mcols[0], mcols[-1]) if mcols else None
    for k, (r, _) in enumerate(mat):
        inside = span and span[0] <= k <= span[1]
        if r is None:
            ax_s.add_patch(Rectangle((k - .5, -.5), 1, nrow,
                                     fc="#DCE3EE" if inside else "#EEF1F6", ec="none", zorder=0))
        elif k in mcols:
            ax_s.add_patch(Rectangle((k - .5, -.5), 1, nrow, fc="#FBE9E9", ec="none", zorder=0))
            ax_l.add_patch(Rectangle((k - .5, 0), 1, 1.02, fc="#FBE9E9", ec="none", zorder=0))
    if span:
        ax_s.add_patch(Rectangle((span[0] - .5, -.5), span[1] - span[0] + 1, nrow,
                                 fc="none", ec="#D9A3A3", lw=.7, zorder=1))
    if anc_rows:
        ax_s.axhline(len(tips) - .5 + .3, color=RULE, lw=.8)
        ax_t.axhline(len(tips) - .5 + .3, color=RULE, lw=.8, xmin=.0, xmax=1)
    ref = {k: c["Lonchura_striata"] for k, (_, c) in enumerate(mat)}
    for g in rows_g:
        y = ypos[g]
        for k, (_, c) in enumerate(mat):
            ch = c[g]
            if g == "Lonchura_striata":
                ax_s.text(k, y, ch, ha="center", va="center", fontsize=7.6,
                          fontproperties=MONO, color=BASE.get(ch, INK3), fontweight="bold")
                continue
            if ch == " ":
                continue
            if ch == ref[k] and ch != "-":
                ax_s.text(k, y, ".", ha="center", va="center", fontsize=7.6,
                          fontproperties=MONO, color="#A9B1BC")
            elif ch == "-":
                # a deletion inside the motif breaks the site just as a substitution
                # does, so it gets the same weight on the page
                if k in mcols:
                    ax_s.add_patch(Rectangle((k - .48, y - .44), .96, .88,
                                             fc="#5A6675", alpha=.22, ec="none", zorder=1))
                    ax_s.text(k, y, "\u2013", ha="center", va="center", fontsize=8.4,
                              color=INK, fontweight="bold")
                else:
                    ax_s.text(k, y, "-", ha="center", va="center", fontsize=7.6,
                              fontproperties=MONO, color="#8C96A3")
            else:
                if k in mcols:
                    ax_s.add_patch(Rectangle((k - .48, y - .44), .96, .88,
                                             fc=BASE.get(ch, INK3), alpha=.16, ec="none",
                                             zorder=1))
                ax_s.text(k, y, ch, ha="center", va="center", fontsize=7.6,
                          fontproperties=MONO, color=BASE.get(ch, INK3), fontweight="bold")

    # logo over the motif columns
    ax_l.set_xlim(-.5, ncol - .5); ax_l.set_ylim(0, 1.05)
    if mcols and len(mcols) == len(p):
        ic = logo(ax_l, p, mcols, 0, 1.0, st["strand"])
        pass
    ax_l.text(-.5, .5, f"{mlab}\n({st['strand']} strand)", ha="right", va="center",
              fontsize=6.8, color=INK2, transform=ax_l.transData, clip_on=False)

    # coordinate ticks under the sequence
    for k in (mcols[0], mcols[-1]) if mcols else ():
        r = mat[k][0]
        ax_s.text(k, nrow - .5 + .55, f"{r + 1:,}", ha="center", va="top", fontsize=6.2,
                  color=INK3, rotation=0, clip_on=False)
    ax_s.text((ncol - 1) / 2, nrow - .5 + 1.35, f"{st['chrom']}  (lonStrDom2)", ha="center",
              va="top", fontsize=6.4, color=INK3, clip_on=False)

    # rows with no aligned base across most of the motif: their "absence" is missing data
    nocov = {g for g in rows_g if not covered.get(g, False)}
    # score column + clade brackets
    ax_r.set_xlim(0, sc_w + rp_w + cl_w)                # inches
    ax_r.text(.03, -1.05, "score", fontsize=6.4, color=INK3, ha="left", va="center")
    b0, b1 = .64, 1.04                                    # genomes-represented bar
    ax_r.text(b0, -1.05, "genomes with site", fontsize=6.4, color=INK3, ha="left", va="center")
    for g in tips:
        G = rep_of.get(g, [g]); y = ypos[g]
        k = sum(1 for x in G if pres.get(x)); n = sum(1 for x in G if pres.get(x) is not None)
        m = len(G) - n; x = b0
        for w, fc, ec in ((k, PRES, PRES), (n - k, "white", ABS), (m, "#EEF0F3", "#EEF0F3")):
            if not w: continue
            dx = (b1 - b0) * w / len(G)
            ax_r.add_patch(Rectangle((x, y - .27), dx, .54, fc=fc, ec=ec, lw=.5, zorder=3)); x += dx
        ax_r.text(b1 + .05, y, f"{k}/{n}", ha="left", va="center", fontsize=6.2,
                  color=INK2 if k else INK3, fontproperties=FontProperties(family="DejaVu Sans Mono"))
    for g in rows_g:
        v = score.get(g); y = ypos[g]
        on = pres.get(g)
        if g in nocov:
            ax_r.text(.094, y, "\u00d7", ha="center", va="center", fontsize=7.5, color=INK3)
            ax_r.text(.204, y, "no seq", ha="left", va="center", fontsize=6.2, color=INK3,
                      style="italic")
            continue
        ax_r.plot([.094], [y], "o", ms=4.6, color=PRES if on else "white",
                  mec=PRES if on else ABS, mew=.9)
        ax_r.text(.204, y, "—" if v is None else f"{v:.1f}", ha="left", va="center",
                  fontsize=6.4, color=INK2 if on else INK3,
                  fontproperties=FontProperties(family="DejaVu Sans Mono"))
    seen = collections.OrderedDict()
    for g in tips:
        seen.setdefault(clade_of.get(g, ""), []).append(ypos[g])
    for c, ys in seen.items():
        x0 = sc_w + rp_w + .02
        ax_r.plot([x0, x0], [min(ys) - .32, max(ys) + .32], color=INK3, lw=.8)
        ax_r.text(x0 + .08, np.mean(ys), c, ha="left", va="center", fontsize=6.6,
                  color=OSCC if c == "oscines" else INK2,
                  fontweight="bold" if c == "oscines" else "normal")

    # title + caption
    oc = [g for g in OSC if pres.get(g) is not None]
    no, nO = sum(pres[g] for g in oc), len(oc)
    out_ = [g for g in pres if g not in INTERNAL and g not in OSC
            and pres[g] is not None]
    nu = sum(pres[g] for g in out_)
    nmiss = sum(1 for g in pres if g not in INTERNAL and pres[g] is None)
    fig.text(L(.12), 1 - B(.16), f"{a.site.split('_', 1)[1] if a.site[:3] in ('RA_', 'G4_') else a.site}"
             f"   ·   {mlab}   ·   {st['cls']} DAR",
             fontsize=9.5, fontweight="bold", color=INK, va="top")
    if a.paper:                                   # short subtitle; the legend carries the rest
        fig.text(L(.12), 1 - B(.38),
                 f"{st['chrom']}:{st['s'] + 1:,}-{st['e']:,}   site present in {no}/{nO} oscines, "
                 f"{nu}/{len(out_)} non-oscines ({cfg['label']})",
                 fontsize=6.8, color=INK3, va="top")
    else:
      fig.text(L(.12), 1 - B(.38),
             f"{st['chrom']}:{st['s'] + 1:,}-{st['e']:,}   site present in {no}/{nO} oscines, "
             f"{nu}/{len(out_)} non-oscines with sequence ({cfg['label']}; {nmiss} without)   "
             f"threshold {TH:.2f}" +
             (f"   \u00b7   site spans {nins} insertion column{'s' if nins > 1 else ''}"
              if nins else "") +
             (f"   ·   {len(unrep)} genomes outside every clade shown" if unrep else "") +
             (f"   ·   each row: its clade's most typical member (\u00b1{a.flank} bp)" if not a.as_listed else ""),
             fontsize=6.8, color=INK3, va="top")

    if a.paper:
        # rescale the finished figure: width to the column, rows to --row-mm; text is
        # remapped to Arial 5-6 pt by role (title 6 bold, sequence and labels 5.5, rest 5)
        from matplotlib.lines import Line2D
        from matplotlib.text import Text
        from matplotlib.patches import Patch
        fx = a.paper / 25.4 / W; fy = a.row_mm / 25.4 / rh
        fig.set_size_inches(W * fx, H * fy)
        for t in fig.findobj(Text):
            if not t.get_text(): continue
            z = t.get_fontsize()
            nz = 6.0 if z >= 9 else 5.5 if z >= 7 else 5.0
            t.set_fontproperties(FontProperties(family="Arial", size=nz, weight=t.get_fontweight(),
                                                style=t.get_fontstyle()))
            bb = t.get_bbox_patch()
            if bb is not None: bb.set_linewidth(bb.get_linewidth() * fx)
        for l in fig.findobj(Line2D):
            l.set_linewidth(l.get_linewidth() * fx); l.set_markersize(l.get_markersize() * fx)
            l.set_markeredgewidth(l.get_markeredgewidth() * fx)
        for p_ in fig.findobj(Patch):
            if p_ is not fig.patch: p_.set_linewidth(p_.get_linewidth() * fx)
        matplotlib.rcParams["svg.fonttype"] = "none"

    stem = a.out or f"{cfg['prefix']}{'p' if a.phylogram else ''}_{a.matrix}_{a.site}" + ("_paper" if a.paper else "")
    coll = cfg["collection"] + ("_phylogram" if a.phylogram else "")
    for ext in (("pdf", "png", "svg") if a.paper else ("pdf", "png")):
        fig.savefig(f"{per_site(coll, stem)}/{stem}.{ext}", dpi=600 if a.paper else 300)
    if a.paper:
        print(f"  paper size {W * fx * 25.4:.1f} x {H * fy * 25.4:.1f} mm")
    print(f"  wrote {per_site(coll, stem)}/{stem}.pdf / .png   {ncol} columns x {len(rows_g)} rows"
          f"   events: {len(events)} ({sum(1 for _, k in events if k == '+')} gains)"
          f"   hidden: {dict(hidden)}")


if __name__ == "__main__":
    main()
