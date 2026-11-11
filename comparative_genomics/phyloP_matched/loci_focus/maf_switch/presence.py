#!/usr/bin/env python
"""Reference-projected sequence per genome, for motif presence calls.

Replaces the concatenation used by site_gain.py, gain_trees.py and the tf_*
scripts, which joined every row of every overlapping block. That has two
failure modes:

  1. Where the alignment holds more than one block over the same reference
     positions (a duplication), the same sequence is joined twice and the
     junction between copies can create a motif hit that exists in no genome.
     The call then depends on block order -- one SATB1 site moved from 61/145
     oscines present to 8/145 on reordering alone.
  2. A genome with no orthologous sequence across the motif, only a few
     aligned bases elsewhere in the window, was scored as lacking the site.

Here each reference position in the window is taken ONCE per genome, from the
covering block in which that genome's best row aligns the most bases across
the window. Columns where the reference has a gap (insertions in the other
genome) are kept only between two consecutive reference positions drawn from
the same block, so an insertion inside a motif is preserved -- indels remain a
mechanism of gain -- but nothing is spliced across blocks.

A genome is COVERED at a site when at least COV of the motif's reference
positions carry a base. Uncovered genomes are missing, not absent.
"""
import bisect, collections

COV = 0.8
REF = "Lonchura_striata"


def load(maf):
    """-> (blocks, starts). Each block: (ref_start, ref_end, ref_coord_per_col,
    {genome: [aligned text, ...]})."""
    blocks = []; cur = []
    def flush():
        ref = [r for r in cur if r[0] == REF]
        if not ref: return
        coord = []; p = ref[0][1]
        for c in ref[0][2]:
            if c == "-": coord.append(None)
            else: coord.append(p); p += 1
        rows = collections.defaultdict(list)
        for g, _, t in cur:
            rows[g].append(t)
        blocks.append((ref[0][1], p, coord, rows))
    for l in open(maf):
        if l.startswith("a"):
            if cur: flush()
            cur = []
        elif l.startswith("s\t"):
            f = l.rstrip("\n").split("\t")
            cur.append((f[1].split(".")[0], int(f[2]), f[6]))
    if cur: flush()
    blocks.sort(key=lambda b: (b[0], b[1]))
    return blocks, [b[0] for b in blocks]


def overlapping(blocks, starts, w0, w1):
    i = bisect.bisect_left(starts, w1)
    out = []
    j = i - 1
    while j >= 0:
        b = blocks[j]
        if b[1] > w0: out.append(b)
        if w0 - b[0] > 5000: break        # blocks are <=1000 columns
        j -= 1
    return out[::-1]


def project(bl, w0, w1, s, e, genomes=None):
    """-> {genome: (seq, n_motif_bases, cols)} for reference window [w0, w1).

    bl     the blocks overlapping the window, from overlapping()

    seq    the genome's degapped sequence, reference-ordered, single copy
    n_motif_bases  motif reference positions [s, e) carrying a base
    cols   list of (ref_pos or None, char) for display; char '-' is a gap,
           ' ' means the genome is absent from the chosen block
    """
    if genomes is not None:
        names = set(genomes)
    else:
        names = set()
        for b in bl: names |= set(b[3])
    # per block: column index of each window reference position, in order
    bpos = []
    for b in bl:
        idx = {c: i for i, c in enumerate(b[2]) if c is not None and w0 <= c < w1}
        bpos.append(idx)
    # reference positions held by more than one block (a duplication in the
    # alignment) are the only place the block choice depends on the genome
    holders = collections.defaultdict(list)
    for k, idx in enumerate(bpos):
        for r in idx: holders[r].append(k)
    dup = any(len(v) > 1 for v in holders.values())
    out = {}
    for g in names:
        # this genome's best row in each block, and how well it aligns there
        choice = []
        for b, idx in zip(bl, bpos):
            rows = b[3].get(g)
            if not rows:
                choice.append((None, -1)); continue
            if len(rows) == 1 and not dup:
                choice.append((rows[0], 0)); continue
            best_t, best_n = None, -1
            for t in rows:
                n = sum(1 for i in idx.values() if i < len(t) and t[i] != "-")
                if n > best_n: best_t, best_n = t, n
            choice.append((best_t, best_n))
        # assign each reference position to one block
        assign = {}
        for r in range(w0, w1):
            hk = holders.get(r)
            if not hk: continue
            cand = [(choice[k][1], -k, k) for k in hk if choice[k][0] is not None]
            assign[r] = max(cand)[2] if cand else None
        cols = []; seq = []
        rs = sorted(assign)
        for n_, r in enumerate(rs):
            k = assign[r]
            if k is None:
                cols.append((r, " ")); continue
            t = choice[k][0]; i = bpos[k][r]
            ch = t[i] if i < len(t) else "-"
            cols.append((r, ch.upper()))
            if ch != "-": seq.append(ch)
            # insertion columns that follow, if the next position is from the same block
            if n_ + 1 < len(rs) and assign.get(rs[n_ + 1]) == k:
                j = bpos[k][rs[n_ + 1]]
                for ii in range(i + 1, j):
                    if bl[k][2][ii] is None:
                        c2 = t[ii] if ii < len(t) else "-"
                        cols.append((None, c2.upper()))
                        if c2 != "-": seq.append(c2)
        nmot = sum(1 for r, ch in cols if r is not None and s <= r < e and ch not in "- ")
        out[g] = ("".join(seq), nmot, cols)
    return out


def covered(nmot, s, e):
    return nmot >= COV * (e - s)
