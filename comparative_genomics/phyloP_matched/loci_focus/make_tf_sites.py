#!/usr/bin/env python
"""Merged motif sites for the 13 SCENIC+-nominated MAFB regulators, MAFB locus.

Same unit as the AR site BEDs: overlapping FIMO calls within one (motif,
element) collapse to one site, and the feature is that group's BEST-SCORING
hit, exactly W bp, so strand and score belong to the interval drawn.

Writes tf_sites_MAFB_lsd2.bed (BED6+) over the 285 DARs in the MAFB +/-1.5 Mb
window, and tfsc_site_windows.bed, the merged +/-30 bp scan windows that the
alignment extraction needs.

Matrices: tf_scenic.meme -- one per regulator, AR included, chosen from the motif
evidence SCENIC+ used for that regulator's cistrome (build_scenic_motifs.py).
NR3C2 shares AR's matrix (both backed by metacluster_185.1) and is skipped, so
each site is called once. The superseded tf13.meme set is in superseded/.
"""
import collections, re

F = "/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus"
A = "/private/warm-archive/colquittlab/assembly/cactus/avian363"
PAD = 30
RGB = {"Glut-RA": "227,26,28", "GABA-4": "44,110,127"}
CL = {"RA_over_Arco": "Glut-RA", "GABA4_over_GABA1": "GABA-4"}

W, n = {}, None
for l in open(f"{F}/tf_scenic.meme"):
    if l.startswith("MOTIF"):
        n = l.split()[1]
    elif l.startswith("letter-probability"):
        W[n] = int(re.search(r"w=\s*(\d+)", l).group(1))

el = {}
for l in open(f"{F}/allelem_cls.bed"):
    x = l.rstrip("\n").split("\t")
    el[x[3]] = (x[0], int(x[1]), int(x[2]), x[4])
start = {}
for l in open(f"{F}/am_dars.bed"):
    x = l.rstrip("\n").split("\t")
    start[x[3]] = int(x[1])
mafb = {l.strip() for l in open(f"{A}/vgp577/lifted/mafb_1mb.txt")}

per = collections.defaultdict(list)
for i, l in enumerate(open(f"{F}/tfsc_hits.txt")):
    if not i or not l.strip():
        continue
    p = l.rstrip("\n").split("\t")
    if len(p) < 7:
        continue
    m, e = p[0], p[2]
    if m in ("NR3C2",): continue          # same matrix as AR in tf_scenic.meme
    if e not in mafb or e not in el or e not in start:
        continue
    s, t = int(p[3]), int(p[4])
    if t - s + 1 != W[m]:
        continue
    per[(m, e)].append((start[e] + s - 1, start[e] + t, p[5], float(p[6])))

rows = []
for (m, e), v in sorted(per.items()):
    ch, _, _, cl = el[e]
    c = CL[cl]
    v.sort()
    groups = []
    for h in v:
        if groups and h[0] < groups[-1][-1][1]:
            groups[-1].append(h)
        else:
            groups.append([h])
    for j, g in enumerate(groups):
        a, b, strand, sc = max(g, key=lambda x: x[3])
        rows.append([ch, a, b, f"{c.split('-')[0]}_{e}_{j + 1}", int(min(1000, sc * 50)),
                     strand, a, b, RGB[c], m, c, e, round(sc, 3)])

rows.sort(key=lambda r: (r[0], r[1]))
with open(f"{F}/tracks/mafb_1p5mb/tf_sites_MAFB_lsd2.bed", "w") as o:
    for r in rows:
        o.write("\t".join(str(x) for x in r) + "\n")

iv = sorted((r[0], max(0, r[1] - PAD), r[2] + PAD) for r in rows)
mg = [list(iv[0])]
for c, a, b in iv[1:]:
    if c == mg[-1][0] and a <= mg[-1][2]:
        mg[-1][2] = max(mg[-1][2], b)
    else:
        mg.append([c, a, b])
with open(f"{F}/maf_switch/data/tfsc_site_windows.bed", "w") as o:
    for c, a, b in mg:
        o.write(f"{c}\t{a}\t{b}\n")

n = collections.Counter(r[9] for r in rows)
print(f"  sites: {len(rows):,} over {len({r[11] for r in rows})} DARs")
for m in sorted(n):
    print(f"    {m:<8} {n[m]:>4}   width {W[m]}")
print(f"  scan windows: {len(mg):,}   bases {sum(b - a for _, a, b in mg):,}")
