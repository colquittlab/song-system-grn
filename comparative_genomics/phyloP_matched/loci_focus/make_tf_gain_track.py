#!/usr/bin/env python
"""Candidate motif-site gains for the 13 other predicted MAFB regulators.

The AR pipeline (make_gain_track.py) run on the SCENIC+-nominated inputs
instead: same +/-30 bp presence call, same gain statistic, same thresholds,
same default span -- the gene-bounded DHX35-TOP1 domain (span.py), with the
MAFB +/-1.5 Mb window reported alongside it and tagged "_1p5mb".

  tf_gain_all_lsd2.bed               domain, every site with coverage
  tf_gain_candidates_lsd2.bed        domain, gain >= 0.50 and qDomain < 0.05
  tf_gain_1p5mb_all_lsd2.bed         +/-1.5 Mb, every site with coverage
  tf_gain_1p5mb_candidates_lsd2.bed  +/-1.5 Mb, gain >= 0.50, qWindow < 0.05

FDR is computed over this 13-motif set, separately from AR's, because the two
were scanned as separate families; the combined correction is reported in
tf_gain_compare.py instead.

Presence is the corrected call (maf_switch/presence.py); tiers as in
make_gain_track.py.

Metacluster caveat: KLF12, STAT1, SOX6, TCF12, FOSL2 and MAF are cistarget
metaclusters covering many factors -- the TCF12 one includes Ar. A site called
under one of those names is a site for that matrix, not evidence about that TF.
"""
import collections, json, subprocess
import numpy as np
from scipy import stats

import sys
sys.path.insert(0, "/private/warm-archive/colquittlab/assembly/cactus/avian363"
                   "/phyloP_matched/loci_focus")
from span import F, DS, DE, elem_coords, in_domain, LABEL, LABEL_WIDE

M = f"{F}/maf_switch"
BB = "/private/groups/colquittlab/opt/ucsc_utilities/bedToBigBed"
RGB = {"Glut-RA": "227,26,28", "GABA-4": "44,110,127"}
GAIN_MIN, Q_MAX = 0.50, 0.05
META = {"MAF", "STAT1", "LIN28B"}   # SCENIC+ motif evidence is another family's motif

el = elem_coords()

bed = {}
for l in open(f"{F}/tracks/mafb_1p5mb/tf_sites_MAFB_lsd2.bed"):
    x = l.rstrip("\n").split("\t")
    bed[(x[9], x[3])] = dict(chrom=x[0], s=int(x[1]), e=int(x[2]), strand=x[5],
                             fimo=float(x[12]))

sites = json.load(open(f"{M}/data/tf_site_gain.json"))
try:
    counts = {(c["motif"], c["name"]): c
              for c in json.load(open(f"{M}/trees363/counts_tf.json"))}
except FileNotFoundError:
    counts = {}

p = np.array([stats.fisher_exact([[g["p_osc"], g["n_osc"] - g["p_osc"]],
                                  [g["p_out"], g["n_out"] - g["p_out"]]])[1] for g in sites])


def bh(pv):
    pv = np.asarray(pv, float)
    o = np.argsort(pv)
    out = np.empty(len(pv))
    out[o] = np.minimum(np.minimum.accumulate(
        (pv[o] * len(pv) / np.arange(1, len(pv) + 1))[::-1])[::-1], 1)
    return out


dom = np.array([in_domain(el[g["elem"]]) for g in sites])
qw = bh(p)
qd = np.full(len(p), np.nan)
qd[dom] = bh(p[dom])


def tier(c):
    return c["tier"] if c else "no-coverage"


def location(c):
    """Principal oscine gain on the avian363 tree, e.g. near-stem:3;osc:94%;nonosc:0;regain:2.
    depth 0 = the stem; 1-4 = the deep backbone nodes holding >= 90% of oscines."""
    if not c: return "none"
    if c["principal"] is None: return f"{c['tier']};nonosc:{c['nonosc_gains']}"
    return (f"{c['tier']}:{c['depth']}" if c["depth"] is not None else c["tier"]) + \
           f";osc:{100 * c['osc_cov']:.0f}%;nonosc:{c['nonosc_gains']};regain:{c['other_osc_gains']}"


rows = []
for g, pv, qdv, qwv, ind in zip(sites, p, qd, qw, dom):
    b = bed[(g["motif"], g["name"])]
    c = counts.get((g["motif"], g["name"]))
    gain = g["f_osc"] - g["f_out"]
    rows.append([b["chrom"], g["s"], g["e"], f'{g["name"]}_{g["motif"]}',
                 max(0, min(1000, int(round(1000 * gain)))), b["strand"],
                 g["s"], g["e"], RGB[g["cls"]],
                 g["motif"], g["cls"], g["elem"],
                 g["p_osc"], g["n_osc"], g["p_out"], g["n_out"],
                 round(g["f_osc"], 4), round(g["f_out"], 4), round(gain, 4),
                 f"{pv:.3g}", ("NA" if np.isnan(qdv) else f"{qdv:.3g}"),
                 f"{qwv:.3g}", int(ind),
                 (c["gains"] if c else -1), (c["losses"] if c else -1), tier(c),
                 round(b["fimo"], 3), int(g["motif"] in META), location(c)])

AS = f'''table tfGain
"Candidate motif site gains for the predicted MAFB regulators (lonStrDom2)"
(
string  chrom;          "Reference sequence chromosome"
uint    chromStart;     "Start position"
uint    chromEnd;       "End position"
string  name;           "Site name and motif"
uint    score;          "1000 x gain (oscine minus outgroup presence)"
char[1] strand;         "Strand of the best-scoring hit"
uint    thickStart;     "Start of thick drawing"
uint    thickEnd;       "End of thick drawing"
uint    reserved;       "Item RGB colour, by DAR class"
string  motif;          "Matrix used to call the site"
string  darClass;       "DAR class: Glut-RA or GABA-4"
string  element;        "Parent DAR"
uint    oscineWith;     "Oscine genomes carrying the site"
uint    oscineN;        "Oscine genomes with alignment coverage"
uint    outWith;        "Non-oscine genomes carrying the site"
uint    outN;           "Non-oscine genomes with alignment coverage"
float   fOscine;        "Fraction of oscines with the site"
float   fOutgroup;      "Fraction of the outgroup with the site"
float   gain;           "fOscine minus fOutgroup"
string  fisherP;        "Fisher exact p, oscine vs outgroup presence"
string  qDomain;        "BH q across the DHX35-TOP1 domain, the default span (NA outside)"
string  qWindow;        "BH q across every site in the MAFB +/-1.5 Mb window"
uint    inDomain;       "1 if the parent DAR overlaps the default span"
int     treeGains;      "Gains on the avian363 tree (-1 = too little coverage)"
int     treeLosses;     "Losses on the avian363 tree (-1 = too little coverage)"
string  tier;           "avian363 placement of the principal oscine gain: stem / near-stem (1-4 backbone nodes, >=90% of oscines) / within / inherited / none / no-coverage"
float   fimoScore;      "FIMO log-odds score of the best hit"
uint    otherFamily;    "1 if SCENIC+'s motif evidence for this TF is another family's motif (MAF, STAT1: AP-1; LIN28B)"
string  gainLocation;   "avian363 principal gain: tier:depth; share of oscines below it; independent non-oscine gains; extra gains inside oscines"
)
'''
open(f"{F}/tracks/tf_gain.as", "w").write(AS)


def cand(r, qcol):
    return r[18] >= GAIN_MIN and r[qcol] != "NA" and float(r[qcol]) < Q_MAX


for stem, keep in (("tracks/mafb_domain/tf_gain_all_lsd2", lambda r: r[22] == 1),
                   ("tracks/mafb_domain/tf_gain_candidates_lsd2", lambda r: r[22] == 1 and cand(r, 20)),
                   ("tracks/mafb_1p5mb/tf_gain_1p5mb_all_lsd2", lambda r: True),
                   ("tracks/mafb_1p5mb/tf_gain_1p5mb_candidates_lsd2", lambda r: cand(r, 21))):
    sel = sorted((r for r in rows if keep(r)), key=lambda r: (r[0], r[1]))
    with open(f"{F}/{stem}.bed", "w") as fh:
        for r in sel:
            fh.write("\t".join(str(x) for x in r) + "\n")
    subprocess.run([BB, f"-type=bed9+{len(rows[0]) - 9}", f"-as={F}/tracks/tf_gain.as", "-tab",
                    f"{F}/{stem}.bed", f"{F}/ar_gen.sizes", f"{F}/{stem}.bb"], check=True)
    n = collections.Counter(r[9] for r in sel)
    print(f"  {stem}.bed  {len(sel)} sites   " +
          " ".join(f"{k}:{v}" for k, v in sorted(n.items())))

print()
print(f"  DEFAULT SPAN -- {LABEL}  chr20:{DS:,}-{DE:,}")
for r in sorted((r for r in rows if r[22] == 1 and cand(r, 20)), key=lambda r: -r[18]):
    print(f"    {r[0]}:{r[1] + 1:,}-{r[2]:,} {r[5]}  {r[3]:<42} {r[10]:<8}"
          f" osc {r[12]:>3}/{r[13]:<3} out {r[14]:>3}/{r[15]:<3}"
          f" gain {r[18]:+.2f}  qDom {r[20]:>9}  {r[25]}")
print()
print(f"  expanded span -- {LABEL_WIDE}: "
      f"{sum(1 for r in rows if cand(r, 21))} candidates of {len(rows)} sites")
