#!/usr/bin/env python
"""Trackhub BEDs for candidate AR-site GAINS at the MAFB locus (lonStrDom2).

Each feature is one merged AR site whose presence/absence was called in every
species of the avian363 alignment by degapping that species' own sequence and
rescanning it with the PWM (maf_switch/site_gain.py). A species "has" the site
if its best hit clears the 99.99th percentile of the matrix's uniform null.

gain = fraction of oscines with the site - fraction of the outgroup.

SPAN. The default is the gene-bounded DHX35-TOP1 domain (see span.py); the
MAFB +/-1.5 Mb window is reported alongside it, tagged "_1p5mb".

  ar_gain_all_lsd2.bed               domain, every site with coverage
  ar_gain_candidates_lsd2.bed        domain, gain >= 0.50 and qDomain < 0.05
  ar_gain_1p5mb_all_lsd2.bed         +/-1.5 Mb, every site with coverage
  ar_gain_1p5mb_candidates_lsd2.bed  +/-1.5 Mb, gain >= 0.50, qWindow < 0.05

Two FDRs, because the set tested is part of the claim: qDomain is BH across
the domain sites alone and selects the default candidates; qWindow is BH
across every site in the +/-1.5 Mb window. Both ride on a Fisher test that
treats genomes as independent, which they are not -- they rank sites, they do
not measure evidence. The tier column, from the 74-way branch reconstruction,
is the phylogeny-aware version of the question.

Presence is the corrected call (maf_switch/presence.py): one reference-
projected sequence per genome, and a genome with bases at fewer than 80% of
the motif's positions is missing, not absent.

tier, from the 74-way branch reconstruction:
  oscine-stem     exactly one gain, on the branch into the oscine MRCA -- the
                  oscine-specific gains, the ones to test
  single-other    exactly one gain, elsewhere (above the oscine stem, within
                  oscines, or in a non-oscine lineage); gainLocation says where
  two-gain        two independent gains
  turnover        more than two
  no-gain         none observed
  no-coverage     too few 74-way genomes with sequence to reconstruct

Colour is DAR class, to match the arefull/jasparAR site tracks; score is
1000 x gain, so browser shading tracks the asymmetry.

usage: make_gain_track.py
"""
import collections, json, subprocess, sys
import numpy as np
from scipy import stats

sys.path.insert(0, "/private/warm-archive/colquittlab/assembly/cactus/avian363"
                   "/phyloP_matched/loci_focus")
from span import F, DS, DE, elem_coords, in_domain, LABEL, LABEL_WIDE

M = f"{F}/maf_switch"
BB = "/private/groups/colquittlab/opt/ucsc_utilities/bedToBigBed"
RGB = {"Glut-RA": "227,26,28", "GABA-4": "44,110,127"}
MN = {"full": "ARE-full", "jaspar": "JASPAR_AR"}
GAIN_MIN, Q_MAX = 0.50, 0.05

el = elem_coords()

bed = {}
for stem, mo in (("tracks/mafb_1p5mb/arefull_sites_MAFB_lsd2.bed", "full"),
                 ("tracks/mafb_1p5mb/jasparAR_sites_MAFB_lsd2.bed", "jaspar")):
    for l in open(f"{F}/{stem}"):
        x = l.rstrip("\n").split("\t")
        bed[(mo, x[3])] = dict(chrom=x[0], s=int(x[1]), e=int(x[2]), strand=x[5],
                               fimo=float(x[11]), ph_osc=float(x[12]),
                               ph_all=float(x[13]))

sites = json.load(open(f"{M}/data/site_gain.json"))
counts = {(c["motif"], c["name"]): c
          for c in json.load(open(f"{M}/trees363/counts_ar.json"))}

p = np.array([stats.fisher_exact([[g["p_osc"], g["n_osc"] - g["p_osc"]],
                                  [g["p_out"], g["n_out"] - g["p_out"]]])[1]
              for g in sites])


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
    assert (b["s"], b["e"]) == (g["s"], g["e"]), g["name"]
    c = counts.get((g["motif"], g["name"]))
    gain = g["f_osc"] - g["f_out"]
    rows.append([b["chrom"], g["s"], g["e"], f'{g["name"]}_{MN[g["motif"]]}',
                 max(0, min(1000, int(round(1000 * gain)))), b["strand"],
                 g["s"], g["e"], RGB[g["cls"]],
                 MN[g["motif"]], g["cls"], g["elem"],
                 g["p_osc"], g["n_osc"], g["p_out"], g["n_out"],
                 round(g["f_osc"], 4), round(g["f_out"], 4), round(gain, 4),
                 f"{pv:.3g}", ("NA" if np.isnan(qdv) else f"{qdv:.3g}"),
                 f"{qwv:.3g}", int(ind),
                 (c["gains"] if c else -1), (c["losses"] if c else -1), tier(c),
                 round(b["fimo"], 3), round(b["ph_osc"], 4), round(b["ph_all"], 4),
                 location(c)])

AS = f'''table arGain
"Candidate androgen-receptor motif site gains at the MAFB locus (lonStrDom2)"
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
string  motif;          "PWM used to call the site"
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
uint    inDomain;       "1 if the parent DAR overlaps chr20:{DS:,}-{DE:,}"
int     treeGains;      "Gains on the avian363 tree (-1 = too little coverage)"
int     treeLosses;     "Losses on the avian363 tree (-1 = too little coverage)"
string  tier;           "avian363 placement of the principal oscine gain: stem / near-stem (1-4 backbone nodes, >=90% of oscines) / within / inherited / none / no-coverage"
float   fimoScore;      "FIMO log-odds score of the best hit"
float   phyloPOscine;   "Mean phyloP across the site, oscine subtree"
float   phyloPWholeTree; "Mean phyloP across the site, whole tree"
string  gainLocation;   "avian363 principal gain: tier:depth; share of oscines below it; independent non-oscine gains; extra gains inside oscines"
)
'''
open(f"{F}/tracks/ar_gain.as", "w").write(AS)


def cand(r, qcol):
    return r[18] >= GAIN_MIN and r[qcol] != "NA" and float(r[qcol]) < Q_MAX


for stem, keep in (("tracks/mafb_domain/ar_gain_all_lsd2", lambda r: r[22] == 1),
                   ("tracks/mafb_domain/ar_gain_candidates_lsd2", lambda r: r[22] == 1 and cand(r, 20)),
                   ("tracks/mafb_1p5mb/ar_gain_1p5mb_all_lsd2", lambda r: True),
                   ("tracks/mafb_1p5mb/ar_gain_1p5mb_candidates_lsd2", lambda r: cand(r, 21))):
    sel = sorted((r for r in rows if keep(r)), key=lambda r: (r[0], r[1]))
    with open(f"{F}/{stem}.bed", "w") as fh:
        for r in sel:
            fh.write("\t".join(str(x) for x in r) + "\n")
    subprocess.run([BB, f"-type=bed9+{len(rows[0]) - 9}", f"-as={F}/tracks/ar_gain.as", "-tab",
                    f"{F}/{stem}.bed", f"{F}/ar_gen.sizes", f"{F}/{stem}.bb"],
                   check=True)
    n = collections.Counter(r[25] for r in sel)
    c = collections.Counter(r[10] for r in sel)
    print(f"  {stem}.bed  {len(sel)} sites   Glut-RA {c['Glut-RA']}  GABA-4 {c['GABA-4']}")
    print("    " + "   ".join(f"{k} {v}" for k, v in sorted(n.items())))

print()
print(f"  DEFAULT SPAN -- {LABEL}  chr20:{DS:,}-{DE:,}")
for r in sorted((r for r in rows if r[22] == 1 and cand(r, 20)), key=lambda r: -r[18]):
    print(f"    {r[0]}:{r[1] + 1:,}-{r[2]:,} {r[5]}  {r[3]:<40} {r[10]:<8}"
          f" osc {r[12]:>3}/{r[13]:<3} out {r[14]:>3}/{r[15]:<3}"
          f" gain {r[18]:+.2f}  qDom {r[20]:>9}  {r[25]}")
print()
print(f"  expanded span -- {LABEL_WIDE}: "
      f"{sum(1 for r in rows if cand(r, 21))} candidates of {len(rows)} sites")
