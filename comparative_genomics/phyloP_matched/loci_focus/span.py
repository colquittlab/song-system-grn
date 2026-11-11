#!/usr/bin/env python
"""The analysis span for the MAFB locus, and membership in it.

DEFAULT SPAN -- the gene-bounded regulatory domain. MAFB sits in a large open
region flanked by two protein-coding genes, DHX35 and TOP1, with two putative
lncRNAs inside it. The domain is the intergenic interval between them:

    chr20:12,551,036-13,547,533      996 kb      (DHX35 end -> TOP1 start)

A DAR belongs to the domain if it OVERLAPS that interval, which is the rule
consolidate_ar.py and the phyloP analyses use; a site belongs if its parent
DAR does. Gene coordinates from the lonStrDom2 stringtie GTF:
DHX35 chr20:12,522,214-12,551,036, MAFB chr20:13,365,864-13,369,534,
TOP1 chr20:13,547,533-13,613,432.

EXPANDED SPAN -- MAFB +/-1.5 Mb, 285 DARs spanning chr20:11,861,215-14,892,299
(vgp577/lifted/mafb_1mb.txt). Reported alongside the default for power, never
in place of it. Outputs for this span carry a "_1p5mb" tag; unsuffixed
outputs are the default span.
"""
F = "/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus"
A = "/private/warm-archive/colquittlab/assembly/cactus/avian363"

DS, DE = 12551036, 13547533
CHROM = "chr20"
GENES = {"DHX35": (12522214, 12551036), "MAFB": (13365864, 13369534),
         "TOP1": (13547533, 13613432)}
LABEL = "DHX35–TOP1 domain"
LABEL_WIDE = "MAFB ±1.5 Mb"


def elem_coords():
    el = {}
    for l in open(f"{F}/allelem_cls.bed"):
        x = l.rstrip("\n").split("\t")
        el[x[3]] = (x[0], int(x[1]), int(x[2]))
    return el


def in_domain(c):
    """c = (chrom, start, end) of the parent DAR."""
    return c[0] == CHROM and c[2] > DS and c[1] < DE


def wide_elements():
    return {l.strip() for l in open(f"{A}/vgp577/lifted/mafb_1mb.txt")}


# ------------------------------------------------------------------ outputs --
# Every result is filed under the span it was computed on. Unsuffixed span
# folders hold results computed on one span only; compare_* folders hold
# outputs that put several spans side by side (one figure can't be split).
import os as _os, re as _re

RESULTS = f"{F}/results"
TRACKS = f"{F}/tracks"
SPAN_DIR = {None: "genome_wide", "gw": "genome_wide",
            "MAFB": "mafb_1p5mb", "1.5Mb": "mafb_1p5mb", "1p5mb": "mafb_1p5mb",
            "domain": "mafb_domain", "EMX2": "emx2_1p5mb"}


def result(name, span="gw"):
    """Path for a result computed on `span` (a SPAN_DIR key or a folder name)."""
    d = f"{RESULTS}/{SPAN_DIR.get(span, span)}"
    _os.makedirs(d, exist_ok=True)
    return f"{d}/{name}"


def track(name, span):
    """Path for a browser track computed on `span`; trackDb/.as files live in TRACKS."""
    d = f"{TRACKS}/{SPAN_DIR.get(span, span)}" if span else TRACKS
    _os.makedirs(d, exist_ok=True)
    return f"{d}/{name}"


_ELC = None


def site_span(name):
    """'mafb_domain' or 'mafb_1p5mb_outside_domain' for any site/element name
    (e.g. RA_Glut-RA.395131_1, jaspar_RA_Glut-HVC-2.394967_1_p5)."""
    global _ELC
    if _ELC is None: _ELC = elem_coords()
    m = _re.search(r"([A-Za-z][A-Za-z0-9-]*\.\d+)", name)
    if not m or m.group(1) not in _ELC:
        raise KeyError(f"no element in site name {name!r}")
    return "mafb_domain" if in_domain(_ELC[m.group(1)]) else "mafb_1p5mb_outside_domain"


def per_site(collection, name):
    """Folder for a per-site figure: results/per_site/<collection>/<site span>/."""
    d = f"{RESULTS}/per_site/{collection}/{site_span(name)}"
    _os.makedirs(d, exist_ok=True)
    return d
