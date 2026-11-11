"""Paired motif-vs-own-enhancer values for three AR motifs, three units.

MOTIF  "full"    HOMER de-novo ARE-full   (motif2, 12 bp)
       "half"    HOMER de-novo ARE-half   (motif4,  8 bp)
       "jaspar"  JASPAR MA0007.3          (16 bp)

UNIT   "instance"  one value per FIMO hit. For a palindromic motif this
                   double-counts: FIMO reports each ARE on both strands.
       "site"      one value per site, measured over the BEST-SCORING hit of
                   that site -- exactly W bp, the same interval the trackhub
                   BED draws. The strand reports of one site collapse into one.
       "element"   one value per DAR, pooling that motif's best hits in it.

Background for every unit is the DAR total minus that motif's own MERGED
coverage -- every base any hit touches is removed, so no motif base leaks into
the background even though the numerator uses only the best hit.
"""
import numpy as np, collections, os
F="/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus"
A="/private/warm-archive/colquittlab/assembly/cactus/avian363"
MOTIF=os.environ.get("AR_MOTIF","jaspar")
UNIT=os.environ.get("AR_UNIT","site")
LABEL={"full":"HOMER ARE-full (12 bp)","half":"HOMER ARE-half (8 bp)",
       "jaspar":"JASPAR MA0007.3 (16 bp)"}

def _tab(path):
    """name -> (size, covered, sum)"""
    d={}
    for l in open(path):
        p=l.split("\t"); d[p[0]]=(int(p[1]),int(p[2]),float(p[3]))
    return d

def _dars(tree):
    return _tab(f"{F}/{tree}_am_dars.tab")

def pairs(tree,motif=None,unit=None,keep=None):
    """-> {(class, locus): array of (motif, rest)}

    keep: optional set of element ids. When given, only those elements are used
    and rows are tagged locus="sub", so group(P,cls,"sub") selects them."""
    motif=motif or MOTIF; unit=unit or UNIT
    dar=_dars(tree)
    sites=_tab(f"{F}/{tree}_ar3_sites.tab")
    # background subtrahend: every base any hit of this motif touches
    etot=collections.defaultdict(lambda:[0,0,0.0]); meta={}
    for k,(sz,cv,sm) in sites.items():
        m,s,loc,e,_=k.split("|")
        if m!=motif: continue
        a=etot[e]; a[0]+=sz; a[1]+=cv; a[2]+=sm; meta[e]=(s,loc)
    # numerator: the best-scoring hit of each site, from the per-position table
    num=collections.defaultdict(lambda: collections.defaultdict(lambda:[0,0,0.0]))
    if unit=="instance":
        for k,(sz,cv,sm) in _tab(f"{F}/{tree}_ar3_hits.tab").items():
            m,s,loc,e,j=k.split("|")
            if m!=motif or cv<0.5*sz: continue
            a=num[e][j]; a[0]+=sz; a[1]+=cv; a[2]+=sm
    else:
        for l in open(f"{F}/{tree}_ar3_pos.tab"):
            p=l.split("\t"); m,s,loc,e,site,k=p[0].split("|")
            if m!=motif: continue
            a=num[e][site]; a[0]+=int(p[1]); a[1]+=int(p[2]); a[2]+=float(p[3])
    out=collections.defaultdict(list)
    for e,d in num.items():
        if e not in dar or e not in etot: continue
        if keep is not None and e not in keep: continue
        dsz,dcov,dsum=dar[e]
        esz,ecov,esum=etot[e]
        rcov=dcov-ecov; rsum=dsum-esum
        if ecov<=0 or rcov<0.5*(dsz-esz) or rcov<50: continue
        rest=rsum/rcov; key=(meta[e][0],"sub") if keep is not None else meta[e]
        good=[(cv,sm) for sz,cv,sm in d.values() if cv>=0.5*sz and cv>0]
        if not good: continue
        if unit=="element":
            tc=sum(c for c,_ in good); ts=sum(v for _,v in good)
            out[key].append((ts/tc,rest))
        else:
            for cv,sm in good: out[key].append((sm/cv,rest))
    return {k:np.array(v) for k,v in out.items()}

def group(P,cls,loc=None):
    v=[x for (c,l),arr in P.items() if c==cls and (loc is None or l==loc) for x in arr]
    return np.array(v) if v else np.empty((0,2))
