#!/usr/bin/env python
"""Per-element motif-minus-rest-of-enhancer values, keeping coordinates.

Same statistic as ar_core (merged-site numerator, element-minus-motif
denominator) but returns one row per element so values can be ranked and
looked up in a browser, and so windows can be slid along the genome.
"""
import collections, numpy as np
F="/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus"
def load(tree,motif):
    dar={}
    for l in open(f"{F}/{tree}_am_dars.tab"):
        p=l.split("\t"); dar[p[0]]=(int(p[1]),int(p[2]),float(p[3]))
    sites=collections.defaultdict(lambda:[0,0,0.0]); meta={}
    for l in open(f"{F}/{tree}_ar3_sites.tab"):
        p=l.split("\t"); m,s,loc,e,_=p[0].split("|")
        if m!=motif: continue
        a=sites[e]; a[0]+=int(p[1]); a[1]+=int(p[2]); a[2]+=float(p[3]); meta[e]=s
    num=collections.defaultdict(lambda: collections.defaultdict(lambda:[0,0,0.0]))
    for l in open(f"{F}/{tree}_ar3_pos.tab"):
        p=l.split("\t"); m,s,loc,e,site,k=p[0].split("|")
        if m!=motif: continue
        a=num[e][site]; a[0]+=int(p[1]); a[1]+=int(p[2]); a[2]+=float(p[3])
    out={}
    for e,d in num.items():
        if e not in dar or e not in sites: continue
        dsz,dcov,dsum=dar[e]; esz,ecov,esum=sites[e]
        rcov=dcov-ecov; rsum=dsum-esum
        if ecov<=0 or rcov<0.5*(dsz-esz) or rcov<50: continue
        good=[(cv,sm) for sz,cv,sm in d.values() if cv>=0.5*sz and cv>0]
        if not good: continue
        tc=sum(c for c,_ in good); ts=sum(v for _,v in good)
        out[e]=dict(motif=ts/tc,rest=rsum/rcov,val=ts/tc-rsum/rcov,
                    nsite=len(good),cls=meta[e])
    return out
