#!/usr/bin/env python
"""Per-motif DoD on both trees for the 182 archetypes + the 12 eRegulon PWMs.

One statistic, one DAR set, one filter for everything: motif phyloP minus the
phyloP of the remainder of the same element, element-aligned across the two
trees so the sampling covariance between the axes is estimable.
AR needs no special handling -- it is archetype metacluster_185.1 (MA0007.3).
"""
import numpy as np, collections, json, sys
F="/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus"
A="/private/warm-archive/colquittlab/assembly/cactus/avian363"
cls={}
for l in open(f"{F}/allelem_cls.bed"):
    f=l.rstrip("\n").split("\t"); cls[f[3]]=f[4]
mafb={l.strip() for l in open(f"{A}/vgp577/lifted/mafb_1mb.txt")}
sys.path.insert(0,F)
from span import elem_coords, in_domain, result
_el=elem_coords(); dom={e for e in mafb if in_domain(_el[e])}   # DHX35-TOP1 domain

def agg(path):
    d=collections.defaultdict(lambda:[0,0,0.0])
    for l in open(path):
        p=l.split("\t"); m,e,_=p[0].split("@")
        a=d[(m,e)]; a[0]+=int(p[1]); a[1]+=int(p[2]); a[2]+=float(p[3])
    return d

def dars(tree):
    d={}
    for l in open(f"{F}/{tree}_am_dars.tab"):
        p=l.split("\t"); d[p[0]]=(int(p[1]),int(p[2]),float(p[3]))
    return d

def diffs(aggd,dar):
    """(motif,element) -> value, with the archetype-pipeline QC"""
    out={}
    for (m,e),(msz,mcov,msum) in aggd.items():
        if e not in dar or e not in cls: continue
        dsz,dcov,dsum=dar[e]
        rcov=dcov-mcov; rsum=dsum-msum
        if mcov<0.5*msz or rcov<0.5*(dsz-msz) or rcov<50: continue
        out[(m,e)]=msum/mcov-rsum/rcov
    return out

from tf_core import SET
SETS={"arch":"am_hits.tab","ereg":f"{SET}_hits.tab"}   # regulators: the SCENIC+ matrix set
DO=dict(subtree222={},wholetree={})
for tree in ("subtree222","wholetree"):
    dar=dars(tree)
    for tag,suf in SETS.items():
        DO[tree][tag]=diffs(agg(f"{F}/{tree}_{suf}"),dar)
        print(f"  {tree:11s} {tag:5s} {len(DO[tree][tag]):>9,} (motif,element) values",flush=True)

rows=[]
for tag in SETS:
    kO=DO["subtree222"][tag]; kW=DO["wholetree"][tag]
    by=collections.defaultdict(lambda: collections.defaultdict(list))   # motif -> class -> [(o,w,isMafb)]
    for (m,e),vo in kO.items():
        vw=kW.get((m,e))
        if vw is None: continue
        c=cls.get(e)
        if c is None: continue
        by[m][c].append((vo,vw,e in mafb,e in dom))
    for m in sorted(by):
        rec={"motif":m,"set":tag}
        for loc,tagl in ((None,"gw"),("MAFB","MAFB"),("domain","domain")):
            per={}
            for c in ("RA_over_Arco","GABA4_over_GABA1"):
                v=by[m].get(c,[])
                if loc=="MAFB": v=[t for t in v if t[2]]
                if loc=="domain": v=[t for t in v if t[3]]
                per[c]=(np.array([t[0] for t in v]),np.array([t[1] for t in v])) if v else None
            a,b=per["RA_over_Arco"],per["GABA4_over_GABA1"]
            if a is None or b is None or len(a[0])<2 or len(b[0])<2:
                rec[tagl]=None; continue
            na,nb=len(a[0]),len(b[0])
            rec[tagl]=dict(na=na,nb=nb,
                y=float(a[0].mean()-b[0].mean()), x=float(a[1].mean()-b[1].mean()),
                vy=float(a[0].var(ddof=1)/na+b[0].var(ddof=1)/nb),
                vx=float(a[1].var(ddof=1)/na+b[1].var(ddof=1)/nb),
                cxy=float(np.cov(a[0],a[1])[0,1]/na+np.cov(b[0],b[1])[0,1]/nb))
        rows.append(rec)
    print(f"  {tag}: {len(by)} motifs",flush=True)
json.dump(rows,open(result("pairs_all.json","compare_mafb_spans"),"w"),indent=1)
print("wrote pairs_all.json")
