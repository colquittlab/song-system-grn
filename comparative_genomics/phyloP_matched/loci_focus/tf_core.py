"""Shared definitions for the SCENIC+ MAFB regulator analysis.

Motif set: SET names the scan whose per-hit phyloP tables are read
({tree}_{SET}_hits.tab). The default, "tfsc", is one matrix per regulator chosen
from the motif evidence SCENIC+ itself used for the regulator's cistrome in the
run that nominated it (build_scenic_motifs.py -> tf_scenic.meme). Every regulator,
AR included, is scanned and measured by the same pipeline. The superseded set,
"tf13", is kept readable for comparison.

In tfsc, AR and NR3C2 are both backed by metacluster_185.1 and get the same
matrix, so NR3C2 is not reported separately (SAME_AS). KLF12 is reported as KLF.
AR* is the HOMER de novo ARE-full motif measured per instance (paired_core): a
labelled reference, not a member of the comparison set.
"""
import numpy as np, collections, sys
F="/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus"
A="/private/warm-archive/colquittlab/assembly/cactus/avian363"
sys.path.insert(0,F)
from paired_core import pairs, group
REG={"KLF":"Glut-RA","NR3C2":"Glut-RA","NFIA":"Glut-RA","LIN28B":"Glut-RA",
     "STAT1":"Glut-RA","FOSL2":"Glut-RA","AR":"Glut-RA",
     "MAF":"GABA-4","SOX6":"GABA-4","TCF12":"GABA-4","ARX":"GABA-4",
     "SATB1":"GABA-4","MEIS2":"GABA-4"}
RENAME={"KLF12":"KLF"}; DROP={"KLF15"}
SET="tfsc"
SAME_AS={"NR3C2":"AR"} if SET=="tfsc" else {}     # identical matrix; reported as AR
_cls={l.split("\t")[3]:l.rstrip("\n").split("\t")[4] for l in open(f"{F}/allelem_cls.bed")}
_mafb={l.strip() for l in open(f"{A}/vgp577/lifted/mafb_1mb.txt")}
from span import elem_coords as _ec, in_domain as _ind
_ELC=_ec(); _dom={e for e in _mafb if _ind(_ELC[e])}     # DHX35-TOP1 domain DARs
def diffs(tree,loc=None):
    """-> {(TF, class): array of per-element (motif - rest of own enhancer)}"""
    dar={}
    for l in open(f"{F}/{tree}_am_dars.tab"):
        p=l.split("\t"); dar[p[0]]=(int(p[1]),int(p[2]),float(p[3]))
    agg=collections.defaultdict(lambda:[0,0,0.0])
    for l in open(f"{F}/{tree}_{SET}_hits.tab"):
        p=l.split("\t"); m,e,_=p[0].split("@")
        if m in DROP or m in SAME_AS: continue
        m=RENAME.get(m,m)
        a=agg[(m,e)]; a[0]+=int(p[1]); a[1]+=int(p[2]); a[2]+=float(p[3])
    out=collections.defaultdict(list)
    for (m,e),(msz,mcov,msum) in agg.items():
        if e not in dar or e not in _cls: continue
        if loc=="MAFB" and e not in _mafb: continue
        if loc=="domain" and e not in _dom: continue
        dsz,dcov,dsum=dar[e]
        rcov=dcov-mcov; rsum=dsum-msum
        if mcov<0.5*msz or rcov<50: continue
        out[(m,_cls[e])].append(msum/mcov-rsum/rcov)
    # AR* -- HOMER ARE-full per instance, kept as a labelled reference. Under the
    # superseded tf13 set this was the only AR measurement and was reported as "AR".
    if loc=="domain":
        P=pairs(tree,keep=_dom); ga=group(P,"ra","sub"); gb=group(P,"g4","sub")
    else:
        P=pairs(tree); ga=group(P,"ra",loc); gb=group(P,"g4",loc)
    ref="AR" if SET=="tf13" else "AR*"
    out[(ref,"RA_over_Arco")]=list(ga[:,0]-ga[:,1])
    out[(ref,"GABA4_over_GABA1")]=list(gb[:,0]-gb[:,1])
    return {k:np.array(v) for k,v in out.items()}
def bh(ps):
    ps=np.array(ps); o=np.argsort(ps); q=ps[o]*len(ps)/np.arange(1,len(ps)+1)
    q=np.minimum.accumulate(q[::-1])[::-1]; r=np.empty(len(ps)); r[o]=np.clip(q,0,1); return r
