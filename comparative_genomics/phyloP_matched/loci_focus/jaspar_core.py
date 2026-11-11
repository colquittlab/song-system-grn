"""JASPAR MA0007.3 (AR) paired motif-vs-own-enhancer values, two units.

UNIT="instance"  one value per FIMO hit. The 'rest' is the element minus every
                 hit of this motif in it, so hits sharing an element share a
                 denominator but keep their own motif phyloP -- which is the
                 point: co-occurring sites can be under different constraint.
UNIT="element"   one value per DAR; the motif side pools all hits in it. This
                 is the unit the 182 archetypes and the 13 SCENIC+ TFs use.

Background is built by subtraction from the DAR total (`{tree}_am_dars.tab`),
the same way the archetype and TF pipelines build theirs -- not gwrest.bed,
whose segments were cut around a different motif.

Hit bases are de-duplicated by genomic position before the element total is
formed, so overlapping hits do not double-count a base.
"""
import numpy as np, collections, os
F="/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus"
A="/private/warm-archive/colquittlab/assembly/cactus/avian363"
UNIT=os.environ.get("AR_UNIT","instance")
MOTIF="JASPAR MA0007.3"

_cls={}
for _l in open(f"{F}/allelem_cls.bed"):
    _f=_l.rstrip("\n").split("\t"); _cls[_f[3]]=_f[4]
_mafb={l.strip() for l in open(f"{A}/vgp577/lifted/mafb_1mb.txt")}
# instance index -> element, in the order jaspos.bed was built (ra then g4, one counter)
_i2e={}; _n=0
for _s in ("ra","g4"):
    for _i,_l in enumerate(open(f"{F}/fimo_jaspar_{_s}.txt")):
        if not _i: continue
        _e=_l.split("\t")[2]
        if _e not in _cls: continue
        _n+=1; _i2e[_n]=_e
# instance index -> its genomic positions, for de-duplication
_ipos=collections.defaultdict(list); _ik=collections.defaultdict(list)
for _l in open(f"{F}/jaspos.bed"):
    _f=_l.rstrip("\n").split("\t"); _p=_f[3].split("|"); _x=int(_p[3])
    _ipos[_x].append((_f[0],int(_f[1]))); _ik[_x].append(int(_p[2]))

def _dars(tree):
    d={}
    for l in open(f"{F}/{tree}_am_dars.tab"):
        p=l.split("\t"); d[p[0]]=(int(p[1]),int(p[2]),float(p[3]))
    return d

def pairs(tree,unit=None):
    """-> {(class, locus): array of (motif, rest)} in the requested unit"""
    unit=unit or UNIT
    dar=_dars(tree)
    inst=collections.defaultdict(lambda:[0,0,0.0])     # idx -> size,cov,sum
    pv={}                                              # (idx,k) -> (cov,sum)
    for l in open(f"{F}/{tree}_jaspos.tab"):
        p=l.split("\t"); _,_,k,idx=p[0].split("|"); idx=int(idx)
        sz,cv,sm=int(p[1]),int(p[2]),float(p[3])
        a=inst[idx]; a[0]+=sz; a[1]+=cv; a[2]+=sm
        pv[(idx,int(k))]=(cv,sm)
    # element totals over hit bases de-duplicated by genomic position
    base=collections.defaultdict(dict)
    for idx,lst in _ipos.items():
        e=_i2e.get(idx)
        if e is None: continue
        for n,(chrom,pos) in enumerate(lst):
            k=_ik[idx][n]
            base[e][(chrom,pos)]=pv.get((idx,k),(0,0.0))
    etot={e:(len(d),sum(c for c,_ in d.values()),sum(v for _,v in d.values()))
          for e,d in base.items()}
    out=collections.defaultdict(list)
    for idx,(msz,mcov,msum) in inst.items():
        e=_i2e.get(idx)
        if e is None or e not in dar or e not in etot: continue
        if mcov<0.5*msz: continue
        dsz,dcov,dsum=dar[e]
        esz,ecov,esum=etot[e]
        rcov=dcov-ecov; rsum=dsum-esum
        if rcov<0.5*(dsz-esz) or rcov<50: continue
        loc="MAFB" if e in _mafb else "other"
        key=(("ra" if _cls[e]=="RA_over_Arco" else "g4"),loc)
        out[key].append((esum/ecov if unit=="element" else msum/mcov, rsum/rcov, e))
    if unit=="element":
        ded=collections.defaultdict(dict)
        for k,v in out.items():
            for m,r,e in v: ded[k][e]=(m,r)
        return {k:np.array([[m,r] for m,r in d.values()]) for k,d in ded.items()}
    return {k:np.array([[m,r] for m,r,_ in v]) for k,v in out.items()}

def group(P,cls,loc=None):
    v=[x for (c,l),arr in P.items() if c==cls and (loc is None or l==loc) for x in arr]
    return np.array(v) if v else np.empty((0,2))
