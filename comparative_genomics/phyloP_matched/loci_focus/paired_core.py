"""Shared loader: paired AR-motif vs own-DAR values per group and tree."""
import numpy as np, collections
F="/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus"
def _load(f):
    d={}
    for l in open(f"{F}/{f}"):
        p=l.split("\t"); d[p[0]]=(int(p[1]),int(p[2]),float(p[5]))
    return d
_coord={l.split("\t")[3]:1 for l in open(f"{F}/allelem_cls.bed")}
_h2e={}; _n=0
for _s in ("ra","g4"):
    for _i,_l in enumerate(open(f"{F}/fimo_motif2_{_s}.txt")):
        if not _i: continue
        _e=_l.split("\t")[2]
        if _e not in _coord: continue
        _n+=1; _h2e[_n]=_e
def pairs(tree,keep=None):
    """-> {(class, locus): array of (motif, rest) per instance}
    keep: optional set of element ids; rows are then tagged locus="sub", so
    group(P,cls,"sub") selects them (e.g. the DHX35-TOP1 domain DARs)."""
    mot=_load(f"{tree}_motif2.tab"); rest=_load(f"{tree}_gwrest.tab")
    rb=collections.defaultdict(list)
    for k,v in rest.items():
        if v[1]/max(v[0],1)>=.5: rb[k.split("|")[1]].append((v[2],v[0]))
    out=collections.defaultdict(list)
    for k,v in mot.items():
        if v[1]/max(v[0],1)<.5: continue
        p=k.split("|"); e=_h2e.get(int(p[3]))
        if not e or e not in rb: continue
        if keep is not None and e not in keep: continue
        w=np.array(rb[e]); out[(p[1],"sub" if keep is not None else p[2])].append((v[2],float(np.average(w[:,0],weights=w[:,1]))))
    return {k:np.array(v) for k,v in out.items()}
def group(P,cls,loc=None):
    v=[x for (c,l),arr in P.items() if c==cls and (loc is None or l==loc) for x in arr]
    return np.array(v) if v else np.empty((0,2))
