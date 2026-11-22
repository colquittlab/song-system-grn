#!/usr/bin/env python
"""Within-element control for the PWM-match-change metric.

The raw delta (oscine consensus minus outgroup consensus, scored with the AR
PWM) is inflated because the sites were called by FIMO on an oscine reference:
the oscine sequence sits at a selected optimum, the outgroup sequence does not.

Control: tile every element that carries an AR site with same-width windows
that do not overlap any AR hit, score each the same way, and subtract the
element's control mean from its AR site. The ascertainment affects only the
real sites, so the paired difference isolates it; comparing that paired
quantity between Glut-RA and GABA-4 is the difference-in-differences.
"""
import json, sys, collections, numpy as np
F="/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus"
D=f"{F}/maf_switch"
sys.path.insert(0,D)
from pwm_delta import LO, B, PWM
MINO,MING=40,60
COMP=str.maketrans("ACGT","TGCA")
OSC={l.strip() for l in open(f"{D}/data/osc363.txt") if l.strip()}
OUT={l.strip() for l in open(f"{D}/data/out363.txt") if l.strip()}
DOM={l.strip() for l in open(f"{D}/data/domain_dars.txt") if l.strip()}
cls={}
for l in open(f"{F}/allelem_cls.bed"):
    g=l.rstrip("\n").split("\t"); cls[g[3]]=(g[0],int(g[1]),int(g[2]),g[4])

def load_cols(path):
    col=collections.defaultdict(dict); block=[]
    def flush(b):
        ref=None
        for sp,st,txt in b:
            if sp.split(".")[0]=="Lonchura_striata": ref=(st,txt); break
        if ref is None: return
        st,rt=ref; p=st; idx=[]
        for i,c in enumerate(rt):
            if c!="-": idx.append((i,p)); p+=1
        for sp,_,txt in b:
            g=sp.split(".")[0]
            for i,rp in idx:
                if i<len(txt):
                    c=txt[i].upper()
                    if c in "ACGT": col[rp].setdefault(g,c)
    for l in open(path):
        if l.startswith("a"): flush(block); block=[]
        elif l.startswith("s\t") or l.startswith("s "):
            f=l.rstrip("\n").split("\t") if "\t" in l else l.split()
            block.append((f[1],int(f[2]),f[6]))
    flush(block)
    return col

def consensus(col,rp,grp):
    d=col.get(rp,{})
    c=collections.Counter(b for s,b in d.items() if s in grp)
    n=sum(c.values())
    return (c.most_common(1)[0][0],n) if n else (None,0)

def score_window(col,start,W,mk):
    """-> (delta, ok) using the better-scoring strand, as FIMO would"""
    so=[];sg=[]
    for k in range(W):
        bo,no=consensus(col,start+k,OSC); bg,ng=consensus(col,start+k,OUT)
        if bo is None or bg is None or no<MINO or ng<MING: return None
        so.append(bo); sg.append(bg)
    def sc(seq):
        f=sum(LO[mk][k][B[seq[k]]] for k in range(W))
        r=sum(LO[mk][k][B[seq[W-1-k].translate(COMP)]] for k in range(W))
        return max(f,r)
    return sc(so)-sc(sg)

if __name__=="__main__":
    col=load_cols(f"{D}/data/ar_elements_363.maf")
    print(f"  reference columns loaded: {len(col):,}")
    hits=collections.defaultdict(list)          # element -> list of (start,end) AR hits, any motif
    for bed in ("tracks/mafb_1p5mb/arefull_sites_MAFB_lsd2.bed","tracks/mafb_1p5mb/jasparAR_sites_MAFB_lsd2.bed"):
        for l in open(f"{F}/{bed}"):
            x=l.split("\t"); hits[x[10]].append((int(x[1]),int(x[2])))
    rows=[]
    for bed,mk in (("tracks/mafb_1p5mb/arefull_sites_MAFB_lsd2.bed","full"),("tracks/mafb_1p5mb/jasparAR_sites_MAFB_lsd2.bed","jaspar")):
        W=len(LO[mk])
        for l in open(f"{F}/{bed}"):
            x=l.split("\t"); e=x[10]; s=int(x[1])
            dsite=score_window(col,s,W,mk)
            if dsite is None: continue
            es,ee=cls[e][1],cls[e][2]
            ctrl=[]
            for st in range(es,ee-W+1,W):                 # non-overlapping tiles
                if any(st<he and st+W>hs for hs,he in hits[e]): continue
                v=score_window(col,st,W,mk)
                if v is not None: ctrl.append(v)
            if len(ctrl)<3: continue
            rows.append(dict(motif=mk,site=x[3],elem=e,cls=("Glut-RA" if x[9]=="Glut-RA" else "GABA-4"),
                             dom=e in DOM,site_delta=dsite,ctrl_mean=float(np.mean(ctrl)),
                             n_ctrl=len(ctrl),paired=dsite-float(np.mean(ctrl))))
    json.dump(rows,open(f"{D}/data/pwm_control.json","w"),indent=1,default=float)
    rng=np.random.default_rng(0)
    def p1(a,n=20000):
        s=rng.choice([-1,1],(n,len(a)))
        return a.mean(),((np.abs((s*a).mean(1))>=abs(a.mean())).sum()+1)/(n+1)
    def perm(a,b,n=20000):
        o=a.mean()-b.mean(); c=np.concatenate([a,b]); k=len(a); h=0
        for _ in range(n):
            p=rng.permutation(c); h+=abs(p[:k].mean()-p[k:].mean())>=abs(o)
        return o,(h+1)/(n+1)
    st=lambda p:"*" if p<.05 else " "
    print(f"  scorable AR sites with >=3 control windows: {len(rows)}"
          f"   median controls per site: {int(np.median([r['n_ctrl'] for r in rows]))}\n")
    print("  AR-site delta, control-window delta, and the paired difference")
    print(f"  {'scope':<20}{'motif':<9}{'class':<9}{'site':>9}{'control':>9}{'paired':>18}  n")
    for scope,sel in (("MAFB +/-1.5 Mb",lambda r:True),("DHX35-TOP1 domain",lambda r:r["dom"])):
        for mk,ml in (("full","ARE-full"),("jaspar","JASPAR")):
            for c in ("Glut-RA","GABA-4"):
                S=[r for r in rows if r["motif"]==mk and r["cls"]==c and sel(r)]
                if len(S)<3: print(f"  {scope:<20}{ml:<9}{c:<9}{'too few':>9}  n={len(S)}"); continue
                a=np.array([r["paired"] for r in S])
                m,p=p1(a)
                print(f"  {scope:<20}{ml:<9}{c:<9}"
                      f"{np.mean([r['site_delta'] for r in S]):>+9.2f}"
                      f"{np.mean([r['ctrl_mean'] for r in S]):>+9.2f}"
                      f"{m:>+11.2f} (p{p:5.3f}){st(p)}  {len(S)}")
            A=[r["paired"] for r in rows if r["motif"]==mk and r["cls"]=="Glut-RA" and sel(r)]
            Bq=[r["paired"] for r in rows if r["motif"]==mk and r["cls"]=="GABA-4" and sel(r)]
            if len(A)>=3 and len(Bq)>=3:
                d,pd=perm(np.array(A),np.array(Bq))
                print(f"  {'':<20}{ml:<9}{'DoD':<9}{'':>18}{d:>+11.2f} (p{pd:5.3f}){st(pd)}"
                      f"  {len(A)}/{len(Bq)}")
        print()
