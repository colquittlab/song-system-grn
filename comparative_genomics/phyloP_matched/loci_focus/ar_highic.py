#!/usr/bin/env python
"""Variant of the paired statistic restricted to the motif's high-IC columns.

Baseline metric (unchanged): for each site, mean phyloP across the motif match
minus mean phyloP across the REST of the same DAR, where "rest" = the whole
element minus every base any hit of that motif touches.

High-IC variant: the motif side uses only columns with IC >= THR bits; the
rest-of-enhancer denominator is identical, so the two variants differ only in
which motif bases enter the numerator.
"""
import sys, collections, numpy as np, pickle
sys.path.insert(0,".")
import ar_core as C
F="/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus"
HM="/private/home/bcolquit/group/motifs/homer/ra_pos_vs_neg/homerResults_manual"
THR=1.0
def homer_pwm(f): return np.array([[float(x) for x in l.split()] for l in open(f)
                                   if not l.startswith(">")])
def meme_pwm(f,W):
    b=open(f).read().split("letter-probability")[1].split("\n")[1:]
    return np.array([[float(x) for x in l.split()] for l in b if len(l.split())==4][:W])
PWM={"full":homer_pwm(f"{HM}/motif2.motif"),"half":homer_pwm(f"{HM}/motif4.motif"),
     "jaspar":meme_pwm(f"{F}/ar_jaspar.meme",16)}
IC={k:np.array([2+(p*np.log2(np.clip(p,1e-9,None))).sum() for p in v]) for k,v in PWM.items()}
HI={k:set(np.where(v>=THR)[0]) for k,v in IC.items()}

def pairs_ic(tree,motif,keep=None,cols="hi"):
    """-> {(class,locus): array of (motif_side, rest)} using chosen columns"""
    dar=C._dars(tree); sites=C._tab(f"{F}/{tree}_ar3_sites.tab")
    etot=collections.defaultdict(lambda:[0,0,0.0]); meta={}
    for k,(sz,cv,sm) in sites.items():
        m,s,loc,e,_=k.split("|")
        if m!=motif: continue
        a=etot[e]; a[0]+=sz; a[1]+=cv; a[2]+=sm; meta[e]=(s,loc)
    sel=HI[motif] if cols=="hi" else (set(range(len(IC[motif])))-HI[motif] if cols=="lo"
                                      else set(range(len(IC[motif]))))
    num=collections.defaultdict(lambda: collections.defaultdict(lambda:[0,0.0]))
    for l in open(f"{F}/{tree}_ar3_pos.tab"):
        p=l.split("\t"); m,s,loc,e,site,k=p[0].split("|")
        if m!=motif or int(k) not in sel or int(p[2])==0: continue
        a=num[e][site]; a[0]+=1; a[1]+=float(p[3])
    out=collections.defaultdict(list)
    for e,d in num.items():
        if e not in dar or e not in etot: continue
        if keep is not None and e not in keep: continue
        dsz,dcov,dsum=dar[e]; esz,ecov,esum=etot[e]
        rcov=dcov-ecov; rsum=dsum-esum
        if ecov<=0 or rcov<0.5*(dsz-esz) or rcov<50: continue
        rest=rsum/rcov
        key=(meta[e][0],"sub") if keep is not None else meta[e]
        for cv,sm in d.values():
            if cv>=max(2,0.5*len(sel)): out[key].append((sm/cv,rest))
    return {k:np.array(v) for k,v in out.items()}

rng=np.random.default_rng(0)
def p1(a,n=20000):
    s=rng.choice([-1,1],(n,len(a)))
    return a.mean(),((np.abs((s*a).mean(1))>=abs(a.mean())).sum()+1)/(n+1)
def perm(a,b,n=20000):
    o=a.mean()-b.mean(); c=np.concatenate([a,b]); k=len(a); h=0
    for _ in range(n):
        p=rng.permutation(c); h+=abs(p[:k].mean()-p[k:].mean())>=abs(o)
    return o,(h+1)/(n+1)
if __name__=="__main__":
    dom=pickle.load(open("/tmp/mafb_domain.pkl","rb"))
    st=lambda p:"*" if p<.05 else " "
    for mk in ("full","half","jaspar"):
        print(f"  {mk:<7} IC per column: "+" ".join(f"{x:.1f}" for x in IC[mk])
              +f"   high-IC cols (>= {THR} bits): {sorted(i+1 for i in HI[mk])}")
    for scope,keep in (("genome-wide",None),("MAFB domain",dom)):
        for cols,cl in (("all","all 12/8/16 columns"),("hi","high-IC columns only"),
                        ("lo","low-IC columns only")):
            print(f"\n=== {scope}, {cl}   (oscine tree)")
            print(f"  {'motif':<9}{'Glut-RA':>20}{'GABA-4':>20}{'DoD':>20}   n")
            for mk,ml in (("full","ARE-full"),("half","ARE-half"),("jaspar","JASPAR")):
                P=pairs_ic("subtree222",mk,keep=keep,cols=cols)
                tag="sub" if keep is not None else None
                a=C.group(P,"ra",tag); b=C.group(P,"g4",tag)
                if len(a)<4 or len(b)<4:
                    print(f"  {ml:<9}{'too few':>20}  n={len(a)}/{len(b)}"); continue
                av=a[:,0]-a[:,1]; bv=b[:,0]-b[:,1]
                ra,pra=p1(av); g4,pg4=p1(bv); d,pd=perm(av,bv)
                print(f"  {ml:<9}{ra:>+11.3f} (p{pra:5.3f}){st(pra)}"
                      f"{g4:>+10.3f} (p{pg4:5.3f}){st(pg4)}{d:>+10.3f} (p{pd:5.3f}){st(pd)}"
                      f"  {len(av)}/{len(bv)}")
