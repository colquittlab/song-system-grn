#!/usr/bin/env python
"""Symmetric ascertainment with the oscine set downsampled to match outgroup
consensus precision.

The 217 outgroup species span far deeper divergence than the 146 oscines, so
their majority-base consensus is a noisier estimate (split-half score sd 1.09
vs 0.69; majority fraction 0.906 vs 0.951). That matters because the size of
the winner's curse in each ascertainment arm scales with the variance of the
quantity selected on, so unequal precision makes the forward and reverse arms
non-comparable.

Matching by split-half calibration: 217 outgroup species behave like ~42
oscines. Oscines are therefore drawn down to 42 per replicate, which equalises
precision by construction instead of reasoning about the direction of the bias.
"""
import sys, collections, json, numpy as np
sys.path.insert(0,"."); sys.path.insert(0,"..")
from ascertain_sym import load_cols, OSC, OUT, cls, dom, sc, MING
from pwm_delta import LO
NOSC=int(sys.argv[1]) if len(sys.argv)>1 else 42
NDRAW=int(sys.argv[2]) if len(sys.argv)>2 else 5
PCT=99.5
MINO=max(8,int(round(40/146*NOSC)))      # same coverage fraction as the full run
col=load_cols("data/domain_363.maf")
OSCl=sorted(OSC); W={mk:len(LO[mk]) for mk in ("full","jaspar")}
print(f"  oscines downsampled to {NOSC} ({NDRAW} draws); coverage floor {MINO}/{NOSC}"
      f" oscine, {MING}/217 outgroup")
# outgroup consensus is fixed across draws
gcons={}
for rp,d in col.items():
    c=collections.Counter(b for s,b in d.items() if s in OUT)
    n=sum(c.values())
    if n>=MING: gcons[rp]=c.most_common(1)[0][0]
rng=np.random.default_rng(0)
acc=collections.defaultdict(list)
for draw in range(NDRAW):
    sub=set(rng.choice(OSCl,NOSC,replace=False))
    ocons={}
    for rp,d in col.items():
        c=collections.Counter(b for s,b in d.items() if s in sub)
        n=sum(c.values())
        if n>=MINO: ocons[rp]=c.most_common(1)[0][0]
    rec=collections.defaultdict(list)
    for e in dom:
        ch,es,ee,cl=cls[e]
        kls="Glut-RA" if cl=="RA_over_Arco" else "GABA-4"
        for mk in ("full","jaspar"):
            w=W[mk]
            for st in range(es,ee-w+1):
                so=[];sg=[];ok=True
                for k in range(w):
                    a=ocons.get(st+k); b=gcons.get(st+k)
                    if a is None or b is None: ok=False; break
                    so.append(a); sg.append(b)
                if ok: rec[(mk,kls)].append((e,sc(so,mk),sc(sg,mk)))
    for mk in ("full","jaspar"):
        allv=[v for kk in ("Glut-RA","GABA-4") for r in rec[(mk,kk)] for v in (r[1],r[2])]
        th=float(np.percentile(allv,PCT))
        for kls in ("Glut-RA","GABA-4"):
            rows=rec[(mk,kls)]
            f=[r[1]-r[2] for r in rows if r[1]>=th]
            v=[r[1]-r[2] for r in rows if r[2]>=th]
            if len(f)<5 or len(v)<5: continue
            acc[(mk,kls)].append((np.mean(f),np.mean(v),(np.mean(f)+np.mean(v))/2,
                                  len(f),len(v),th))
    print(f"    draw {draw+1}/{NDRAW} done",flush=True)
print(f"\n  delta = logodds(oscine consensus) − logodds(outgroup consensus)")
print(f"  mean over draws, with the across-draw range\n")
print(f"  {'motif':<8}{'class':<9}{'fwd':>16}{'rev':>16}{'midpoint':>20}{'n f/r':>12}")
out={}
for mk in ("full","jaspar"):
    for kls in ("Glut-RA","GABA-4"):
        v=acc.get((mk,kls))
        if not v: print(f"  {mk:<8}{kls:<9}{'too few':>16}"); continue
        a=np.array(v)
        f,r,m=a[:,0].mean(),a[:,1].mean(),a[:,2].mean()
        print(f"  {mk:<8}{kls:<9}{f:>+9.2f} [{a[:,0].min():+.1f},{a[:,0].max():+.1f}]"
              f"{r:>+9.2f} [{a[:,1].min():+.1f},{a[:,1].max():+.1f}]"
              f"{m:>+11.2f} [{a[:,2].min():+.2f},{a[:,2].max():+.2f}]"
              f"{a[:,3].mean():.0f}/{a[:,4].mean():.0f}")
        out[f"{mk}|{kls}"]=dict(fwd=float(f),rev=float(r),mid=float(m),
                                mid_range=[float(a[:,2].min()),float(a[:,2].max())],
                                n_fwd=float(a[:,3].mean()),n_rev=float(a[:,4].mean()),
                                nosc=NOSC,ndraw=NDRAW)
print()
for mk in ("full","jaspar"):
    A=out.get(f"{mk}|Glut-RA"); Bq=out.get(f"{mk}|GABA-4")
    if A and Bq: print(f"  {mk:<8} corrected Glut-RA − GABA-4 = {A['mid']-Bq['mid']:+.2f}")
json.dump(out,open("data/ascertain_down.json","w"),indent=1)
