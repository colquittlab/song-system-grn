#!/usr/bin/env python
"""Symmetric ascertainment: how much of the oscine PWM advantage is winner's curse?

The published sites were called by FIMO on lonStrDom2, an oscine. Any such
position sits at a local maximum for the oscine sequence, so the outgroup can
only regress downward whether or not anything evolved.

Here the same elements are scanned three ways with the same PWM and threshold,
and the same statistic computed each time:

    delta = logodds(oscine consensus) - logodds(outgroup consensus)

  ascertained on the oscine consensus  -> forward,  biased toward +
  ascertained on the outgroup consensus -> reverse, biased toward -
  ascertained on neither (all windows)  -> unselected baseline

If the effect is pure ascertainment, forward and reverse are equal and
opposite. The midpoint (forward + reverse)/2 is the bias-corrected estimate;
the asymmetry is what is left over.
"""
import sys, json, collections, numpy as np
sys.path.insert(0,"."); sys.path.insert(0,"..")
from pwm_delta import LO, B
F="/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus"
D="."
MINO,MING=40,60
THRESH=float(sys.argv[1]) if len(sys.argv)>1 else 6.0   # log-odds cut for "a site"
COMP=str.maketrans("ACGT","TGCA")
OSC={l.strip() for l in open(f"{D}/data/osc363.txt") if l.strip()}
OUT={l.strip() for l in open(f"{D}/data/out363.txt") if l.strip()}
cls={}
for l in open(f"{F}/allelem_cls.bed"):
    g=l.rstrip("\n").split("\t"); cls[g[3]]=(g[0],int(g[1]),int(g[2]),g[4])
dom=[l.strip() for l in open(f"{D}/data/domain_dars.txt") if l.strip()]

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

def maj(col,rp,grp):
    d=col.get(rp,{})
    c=collections.Counter(b for s,b in d.items() if s in grp)
    n=sum(c.values())
    return (c.most_common(1)[0][0],n) if n else (None,0)

def sc(seq,mk):
    W=len(LO[mk])
    f=sum(LO[mk][k][B[seq[k]]] for k in range(W))
    r=sum(LO[mk][k][B[seq[W-1-k].translate(COMP)]] for k in range(W))
    return max(f,r)

if __name__=="__main__":
    col=load_cols(f"{D}/data/domain_363.maf")
    print(f"  columns loaded: {len(col):,}   threshold: log-odds >= {THRESH}")
    res=collections.defaultdict(list)
    for e in dom:
        ch,es,ee,cl = cls[e]
        c="Glut-RA" if cl=="RA_over_Arco" else "GABA-4"
        for mk in ("full","jaspar"):
            W=len(LO[mk])
            for st in range(es,ee-W+1):
                so=[];sg=[];ok=True
                for k in range(W):
                    bo,no=maj(col,st+k,OSC); bg,ng=maj(col,st+k,OUT)
                    if bo is None or bg is None or no<MINO or ng<MING: ok=False; break
                    so.append(bo); sg.append(bg)
                if not ok: continue
                a=sc(so,mk); b=sc(sg,mk)
                res[(mk,c,"all")].append(a-b)
                if a>=THRESH: res[(mk,c,"fwd")].append(a-b)
                if b>=THRESH: res[(mk,c,"rev")].append(a-b)
    rng=np.random.default_rng(0)
    def ci(v,n=4000):
        b=[np.mean(rng.choice(v,len(v))) for _ in range(n)]
        return float(np.mean(v)),float(np.percentile(b,2.5)),float(np.percentile(b,97.5))
    print(f"\n  delta = logodds(oscine consensus) - logodds(outgroup consensus)\n")
    print(f"  {'motif':<9}{'class':<9}{'ascertained on':<22}{'mean delta':>22}{'n windows':>12}")
    summ={}
    for mk in ("full","jaspar"):
        for c in ("Glut-RA","GABA-4"):
            for tag,lab in (("fwd","oscine consensus"),("rev","outgroup consensus"),
                            ("all","nothing (all windows)")):
                v=res.get((mk,c,tag),[])
                if len(v)<5: print(f"  {mk:<9}{c:<9}{lab:<22}{'too few':>22}{len(v):>12}"); continue
                m,lo,hi=ci(np.array(v))
                summ[(mk,c,tag)]=(m,len(v))
                print(f"  {mk:<9}{c:<9}{lab:<22}{m:>+10.2f} [{lo:+.2f},{hi:+.2f}]{len(v):>12,}")
            f=summ.get((mk,c,"fwd")); r=summ.get((mk,c,"rev"))
            if f and r:
                print(f"  {'':<9}{'':<9}{'-> midpoint (corrected)':<22}{(f[0]+r[0])/2:>+10.2f}"
                      f"{'':>12}  asymmetry {f[0]+r[0]:+.2f}")
            print()
    json.dump({f"{k[0]}|{k[1]}|{k[2]}":v for k,v in summ.items()},
              open(f"{D}/data/ascertain_sym.json","w"),indent=1)
