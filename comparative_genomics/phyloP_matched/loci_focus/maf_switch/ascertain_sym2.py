#!/usr/bin/env python
"""Symmetric ascertainment with motif-matched thresholds.

A fixed log-odds cut is not equally stringent across motifs: JASPAR is 16 bp
with 7 sub-1-bit columns, so its whole score distribution sits lower and a
shared cut slices further into its tail. Each motif's threshold is therefore
set at the same PERCENTILE of its own pooled window-score distribution.

delta = logodds(oscine consensus) - logodds(outgroup consensus), and the
midpoint of the oscine-ascertained and outgroup-ascertained means is the
winner's-curse-corrected estimate. CIs bootstrap over ELEMENTS, not windows,
since windows within an element overlap heavily.
"""
import sys, json, collections, numpy as np
sys.path.insert(0,"."); sys.path.insert(0,"..")
from ascertain_sym import load_cols, maj, sc, OSC, OUT, cls, dom, MINO, MING
from pwm_delta import LO
PCT=float(sys.argv[1]) if len(sys.argv)>1 else 99.5
NB=4000
col=load_cols("data/domain_363.maf")
W={mk:len(LO[mk]) for mk in ("full","jaspar")}
rec=collections.defaultdict(list)          # (motif,class) -> (elem, osc_score, out_score)
for e in dom:
    ch,es,ee,cl=cls[e]
    c="Glut-RA" if cl=="RA_over_Arco" else "GABA-4"
    for mk in ("full","jaspar"):
        for st in range(es,ee-W[mk]+1):
            so=[];sg=[];ok=True
            for k in range(W[mk]):
                bo,no=maj(col,st+k,OSC); bg,ng=maj(col,st+k,OUT)
                if bo is None or bg is None or no<MINO or ng<MING: ok=False; break
                so.append(bo); sg.append(bg)
            if ok: rec[(mk,c)].append((e,sc(so,mk),sc(sg,mk)))
TH={}
for mk in ("full","jaspar"):
    allv=[v for c in ("Glut-RA","GABA-4") for r in rec[(mk,c)] for v in (r[1],r[2])]
    TH[mk]=float(np.percentile(allv,PCT))
print(f"  matched cuts at the {PCT}th percentile of each motif's own window scores:")
for mk in ("full","jaspar"): print(f"    {mk:<8}{TH[mk]:+.2f} log-odds")
rng=np.random.default_rng(0)
def boot(rows,mk):
    """-> (fwd, rev, mid) point estimates and element-bootstrap CIs"""
    els=sorted({r[0] for r in rows})
    by=collections.defaultdict(list)
    for r in rows: by[r[0]].append(r)
    def stat(sample):
        f=[r[1]-r[2] for e in sample for r in by[e] if r[1]>=TH[mk]]
        v=[r[1]-r[2] for e in sample for r in by[e] if r[2]>=TH[mk]]
        if not f or not v: return None
        return np.mean(f),np.mean(v),(np.mean(f)+np.mean(v))/2
    obs=stat(els)
    if obs is None: return None
    bs=[]
    for _ in range(NB):
        s=list(rng.choice(els,len(els),replace=True))
        t=stat(s)
        if t: bs.append(t)
    bs=np.array(bs)
    return obs,np.percentile(bs,[2.5,97.5],axis=0),len(els)
print(f"\n  delta = logodds(oscine consensus) − logodds(outgroup consensus)")
print(f"  CIs bootstrap over elements ({NB} resamples)\n")
print(f"  {'motif':<8}{'class':<9}{'ascert. oscine':>22}{'ascert. outgroup':>22}"
      f"{'midpoint':>22}{'n win':>12}")
out={}
for mk in ("full","jaspar"):
    for c in ("Glut-RA","GABA-4"):
        rows=rec[(mk,c)]
        nf=sum(1 for r in rows if r[1]>=TH[mk]); nv=sum(1 for r in rows if r[2]>=TH[mk])
        res=boot(rows,mk)
        if res is None:
            print(f"  {mk:<8}{c:<9}{'too few':>22}"); continue
        (f,v,m),ci,ne=res
        print(f"  {mk:<8}{c:<9}{f:>+10.2f} [{ci[0][0]:+.1f},{ci[1][0]:+.1f}]"
              f"{v:>+10.2f} [{ci[0][1]:+.1f},{ci[1][1]:+.1f}]"
              f"{m:>+10.2f} [{ci[0][2]:+.1f},{ci[1][2]:+.1f}]"
              f"{nf}/{nv:>7}")
        out[f"{mk}|{c}"]=dict(fwd=f,rev=v,mid=m,ci_mid=[ci[0][2],ci[1][2]],
                              n_fwd=nf,n_rev=nv,n_elem=ne,thresh=TH[mk])
print()
for mk in ("full","jaspar"):
    a=out.get(f"{mk}|Glut-RA"); b=out.get(f"{mk}|GABA-4")
    if a and b:
        print(f"  {mk:<8} corrected Glut-RA − GABA-4 = {a['mid']-b['mid']:+.2f} log-odds")
json.dump(out,open("data/ascertain_sym2.json","w"),indent=1)
