#!/usr/bin/env python
"""Motif scoring across the SCENIC+ nominated MAFB regulators.

Same statistic as the AR motif-scoring analysis: delta = logodds(oscine
consensus) - logodds(outgroup consensus), with symmetric ascertainment
(forward = windows where the oscine consensus scores high, reverse = where the
outgroup does) and the midpoint as the winner's-curse-corrected value.
Thresholds are per-motif, at a matched percentile of each motif's own score
distribution, since a fixed cut is not equally stringent across widths.

Run inside the DHX35-TOP1 domain, where the alignment is extracted.
"""
import sys, collections, json, numpy as np
sys.path.insert(0,"."); sys.path.insert(0,"..")
from ascertain_sym import load_cols, OSC, OUT, cls, dom, MINO, MING
F="/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus"
HM="/private/home/bcolquit/group/motifs/homer/ra_pos_vs_neg/homerResults_manual"
PC=0.001; PCT=99.5; NB=2000
COMP=str.maketrans("ACGT","TGCA"); B={"A":0,"C":1,"G":2,"T":3}
REG={"KLF":"Glut-RA","NR3C2":"Glut-RA","NFIA":"Glut-RA","LIN28B":"Glut-RA",
     "STAT1":"Glut-RA","FOSL2":"Glut-RA","AR":"Glut-RA",
     "MAF":"GABA-4","SOX6":"GABA-4","TCF12":"GABA-4","ARX":"GABA-4",
     "SATB1":"GABA-4","MEIS2":"GABA-4"}
def meme_all(path):
    """-> {name: prob matrix} for every MOTIF block"""
    out={}; cur=None; rows=[]; take=False
    for l in open(path):
        if l.startswith("MOTIF"):
            if cur and rows: out[cur]=np.array(rows)
            cur=l.split()[1]; rows=[]; take=False
        elif l.startswith("letter-probability"): take=True
        elif take:
            p=l.split()
            if len(p)==4:
                try: rows.append([float(x) for x in p])
                except ValueError: take=False
            elif rows: take=False
    if cur and rows: out[cur]=np.array(rows)
    return out
def lo(p): return np.log2((p+PC)/(1+4*PC)/.25)
PWM={k:lo(v) for k,v in meme_all(f"{F}/tf13.meme").items()}
PWM={("KLF" if k=="KLF12" else k):v for k,v in PWM.items() if k!="KLF15"}
PWM["AR"]=lo(np.array([[float(x) for x in l.split()] for l in
       open(f"{F}/ar_jaspar.meme").read().split("letter-probability")[1].split("\n")[1:]
       if len(l.split())==4][:16]))
PWM["ARE-full(HOMER)"]=lo(np.array([[float(x) for x in l.split()]
       for l in open(f"{HM}/motif2.motif") if not l.startswith(">")]))
col=load_cols("data/domain_363.maf")
ocons={};gcons={}
for rp,d in col.items():
    c=collections.Counter(b for s,b in d.items() if s in OSC)
    if sum(c.values())>=MINO: ocons[rp]=c.most_common(1)[0][0]
    c=collections.Counter(b for s,b in d.items() if s in OUT)
    if sum(c.values())>=MING: gcons[rp]=c.most_common(1)[0][0]
def sc(seq,m):
    W=len(m)
    f=sum(m[k][B[seq[k]]] for k in range(W))
    r=sum(m[k][B[seq[W-1-k].translate(COMP)]] for k in range(W))
    return max(f,r)
rng=np.random.default_rng(0)
print(f"  {len(PWM)} matrices; consensus columns {len(ocons):,} oscine / {len(gcons):,} outgroup\n")
print(f"  {'TF':<17}{'regulon':<9}{'W':>3}{'RA mid':>9}{'G4 mid':>9}{'RA-G4':>9}"
      f"{'95% CI':>18}{'n RA f/r':>11}")
res={}
for nm,m in sorted(PWM.items(),key=lambda kv:kv[0]):
    W=len(m); rec=collections.defaultdict(list)
    for e in dom:
        ch,es,ee,cl=cls[e]
        kls="Glut-RA" if cl=="RA_over_Arco" else "GABA-4"
        for st in range(es,ee-W+1):
            so=[];sg=[];ok=True
            for k in range(W):
                a=ocons.get(st+k); b=gcons.get(st+k)
                if a is None or b is None: ok=False; break
                so.append(a); sg.append(b)
            if ok: rec[kls].append((e,sc(so,m),sc(sg,m)))
    allv=[v for kk in rec for r in rec[kk] for v in (r[1],r[2])]
    if not allv: continue
    th=float(np.percentile(allv,PCT))
    def mid(rows,els=None):
        by=collections.defaultdict(list)
        for r in rows: by[r[0]].append(r)
        src=els if els is not None else list(by)
        f=[r[1]-r[2] for e in src for r in by[e] if r[1]>=th]
        v=[r[1]-r[2] for e in src for r in by[e] if r[2]>=th]
        if len(f)<4 or len(v)<4: return None
        return (np.mean(f)+np.mean(v))/2,len(f),len(v)
    A=mid(rec["Glut-RA"]); Bv=mid(rec["GABA-4"])
    if not A or not Bv:
        print(f"  {nm:<17}{REG.get(nm,'-'):<9}{W:>3}{'too few windows':>27}"); continue
    eA=sorted({r[0] for r in rec["Glut-RA"]}); eB=sorted({r[0] for r in rec["GABA-4"]})
    bs=[]
    for _ in range(NB):
        x=mid(rec["Glut-RA"],list(rng.choice(eA,len(eA),replace=True)))
        y=mid(rec["GABA-4"],list(rng.choice(eB,len(eB),replace=True)))
        if x and y: bs.append(x[0]-y[0])
    lo_,hi_=np.percentile(bs,[2.5,97.5]) if bs else (np.nan,np.nan)
    star="*" if lo_>0 or hi_<0 else " "
    print(f"  {nm:<17}{REG.get(nm,'-'):<9}{W:>3}{A[0]:>+9.2f}{Bv[0]:>+9.2f}"
          f"{A[0]-Bv[0]:>+9.2f}{f'[{lo_:+.2f},{hi_:+.2f}]':>18}{f'{A[1]}/{A[2]}':>11}{star}")
    res[nm]=dict(reg=REG.get(nm,"-"),W=W,ra=float(A[0]),g4=float(Bv[0]),
                 diff=float(A[0]-Bv[0]),ci=[float(lo_),float(hi_)],
                 n=[A[1],A[2],Bv[1],Bv[2]])
json.dump(res,open("data/ascertain_tf.json","w"),indent=1)
