#!/usr/bin/env python
"""Does the Glut-RA midpoint follow the motif, or the way the motif was found?

ARE-full was discovered by HOMER de novo on Glut-RA versus Glut-Arco-1, so the
matrix is fit to discriminate oscine Glut-RA sequence. Symmetric ascertainment
cannot correct that: the bias is in the PWM, not in the window selection. If
discovery bias is the explanation, every de novo motif from the same scan
should show the same positive Glut-RA midpoint. If only the AR motifs do, it
is not the explanation.

Curated JASPAR AR is carried through as the external reference.
"""
import sys, collections, json, numpy as np
sys.path.insert(0,"."); sys.path.insert(0,"..")
from ascertain_sym import load_cols, OSC, OUT, cls, dom, MINO, MING
HM="/private/home/bcolquit/group/motifs/homer/ra_pos_vs_neg/homerResults_manual"
F="/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus"
PC=0.001; PCT=99.5; COMP=str.maketrans("ACGT","TGCA"); B={"A":0,"C":1,"G":2,"T":3}
def homer(f):
    p=np.array([[float(x) for x in l.split()] for l in open(f) if not l.startswith(">")])
    return np.log2((p+PC)/(1+4*PC)/.25)
def meme(f,W):
    b=open(f).read().split("letter-probability")[1].split("\n")[1:]
    p=np.array([[float(x) for x in l.split()] for l in b if len(l.split())==4][:W])
    return np.log2((p+PC)/(1+4*PC)/.25)
NAMES={"motif1":"AP-1 (de novo)","motif2":"ARE-full (de novo)","motif3":"NFI (de novo)",
       "motif4":"ARE-half (de novo)","motif7":"bHLH (de novo)","motif8":"MAFB (de novo)",
       "motif12":"CEBP (de novo)","motif19":"KLF (de novo)"}
LOm={k:homer(f"{HM}/{k}.motif") for k in NAMES}
LOm["jaspar"]=meme(f"{F}/ar_jaspar.meme",16); NAMES["jaspar"]="JASPAR AR (curated)"
col=load_cols("data/domain_363.maf")
# consensus per column, once
ocons={};gcons={}
for rp,d in col.items():
    c=collections.Counter(b for s,b in d.items() if s in OSC)
    if sum(c.values())>=MINO: ocons[rp]=c.most_common(1)[0][0]
    c=collections.Counter(b for s,b in d.items() if s in OUT)
    if sum(c.values())>=MING: gcons[rp]=c.most_common(1)[0][0]
print(f"  consensus columns: oscine {len(ocons):,}  outgroup {len(gcons):,}\n")
def score(seq,lo):
    W=len(lo)
    f=sum(lo[k][B[seq[k]]] for k in range(W))
    r=sum(lo[k][B[seq[W-1-k].translate(COMP)]] for k in range(W))
    return max(f,r)
res={}
print(f"  {'motif':<22}{'W':>3}{'Glut-RA mid':>14}{'GABA-4 mid':>13}{'RA-G4':>9}{'n f/r RA':>12}")
for mk,nm in NAMES.items():
    lo=LOm[mk]; W=len(lo); rec=collections.defaultdict(list)
    for e in dom:
        ch,es,ee,cl=cls[e]
        kls="Glut-RA" if cl=="RA_over_Arco" else "GABA-4"
        for st in range(es,ee-W+1):
            so=[];sg=[];ok=True
            for k in range(W):
                a=ocons.get(st+k); b=gcons.get(st+k)
                if a is None or b is None: ok=False; break
                so.append(a); sg.append(b)
            if ok: rec[kls].append((score(so,lo),score(sg,lo)))
    allv=[v for kk in rec for r in rec[kk] for v in r]
    if not allv: continue
    th=float(np.percentile(allv,PCT)); mids={}
    for kls in ("Glut-RA","GABA-4"):
        rows=rec[kls]
        f=[a-b for a,b in rows if a>=th]; v=[a-b for a,b in rows if b>=th]
        mids[kls]=((np.mean(f)+np.mean(v))/2,len(f),len(v)) if len(f)>=5 and len(v)>=5 else None
    a=mids["Glut-RA"]; b=mids["GABA-4"]
    if a and b:
        print(f"  {nm:<22}{W:>3}{a[0]:>+14.2f}{b[0]:>+13.2f}{a[0]-b[0]:>+9.2f}"
              f"{f'{a[1]}/{a[2]}':>12}")
        res[mk]=dict(name=nm,W=W,ra=float(a[0]),g4=float(b[0]),diff=float(a[0]-b[0]),
                     n=[a[1],a[2],b[1],b[2]])
    else:
        print(f"  {nm:<22}{W:>3}{'too few windows':>14}")
json.dump(res,open("data/ascertain_panel.json","w"),indent=1)
