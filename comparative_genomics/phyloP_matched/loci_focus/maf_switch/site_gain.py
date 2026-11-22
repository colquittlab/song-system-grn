#!/usr/bin/env python
"""Per-lineage presence of an AR site: candidate gains at the MAFB locus.

The question is whether a scoring AR site exists in each species at the
orthologous position, not whether individual bases are conserved. That makes
indels a mechanism of gain rather than a confound: a deletion that brings two
half-sites into register creates a site just as a substitution does.

For each AR site called on lonStrDom2, the aligned sequence of every species
across a padded window is degapped back to that species' own sequence, scanned
with the PWM on both strands, and scored present if its best hit clears the
same threshold used to call the reference site. A site present in most oscines
and absent from most outgroup species is a candidate gain.
"""
import sys, collections, json, numpy as np
sys.path.insert(0,"."); sys.path.insert(0,"..")
F="/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus"
HM="/private/home/bcolquit/group/motifs/homer/ra_pos_vs_neg/homerResults_manual"
PC=0.001; PAD=30
COMP=str.maketrans("ACGTacgt","TGCAtgca"); B={"A":0,"C":1,"G":2,"T":3}
OSC={l.strip() for l in open("data/osc363.txt") if l.strip()}
OUT={l.strip() for l in open("data/out363.txt") if l.strip()}
def homer_lo(f):
    p=np.array([[float(x) for x in l.split()] for l in open(f) if not l.startswith(">")])
    return np.log2((p+PC)/(1+4*PC)/.25)
def meme_lo(f,W):
    b=open(f).read().split("letter-probability")[1].split("\n")[1:]
    p=np.array([[float(x) for x in l.split()] for l in b if len(l.split())==4][:W])
    return np.log2((p+PC)/(1+4*PC)/.25)
LO={"full":homer_lo(f"{HM}/motif2.motif"),"jaspar":meme_lo(f"{F}/ar_jaspar.meme",16)}
# FIMO p<1e-4 corresponds to roughly the 99.99th percentile of the score null;
# calibrate each matrix on its own uniform-background distribution
def thresh(lo,n=200000,q=99.99):
    rng=np.random.default_rng(0); W=len(lo)
    s=rng.integers(0,4,(n,W))
    v=np.array([lo[k][s[:,k]] for k in range(W)]).sum(0)
    return float(np.percentile(v,q))
TH={k:thresh(v) for k,v in LO.items()}
print(f"  thresholds (99.99th pct of the uniform null): "
      f"ARE-full {TH['full']:+.2f}, JASPAR {TH['jaspar']:+.2f}")
def best(seq,lo):
    W=len(lo); s=seq.upper()
    if len(s)<W: return -99.0
    bs=-99.0
    for i in range(len(s)-W+1):
        sub=s[i:i+W]
        if any(c not in B for c in sub): continue
        f=sum(lo[k][B[sub[k]]] for k in range(W))
        r=sum(lo[W-1-k][B[sub[k].translate(COMP)]] for k in range(W))
        bs=max(bs,f,r)
    return bs
# load the gapped MAF, keyed by reference window
blocks=[];cur=[]
for l in open("data/ar_windows_gap.maf"):
    if l.startswith("a"):
        if cur: blocks.append(cur); 
        cur=[]
    elif l.startswith("s\t"):
        f=l.rstrip("\n").split("\t"); cur.append((f[1],int(f[2]),int(f[3]),f[6]))
if cur: blocks.append(cur)
byref=collections.defaultdict(list)
for b in blocks:
    ref=[r for r in b if r[0].startswith("Lonchura_striata.")]
    if ref: byref[(ref[0][1],ref[0][1]+ref[0][2])].append(b)
print(f"  alignment blocks: {len(blocks):,} over {len(byref)} reference spans")
sites=[]
for bed,mo in ((f"{F}/tracks/mafb_1p5mb/arefull_sites_MAFB_lsd2.bed","full"),
               (f"{F}/tracks/mafb_1p5mb/jasparAR_sites_MAFB_lsd2.bed","jaspar")):
    for l in open(bed):
        x=l.rstrip("\n").split("\t")
        sites.append(dict(motif=mo,name=x[3],s=int(x[1]),e=int(x[2]),
                          cls=x[9],elem=x[10]))
res=[]
for st in sites:
    lo=LO[st["motif"]]; w0,w1=st["s"]-PAD,st["e"]+PAD
    seq=collections.defaultdict(str)
    for (bs,be),bl in byref.items():
        if be<=w0 or bs>=w1: continue
        for b in bl:
            ref=[r for r in b if r[0].startswith("Lonchura_striata.")][0]
            coord=[];p=ref[1]
            for c in ref[3]:
                if c=="-": coord.append(None)
                else: coord.append(p); p+=1
            keep=[i for i,c in enumerate(coord) if c is None or w0<=c<w1]
            if not keep: continue
            for sp,_,_,t in b:
                g=sp.split(".")[0]
                if g.startswith("birdAnc"): continue
                seq[g]+="".join(t[i] for i in keep if i<len(t)).replace("-","")
    o=[g for g in seq if g in OSC]; u=[g for g in seq if g in OUT]
    if len(o)<40 or len(u)<60: continue
    po=sum(1 for g in o if best(seq[g],lo)>=TH[st["motif"]])
    pu=sum(1 for g in u if best(seq[g],lo)>=TH[st["motif"]])
    res.append(dict(**{k:st[k] for k in ("motif","name","cls","elem","s","e")},
                    n_osc=len(o),n_out=len(u),p_osc=po,p_out=pu,
                    f_osc=po/len(o),f_out=pu/len(u)))
json.dump(res,open("data/site_gain.json","w"),indent=2)
print(f"  sites scored: {len(res)}")
