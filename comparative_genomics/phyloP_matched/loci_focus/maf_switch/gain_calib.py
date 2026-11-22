#!/usr/bin/env python
"""Reverse-ascertainment calibration for the site-gain statistic, with the
comparison sets matched on phylogenetic distance.

Sites were called on lonStrDom2, an oscine, so oscines are both closer to the
reference and guaranteed to carry the site. The gain score (fraction of oscines
with a site minus fraction of outgroup) is therefore inflated before any
biology.

Calibration: tile the domain, and let EACH of several genomes in turn play the
role lonStrDom2 played -- a window counts for that genome if its own sequence
carries a scoring site, and the gain score is then computed exactly as before.
Ascertaining on an oscine should give a positive gain, on an outgroup genome a
negative one. If the asymmetry is pure ascertainment the two are equal and
opposite; the midpoint is the corrected estimate.
"""
import sys, collections, json, numpy as np
F="/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus"
HM="/private/home/bcolquit/group/motifs/homer/ra_pos_vs_neg/homerResults_manual"
PC=0.001; WIN=76; STEP=38
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
def thresh(lo,n=200000,q=99.99):
    rng=np.random.default_rng(0); W=len(lo)
    s=rng.integers(0,4,(n,W))
    return float(np.percentile(np.array([lo[k][s[:,k]] for k in range(W)]).sum(0),q))
TH={k:thresh(v) for k,v in LO.items()}
def best(seq,lo):
    W=len(lo); s=seq.upper()
    if len(s)<W: return -99.0
    bs=-99.0
    for i in range(len(s)-W+1):
        sub=s[i:i+W]
        if any(c not in B for c in sub): continue
        f=sum(lo[k][B[sub[k]]] for k in range(W))
        r=sum(lo[W-1-k][B[sub[k].translate(COMP)]] for k in range(W))
        if f>bs: bs=f
        if r>bs: bs=r
    return bs
KNEAR=40          # species counted as "near" the ascertainment genome
def patristic(tree_path):
    """-> {leaf: {leaf: distance}} via root paths; good enough for ranking"""
    import re
    t=open(tree_path).read().strip().rstrip(";")
    pos=[0]
    def node():
        ch=[]
        if t[pos[0]]=="(":
            pos[0]+=1
            while True:
                ch.append(node())
                if pos[0]<len(t) and t[pos[0]]==",": pos[0]+=1
                else: break
            pos[0]+=1
        m=re.match(r"[^,()\[\]:;]*",t[pos[0]:]); nm=m.group(0); pos[0]+=len(nm)
        bl=0.0
        if pos[0]<len(t) and t[pos[0]]==":":
            m=re.match(r":([-\d.eE+]*)",t[pos[0]:]); bl=float(m.group(1) or 0)
            pos[0]+=len(m.group(0))
        return (nm,bl,ch)
    root=node()
    paths={}
    def walk(n,anc,d):
        nm,bl,ch=n; a=anc+[(id(n),d)]
        if not ch: paths[nm]=a
        for c in ch: walk(c,a,d+c[1])
    walk(root,[],0.0)
    return paths
PATHS=patristic("/private/warm-archive/colquittlab/assembly/cactus/avian363/363-avian-2020_lsd_tree.nw")
def dist(a,b):
    pa,pb=PATHS.get(a),PATHS.get(b)
    if not pa or not pb: return None
    i=0
    while i<len(pa) and i<len(pb) and pa[i][0]==pb[i][0]: i+=1
    anc=pa[i-1][1]
    return (pa[-1][1]-anc)+(pb[-1][1]-anc)

ASC_OSC=["Lonchura_striata","Taeniopygia_guttata","Corvus_brachyrhynchos",
         "Hirundo_rustica","Melospiza_melodia","Parus_major"]
ASC_OUT=["Gallus_gallus","Columba_livia","Falco_peregrinus","Melopsittacus_undulatus",
         "Acanthisitta_chloris","Anas_platyrhynchos","Apus_apus","Struthio_camelus"]
if __name__=="__main__":
    blocks=[];cur=[]
    for l in open("data/domain_gap.maf"):
        if l.startswith("a"):
            if cur: blocks.append(cur)
            cur=[]
        elif l.startswith("s\t"):
            f=l.rstrip("\n").split("\t"); cur.append((f[1],int(f[2]),int(f[3]),f[6]))
    if cur: blocks.append(cur)
    print(f"  blocks: {len(blocks):,}")
    # per reference window, each genome's own (degapped) sequence
    wins=collections.defaultdict(lambda: collections.defaultdict(str))
    for b in blocks:
        ref=[r for r in b if r[0].startswith("Lonchura_striata.")]
        if not ref: continue
        st=ref[0][1]; coord=[];p=st
        for c in ref[0][3]:
            if c=="-": coord.append(None)
            else: coord.append(p); p+=1
        real=[c for c in coord if c is not None]
        if not real: continue
        w0=(real[0]//STEP)*STEP
        for w in range(w0, real[-1]+1, STEP):
            keep=[i for i,c in enumerate(coord) if c is None or w<=c<w+WIN]
            if not keep: continue
            for sp,_,_,t in b:
                g=sp.split(".")[0]
                if g.startswith("birdAnc"): continue
                wins[w][g]+="".join(t[i] for i in keep if i<len(t)).replace("-","")
    print(f"  windows: {len(wins):,}")
    res=collections.defaultdict(list)
    for mo,lo in LO.items():
        sc={}
        for w,d in wins.items():
            if len(d)<200: continue
            sc[w]={g:best(s,lo) for g,s in d.items() if len(s)>=len(lo)}
        for G in ASC_OSC+ASC_OUT:
            gains=[]
            for w,d in sc.items():
                if d.get(G,-99)<TH[mo]: continue
                pres=[g for g in d]
                if len(pres)<200: continue
                dd=[(dist(G,g),g) for g in pres if g!=G]
                dd=[(x,g) for x,g in dd if x is not None]
                if len(dd)<KNEAR+60: continue
                dd.sort()
                near=[g for _,g in dd[:KNEAR]]; far=[g for _,g in dd[KNEAR:]]
                fo=sum(1 for g in near if d[g]>=TH[mo])/len(near)
                fu=sum(1 for g in far if d[g]>=TH[mo])/len(far)
                gains.append(fo-fu)
            if len(gains)>=10:
                res[(mo,G)]=gains
    print(f"\n  gain = fraction of oscines with a site − fraction of outgroup")
    print(f"  {'motif':<9}{'ascertained on':<26}{'clade':<10}{'windows':>9}{'mean gain':>12}")
    summ=collections.defaultdict(list)
    for (mo,G),v in sorted(res.items()):
        cl="oscine" if G in OSC else "outgroup"
        print(f"  {mo:<9}{G:<26}{cl:<10}{len(v):>9}{np.mean(v):>+12.3f}")
        summ[(mo,cl)].append(np.mean(v))
    print()
    for mo in LO:
        a=summ.get((mo,"oscine")); b=summ.get((mo,"outgroup"))
        if a and b:
            print(f"  {mo:<9} oscine-ascertained {np.mean(a):+.3f}   "
                  f"outgroup-ascertained {np.mean(b):+.3f}   "
                  f"midpoint {(np.mean(a)+np.mean(b))/2:+.3f}")
    json.dump({f"{k[0]}|{k[1]}":float(np.mean(v)) for k,v in res.items()},
              open("data/gain_calib.json","w"),indent=1)
