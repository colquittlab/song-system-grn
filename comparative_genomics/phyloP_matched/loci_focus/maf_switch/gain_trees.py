#!/usr/bin/env python
"""Per-site AR gain/loss mapped onto the 74-way tree.

Presence of a site is read directly from each node's own sequence -- cactus
stores a reconstructed sequence for every ancestral node, so a gain or loss is
observed on a branch rather than inferred from the tips. Each node's aligned
characters across the site window are degapped back to that node's own
sequence and scanned with the PWM on both strands; the node has the site if its
best hit clears the calling threshold. Branches where presence changes are the
gain and loss events.
"""
import sys, collections, json, os, numpy as np
sys.path.insert(0,".")
from treelib import parse_newick, leaves
import sys as _s; _s.path.insert(0, "/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus"); from span import per_site  # results/per_site/
import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams.update({"font.family":"sans-serif",
    "font.sans-serif":["Arial","Liberation Sans","DejaVu Sans"],
    "pdf.fonttype":42,"ps.fonttype":42,"font.size":5,"axes.linewidth":.4})
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
F="/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus"
HM="/private/home/bcolquit/group/motifs/homer/ra_pos_vs_neg/homerResults_manual"
PC=0.001; PAD=30
COMP=str.maketrans("ACGTacgt","TGCAtgca"); B={"A":0,"C":1,"G":2,"T":3}
INK="#171B22"; INK3="#6E7887"; RULE="#DDE1E7"
PRES="#C8161A"; ABS="#AEB6C2"; GAIN="#C7521B"; LOSS="#2C6E7F"
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
        bs=max(bs,f,r)
    return bs
def count_events(pres,tree):
    root=parse_newick(tree) if isinstance(tree,str) else tree
    g=l=0
    def rec(n):
        nonlocal g,l
        for c in n["ch"]:
            a,b=pres.get(n["name"]),pres.get(c["name"])
            if a is not None and b is not None and a!=b:
                if b: g+=1
                else: l+=1
            rec(c)
    rec(root); return g,l

def draw(pres,tree,stem,osc,title,collapse_all=False):
    """collapse_all: also collapse uniform oscine clades (needed on avian363, 144 oscines);
    collapsed clades are then labelled with their genus range"""
    root=parse_newick(tree) if isinstance(tree,str) else tree
    has=lambda n: any(x["name"] in osc for x in leaves(n))
    def nodes(n):
        yield n
        for c in n["ch"]: yield from nodes(c)
    def mark(n):
        # collapse only if EVERY node in the clade shares one state -- tips alone can
        # agree while an internal branch carries a gain and a later loss
        st={pres.get(x["name"]) for x in nodes(n)}
        n["collapse"]=(collapse_all or not has(n)) and len(st-{None})==1 and len(leaves(n))>1
        if not n["collapse"]:
            for c in n["ch"]: mark(c)
    mark(root)
    def rot(n):
        if not n["ch"]: return
        for c in n["ch"]: rot(c)
        b=[c for c in n["ch"] if has(c)]
        if len(b)==1: n["ch"].remove(b[0]); n["ch"].append(b[0])
    rot(root)
    def setx(n,x):
        n["x"]=x+n["bl"]
        for c in n["ch"]: setx(c,n["x"])
    setx(root,0.0)
    order=[]
    def lay(n):
        if n["collapse"] or not n["ch"]: n["y"]=len(order); order.append(n); return
        for c in n["ch"]: lay(c)
        n["y"]=float(np.mean([c["y"] for c in n["ch"]]))
    lay(root); N=len(order)
    fig,ax=plt.subplots(figsize=(4.6,max(1.9,N*0.108+0.62)))
    fig.patch.set_facecolor("white"); ev=[0,0]
    def rec(n):
        for c in (n["ch"] if not n["collapse"] else []):
            a,b2=pres.get(n["name"]),pres.get(c["name"])
            ch=(a is not None and b2 is not None and a!=b2)
            col=(GAIN if b2 else LOSS) if ch else RULE
            ax.plot([n["x"]]*2,[n["y"],c["y"]],color=RULE,lw=.55,zorder=1)
            ax.plot([n["x"],c["x"]],[c["y"]]*2,color=col,lw=1.6 if ch else .55,zorder=2,
                    solid_capstyle="round")
            if ch:
                ev[0 if b2 else 1]+=1
                ax.annotate("gain" if b2 else "loss",xy=((n["x"]+c["x"])/2,c["y"]),
                            xytext=(0,3.0),textcoords="offset points",ha="center",
                            va="bottom",fontsize=4.3,color=col,fontweight="700")
            rec(c)
        p=pres.get(n["name"])
        if n["collapse"]:
            L=leaves(n); xm=max(x["x"] for x in L)
            ax.fill([n["x"],xm,xm],[n["y"],n["y"]-.36,n["y"]+.36],
                    color=PRES if p else ABS,alpha=.30,lw=0,zorder=1)
            if collapse_all:
                g=[x["name"].split("_")[0] for x in L]; ao=all(x["name"] in osc for x in L)
                ax.text(xm+.012,n["y"],f"  {g[0]} \u2026 {g[-1]}  ({len(L)})",va="center",ha="left",
                        fontsize=4.2,style="italic",color="#8A1416" if ao else INK3,
                        fontweight="700" if ao else "400")
            else:
                ax.text(xm+.012,n["y"],f"  {len(L)} taxa",va="center",ha="left",
                        fontsize=4.2,color=INK3,style="italic")
        if p is not None:
            ax.scatter([n["x"]],[n["y"]],s=12 if not n["ch"] else 8,
                       color=PRES if p else ABS,zorder=4,linewidths=0,
                       marker="o" if n["ch"] else "s")
        if not n["ch"]:
            ax.text(n["x"]+.012,n["y"],n["name"].replace("_"," "),va="center",ha="left",
                    fontsize=4.2,style="italic",
                    color="#8A1416" if n["name"] in osc else INK,
                    fontweight="700" if n["name"] in osc else "400")
    rec(root)
    ax.set_ylim(-.9,N-.1); ax.invert_yaxis()
    ax.set_xlim(-.02,max(n["x"] for n in order)*1.40); ax.axis("off")
    ax.set_title(title,loc="left",fontsize=6.5,color=INK,fontweight="600",pad=3)
    ax.legend(handles=[Line2D([],[],marker="s",color=PRES,ls="",ms=3,label="site present"),
                       Line2D([],[],marker="s",color=ABS,ls="",ms=3,label="absent"),
                       Line2D([],[],color=GAIN,lw=1.5,label="gain"),
                       Line2D([],[],color=LOSS,lw=1.5,label="loss"),
                       Line2D([],[],color="none",label="bold = oscine")],
              frameon=False,fontsize=4.5,loc="lower left",ncol=5,
              bbox_to_anchor=(0,-.045),handletextpad=.3,columnspacing=.7)
    fig.subplots_adjust(left=.01,right=.99,top=.95,bottom=.05)
    for e in ("png","pdf"):
        fig.savefig(f"{stem}.{e}",dpi=240,facecolor="white",bbox_inches="tight",pad_inches=.015)
    plt.close(fig)
    return N,ev[0],ev[1]

if __name__=="__main__":
    OSC={l.strip() for l in open("data74/osc74.txt") if l.strip()}
    tree=open("data74/tree74.nw").read().strip()
    # load the gapped 74-way, keeping ancestral rows
    blocks=[];cur=[]
    for l in open("data74/ar_windows_gap_74.maf"):
        if l.startswith("# hal") and not os.path.exists("data74/tree74_gap.nw"):
            open("data74/tree74_gap.nw","w").write(l[6:].strip()+"\n")
        if l.startswith("a"):
            if cur: blocks.append(cur)
            cur=[]
        elif l.startswith("s\t"):
            f=l.rstrip("\n").split("\t"); cur.append((f[1],int(f[2]),int(f[3]),f[6]))
    if cur: blocks.append(cur)
    print(f"  blocks {len(blocks):,}")
    sites=[]
    for bed,mo in ((f"{F}/tracks/mafb_1p5mb/arefull_sites_MAFB_lsd2.bed","full"),
                   (f"{F}/tracks/mafb_1p5mb/jasparAR_sites_MAFB_lsd2.bed","jaspar")):
        for l in open(bed):
            x=l.rstrip("\n").split("\t")
            sites.append(dict(motif=mo,name=x[3],s=int(x[1]),e=int(x[2]),
                              cls=x[9],elem=x[10]))
    man=[]; counts=[]
    for st in sites:
        lo=LO[st["motif"]]; w0,w1=st["s"]-PAD,st["e"]+PAD
        seq=collections.defaultdict(str)
        for b in blocks:
            ref=[r for r in b if r[0].startswith("Lonchura_striata.")]
            if not ref: continue
            coord=[];p=ref[0][1]
            for c in ref[0][3]:
                if c=="-": coord.append(None)
                else: coord.append(p); p+=1
            real=[c for c in coord if c is not None]
            if not real or real[-1]<w0 or real[0]>=w1: continue
            keep=[i for i,c in enumerate(coord) if c is None or w0<=c<w1]
            if not keep: continue
            for sp,_,_,t in b:
                seq[sp.split(".")[0]]+="".join(t[i] for i in keep if i<len(t)).replace("-","")
        pres={g:(best(s,lo)>=TH[st["motif"]]) for g,s in seq.items() if len(s)>=len(lo)}
        ext=[g for g in pres if not g.startswith("Anc")]
        o=[g for g in ext if g in OSC]; u=[g for g in ext if g not in OSC]
        if len(o)<8 or len(u)<20: continue
        fo=sum(pres[g] for g in o)/len(o); fu=sum(pres[g] for g in u)/len(u)
        # branch events for EVERY covered site, so the dropped ones stay on record
        ng_,nl_=count_events(pres,tree)
        counts.append(dict(**{k:st[k] for k in ("motif","name","cls","elem")},
                           f_osc=round(fo,3),f_out=round(fu,3),gains=ng_,losses=nl_,
                           n_osc=len(o),n_out=len(u),
                           pass_freq=bool(fo>=0.6 and fu<=0.35)))
        if not (fo>=0.6 and fu<=0.35): continue          # candidate gains only
        if ng_>2: 
            print(f"    skip {st['motif']}_{st['name']:<34} {ng_} gains (turnover, not a gain)")
            continue
        nm=f"{st['motif']}_{st['name']}"
        N,ng,nl=draw(pres,tree,f"{per_site('gain_trees_ar', nm)}/{nm}",OSC,
                     f"{st['name']}  ({'ARE-full' if st['motif']=='full' else 'JASPAR AR'})"
                     f"  chr20:{st['s']+1:,}-{st['e']:,}   oscine {sum(pres[g] for g in o)}/{len(o)}"
                     f"  outgroup {sum(pres[g] for g in u)}/{len(u)}")
        man.append(dict(file=nm+".pdf",**{k:st[k] for k in ("motif","name","cls","elem")},
                        f_osc=round(fo,3),f_out=round(fu,3),tips=N,gains=ng,losses=nl))
        print(f"    {nm:<44} osc {fo:.2f} out {fu:.2f}   {ng} gains, {nl} losses")
    json.dump(man,open("trees_gain/manifest.json","w"),indent=1)
    json.dump(counts,open("trees_gain/counts.json","w"),indent=1)
    print(f"  trees drawn: {len(man)}")
