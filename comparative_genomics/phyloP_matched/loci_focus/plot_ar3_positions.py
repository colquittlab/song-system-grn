#!/usr/bin/env python
"""Per-position phyloP across the three AR site types.

One representative hit per merged site, so the two strand reports of a
palindromic ARE are not folded onto each other. Position 1 is always the
motif's first column.
"""
import numpy as np, collections, matplotlib
matplotlib.use("Agg")
matplotlib.rcParams.update({"font.family":"sans-serif",
    "font.sans-serif":["Arial","Liberation Sans","DejaVu Sans"],
    "mathtext.fontset":"custom","mathtext.rm":"Arial",
    "pdf.fonttype":42,"ps.fonttype":42,
    "font.size":5,"axes.labelsize":6,"axes.titlesize":7,
    "xtick.labelsize":5,"ytick.labelsize":5,"legend.fontsize":5,
    "axes.linewidth":.4,"xtick.major.width":.4,"ytick.major.width":.4,
    "xtick.major.size":1.5,"ytick.major.size":1.5,"xtick.major.pad":1,
    "ytick.major.pad":1,"axes.labelpad":1,"axes.titlepad":2})
import matplotlib.pyplot as plt
import pandas as pd, logomaker
F="/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus"
HM="/private/home/bcolquit/group/motifs/homer/ra_pos_vs_neg/homerResults_manual"
rng=np.random.default_rng(0)
INK="#171B22"; INK3="#6E7887"; RULE="#DDE1E7"; RA="#e31a1c"; G4="#2C6E7F"
MOT=[("full","ARE-full",12),("half","ARE-half",8),("jaspar","JASPAR AR",16)]
def homer_pwm(f):
    return np.array([[float(x) for x in l.split()] for l in open(f) if not l.startswith(">")])
def meme_pwm(f,W):
    blk=open(f).read().split("letter-probability")[1].split("\n")[1:]
    return np.array([[float(x) for x in l.split()] for l in blk if len(l.split())==4][:W])
PWM={"full":homer_pwm(f"{HM}/motif2.motif"),"half":homer_pwm(f"{HM}/motif4.motif"),
     "jaspar":meme_pwm(f"{F}/ar_jaspar.meme",16)}
# load per-position values
V={}
for tree in ("subtree222","wholetree"):
    d=collections.defaultdict(list)
    for l in open(f"{F}/{tree}_ar3_pos.tab"):
        p=l.split("\t")
        if int(p[2])==0: continue
        m,s,loc,e,site,k=p[0].split("|")
        d[(m,s,loc,int(k))].append(float(p[5]))
    V[tree]=d
def series(tree,m,s,W,loc=None):
    mu=[];lo=[];hi=[];n=0
    for k in range(W):
        v=np.array([x for (mm,ss,ll,kk),vals in V[tree].items()
                    if mm==m and ss==s and kk==k and (loc is None or ll==loc) for x in vals])
        n=len(v); mu.append(v.mean())
        b=[np.mean(rng.choice(v,len(v))) for _ in range(2000)]
        a,c=np.percentile(b,[2.5,97.5]); lo.append(a); hi.append(c)
    return np.array(mu),np.array(lo),np.array(hi),n
fig,axes=plt.subplots(4,3,figsize=(5.3,3.72),
                      gridspec_kw={"height_ratios":[1,1.55,1.55,1.55],
                                   "width_ratios":[12,8,16],"hspace":.58,"wspace":.26})
fig.patch.set_facecolor("white")
PAN="ABC"
for c,(mk,ml,W) in enumerate(MOT):
    ic=np.array([2+(p*np.log2(np.clip(p,1e-9,None))).sum() for p in PWM[mk]])
    ax=axes[0][c]
    logomaker.Logo(pd.DataFrame(PWM[mk]*ic[:,None],columns=list("ACGT")),
                   ax=ax,color_scheme="classic",show_spines=False,font_name="Arial")
    ax.set_ylim(0,2); ax.set_yticks([0,1,2])
    if c==0: ax.set_ylabel("bits",fontsize=6,color=INK)
    ax.set_title(f"{PAN[c]}   {ml} ({W} bp)",loc="left",fontsize=7,color=INK,fontweight="600",pad=2)
    ax.tick_params(colors=INK3,labelsize=5); ax.set_xticks([])
    for r,(tree,tl,loc) in enumerate((("subtree222","oscine",None),
                                      ("wholetree","whole tree",None),
                                      ("subtree222","oscine, MAFB","MAFB"))):
        ax=axes[r+1][c]
        for s,col,lab in (("ra",RA,"Glut-RA"),("g4",G4,"GABA-4")):
            mu,lo,hi,n=series(tree,mk,s,W,loc); x=np.arange(W)
            ax.fill_between(x,lo,hi,color=col,alpha=.16,linewidth=0)
            ax.plot(x,mu,color=col,lw=.9,marker="o",ms=1.9,label=f"{lab} ({n:,})")
        ax.axhline(0,color=RULE,lw=.5)
        if c==0: ax.set_ylabel(f"{tl} phyloP",fontsize=6,color=INK)
        ax.set_xticks(range(W)); ax.set_xticklabels(range(1,W+1),fontsize=4.3)
        if r==2: ax.set_xlabel("motif position",fontsize=6,color=INK)
        ax.tick_params(colors=INK3,labelsize=5)
        ax.margins(x=.04,y=.10); ax.locator_params(axis="y",nbins=4)
        ax.legend(frameon=False,fontsize=4.3,loc="upper center",ncol=2,
                  bbox_to_anchor=(.5,1.30),columnspacing=.7,handletextpad=.25,
                  handlelength=1.0)
        for sp in ("top","right"): ax.spines[sp].set_visible(False)
        for sp in ("left","bottom"): ax.spines[sp].set_color(RULE)
        if loc=="MAFB":
            for sp in ("left","bottom"): ax.spines[sp].set_color("#C8A3A3")
fig.subplots_adjust(left=.075,right=.99,top=.945,bottom=.105)
fig.text(.075,.004,"Bottom row: MAFB locus only \u2014 13\u201352 sites per curve. "
         "Shown for comparability with the genome-wide rows; not a basis for per-position inference.",
         fontsize=4.4,color="#9A6B6B",va="bottom",ha="left",style="italic")
for e in ("png","pdf"):
    fig.savefig(f"{F}/results/compare_gw_mafb1p5mb/fig_ar3_positions.{e}",dpi=200,facecolor="white",
                bbox_inches="tight",pad_inches=.01)
print("  wrote fig_ar3_positions")
# high-IC vs low-IC contrast, per motif and class
def pv(d,n=20000):
    b=np.array([np.mean(rng.choice(d,len(d))) for _ in range(n)])
    return max(2*min((b<=0).mean(),(b>=0).mean()),1/n)
print("\n  high-IC minus low-IC columns (IC >= 1.0 bits), oscine tree:\n")
for mk,ml,W in MOT:
    ic=np.array([2+(p*np.log2(np.clip(p,1e-9,None))).sum() for p in PWM[mk]])
    HI=ic>=1.0
    if HI.all() or not HI.any():
        print(f"    {ml:<11} no IC split (all columns {'high' if HI.all() else 'low'})"); continue
    inst=collections.defaultdict(dict)
    for l in open(f"{F}/subtree222_ar3_pos.tab"):
        p=l.split("\t")
        if int(p[2])==0: continue
        m,s,loc,e,site,k=p[0].split("|")
        if m!=mk: continue
        inst[(s,e,site)][int(k)]=float(p[5])
    for s,lab in (("ra","Glut-RA"),("g4","GABA-4")):
        M=np.array([[v[k] for k in range(W)] for (ss,e,si),v in inst.items()
                    if ss==s and len(v)==W])
        d=M[:,HI].mean(1)-M[:,~HI].mean(1)
        print(f"    {ml:<11} {lab:<8} high {M[:,HI].mean():+.3f}  low {M[:,~HI].mean():+.3f}"
              f"  diff {d.mean():+.3f}  p={pv(d):.4f}  n={len(M):,}  ({HI.sum()}/{W} cols high)")
