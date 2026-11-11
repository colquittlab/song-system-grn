#!/usr/bin/env python
"""AR motif vs the enhancer it sits in, as levels and the gap between them.

Shown as a difference rather than a ratio: phyloP crosses zero in several of
these groups, so a ratio would be undefined in 8 of 12 cells and dominated by
the denominator in the rest.
"""
import numpy as np, sys, os, matplotlib
matplotlib.use("Agg")
matplotlib.rcParams.update({"font.family":"sans-serif",
    "font.sans-serif":["Arial","Liberation Sans","DejaVu Sans"],
    "mathtext.fontset":"custom","mathtext.rm":"Arial",
    "mathtext.it":"Arial:italic","mathtext.bf":"Arial:bold",
    "pdf.fonttype":42,"ps.fonttype":42,
    "font.size":5,"axes.labelsize":6,"axes.titlesize":7,
    "xtick.labelsize":5,"ytick.labelsize":5,"legend.fontsize":5,
    "axes.linewidth":.4,"xtick.major.width":.4,"ytick.major.width":.4,
    "xtick.major.size":1.5,"ytick.major.size":1.5,"xtick.major.pad":1,
    "ytick.major.pad":1,"axes.labelpad":1,"axes.titlepad":2,
    "legend.handletextpad":.4,"legend.borderpad":.2,"legend.labelspacing":.3})
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
sys.path.insert(0,"/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus")
from jaspar_core import pairs, group
F="/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus"
U=os.environ.get("AR_UNIT","instance")
rng=np.random.default_rng(0)
INK="#171B22"; INK3="#6E7887"; RULE="#DDE1E7"; RA="#e31a1c"; G4="#2C6E7F"
def ci(v,n=4000):
    b=np.array([np.mean(rng.choice(v,len(v))) for _ in range(n)]); return np.percentile(b,[2.5,97.5])
def pv(d,n=20000):
    b=np.array([np.mean(rng.choice(d,len(d))) for _ in range(n)])
    return max(2*min((b<=0).mean(),(b>=0).mean()),1/n)
ROWS=[("Glut-RA, genome-wide","ra",None,RA),("GABA-4, genome-wide","g4",None,G4),
      ("Glut-RA at MAFB","ra","MAFB",RA),("GABA-4 at MAFB","g4","MAFB",G4),
]
P={"subtree222":pairs("subtree222"),"wholetree":pairs("wholetree")}
fig,axes=plt.subplots(1,2,figsize=(3.51,1.40),sharey=True)
fig.patch.set_facecolor("white")
for ax,(tree,ttl) in zip(axes,(("subtree222",f"A   oscine subtree"),
                               ("wholetree","B   whole tree (pan-avian)"))):
    ys=[];LABS=[];XS=[]
    for i,(lab,c,loc,col) in enumerate(ROWS):
        g=group(P[tree],c,loc); y=len(ROWS)-1-i; ys.append(y)
        LABS.append(f"{lab}   (n={len(g):,})")
        m,r=g[:,0].mean(),g[:,1].mean(); d=g[:,0]-g[:,1]; p=pv(d)
        ml,mh=ci(g[:,0]); rl,rh=ci(g[:,1]); XS+=[ml,mh,rl,rh,0.0]
        # each mean sits on its own interval: rest below the row, motif above
        yr,ym=y-.16,y+.16
        # arrow runs horizontally down the row centre, between the two point rows
        ax.annotate("",xy=(m,y),xytext=(r,y),
                    arrowprops=dict(arrowstyle="-|>,head_width=.17,head_length=.38",
                                    color=col,lw=0.85,alpha=.85 if p<.05 else .35,
                                    shrinkA=0,shrinkB=0))
        ax.plot([rl,rh],[yr]*2,color=INK3,lw=0.65,alpha=.65,solid_capstyle="round")
        ax.plot([ml,mh],[ym]*2,color=col,lw=0.65,alpha=.65,solid_capstyle="round")
        ax.scatter([r],[yr],s=8.2,facecolor="white",edgecolor=INK3,linewidths=.7,zorder=4)
        ax.scatter([m],[ym],s=12.1,color=col,zorder=5,linewidths=0)
        ax.text(.985,y,f"{d.mean():+.3f}{'*' if p<.05 else ''}",va="center",ha="right",
                fontsize=5,color=col,fontweight="600" if p<.05 else "400",
                transform=ax.get_yaxis_transform())
    ax.axvline(0,color=INK,lw=0.50,zorder=1)
    ax.axhline(1.5,color=RULE,lw=0.50,ls=(0,(3,3)),zorder=1)
    ax.set_yticks(ys); ax.set_yticklabels(LABS,fontsize=5,color=INK)
    ax.set_xlabel("phyloP",fontsize=6,color=INK)
    ax.set_title(ttl,loc="left",fontsize=7,color=INK,fontweight="600",pad=10)
    ax.tick_params(colors=INK3,labelsize=5)
    ax.margins(x=.03,y=.05)
    ax.locator_params(axis="x",nbins=5)
    for sp in ("top","right"): ax.spines[sp].set_visible(False)
    for sp in ("left","bottom"): ax.spines[sp].set_color(RULE)
    lo,hi=min(XS),max(XS); span=hi-lo
    ax.set_xlim(lo-.07*span, hi+.34*span)      # headroom for the right-aligned gap label
fig.legend(handles=[Line2D([],[],marker="o",mfc="white",mec=INK3,mew=1.5,color="none",
                           ms=4.0,label="rest of the enhancer"),
                    Line2D([],[],marker="o",color=INK3,ls="",ms=4.0,label="the AR motif"),
                    Line2D([],[],color=INK3,lw=1.10,label="arrow = the gap, solid if p < 0.05")],
           frameon=False,fontsize=5,loc="lower center",ncol=3,bbox_to_anchor=(.55,-.10))
fig.subplots_adjust(left=.215,right=.975,top=.88,bottom=.30,wspace=.16)
fig.savefig(f"{F}/results/compare_gw_mafb1p5mb/fig_ar_fold_{U}.png",dpi=200,facecolor="white",bbox_inches="tight",pad_inches=.01)
fig.savefig(f"{F}/results/compare_gw_mafb1p5mb/fig_ar_fold_{U}.pdf",dpi=200,facecolor="white",bbox_inches="tight",pad_inches=.01)
print(f"  wrote fig_ar_fold_{U}.png/.pdf")
for tree in ("subtree222","wholetree"):
    for lab,c,loc,_ in ROWS:
        g=group(P[tree],c,loc)
        print(f"  {tree:<11}{lab:<24} rest {g[:,1].mean():+.3f} -> motif {g[:,0].mean():+.3f}"
              f"   gap {(g[:,0]-g[:,1]).mean():+.3f}")
