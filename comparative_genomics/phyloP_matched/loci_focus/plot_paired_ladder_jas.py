#!/usr/bin/env python
"""Paired AR-vs-enhancer relationships in the oscine-ladder style."""
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
sys.path.insert(0,"/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus")
from jaspar_core import pairs, group
F="/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus"
U=os.environ.get("AR_UNIT","instance")
rng=np.random.default_rng(0)
INK="#171B22"; INK3="#6E7887"; RULE="#DDE1E7"; RA="#e31a1c"; G4="#2C6E7F"
def ci(d,n=4000):
    b=np.array([np.mean(rng.choice(d,len(d))) for _ in range(n)]); return np.percentile(b,[2.5,97.5])
def pv(d,n=20000):
    b=np.array([np.mean(rng.choice(d,len(d))) for _ in range(n)])
    return max(2*min((b<=0).mean(),(b>=0).mean()),1/n)
PO,PW=pairs("subtree222"),pairs("wholetree")
ROWS=[("Glut-RA, genome-wide","ra",None,RA),("GABA-4, genome-wide","g4",None,G4),
      ("Glut-RA at MAFB","ra","MAFB",RA),("GABA-4 at MAFB","g4","MAFB",G4),
]
fig,axes=plt.subplots(1,2,figsize=(3.56,1.40),gridspec_kw={"width_ratios":[1.25,1]})
fig.patch.set_facecolor("white")
ax=axes[0]; ys=[]; LABS=[]; XS=[0.0]
for i,(lab,c,loc,col) in enumerate(ROWS):
    g=group(PO,c,loc); d=g[:,0]-g[:,1]
    lo,hi=ci(d); y=len(ROWS)-1-i; ys.append(y); XS+=[lo,hi]
    p=pv(d); LABS.append(f"{lab}   (n={len(d):,})")
    ax.plot([lo,hi],[y,y],color=col,lw=1.25,solid_capstyle="round",zorder=2)
    ax.scatter([d.mean()],[y],s=14.9,color=col,zorder=3,linewidths=0)
    ax.text(hi+.012,y,f"{d.mean():+.3f}{'*' if p<.05 else ''}",va="center",ha="left",
            fontsize=5,color=col,fontweight="600" if p<.05 else "400")
ax.axvline(0,color=INK,lw=0.50,zorder=1)
ax.axhline(1.5,color=RULE,lw=0.50,ls=(0,(3,3)),zorder=1)
ax.set_yticks(ys); ax.set_yticklabels(LABS,fontsize=5,color=INK)
ax.set_xlabel("motif − enhancer (oscine)",fontsize=5,color=INK)
ax.set_title("A   motif vs its enhancer",loc="left",
             fontsize=7,color=INK,fontweight="600",pad=10)
_l,_h=min(XS),max(XS); _s=_h-_l
ax.set_xlim(_l-.07*_s,_h+.30*_s)   # headroom for the right-aligned label
ax=axes[1]
PTS=[("Glut-RA genome-wide","ra",None,RA,"o",(9,10),"left"),
     ("GABA-4 genome-wide","g4",None,G4,"o",(-9,-20),"right"),
     ("Glut-RA at MAFB","ra","MAFB",RA,"D",(-6,16),"center"),
     ("GABA-4 at MAFB","g4","MAFB",G4,"D",(12,-16),"left"),
]
BX=[0.0]; BY=[0.0]
for lab,c,loc,col,mk,off,ha in PTS:
    gw=group(PW,c,loc); go=group(PO,c,loc)
    x=(gw[:,0]-gw[:,1]); y=(go[:,0]-go[:,1])
    xl,xh=ci(x); yl,yh=ci(y); BX+=[xl,xh]; BY+=[yl,yh]
    ax.plot([xl,xh],[y.mean()]*2,color=col,lw=0.70,alpha=.5,zorder=2)
    ax.plot([x.mean()]*2,[yl,yh],color=col,lw=0.70,alpha=.5,zorder=2)
    ax.scatter([x.mean()],[y.mean()],s=18.7 if mk!="o" else 14.9,color=col,marker=mk,zorder=3,linewidths=0)

ax.axhline(0,color=RULE,lw=0.50,zorder=1); ax.axvline(0,color=RULE,lw=0.50,zorder=1)
ax.set_xlabel("whole-tree phyloP",fontsize=5,color=INK)
ax.set_ylabel("oscine phyloP",fontsize=6,color=INK)
ax.set_title("B   two trees",loc="left",fontsize=7,color=INK,fontweight="600",pad=10)
_bl,_bh=min(BX),max(BX); _bs=_bh-_bl
_cl,_ch=min(BY),max(BY); _cs=_ch-_cl
ax.set_xlim(_bl-.08*_bs,_bh+.08*_bs); ax.set_ylim(_cl-.10*_cs,_ch+.30*_cs)
from matplotlib.lines import Line2D
ax.legend(handles=[Line2D([],[],marker="o",color=INK3,ls="",ms=4.0,label="genome-wide"),
                   Line2D([],[],marker="D",color=INK3,ls="",ms=4.0,label="at MAFB"),
                   Line2D([],[],marker="o",color=RA,ls="",ms=4.5,label="Glut-RA DARs"),
                   Line2D([],[],marker="o",color=G4,ls="",ms=4.5,label="GABA-4 DARs")],
          frameon=False,fontsize=5,loc="upper left",ncol=1,handletextpad=.3,
          labelspacing=.25,borderaxespad=.1)
for ax in axes:
    ax.margins(x=.03,y=.05)
    ax.locator_params(axis="x",nbins=5)
    for sp in ("top","right"): ax.spines[sp].set_visible(False)
    for sp in ("left","bottom"): ax.spines[sp].set_color(RULE)
    ax.tick_params(colors=INK3,labelsize=5)
fig.subplots_adjust(left=.20,right=.985,top=.90,bottom=.185,wspace=.16)
fig.savefig(f"{F}/results/compare_gw_mafb1p5mb/fig_ar_paired_ladder_{U}.png",dpi=200,facecolor="white",bbox_inches="tight",pad_inches=.01)
fig.savefig(f"{F}/results/compare_gw_mafb1p5mb/fig_ar_paired_ladder_{U}.pdf",dpi=200,facecolor="white",bbox_inches="tight",pad_inches=.01)
print(f"  wrote fig_ar_paired_ladder_{U}.png/.pdf")
