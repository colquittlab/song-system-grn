#!/usr/bin/env python
"""AR-site constraint on the oscine lineage, full and half response elements."""
import numpy as np, matplotlib
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
F="/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus"
rng=np.random.default_rng(0)
INK="#171B22"; INK2="#454E5C"; INK3="#6E7887"; RULE="#DDE1E7"
RA="#e31a1c"; G4="#2C6E7F"; GREY="#9AA3AF"
def load(t,f):
    d={}
    for l in open(f"{F}/{t}_{f}.tab"):
        p=l.split("\t"); d[p[0]]=(int(p[1]),int(p[2]),float(p[5]))
    return d
def grab(f,cls=None,loc=None,tree="subtree222"):
    D=load(tree,f); out=[]
    for k,v in D.items():
        if v[1]/max(v[0],1)<.5: continue
        p=k.split("|")
        if cls and p[1]!=cls: continue
        if loc and p[2]!=loc: continue
        out.append(v[2])
    return np.array(out)
def ci(v,n=4000):
    return np.percentile([np.mean(rng.choice(v,len(v))) for _ in range(n)],[2.5,97.5])
# label, file, class, locus, colour, marker  (o = full site / other, ^ = half site)
ROWS=[
 ("non-motif bases","nm2",None,None,GREY,"o"),
 ("all motifs","m2",None,None,GREY,"o"),
 ("AR, GABA-4","motif2","g4",None,G4,"o"),
 ("AR, Glut-RA","motif2","ra",None,RA,"o"),
 ("AR at MAFB, GABA-4","motif2","g4","MAFB",G4,"o"),
 ("AR at MAFB, Glut-RA","motif2","ra","MAFB",RA,"o"),
]
fig,axes=plt.subplots(1,2,figsize=(3.51,1.32),gridspec_kw={"width_ratios":[1.35,1]})
fig.patch.set_facecolor("white")

ax=axes[0]; ys=[]; LABS=[]
for i,(lab,f,c,loc,col,mk) in enumerate(ROWS):
    v=grab(f,c,loc); lo,hi=ci(v); y=len(ROWS)-1-i; ys.append(y)
    LABS.append(f"{lab}   (n={len(v):,})")
    ax.plot([lo,hi],[y,y],color=col,lw=1.20,solid_capstyle="round",zorder=2)
    ax.scatter([v.mean()],[y],s=16.0 if mk=="^" else 13.8,color=col,marker=mk,zorder=3,linewidths=0)
    ax.text(hi+.028,y,f"{v.mean():+.2f}",va="center",ha="left",fontsize=5,
            color=col,fontweight="600")
ax.axvline(0,color=INK,lw=0.50,zorder=1)
ax.axhline(1.5,color=RULE,lw=0.50,ls=(0,(3,3)),zorder=1)
ax.text(-.12,1.58,"gw",ha="left",va="bottom",fontsize=5,color=INK3,style="italic")
ax.text(-.12,1.42,"MAFB",ha="left",va="top",fontsize=5,color=INK3,style="italic")
ax.set_yticks(ys); ax.set_yticklabels(LABS,fontsize=5,color=INK)
ax.set_xlabel("oscine phyloP",fontsize=6,color=INK)
ax.set_title("A   constraint ladder",loc="left",
             fontsize=7,color=INK,fontweight="600",pad=12)
ax.set_xlim(-.14,.72); ax.set_ylim(-.6,len(ROWS)-.4)

ax=axes[1]
PTS=[("all motifs","m2",None,None,GREY,"o",(7,-1)),
     ("RA gw","motif2","ra",None,RA,"o",(1,-8)),
     ("G4 gw","motif2","g4",None,G4,"o",(5,2)),
     ("RA MAFB","motif2","ra","MAFB",RA,"D",(5,-1)),
     ("G4 MAFB","motif2","g4","MAFB",G4,"D",(-4,5))]
for lab,f,c,loc,col,mk,off in PTS:
    xs=grab(f,c,loc,tree="wholetree"); yv=grab(f,c,loc,tree="subtree222")
    xl,xh=ci(xs); yl,yh=ci(yv)
    ax.plot([xl,xh],[yv.mean()]*2,color=col,lw=0.75,alpha=.5,zorder=2)
    ax.plot([xs.mean()]*2,[yl,yh],color=col,lw=0.75,alpha=.5,zorder=2)
    ax.scatter([xs.mean()],[yv.mean()],s=19.2 if mk in "D^" else 14.9,color=col,marker=mk,
               zorder=3,linewidths=0)
    ax.annotate(lab,(xs.mean(),yv.mean()),textcoords="offset points",xytext=off,
                fontsize=5,color=col,fontweight="600" if "MAFB" in lab else "400")
ax.axhline(0,color=RULE,lw=0.50,zorder=1); ax.axvline(0,color=RULE,lw=0.50,zorder=1)
ax.set_xlabel("whole-tree phyloP",fontsize=5,color=INK)
ax.set_ylabel("oscine phyloP",fontsize=6,color=INK)
ax.set_title("B   two trees",loc="left",
             fontsize=7,color=INK,fontweight="600",pad=12)
ax.set_xlim(-.04,.78); ax.set_ylim(.14,.52)
for ax in axes:
    ax.margins(x=.03,y=.05)
    ax.locator_params(axis="x",nbins=5)
    for sp in ("top","right"): ax.spines[sp].set_visible(False)
    for sp in ("left","bottom"): ax.spines[sp].set_color(RULE)
    ax.tick_params(colors=INK3,labelsize=5)
fig.tight_layout(rect=[0,.075,1,1])
fig.savefig(f"{F}/results/compare_gw_mafb1p5mb/fig_ar_oscine.png",dpi=200,facecolor="white",bbox_inches="tight",pad_inches=.01)
fig.savefig(f"{F}/results/compare_gw_mafb1p5mb/fig_ar_oscine.pdf",dpi=200,facecolor="white",bbox_inches="tight",pad_inches=.01)
print("  wrote fig_ar_oscine.png")
for lab,f,c,loc,_,_ in ROWS:
    v=grab(f,c,loc); w=grab(f,c,loc,tree="wholetree")
    print(f"    {lab:<32} n={len(v):>5}  oscine {v.mean():+.3f}  whole {w.mean():+.3f}")
