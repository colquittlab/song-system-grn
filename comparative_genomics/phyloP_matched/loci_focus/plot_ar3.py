#!/usr/bin/env python
"""Three AR motifs in parallel: levels at the merged-site unit, and the
difference-of-differences across all three units."""
import numpy as np, sys, matplotlib
matplotlib.use("Agg")
matplotlib.rcParams.update({"font.family":"sans-serif",
    "font.sans-serif":["Arial","Liberation Sans","DejaVu Sans"],
    "mathtext.fontset":"custom","mathtext.rm":"Arial",
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
sys.path.insert(0,F)
import ar_core as C
rng=np.random.default_rng(0)
INK="#171B22"; INK3="#6E7887"; RULE="#DDE1E7"; RA="#e31a1c"; G4="#2C6E7F"
UC={"instance":"#B9C0CB","site":"#6E7887","element":"#2B313B"}
MOT=[("full","ARE-full"),("half","ARE-half"),("jaspar","JASPAR AR")]
def ci(v,n=4000):
    b=np.array([np.mean(rng.choice(v,len(v))) for _ in range(n)]); return np.percentile(b,[2.5,97.5])
def pv(d,n=20000):
    b=np.array([np.mean(rng.choice(d,len(d))) for _ in range(n)])
    return max(2*min((b<=0).mean(),(b>=0).mean()),1/n)

# ---------- figure 1: levels, merged-site unit ----------
ROWS=[]
for mk,ml in MOT:
    for lab,c,loc,col in (("Glut-RA, gw","ra",None,RA),("GABA-4, gw","g4",None,G4),
                          ("Glut-RA, MAFB","ra","MAFB",RA),("GABA-4, MAFB","g4","MAFB",G4)):
        ROWS.append((mk,ml,lab,c,loc,col))
P={t:{m:C.pairs(t,m,"site") for m,_ in MOT} for t in ("subtree222","wholetree")}
fig,axes=plt.subplots(1,2,figsize=(4.05,2.62),sharey=True); fig.patch.set_facecolor("white")
for ax,(tree,ttl) in zip(axes,(("subtree222","A   oscine subtree"),
                               ("wholetree","B   whole tree (pan-avian)"))):
    ys=[];LABS=[];XS=[0.0]
    for i,(mk,ml,lab,c,loc,col) in enumerate(ROWS):
        g=C.group(P[tree][mk],c,loc); y=len(ROWS)-1-i; ys.append(y)
        LABS.append(f"{lab}  (n={len(g):,})")
        m,r=g[:,0].mean(),g[:,1].mean(); d=g[:,0]-g[:,1]; p=pv(d)
        ml_,mh=ci(g[:,0]); rl,rh=ci(g[:,1]); XS+=[ml_,mh,rl,rh]
        ax.annotate("",xy=(m,y),xytext=(r,y),
                    arrowprops=dict(arrowstyle="-|>,head_width=.17,head_length=.38",
                                    color=col,lw=.85,alpha=.85 if p<.05 else .35,shrinkA=0,shrinkB=0))
        ax.plot([rl,rh],[y-.17]*2,color=INK3,lw=.65,alpha=.65,solid_capstyle="round")
        ax.plot([ml_,mh],[y+.17]*2,color=col,lw=.65,alpha=.65,solid_capstyle="round")
        ax.scatter([r],[y-.17],s=8.2,facecolor="white",edgecolor=INK3,linewidths=.7,zorder=4)
        ax.scatter([m],[y+.17],s=12.1,color=col,zorder=5,linewidths=0)
        ax.text(.985,y,f"{d.mean():+.3f}{'*' if p<.05 else ''}",va="center",ha="right",
                fontsize=5,color=col,fontweight="600" if p<.05 else "400",
                transform=ax.get_yaxis_transform())
    ax.axvline(0,color=INK,lw=.5,zorder=1)
    for b in (3.5,7.5): ax.axhline(b,color=RULE,lw=.5,ls=(0,(3,3)),zorder=1)
    ax.set_yticks(ys); ax.set_yticklabels(LABS,fontsize=5,color=INK)
    ax.set_xlabel("phyloP",fontsize=6,color=INK)
    ax.set_title(ttl,loc="left",fontsize=7,color=INK,fontweight="600",pad=3)
    ax.tick_params(colors=INK3,labelsize=5); ax.margins(y=.025)
    lo,hi=min(XS),max(XS); s=hi-lo; ax.set_xlim(lo-.07*s,hi+.34*s)
    ax.locator_params(axis="x",nbins=4)
    for sp in ("top","right"): ax.spines[sp].set_visible(False)
    for sp in ("left","bottom"): ax.spines[sp].set_color(RULE)
for j,(mk,ml) in enumerate(MOT):
    axes[0].text(-.72,11-4*j-1.5,ml,transform=axes[0].get_yaxis_transform(),rotation=90,
                 va="center",ha="center",fontsize=5.5,color=INK,fontweight="700")
fig.legend(handles=[Line2D([],[],marker="o",mfc="white",mec=INK3,mew=1.5,color="none",ms=4,
                           label="rest of the enhancer"),
                    Line2D([],[],marker="o",color=INK3,ls="",ms=4,label="the motif"),
                    Line2D([],[],color=INK3,lw=1.1,label="arrow = gap, solid if p < 0.05")],
           frameon=False,fontsize=5,loc="lower center",ncol=3,bbox_to_anchor=(.55,-.055))
fig.subplots_adjust(left=.345,right=.975,top=.945,bottom=.145,wspace=.14)
for e in ("png","pdf"):
    fig.savefig(f"{F}/results/compare_gw_mafb1p5mb/fig_ar3_fold.{e}",dpi=200,facecolor="white",bbox_inches="tight",pad_inches=.01)
print("  wrote fig_ar3_fold")

# ---------- figure 2: DoD across all three units ----------
SC=[(None,"genome-wide"),("MAFB","at MAFB")]
fig,axes=plt.subplots(1,2,figsize=(4.05,1.78),sharey=True); fig.patch.set_facecolor("white")
rows=[(mk,ml,loc,sl) for mk,ml in MOT for loc,sl in SC]
for ax,(tree,ttl) in zip(axes,(("subtree222","A   oscine subtree"),
                               ("wholetree","B   whole tree (pan-avian)"))):
    ys=[];LABS=[];XS=[0.0]
    for i,(mk,ml,loc,sl) in enumerate(rows):
        y=len(rows)-1-i; ys.append(y); LABS.append(f"{ml}, {sl}")
        for j,u in enumerate(("instance","site","element")):
            Pu=C.pairs(tree,mk,u)
            ga=C.group(Pu,"ra",loc); gb=C.group(Pu,"g4",loc)
            if len(ga)<5 or len(gb)<5: continue
            a=ga[:,0]-ga[:,1]; b=gb[:,0]-gb[:,1]; d=a.mean()-b.mean()
            bs=np.array([np.mean(rng.choice(a,len(a)))-np.mean(rng.choice(b,len(b)))
                         for _ in range(3000)])
            lo,hi=np.percentile(bs,[2.5,97.5]); XS+=[lo,hi]
            yy=y+(1-j)*.24
            ax.plot([lo,hi],[yy]*2,color=UC[u],lw=.7,solid_capstyle="round",zorder=2)
            ax.scatter([d],[yy],s=13 if u=="site" else 8.5,color=UC[u],zorder=3,linewidths=0,
                       marker="D" if u=="site" else "o")
    ax.axvline(0,color=INK,lw=.5,zorder=1)
    for b in (1.5,3.5): ax.axhline(b,color=RULE,lw=.5,ls=(0,(3,3)),zorder=1)
    ax.set_yticks(ys); ax.set_yticklabels(LABS,fontsize=5,color=INK)
    ax.set_xlabel("Glut-RA − GABA-4  (motif − enhancer)",fontsize=6,color=INK)
    ax.set_title(ttl,loc="left",fontsize=7,color=INK,fontweight="600",pad=3)
    ax.tick_params(colors=INK3,labelsize=5); ax.margins(y=.06)
    lo,hi=min(XS),max(XS); s=hi-lo; ax.set_xlim(lo-.08*s,hi+.08*s)
    ax.locator_params(axis="x",nbins=4)
    for sp in ("top","right"): ax.spines[sp].set_visible(False)
    for sp in ("left","bottom"): ax.spines[sp].set_color(RULE)
axes[0].legend(handles=[Line2D([],[],marker="o",color=UC["instance"],ls="",ms=3,label="instance"),
                        Line2D([],[],marker="D",color=UC["site"],ls="",ms=3.3,label="merged site"),
                        Line2D([],[],marker="o",color=UC["element"],ls="",ms=3,label="element")],
               frameon=True,framealpha=.9,facecolor="white",edgecolor="none",
               fontsize=4.6,loc="lower left",borderaxespad=.12).get_frame().set_linewidth(0)
fig.subplots_adjust(left=.235,right=.985,top=.90,bottom=.185,wspace=.10)
for e in ("png","pdf"):
    fig.savefig(f"{F}/results/compare_gw_mafb1p5mb/fig_ar3_dod.{e}",dpi=200,facecolor="white",bbox_inches="tight",pad_inches=.01)
print("  wrote fig_ar3_dod")
