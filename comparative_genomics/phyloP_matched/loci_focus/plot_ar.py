#!/usr/bin/env python
"""AR-site constraint at MAFB: full and half response elements, Glut-RA vs GABA-4."""
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
from matplotlib.patches import Patch
F="/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus"
rng=np.random.default_rng(0)
INK="#171B22"; INK2="#454E5C"; INK3="#6E7887"; RULE="#DDE1E7"
RA="#e31a1c"; G4="#2C6E7F"
def load(f):
    d={}
    for l in open(f):
        p=l.split("\t"); d[p[0]]=(int(p[1]),int(p[2]),float(p[5]))
    return d
def sel(f,s,loc=None):
    D=load(f)
    return np.array([v[2] for k,v in D.items() if v[1]/max(v[0],1)>=.5
                     and k.split("|")[1]==s and (loc is None or k.split("|")[2]==loc)])
def ci(v,n=4000):
    return np.percentile([np.mean(rng.choice(v,len(v))) for _ in range(n)],[2.5,97.5])
def perm(a,b,n=20000):
    o=a.mean()-b.mean(); pool=np.concatenate([a,b]); k=len(a)
    d=np.array([(lambda p:p[:k].mean()-p[k:].mean())(rng.permutation(pool)) for _ in range(n)])
    return max((np.abs(d)>=abs(o)).mean(),1/n)
PW=[("motif2","HOMER ARE-full"),("ar","JASPAR MA0007.3")]
fig,axes=plt.subplots(1,2,figsize=(3.19,1.52),gridspec_kw={"width_ratios":[1.5,1]})
fig.patch.set_facecolor("white")

ax=axes[0]; xs=[]; labs=[]
for i,(m,nm) in enumerate(PW):
    a=sel(f"{F}/subtree222_{m}.tab","ra","MAFB"); b=sel(f"{F}/subtree222_{m}.tab","g4","MAFB")
    for j,(v,c,lab) in enumerate(((a,RA,"Glut-RA"),(b,G4,"GABA-4"))):
        x=i*2.6+j
        ax.scatter(x+rng.normal(0,.075,len(v)),v,s=2.8,color=c,alpha=.40,linewidths=0,zorder=2)
        lo,hi=ci(v)
        ax.plot([x,x],[lo,hi],color=c,lw=1.20,zorder=3,solid_capstyle="round")
        ax.plot([x-.27,x+.27],[v.mean()]*2,color=c,lw=1.40,zorder=4,solid_capstyle="round")
        ax.text(x,hi+.12,f"{v.mean():+.2f}",ha="center",va="bottom",fontsize=5,
                color=c,fontweight="600")
        xs.append(x); labs.append(f"{lab}\nn={len(v)}")
    p=perm(a,b); y=max(a.max(),b.max())*.98
    ax.plot([i*2.6,i*2.6+1],[y,y],color=INK3,lw=0.50)
    ax.text(i*2.6+.5,y+.03,f"p = {p:.2f}",ha="center",va="bottom",fontsize=5,color=INK2)
    ax.text(i*2.6+.5,1.035,nm,ha="center",va="bottom",fontsize=5,color=INK,
            fontweight="600",transform=ax.get_xaxis_transform())
ax.axhline(0,color=RULE,lw=0.50,zorder=1)
ax.set_xticks(xs); ax.set_xticklabels(labs,fontsize=5,color=INK2)
ax.set_ylabel("oscine phyloP",fontsize=5,color=INK)
ax.set_title("A   AR full sites at the MAFB locus",loc="left",
             fontsize=7,color=INK,fontweight="600",pad=30)
ax.set_ylim(-1.85,2.05); ax.set_xlim(-.62,3.62)

ax=axes[1]
GR=[("genome-wide","ar",None),("MAFB locus","ar","MAFB")]
w=.32
for i,(lab,m,loc) in enumerate(GR):
    for j,(s,c) in enumerate((("ra",RA),("g4",G4))):
        v=sel(f"{F}/subtree222_{m}.tab",s,loc); lo,hi=ci(v); x=i+(j-.5)*w
        ax.plot([x,x],[lo,hi],color=c,lw=1.30,zorder=2,solid_capstyle="round")
        ax.scatter([x],[v.mean()],s=13.8,color=c,zorder=3,linewidths=0)
        ax.text(x+.055,v.mean(),f"{v.mean():+.2f}",ha="left",va="center",fontsize=5,
                color=c,fontweight="600")
        ax.text(x,-.20,f"n={len(v)}",ha="center",va="top",fontsize=5,color=INK3,
                transform=ax.get_xaxis_transform())
ax.axhline(0,color=INK,lw=0.50,zorder=1)
ax.set_xticks(range(len(GR))); ax.set_xticklabels([g[0] for g in GR],fontsize=5,color=INK)
ax.set_ylabel("oscine phyloP",fontsize=6,color=INK)
ax.set_title("B   the elevation is local to MAFB",loc="left",fontsize=7,
             color=INK,fontweight="600",pad=12)
ax.legend(handles=[Patch(facecolor=RA,label="Glut-RA DARs"),
                   Patch(facecolor=G4,label="GABA-4 DARs")],
          frameon=False,fontsize=5,loc="lower left")
ax.set_ylim(-.22,.80); ax.set_xlim(-.40,1.40)
for ax in axes:
    ax.margins(x=.03,y=.05)
    ax.locator_params(nbins=5)
    for sp in ("top","right"): ax.spines[sp].set_visible(False)
    for sp in ("left","bottom"): ax.spines[sp].set_color(RULE)
    ax.tick_params(colors=INK3,labelsize=5)
fig.tight_layout(rect=[0,.145,1,1])
fig.savefig(f"{F}/results/compare_gw_mafb1p5mb/fig_ar_mafb.png",dpi=200,facecolor="white",bbox_inches="tight",pad_inches=.01)
fig.savefig(f"{F}/results/compare_gw_mafb1p5mb/fig_ar_mafb.pdf",dpi=200,facecolor="white",bbox_inches="tight",pad_inches=.01)
print("  wrote fig_ar_mafb.png")
for m,nm in PW:
    a=sel(f"{F}/subtree222_{m}.tab","ra","MAFB"); b=sel(f"{F}/subtree222_{m}.tab","g4","MAFB")
    print(f"    {nm:<18} RA {a.mean():+.3f} (n={len(a)})  G4 {b.mean():+.3f} (n={len(b)})  p={perm(a,b):.3f}")
