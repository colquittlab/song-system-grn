#!/usr/bin/env python
"""Per-position phyloP across the 16 bp JASPAR MA0007.3 AR motif."""
import numpy as np, collections, re, matplotlib
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
import pandas as pd, logomaker
F="/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus"
rng=np.random.default_rng(0)
INK="#171B22"; INK3="#6E7887"; RULE="#DDE1E7"; RA="#e31a1c"; G4="#2C6E7F"
W=16
txt=open(f"{F}/ar_jaspar.meme").read()
blk=txt.split("letter-probability")[1].split("\n")[1:]
pwm=np.array([[float(x) for x in l.split()] for l in blk if len(l.split())==4][:W])
ic=np.array([2+(p*np.log2(np.clip(p,1e-9,None))).sum() for p in pwm])
info=pd.DataFrame(pwm*ic[:,None],columns=list("ACGT"))
d=collections.defaultdict(list)
for l in open(f"{F}/subtree222_jaspos.tab"):
    p=l.split("\t")
    if int(p[2])==0: continue
    s,loc,k,_=p[0].split("|"); d[(s,loc,int(k))].append(float(p[5]))
def series(s,loc=None):
    m=[];lo=[];hi=[];nn=0
    for k in range(W):
        v=np.array([x for (ss,ll,kk),vals in d.items()
                    if ss==s and kk==k and (loc is None or ll==loc) for x in vals])
        nn=len(v); m.append(v.mean())
        b=[np.mean(rng.choice(v,len(v))) for _ in range(2000)]
        a,c=np.percentile(b,[2.5,97.5]); lo.append(a); hi.append(c)
    return np.array(m),np.array(lo),np.array(hi),nn
fig,axes=plt.subplots(3,1,figsize=(2.86,2.36),sharex=True,
                      gridspec_kw={"height_ratios":[1,1.5,1.5]})
fig.patch.set_facecolor("white")
ax=axes[0]
logomaker.Logo(info,ax=ax,color_scheme="classic",show_spines=False,font_name="Arial")
ax.set_ylabel("bits",fontsize=6,color=INK); ax.set_yticks([0,1,2])
ax.set_title("A   JASPAR MA0007.3 (AR), 16 bp",loc="left",fontsize=5,color=INK,fontweight="600",pad=8)
ax.tick_params(colors=INK3,labelsize=5)
for i,(ttl,loc) in enumerate((("B   per-position oscine phyloP, genome-wide",None),
                              ("C   per-position oscine phyloP at MAFB","MAFB"))):
    ax=axes[i+1]
    for s,col,lab in (("ra",RA,"Glut-RA DARs"),("g4",G4,"GABA-4 DARs")):
        m,lo,hi,n=series(s,loc); x=np.arange(W)
        ax.fill_between(x,lo,hi,color=col,alpha=.16,linewidth=0)
        ax.plot(x,m,color=col,lw=1.00,marker="o",ms=2.2,label=f"{lab}  (n={n})")
    ax.axhline(0,color=RULE,lw=0.50)
    ax.set_ylabel("oscine phyloP",fontsize=6,color=INK)
    ax.set_title(ttl,loc="left",fontsize=7,color=INK,fontweight="600",pad=8)
    ax.legend(frameon=False,fontsize=5,loc="upper center",ncol=2,
              bbox_to_anchor=(.5,1.30),columnspacing=1.0,handletextpad=.3)
    ax.tick_params(colors=INK3,labelsize=5)
    ax.margins(x=.03,y=.05)
    ax.locator_params(nbins=5)
    for sp in ("top","right"): ax.spines[sp].set_visible(False)
    for sp in ("left","bottom"): ax.spines[sp].set_color(RULE)
axes[2].set_xlabel("motif position",fontsize=6,color=INK)
axes[2].set_xticks(range(W)); axes[2].set_xticklabels(range(1,W+1))
HI=ic>=1.4
stat=[]
inst=collections.defaultdict(dict)
for l in open(f"{F}/subtree222_jaspos.tab"):
    p=l.split("\t")
    if int(p[2])==0: continue
    s,loc,k,idx=p[0].split("|"); inst[(s,loc,idx)][int(k)]=float(p[5])
for s,loc,lab in (("ra",None,"Glut-RA gw"),("g4",None,"GABA-4 gw"),
                  ("ra","MAFB","Glut-RA MAFB"),("g4","MAFB","GABA-4 MAFB")):
    V=np.array([[v[k] for k in range(W)] for (ss,ll,ix),v in inst.items()
                if ss==s and (loc is None or ll==loc) and len(v)==W])
    if len(V)<10: continue
    dd=V[:,HI].mean(1)-V[:,~HI].mean(1)
    b=np.array([np.mean(rng.choice(dd,len(dd))) for _ in range(20000)])
    p=max(2*min((b<=0).mean(),(b>=0).mean()),1/20000)
    stat.append(f"{lab} {dd.mean():+.3f} (p={p:.2f})")
    print(f"  {lab:<15} high-IC {V[:,HI].mean():+.3f}  low-IC {V[:,~HI].mean():+.3f}  "
          f"diff {dd.mean():+.3f}  p={p:.4f}  n={len(V)}")
fig.subplots_adjust(left=.105,right=.985,top=.945,bottom=.125,hspace=.62)
fig.savefig(f"{F}/results/compare_gw_mafb1p5mb/fig_jaspar_positions.png",dpi=200,facecolor="white",bbox_inches="tight",pad_inches=.01)
fig.savefig(f"{F}/results/compare_gw_mafb1p5mb/fig_jaspar_positions.pdf",dpi=200,facecolor="white",bbox_inches="tight",pad_inches=.01)
print("  wrote fig_jaspar_positions.png")
print(f"  IC: {np.round(ic,2)}")
