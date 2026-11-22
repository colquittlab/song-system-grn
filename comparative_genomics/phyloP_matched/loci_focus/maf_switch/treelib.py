#!/usr/bin/env python
"""Map one alignment column onto a species tree, leaves and internal nodes alike.

Internal-node states are read from the alignment, not inferred here: cactus
stores a reconstructed sequence for every ancestral node, so an Anc/birdAnc row
in the MAF is that node's base. Clades that are monomorphic and contain no
oscine are collapsed to a wedge, so the substitutions are what the eye lands on.
Node order is rotated so the oscine clade sits at the bottom.
"""
import re, numpy as np, matplotlib
matplotlib.use("Agg")
matplotlib.rcParams.update({"font.family":"sans-serif",
    "font.sans-serif":["Arial","Liberation Sans","DejaVu Sans"],
    "pdf.fonttype":42,"ps.fonttype":42,
    "font.size":5,"axes.labelsize":6,"axes.titlesize":7,"axes.linewidth":.4})
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
INK="#171B22"; INK3="#6E7887"; RULE="#DDE1E7"
BASE={"A":"#2E8B57","C":"#e31a1c","G":"#B8860B","T":"#2C6E7F"}
GAP="#C9CED6"

def parse_newick(s):
    s=s.strip().rstrip(";"); pos=[0]
    def node():
        ch=[]
        if s[pos[0]]=="(":
            pos[0]+=1
            while True:
                ch.append(node())
                if pos[0]<len(s) and s[pos[0]]==",": pos[0]+=1
                else: break
            pos[0]+=1
        m=re.match(r"[^,()\[\]:;]*",s[pos[0]:]); nm=m.group(0); pos[0]+=len(nm)
        bl=0.0
        if pos[0]<len(s) and s[pos[0]]==":":
            m=re.match(r":([-\d.eE+]*)",s[pos[0]:]); bl=float(m.group(1) or 0)
            pos[0]+=len(m.group(0))
        return dict(name=nm,bl=bl,ch=ch)
    return node()

def leaves(n): return [n] if not n["ch"] else [x for c in n["ch"] for x in leaves(c)]

def draw_tree(col,tree,stem,osc,title="",width=4.3,row=0.105):
    """col {genome: base}; tree newick str; osc set of oscine tip names."""
    root=parse_newick(tree) if isinstance(tree,str) else tree
    has_osc=lambda n: any(x["name"] in osc for x in leaves(n))
    def states(n): return {col.get(x["name"]) for x in leaves(n)}|{col.get(n["name"])}
    def mark(n):
        st=states(n)-{None}
        n["collapse"]=(not has_osc(n)) and len(st)==1 and len(leaves(n))>1
        if not n["collapse"]:
            for c in n["ch"]: mark(c)
    mark(root)
    def rotate(n):                      # oscine-bearing child last -> drawn at the bottom
        if not n["ch"]: return
        for c in n["ch"]: rotate(c)
        bear=[c for c in n["ch"] if has_osc(c)]
        if len(bear)==1:
            n["ch"].remove(bear[0]); n["ch"].append(bear[0])
    rotate(root)
    def setx(n,x):
        n["x"]=x+n["bl"]
        for c in n["ch"]: setx(c,n["x"])
    setx(root,0.0)
    order=[]
    def layout(n):
        if n["collapse"] or not n["ch"]:
            n["y"]=len(order); order.append(n); return
        for c in n["ch"]: layout(c)
        n["y"]=float(np.mean([c["y"] for c in n["ch"]]))
    layout(root)
    N=len(order)
    fig,ax=plt.subplots(figsize=(width,max(1.9,N*row+0.60)))
    fig.patch.set_facecolor("white")
    nsub=[0]
    def draw(n):
        for c in (n["ch"] if not n["collapse"] else []):
            sp,sc=col.get(n["name"]),col.get(c["name"])
            chg=bool(sp and sc and sp!=sc and sp in BASE and sc in BASE)
            ax.plot([n["x"]]*2,[n["y"],c["y"]],color=RULE,lw=.55,zorder=1,solid_capstyle="round")
            ax.plot([n["x"],c["x"]],[c["y"]]*2,color="#C7521B" if chg else RULE,
                    lw=1.5 if chg else .55,zorder=2,solid_capstyle="round")
            if chg:
                nsub[0]+=1
                ax.annotate(f"{sp}→{sc}",xy=((n["x"]+c["x"])/2,c["y"]),
                            xytext=(0,3.1),textcoords="offset points",ha="center",
                            va="bottom",fontsize=4.4,color="#8A3A12",fontweight="700")
            draw(c)
        b=col.get(n["name"])
        if n["collapse"]:
            L=leaves(n); xm=max(x["x"] for x in L)
            ax.fill([n["x"],xm,xm],[n["y"],n["y"]-.36,n["y"]+.36],
                    color=BASE.get(b,GAP),alpha=.30,lw=0,zorder=1)
            ax.text(xm+.012,n["y"],f"  {len(L)} taxa",va="center",ha="left",
                    fontsize=4.2,color=INK3,style="italic")
        if b:
            ax.scatter([n["x"]],[n["y"]],s=11 if not n["ch"] else 7.5,
                       color=BASE.get(b,GAP),zorder=4,linewidths=0,
                       marker="o" if n["ch"] else "s")
        if not n["ch"]:
            ax.text(n["x"]+.012,n["y"],n["name"].replace("_"," "),va="center",ha="left",
                    fontsize=4.2,style="italic",
                    color="#8A1416" if n["name"] in osc else INK,
                    fontweight="700" if n["name"] in osc else "400")
    draw(root)
    ax.set_ylim(-.9,N-.1); ax.invert_yaxis()
    ax.set_xlim(-.02,max(n["x"] for n in order)*1.37); ax.axis("off")
    if title: ax.set_title(title,loc="left",fontsize=6.5,color=INK,fontweight="600",pad=3)
    seen=[b for b in "ACGT" if b in set(col.values())]
    ax.legend(handles=[Line2D([],[],marker="s",color=BASE[b],ls="",ms=3,label=b) for b in seen]
              +[Line2D([],[],color="#C7521B",lw=1.4,label="substitution"),
                Line2D([],[],color="none",label="bold = oscine (Passeri)")],
              frameon=False,fontsize=4.6,loc="lower left",ncol=len(seen)+2,
              bbox_to_anchor=(0,-.045),handletextpad=.3,columnspacing=.8)
    fig.subplots_adjust(left=.01,right=.99,top=.955,bottom=.045)
    for e in ("png","pdf"):
        fig.savefig(f"{stem}.{e}",dpi=240,facecolor="white",bbox_inches="tight",pad_inches=.015)
    plt.close(fig)
    return N,sum(1 for n in order if n["collapse"]),nsub[0]
