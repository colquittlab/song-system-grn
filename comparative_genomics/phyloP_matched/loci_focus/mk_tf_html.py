import json
import sys as _s; _s.path.insert(0, "/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus"); from span import result  # results/<span>/
def f(x,d=3): 
    s=("%+."+str(d)+"f")%x
    return s.replace("-","&minus;").replace("+","+")
def pf(p):
    if p<1e-4: return '<td class="n sig">1e-4</td>'
    c="n sig" if p<.05 else ("n" if p<.1 else "n ns")
    return '<td class="%s">%s</td>'%(c, ("%.4f"%p).rstrip("0") if p<.01 else "%.2f"%p)
def qf(q):
    c="n sig" if q<.05 else ("n" if q<.1 else "n ns")
    return '<td class="%s">%s</td>'%(c, "1e-4" if q<1e-4 else (("%.4f"%q).rstrip("0") if q<.01 else "%.2f"%q))
O=json.load(open(result("tf_tab_subtree222_gw.json", "gw"))); W={r["tf"]:r for r in json.load(open(result("tf_tab_wholetree_gw.json", "gw")))}
M=json.load(open(result("tf_tab_subtree222_MAFB.json", "MAFB")))

# TABLE 1: DoD, oscine vs whole tree
rows=sorted(O,key=lambda r:-r["dod"]); out=[]
for r in rows:
    w=W[r["tf"]]; hl=' class="hl"' if r["tf"]=="AR" else ""
    out.append('      <tr%s><td>%s</td><td>%s</td><td class="n">%s</td>%s\n'
               '        <td class="n">%s</td>%s<td class="n">%s / %s</td></tr>'%(
        hl,r["tf"],r["reg"],f(r["dod"]),qf(r["q_p"]),f(w["dod"]),qf(w["q_p"]),
        format(r["na"],","),format(r["nb"],",")))
open("results/genome_wide/t_dod.html","w").write("\n".join(out))

# TABLE 2: motif minus own enhancer, Glut-RA, oscine vs whole tree
rows=sorted(O,key=lambda r:-r["ra"]); out=[]
for r in rows:
    w=W[r["tf"]]; hl=' class="hl"' if r["tf"]=="AR" else ""
    out.append('      <tr%s><td>%s</td><td>%s</td><td class="n">%s</td>%s\n'
               '        <td class="n">%s</td>%s</tr>'%(
        hl,r["tf"],r["reg"],f(r["ra"]),qf(r["q_pra"]),f(w["ra"]),qf(w["q_pra"])))
open("results/genome_wide/t_ra.html","w").write("\n".join(out))

# TABLE 3: MAFB oscine
rows=sorted(M,key=lambda r:-r["ra"]); out=[]
for r in rows:
    hl=' class="hl"' if r["tf"]=="AR" else ""
    out.append('      <tr%s><td>%s</td><td>%s</td><td class="n">%s</td>%s\n'
               '        <td class="n">%s</td>%s<td class="n">%s</td>\n'
               '        %s<td class="n">%d / %d</td></tr>'%(
        hl,r["tf"],r["reg"],f(r["ra"]),pf(r["pra"]),f(r["g4"]),pf(r["pg4"]),f(r["dod"]),
        pf(r["p"]),r["na"],r["nb"]))
open("results/mafb_1p5mb/t_mafb.html","w").write("\n".join(out))
print("ok")
