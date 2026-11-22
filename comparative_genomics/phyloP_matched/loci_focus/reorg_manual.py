"""Hand edits for reorg_outputs.py: references the static prefix rewrite can't
handle. (script relative to loci_focus, old, new); applied after the static
rewrites, and each old string must occur exactly once."""
LF = "/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus"
IMP = f'import sys as _s; _s.path.insert(0, "{LF}"); from span import result  # results/<span>/\n'
IMPS = f'import sys as _s; _s.path.insert(0, "{LF}"); from span import per_site  # results/per_site/\n'

MANUAL = [
    # ---- outputs whose span is chosen at run time -> span.result()
    ("analyse_loci_motif.py", 'open(f"{F}/locimotif_{LOC}.json","w")',
     'open(result(f"locimotif_{LOC}.json", LOC),"w")'),
    ("analyse_loci_motif.py", "import numpy as np, collections, json, os\n",
     "import numpy as np, collections, json, os\n" + IMP),
    ("consolidate_ar2.py", 'json.load(open(f"tf_tab_{k}.json"))',
     'json.load(open(result(f"tf_tab_{k}.json", k.rsplit("_", 1)[1])))'),
    ("consolidate_ar2.py", 'json.load(open(f"allmotif_outlier_{k}.json"))',
     'json.load(open(result(f"allmotif_outlier_{k}.json", k)))'),
    ("consolidate_ar2.py", "import json, collections, sys, numpy as np\n",
     "import json, collections, sys, numpy as np\n" + IMP),
    ("tf_outlier.py", "open(f\"tf_outlier_{loc or 'gw'}.json\",\"w\")",
     "open(result(f\"tf_outlier_{loc or 'gw'}.json\", loc),\"w\")"),
    ("tf_outlier.py", "import numpy as np, sys, json\n", "import numpy as np, sys, json\n" + IMP),
    ("tf_outlier_perm.py", "open(f\"tf_outlier_perm_{loc or 'gw'}.json\",\"w\")",
     "open(result(f\"tf_outlier_perm_{loc or 'gw'}.json\", loc),\"w\")"),
    ("tf_outlier_perm.py", "import numpy as np, sys, json\n", "import numpy as np, sys, json\n" + IMP),
    ("tf_tables.py", "open(f\"tf_tab_{tree}_{loc or 'gw'}.json\",\"w\")",
     "open(result(f\"tf_tab_{tree}_{loc or 'gw'}.json\", loc),\"w\")"),
    ("tf_tables.py", "import numpy as np, sys, json\n", "import numpy as np, sys, json\n" + IMP),
    ("plot_tf_main.py", "json.load(open(f\"{F}/tf_tab_{t}_{l or 'gw'}.json\"))",
     "json.load(open(result(f\"tf_tab_{t}_{l or 'gw'}.json\", l)))"),
    ("plot_tf_main.py", "import numpy as np, sys, json, matplotlib\n",
     "import numpy as np, sys, json, matplotlib\n" + IMP),
    ("mk_tf_html.py", 'json.load(open("tf_tab_subtree222_gw.json"))',
     'json.load(open(result("tf_tab_subtree222_gw.json", "gw")))'),
    ("mk_tf_html.py", 'json.load(open("tf_tab_wholetree_gw.json"))',
     'json.load(open(result("tf_tab_wholetree_gw.json", "gw")))'),
    ("mk_tf_html.py", 'json.load(open("tf_tab_subtree222_MAFB.json"))',
     'json.load(open(result("tf_tab_subtree222_MAFB.json", "MAFB")))'),
    ("mk_tf_html.py", "import json\n", "import json\n" + IMP),
    ("plot_allmotif_outlier.py", 'json.load(open(f"{F}/allmotif_outlier_{tagl}.json"))',
     'json.load(open(result(f"allmotif_outlier_{tagl}.json", tagl)))'),
    ("plot_allmotif_outlier.py", "import numpy as np, json, sys, matplotlib\n",
     "import numpy as np, json, sys, matplotlib\n" + IMP),
    ("allmotif_outlier.py", 'open(f"{F}/allmotif_outlier_{tagl}.json","w")',
     'open(result(f"allmotif_outlier_{tagl}.json", tagl),"w")'),
    ("allmotif_outlier.py", "import numpy as np, json, sys\n", "import numpy as np, json, sys\n" + IMP),
    ("plot_tf_dod_scatter.py", 'json.load(open(f"{F}/tf_tab_subtree222_{tag}.json"))',
     'json.load(open(result(f"tf_tab_subtree222_{tag}.json", tag)))'),
    ("plot_tf_dod_scatter.py", 'json.load(open(f"{F}/tf_tab_wholetree_{tag}.json"))',
     'json.load(open(result(f"tf_tab_wholetree_{tag}.json", tag)))'),
    ("plot_tf_dod_scatter.py", "import numpy as np, json, matplotlib\n",
     "import numpy as np, json, matplotlib\n" + IMP),
    # ---- dynamic track names
    ("make_ar_track.py",
     'for sub,fn in ((None,f"{stem}_sites_lsd2.bed"),("MAFB",f"{stem}_sites_MAFB_lsd2.bed")):',
     'for sub,fn in ((None,f"tracks/genome_wide/{stem}_sites_lsd2.bed"),'
     '("MAFB",f"tracks/mafb_1p5mb/{stem}_sites_MAFB_lsd2.bed")):'),
    ("make_ar_track.py", 'with open(f"{F}/{stem}.as","w") as o:', 'with open(f"{F}/tracks/{stem}.as","w") as o:'),
    # ---- per-site figures -> results/per_site/<collection>/<site span>/
    ("maf_switch/gain_trees.py", 'draw(pres,tree,f"trees_gain/{nm}",OSC,',
     "draw(pres,tree,f\"{per_site('gain_trees_ar', nm)}/{nm}\",OSC,"),
    ("maf_switch/gain_trees.py", "from treelib import parse_newick, leaves\n",
     "from treelib import parse_newick, leaves\n" + IMPS),
    ("maf_switch/tf_gain_trees.py", 'draw(pres, tree, f"trees_tf/{nm}", OSC,',
     "draw(pres, tree, f\"{per_site('gain_trees_tf', nm)}/{nm}\", OSC,"),
    ("maf_switch/tf_gain_trees.py", "from gain_trees import count_events, draw\n",
     "from gain_trees import count_events, draw\n" + IMPS),
    ("maf_switch/batch_trees.py", 'stem=f"{D}/trees/{nm}"', "stem=f\"{per_site('switch_trees_363', nm)}/{nm}\""),
    ("maf_switch/batch_trees.py", "from treelib import draw_tree\n", "from treelib import draw_tree\n" + IMPS),
    ("maf_switch/batch_trees74.py", 'stem=f"{D}/trees74/{nm}"', "stem=f\"{per_site('switch_trees_74', nm)}/{nm}\""),
    ("maf_switch/batch_trees74.py", "from treelib import draw_tree\n", "from treelib import draw_tree\n" + IMPS),
    ("maf_switch/plot_site_alignment.py", 'fig.savefig(f"site_aln/{stem}.{ext}", dpi=300)',
     "fig.savefig(f\"{per_site('site_alignments', stem)}/{stem}.{ext}\", dpi=300)"),
    ("maf_switch/plot_site_alignment.py", 'print(f"  wrote site_aln/{stem}.pdf / .png',
     "print(f\"  wrote {per_site('site_alignments', stem)}/{stem}.pdf / .png"),
    ("maf_switch/plot_site_alignment.py", "import presence as PR\n", "import presence as PR\n" + IMPS),
    ("maf_switch/rebuild_gain.py",
     '        for f in glob.glob(f"{tdir}/*.pdf") + glob.glob(f"{tdir}/*.png"):\n'
     '            os.remove(f)                      # no stale trees from the old call\n',
     "        # trees are written to results/per_site/; the manifest says which are current\n"),
    ("maf_switch/rebuild_gain.py", 'N, g2, l2 = draw(p7, tree_str, f"{tdir}/{stem}", OSC74,',
     'N, g2, l2 = draw(p7, tree_str, f"{per_site(COLL[fam], stem)}/{stem}", OSC74,'),
    ("maf_switch/rebuild_gain.py", 'STEM = "Anc54"\n',
     'STEM = "Anc54"\nCOLL = {"AR": "gain_trees_ar", "TF": "gain_trees_tf"}\n'),
    ("maf_switch/rebuild_gain.py", "import presence as PR\n", "import presence as PR\n" + IMPS),
    # ---- the report's embedded per-site figures (both sites lie in the domain)
    ("build_ar_report.py", 'img(f"{M}/site_aln/aln_jaspar_RA_Glut-HVC-2.394967_1.png",',
     'img(f"{F}/results/per_site/site_alignments/mafb_domain/aln_jaspar_RA_Glut-HVC-2.394967_1.png",'),
    ("build_ar_report.py", 'img(f"{M}/trees_gain/jaspar_RA_Glut-HVC-2.394967_1.png",',
     'img(f"{F}/results/per_site/gain_trees_ar/mafb_domain/jaspar_RA_Glut-HVC-2.394967_1.png",'),
    ("build_ar_report.py", 'img(f"{M}/trees/jaspar_RA_Glut-RA.395131_1_p5.png",',
     'img(f"{F}/results/per_site/switch_trees_363/mafb_domain/jaspar_RA_Glut-RA.395131_1_p5.png",'),
]

# Applied by hand after the run: the static pattern required "." or "_" after a
# family root ("ar_report_part1" has a digit), and bare names inside maf_switch
# scripts were not rewritten (pwm_control.py joins them to {F} later).
POST_FIXES = [
    ("consolidate_ar.py", 'open("ar_report_part1.json","w")', 'open("results/compare_mafb_spans/ar_report_part1.json","w")'),
    ("consolidate_ar2.py", 'open("ar_report_part23.json","w")', 'open("results/compare_mafb_spans/ar_report_part23.json","w")'),
    ("maf_switch/pwm_control.py", '"arefull_sites_MAFB_lsd2.bed"', '"tracks/mafb_1p5mb/arefull_sites_MAFB_lsd2.bed"'),
    ("maf_switch/pwm_control.py", '"jasparAR_sites_MAFB_lsd2.bed"', '"tracks/mafb_1p5mb/jasparAR_sites_MAFB_lsd2.bed"'),
]
POST_FIXES.append(("maf_switch/plot_site_alignment.py",
    "f\"{F}/{'arefull' if key == 'full' else 'jasparAR'}_sites_MAFB_lsd2.bed\"",
    "f\"{F}/tracks/mafb_1p5mb/{'arefull' if key == 'full' else 'jasparAR'}_sites_MAFB_lsd2.bed\""))
