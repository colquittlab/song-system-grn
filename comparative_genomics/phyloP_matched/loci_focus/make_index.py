#!/usr/bin/env python
"""Write INDEX.md: every result in loci_focus, the span it was computed on,
what it shows, and the script that makes it. Rerun after adding outputs."""
import collections, glob, os, re

F = "/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus"

# family root -> (what it shows, producing script)
DESC = [
    ("fig_motifscan", "Where AR sits among 182 motif archetypes", "plot_motifscan.py"),
    ("fig_tf13_domain", "SCENIC+ regulators, motif minus rest of enhancer: genome-wide (A) vs DHX35-TOP1 domain (B), oscine", "plot_tf.py"),
    ("fig_tf_domain_main", "Regulators ranked, motif minus rest of enhancer, in the domain; oscine (A), whole tree (B)", "plot_tf_main.py"),
    ("fig_tf_dod_scatter_domain", "Whole-tree vs oscine DoD for the regulators: genome-wide (A) vs domain (B)", "plot_tf_dod_scatter.py"),
    ("scenic_motif_evidence", "Per regulator: every enriched motif SCENIC+ annotated to it in the nominating run (config4), with NES", "inline (see build_scenic_motifs.py)"),
    ("fig_motif_sensitivity", "Every statistic for every defensible matrix per regulator (SCENIC+ cluster members, HOMER known, HOMER de novo)", "sens_summary.py"),
    ("motif_panel", "The 330-matrix sensitivity panel: regulator, source, relation, consensus", "build_motif_panel.py"),
    ("sens_phylop", "Per matrix: Glut-RA/GABA-4 DoD on both trees in each span, outlier z genome-wide", "sens_phylop.py"),
    ("sens_gain_sites", "Per site and matrix: presence counts, gain, 74-way tier, at the locus", "sens_gain.py"),
    ("sens_gain", "Per matrix: site gain, candidates and oscine-stem gains in the domain and +/-1.5 Mb", "sens_gain.py"),
    ("sens_summary", "Per regulator and statistic: medoid, spread, sign agreement, robustness call; AR rank under every AR matrix", "sens_summary.py"),
    ("tf_scenic_choice", "The matrix chosen for each regulator: cluster, NES, medoid member, consensus", "build_scenic_motifs.py"),
    ("fig_tf13_mafb", "SCENIC+ regulators, motif minus rest of enhancer: genome-wide (A) vs at MAFB (B), oscine", "plot_tf13_mafb.py / plot_tf.py"),
    ("fig_tf13", "SCENIC+ regulators, Glut-RA vs GABA-4 asymmetry, oscine vs whole tree", "plot_tf.py / plot_tf13.py"),
    ("fig_tf_main", "Regulators ranked, motif minus rest of enhancer, oscine (A) and whole tree (B); no span on the figure", "plot_tf_main.py (earlier version)"),
    ("fig_tf_mafb_main", "Regulators ranked, motif minus rest of enhancer, at MAFB; oscine (A), whole tree (B)", "plot_tf_main.py"),
    ("fig_tf_dod_scatter", "Whole-tree vs oscine difference-of-differences for the regulators", "plot_tf_dod_scatter.py"),
    ("fig_tf_gain_rank", "Motif site gains per regulator vs pooled expectation; candidate and oscine-stem gains", "plot_tf_gain_rank.py"),
    ("fig_tf_gain", "Mean oscine site gain per regulator matrix, element-bootstrap CIs", "tf_gain_compare.py"),
    ("fig_allmotif_outlier_domain", "Regulators against the archetype DoD trend: genome-wide (A), domain at >=8 (B) and >=4 (C) elements per class", "plot_allmotif_outlier.py"),
    ("fig_allmotif_outlier", "Regulators against the 182-archetype whole-tree/oscine DoD trend; AR = MA0007.3 archetype, AR* = HOMER ARE-full", "plot_allmotif_outlier.py"),
    ("fig_ar3_dod", "Three AR motifs: difference-of-differences across units", "plot_ar3.py"),
    ("fig_ar3_fold", "Three AR motifs: motif vs enhancer levels, merged-site unit", "plot_ar3.py"),
    ("fig_ar3_positions", "Per-position phyloP for the three AR site types", "plot_ar3_positions.py"),
    ("fig_ar_fold", "AR (ARE-full) motif vs its enhancer, levels and gap", "plot_ar_fold.py"),
    ("fig_ar_heatmap", "Per-instance phyloP across the JASPAR AR motif, one row per site", "plot_ar_heatmap.py"),
    ("fig_ar_mafb", "AR full/half sites at MAFB, Glut-RA vs GABA-4; elevation vs genome-wide", "plot_ar.py"),
    ("fig_ar_oscine", "AR-site constraint ladder on the oscine lineage, two trees", "plot_ar_oscine.py"),
    ("fig_ar_paired_ladder", "Paired AR-vs-enhancer, ladder style", "plot_paired_ladder.py / plot_paired_ladder_jas.py"),
    ("fig_ar_paired", "AR motif minus remainder of its own enhancer, both trees", "plot_paired.py / plot_paired_jas.py"),
    ("fig_ar_positions", "Per-position phyloP across the ARE-full site with logo", "plot_ar_positions.py"),
    ("fig_jaspar_positions", "Per-position phyloP across the JASPAR MA0007.3 site", "plot_jaspar_positions.py"),
    ("fig_maf_switches", "Per-site oscine/outgroup base switches across AR motif sites", "plot_maf_switches.py"),
    ("fig_maf_tree", "One alignment column mapped onto the species tree; the site is not recorded in the file", "plot_maf_tree.py"),
    ("fig_maskcurve", "Can motif masking produce a neutral reference class inside enhancers? MAFB, EMX2, other", "plot_maskcurve.py"),
    ("fig_switches363_domain", "Oscine/outgroup switches inside the domain, avian363", "maf_switch/plot_switches_domain.py"),
    ("fig_switches74_domain", "Oscine/outgroup switches inside the domain, 74-way", "maf_switch/plot_switches_domain.py"),
    ("fig_switches363", "Oscine/outgroup switches at all AR sites, avian363", "maf_switch/plot_switches363.py"),
    ("fig_score_dists", "Motif-score distributions, oscine vs outgroup consensus, by DAR class", "maf_switch/plot_score_dists.py"),
    ("allmotif_subtree222", "182-archetype paired statistic, oscine subtree", "analyse_allmotif.py"),
    ("allmotif_outlier_domain_min4", "Outlier tests in the domain at >=4 elements per class (sensitivity; admits the AR archetype)", "allmotif_outlier.py"),
    ("allmotif_outlier_", "Outlier tests for the regulators on the archetype DoD trend", "allmotif_outlier.py"),
    ("tf_outlier_perm_", "Westfall-Young max-statistic outlier test on the DoD line", "tf_outlier_perm.py"),
    ("tf_outlier_", "Deleted-residual and meta-regression outlier tests on the DoD line", "tf_outlier.py"),
    ("tf_tab_", "Per-regulator table (motif minus rest, by class)", "tf_tables.py"),
    ("t_dod", "HTML table: DoD, oscine vs whole tree (patched in place by patch_tf.py)", "mk_tf_html.py"),
    ("t_ra", "HTML table: motif minus own enhancer, Glut-RA, oscine vs whole tree", "mk_tf_html.py"),
    ("t_mafb", "HTML table: at MAFB, oscine", "mk_tf_html.py"),
    ("ar3_stats_merged_prev", "Previous version of ar3_stats.json (kept for reference)", "ar3_stats.py (earlier run)"),
    ("ar3_stats", "Three AR motifs x units x scopes: Glut-RA/GABA-4 values and p", "ar3_stats.py"),
    ("pairs_all", "Per-motif DoD on both trees, 182 archetypes + 13 regulators; gw, MAFB +/-1.5 Mb and domain sets", "build_pairs_all.py"),
    ("tf13_results", "Paired statistic for the 13 regulators", "analyse_tf13.py"),
    ("tf14_results", "Paired statistic for 13 regulators + AR (version)", "earlier analysis"),
    ("tf_final", "Regulator results with q values (feeds fig_tf13)", "earlier analysis"),
    ("tf_gain_compare", "Per-matrix gain, tiers and AR-vs-others tests, both spans", "tf_gain_compare.py"),
    ("tf_gain_rank", "Observed/expected gain counts per matrix", "plot_tf_gain_rank.py"),
    ("ar_report_part1", "AR report data, Part 1: genome-wide, domain, gene +/-500 kb, +/-1 Mb, +/-1.5 Mb", "consolidate_ar.py"),
    ("ar_report_part23", "AR report data, Parts 2-3: +/-1.5 Mb and domain, plus outlier sets", "consolidate_ar2.py"),
    ("slide_summary", "Sliding-window AR effect along chr20", "inline analysis"),
    ("slide_winsize", "AR effect vs window size, 0.5-10 Mb", "inline analysis"),
    ("slide_", "Sliding-window arrays", "inline analysis"),
    ("locus_contrast", "phyloP matched-set contrast: genome-wide, MAFB, EMX2 (labelled '+/-1Mb' in the file = gene +/-1.5 Mb)", "locus_contrast.py"),
    ("locimotif_", "192-archetype paired statistic at the locus, with genome-wide reference", "analyse_loci_motif.py"),
    ("maf_ar_mafb", "Per-site, per-column base composition across species at AR sites", "maf_sites.py"),
    ("maf_target_col", "Bases at one alignment column (input to fig_maf_tree)", "plot_maf_tree.py"),
    ("arefull_sites_MAFB_lsd2", "ARE-full sites in Glut-RA/GABA-4 DARs at MAFB", "make_ar_track.py"),
    ("jasparAR_sites_MAFB_lsd2", "JASPAR AR sites in Glut-RA/GABA-4 DARs at MAFB", "make_ar_track.py"),
    ("tf_sites_MAFB_lsd2", "Sites for the 13 regulators' matrices at MAFB", "make_tf_sites.py"),
    ("notable_sites_lsd2", "Hand-picked DARs and ARE-full sites", "earlier analysis"),
    ("arefull_sites_lsd2", "ARE-full sites in all Glut-RA/GABA-4 DARs", "make_ar_track.py"),
    ("jasparAR_sites_lsd2", "JASPAR AR sites in all Glut-RA/GABA-4 DARs", "make_ar_track.py"),
    ("ar_gain_1p5mb_candidates", "AR site-gain candidates, gain >= 0.50, q < 0.05", "make_gain_track.py"),
    ("ar_gain_1p5mb_all", "All AR sites with coverage, scored by gain", "make_gain_track.py"),
    ("tf_gain_1p5mb_candidates", "Regulator site-gain candidates", "make_tf_gain_track.py"),
    ("tf_gain_1p5mb_all", "All regulator sites with coverage", "make_tf_gain_track.py"),
    ("ar_gain_candidates", "AR site-gain candidates, gain >= 0.50, within-domain q < 0.05", "make_gain_track.py"),
    ("ar_gain_all", "All AR sites with coverage, scored by gain", "make_gain_track.py"),
    ("tf_gain_candidates", "Regulator site-gain candidates", "make_tf_gain_track.py"),
    ("tf_gain_all", "All regulator sites with coverage", "make_tf_gain_track.py"),
    ("trackDb_", "UCSC trackDb stanza (bigDataUrl paths are relative to tracks/)", "track builders"),
    ("arefull.as", "autoSql", "make_ar_track.py"), ("jasparAR.as", "autoSql", "make_ar_track.py"),
    ("ar_gain.as", "autoSql", "make_gain_track.py"), ("tf_gain.as", "autoSql", "make_tf_gain_track.py"),
]
FOLDER = {
    "genome_wide": "All Glut-RA and GABA-4 DARs genome-wide.",
    "mafb_1p5mb": "MAFB gene +/-1.5 Mb (285 DARs, `vgp577/lifted/mafb_1mb.txt`). "
                  "**Every \"MAFB\" / \"at MAFB\" label in these figures and every `_MAFB` filename means this span, not the domain.**",
    "mafb_domain": "The default span: the gene-bounded DHX35-TOP1 domain, chr20:12,551,036-13,547,533 (118 DARs; a DAR counts if it overlaps).",
    "emx2_1p5mb": "EMX2 gene +/-1.5 Mb (`vgp577/lifted/emx2_1mb.txt`).",
    "compare_gw_mafb1p5mb": "Genome-wide and MAFB +/-1.5 Mb side by side in one figure or table (usually panel/row A genome-wide, B at MAFB).",
    "compare_domain_mafb1p5mb": "The domain and MAFB +/-1.5 Mb side by side (domain = left panel / first table).",
    "compare_mafb_spans": "Several MAFB spans: domain, gene +/-500 kb, +/-1 Mb, +/-1.5 Mb, genome-wide; and the sliding/size-sweep windows.",
    "compare_gw_mafb_emx2": "Genome-wide, MAFB +/-1.5 Mb and EMX2 +/-1.5 Mb.",
    "compare_gw_domain": "Genome-wide and the DHX35-TOP1 domain side by side.",
    "regulator_motifs": "How each MAFB regulator's matrix was chosen from SCENIC+'s own motif evidence (not a span).",
}
PER_SITE = {
    "site_alignments": ("Compact alignment + tree per candidate site, 74-way", "maf_switch/plot_site_alignment.py"),
    "site_alignments_363": ("Compact alignment + tree per candidate site, avian363 (basal oscines and both suboscine radiations)", "maf_switch/plot_site_alignment.py --aln 363"),
    "site_alignments_phylogram": ("As site_alignments, tree scaled by branch length (74-way: relative time)", "maf_switch/plot_site_alignment.py --phylogram"),
    "site_alignments_363_phylogram": ("As site_alignments_363, tree scaled by branch length (substitutions/site; dotted leaders to labels)", "maf_switch/plot_site_alignment.py --aln 363 --phylogram"),
    "gain_trees_ar": ("AR site gain/loss on the 74-way tree", "maf_switch/rebuild_gain.py"),
    "gain_trees_tf": ("Regulator site gain/loss on the 74-way tree", "maf_switch/rebuild_gain.py"),
    "gain_trees_tf_363": ("Regulator site gain/loss on the avian363 tree (oscine-enriched or stem/near-stem sites; uniform clades collapsed)", "maf_switch/gain_trees363.py"),
    "switch_trees_363": ("One oscine/outgroup switch column on the avian363 tree", "maf_switch/batch_trees.py"),
    "switch_trees_74": ("One switch column on the 74-way tree", "maf_switch/batch_trees74.py"),
}


def describe(n):
    for root, what, src in DESC:
        if n.startswith(root): return what, src
    return "", ""


def main():
    by = collections.defaultdict(set)
    for p in glob.glob(f"{F}/results/*/*") + glob.glob(f"{F}/tracks/*/*") + glob.glob(f"{F}/tracks/*"):
        if os.path.isfile(p) and not p.endswith("MOVES.tsv"):
            rel = os.path.relpath(p, F); d = os.path.dirname(rel)
            by[d].add(re.sub(r"\.(png|pdf)$", ".{png,pdf}", os.path.basename(p)))
    L = ["# loci_focus results index", "",
         "Every result, filed by the span it was computed on. Generated by `make_index.py`;",
         "moves are logged in `results/MOVES.tsv` (old path → new path).", "",
         "## Spans", ""]
    for k, v in FOLDER.items(): L.append(f"- **`{k}/`** — {v}")
    L += ["", "Outputs are routed by `span.py`: `result(name, span)`, `track(name, span)` and "
          "`per_site(collection, site)` -- new scripts should use them rather than writing to the top level."]
    L += ["", "Notes:", "",
          "- The ±1.5 Mb set is the gene ±500 kb locus expanded by 1 Mb. `locus_contrast.py` labels it "
          "\"±1Mb\" (±1 Mb of the locus); it is the same 285-DAR set.",
          "- Site-gain results (`tf_gain_*`, `ar_gain_*`, gain trees, `fig_tf_gain*`) use the corrected presence call "
          "(`maf_switch/presence.py`, 2026-11-08). Gain placement (the `tier` column, the ranking figure's second "
          "measure, the rate test) is on the avian363 tree (`maf_switch/tiers363.py`): stem, near-stem (the four "
          "backbone nodes still holding >=90% of oscines), within, inherited. The 74-way gain trees in "
          "`per_site/gain_trees_*` lack basal oscines, so their 'oscine stem' is the core-oscine stem.",
          "- Not moved: pipeline inputs and intermediates (DAR sets, FIMO hits, per-element `.tab` files, FASTA, MAFs, "
          "`maf_switch/data*`, tree manifests). Scripts read them back, and most span several regions.",
          "- `scripts_before_reorg/` holds every script as it was before the reorganization.",
          "- Regulator comparisons use `tf_scenic.meme`: one matrix per regulator, AR included, chosen from SCENIC+'s motif "
          "evidence for its cistrome (`build_scenic_motifs.py`). NR3C2 shares AR's matrix and is reported as AR. "
          "The superseded `tf13.meme`-based results are in `superseded/tf13_motifs_2026-11-09/`. "
          "`maf_switch/ascertain_tf*.py` still use the old set and feed nothing downstream.", ""]
    for d in sorted(k for k in by if k.startswith("results/")) + sorted(k for k in by if k.startswith("tracks")):
        L += [f"## `{d}/`", "", "| file | shows | made by |", "|---|---|---|"]
        for n in sorted(by[d]):
            what, src = describe(n)
            L.append(f"| `{n}` | {what} | `{src}` |" if src else f"| `{n}` | {what} | |")
        L.append("")
    L += ["## `results/per_site/`", "",
          "One figure per site, split by where the site lies: `mafb_domain/` (inside DHX35–TOP1) or "
          "`mafb_1p5mb_outside_domain/` (in the ±1.5 Mb set, beyond the domain).", "",
          "| collection | shows | made by | in domain | outside domain |", "|---|---|---|---|---|"]
    for c, (what, src) in PER_SITE.items():
        n = {s: len(glob.glob(f"{F}/results/per_site/{c}/{s}/*.png"))
             for s in ("mafb_domain", "mafb_1p5mb_outside_domain")}
        L.append(f"| `{c}/` | {what} | `{src}` | {n['mafb_domain']} | {n['mafb_1p5mb_outside_domain']} |")
    open(f"{F}/INDEX.md", "w").write("\n".join(L) + "\n")
    print(f"  wrote INDEX.md: {sum(len(v) for v in by.values())} entries in {len(by)} folders")


if __name__ == "__main__":
    main()
