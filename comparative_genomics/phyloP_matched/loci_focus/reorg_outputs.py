#!/usr/bin/env python
"""File every result in loci_focus under the span it was computed on.

Moves only; nothing is deleted. Each move is logged to results/MOVES.tsv
(old path, new path). Script references to moved files are rewritten so every
writer and reader keeps working:

  * a file family that moves as a whole (fig_ar_fold*, tf_gain_*, ...) gets a
    path prefix inserted wherever a script names it;
  * outputs whose span is chosen at run time (tf_outlier_{gw,MAFB}.json and
    friends) are routed through span.result() by hand -- see MANUAL below;
  * per-site figures go to results/per_site/<collection>/<site span>/, and
    their writers route through span.per_site().

Pipeline inputs and intermediates (DAR sets, FIMO hits, per-element .tab
files, FASTA, MAFs, maf_switch/data*) stay where they are: they are read back
by the scripts and are not span-specific results.

usage: reorg_outputs.py [--apply]      (default is a dry run)
"""
import glob, os, re, shutil, sys

F = "/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus"
M = f"{F}/maf_switch"
APPLY = "--apply" in sys.argv
sys.path.insert(0, F)
from span import site_span
from reorg_manual import MANUAL

# ---------------------------------------------------------------- families --
# (family root, destination relative to loci_focus). Longest roots first, so
# fig_tf13_mafb is matched before fig_tf13.
R = "results"
FAMILIES = [
    # one span only
    ("fig_motifscan", f"{R}/genome_wide"), ("fig_tf_main", f"{R}/genome_wide"),
    ("allmotif_subtree222", f"{R}/genome_wide"),
    ("t_dod.html", f"{R}/genome_wide"), ("t_ra.html", f"{R}/genome_wide"),
    ("fig_maf_switches", f"{R}/mafb_1p5mb"), ("fig_maf_tree", f"{R}/mafb_1p5mb"),
    ("maf_ar_mafb", f"{R}/mafb_1p5mb"), ("maf_target_col", f"{R}/mafb_1p5mb"),
    ("fig_tf_mafb_main", f"{R}/mafb_1p5mb"), ("t_mafb.html", f"{R}/mafb_1p5mb"),
    # side by side
    ("fig_tf13_mafb", f"{R}/compare_gw_mafb1p5mb"),
    ("fig_tf13", f"{R}/genome_wide"),
    ("fig_allmotif_outlier", f"{R}/compare_gw_mafb1p5mb"),
    ("fig_ar3_", f"{R}/compare_gw_mafb1p5mb"), ("fig_ar_", f"{R}/compare_gw_mafb1p5mb"),
    ("fig_jaspar_positions", f"{R}/compare_gw_mafb1p5mb"),
    ("fig_tf_dod_scatter", f"{R}/compare_gw_mafb1p5mb"),
    ("ar3_stats", f"{R}/compare_gw_mafb1p5mb"), ("pairs_all", f"{R}/compare_gw_mafb1p5mb"),
    ("tf13_results", f"{R}/compare_gw_mafb1p5mb"), ("tf14_results", f"{R}/compare_gw_mafb1p5mb"),
    ("tf_final", f"{R}/compare_gw_mafb1p5mb"),
    ("fig_tf_gain", f"{R}/compare_domain_mafb1p5mb"),
    ("tf_gain_compare", f"{R}/compare_domain_mafb1p5mb"),
    ("tf_gain_rank", f"{R}/compare_domain_mafb1p5mb"),
    ("ar_report_part", f"{R}/compare_mafb_spans"), ("slide_", f"{R}/compare_mafb_spans"),
    ("locus_contrast", f"{R}/compare_gw_mafb_emx2"), ("fig_maskcurve", f"{R}/compare_gw_mafb_emx2"),
    # tracks
    ("arefull_sites_MAFB_lsd2", "tracks/mafb_1p5mb"), ("jasparAR_sites_MAFB_lsd2", "tracks/mafb_1p5mb"),
    ("tf_sites_MAFB_lsd2", "tracks/mafb_1p5mb"), ("notable_sites_lsd2", "tracks/mafb_1p5mb"),
    ("arefull_sites_lsd2", "tracks/genome_wide"), ("jasparAR_sites_lsd2", "tracks/genome_wide"),
    ("ar_gain_1p5mb_", "tracks/mafb_1p5mb"), ("tf_gain_1p5mb_", "tracks/mafb_1p5mb"),
    ("ar_gain_all_lsd2", "tracks/mafb_domain"), ("ar_gain_candidates_lsd2", "tracks/mafb_domain"),
    ("tf_gain_all_lsd2", "tracks/mafb_domain"), ("tf_gain_candidates_lsd2", "tracks/mafb_domain"),
    ("trackDb_", "tracks"), ("arefull.as", "tracks"), ("jasparAR.as", "tracks"),
    ("ar_gain.as", "tracks"), ("tf_gain.as", "tracks"),
]
# split by span at run time: routed by hand through span.result()
SPLIT = [(r"^allmotif_outlier_(gw|MAFB)\.json$",), (r"^tf_outlier_(gw|MAFB)\.json$",),
         (r"^tf_outlier_perm_(gw|MAFB)\.json$",), (r"^tf_tab_\w+?_(gw|MAFB)\.json$",),
         (r"^locimotif_(MAFB|EMX2)\.json$",)]
SPLIT_DIR = {"gw": "genome_wide", "MAFB": "mafb_1p5mb", "EMX2": "emx2_1p5mb"}
# maf_switch figures (named relative to maf_switch/)
MS_FAMILIES = [("fig_switches363_domain", f"{R}/mafb_domain"), ("fig_switches74_domain", f"{R}/mafb_domain"),
               ("fig_score_dists", f"{R}/mafb_domain"), ("fig_switches363", f"{R}/mafb_1p5mb")]
# per-site collections: (source dir under maf_switch, collection name)
PER_SITE = [("site_aln", "site_alignments"), ("trees", "switch_trees_363"),
            ("trees74", "switch_trees_74"), ("trees_gain", "gain_trees_ar"),
            ("trees_tf", "gain_trees_tf")]


def family_of(name, fams):
    for root, dst in fams:
        if name.startswith(root) or name == root:
            return root, dst
    return None


def plan():
    moves = []
    for p in sorted(glob.glob(f"{F}/*")):
        n = os.path.basename(p)
        if not os.path.isfile(p) or n.endswith((".py", ".sh")): continue
        sp = None
        for (rx,) in SPLIT:
            m = re.match(rx, n)
            if m: sp = SPLIT_DIR[m.group(1)]; break
        if sp:
            moves.append((p, f"{F}/{R}/{sp}/{n}")); continue
        f = family_of(n, FAMILIES)
        if f: moves.append((p, f"{F}/{f[1]}/{n}"))
    for p in sorted(glob.glob(f"{M}/fig_*")):
        n = os.path.basename(p); f = family_of(n, MS_FAMILIES)
        if f: moves.append((p, f"{F}/{f[1]}/{n}"))
    for src, coll in PER_SITE:
        for p in sorted(glob.glob(f"{M}/{src}/*.png") + glob.glob(f"{M}/{src}/*.pdf")):
            n = os.path.basename(p)
            moves.append((p, f"{F}/{R}/per_site/{coll}/{site_span(n)}/{n}"))
    return moves


# ------------------------------------------------------------- rewriting --
PRE = r"(\{F\}/|\{P\}/loci_focus/|\{F\}\"\+\"/|(?<=[\"'(]))"


def rewrite_static(src, in_maf_switch):
    """Insert destination prefixes for whole-family moves."""
    out = src
    for root, dst in FAMILIES:
        # {F}/root  or  "root (bare, loci_focus-relative)
        pat = re.compile(r"(\{F\}/|\{P\}/loci_focus/)(?!results/|tracks/)" + re.escape(root))
        out = pat.sub(lambda m: m.group(1) + dst + "/" + root, out)
        if not in_maf_switch:
            # a quoted string that starts with the family name (roots ending in "_"
            # are prefixes, so no terminator is required after them)
            tail = r"" if root.endswith("_") else r"(?=[._{\"'])"
            pat = re.compile(r"(?<=[\"'])(?!results/|tracks/)" + re.escape(root) + tail)
            out = pat.sub(dst + "/" + root, out)
    for root, dst in MS_FAMILIES:
        pat = re.compile(r"(\{D\}/|\{M\}/)" + re.escape(root))
        out = pat.sub(lambda m: m.group(1) + "../" + dst + "/" + root, out)
        if in_maf_switch:
            pat = re.compile(r"(?<=[\"'])" + re.escape(root) + r"(?=[._{\"'])")
            out = pat.sub("../" + dst + "/" + root, out)
    return out


def main():
    moves = plan()
    dests = {}
    for a, b in moves:
        dests.setdefault(os.path.dirname(b).replace(F + "/", ""), []).append(os.path.basename(a))
    print(f"  {len(moves)} files to move")
    for d in sorted(dests):
        print(f"    {d:58s} {len(dests[d]):4d}   e.g. {', '.join(sorted(dests[d])[:3])}")
    scripts = sorted(glob.glob(f"{F}/*.py") + glob.glob(f"{F}/*.sh") + glob.glob(f"{M}/*.py"))
    changed = {}
    for s in scripts:
        if os.path.basename(s) in ("reorg_outputs.py", "reorg_manual.py"): continue
        t = open(s).read()
        u = rewrite_static(t, s.startswith(M + "/"))
        if u != t: changed[s] = (t, u)
    print(f"\n  {len(changed)} scripts get static path rewrites")
    bad = 0
    for rel, old, new in MANUAL:
        s = f"{F}/{rel}"
        t, u = changed.get(s, (open(s).read(), open(s).read()))
        n = u.count(old)
        if n != 1:
            print(f"  MANUAL edit does not match exactly once ({n}x): {rel}: {old[:70]!r}"); bad += 1
            continue
        changed[s] = (t, u.replace(old, new))
    print(f"  {len(MANUAL)} hand edits, {bad} failed to match; {len(changed)} scripts changed in total")
    if bad: raise SystemExit("fix the hand edits before applying")
    if not APPLY:
        for s, (t, u) in list(changed.items()):
            for a, b in zip(t.splitlines(), u.splitlines()):
                if a != b:
                    print(f"    {os.path.relpath(s, F)}:  {a.strip()[:90]}\n{'':6}-> {b.strip()[:110]}")
        return
    log = open(f"{F}/{R}/MOVES.tsv", "a") if os.path.isdir(f"{F}/{R}") else None
    os.makedirs(f"{F}/{R}", exist_ok=True)
    log = log or open(f"{F}/{R}/MOVES.tsv", "a")
    for a, b in moves:
        os.makedirs(os.path.dirname(b), exist_ok=True)
        if os.path.exists(b): raise SystemExit(f"refusing to overwrite {b}")
        shutil.move(a, b)
        log.write(f"{os.path.relpath(a, F)}\t{os.path.relpath(b, F)}\n")
    log.close()
    for s, (t, u) in changed.items():
        open(s, "w").write(u)
    print(f"  moved {len(moves)} files, rewrote {len(changed)} scripts; log in {R}/MOVES.tsv")


if __name__ == "__main__":
    main()
