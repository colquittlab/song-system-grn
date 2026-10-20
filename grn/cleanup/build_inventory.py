#!/usr/bin/env python3
"""Inventory of the GRN files tracked in a git ref (default origin/main): who added them, when, what they are, and a provisional
category: PUBLICATION (candidate to stay in main), DEVELOPMENT (candidate to move to grn-dev) or DECIDE (needs the repo owner's call).
The categories come from the rules in classify() below (path patterns plus evidence from make_configs.py, the commit history, README.md and
config/paths.R); they are provisional and meant to be edited. Writes inventory.csv next to this script and prints a dependency scan
(PUBLICATION files that mention DEVELOPMENT files).

    python build_inventory.py [ref=origin/main]
"""
import os, sys
REF = sys.argv[1] if len(sys.argv) > 1 else "origin/main"
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = HERE + "/"

# ---- 1. collect: size and first/last commit of every file
import subprocess, csv, os, collections
def git(*a): return subprocess.run(["git", *a], capture_output=True, text=True).stdout
files = []
for l in git("ls-tree", "-r", "-l", REF, "grn").splitlines():
    meta, path = l.split("\t"); files.append((path, int(meta.split()[3])))
rows = []
for path, size in files:
    last = git("log", "-1", "--format=%ad|%h|%s", "--date=short", REF, "--", path).strip().split("|", 2)
    first = git("log", "--diff-filter=A", "--format=%ad|%h|%s", "--date=short", REF, "--", path).strip().splitlines()
    f = first[-1].split("|", 2) if first else ["", "", ""]
    rows.append(dict(path=path, bytes=size, first_date=f[0], first_commit=f[1], first_subject=f[2], last_date=last[0], last_commit=last[1], last_subject=last[2]))
out = OUT + "files_raw.csv"
with open(out, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
print(len(rows), "files ->", out)
# history narrative: commits touching grn/ on main
print()
pass

# ---- 2. classify
import csv, os, re, subprocess, collections

T = OUT
rows = list(csv.DictReader(open(T + "files_raw.csv")))

CONFIG_PURPOSE = [  # from grn/ra-arco-hvc-nc/scenicplus/make_configs.py comments
    (range(1, 12), "sweep of one-at-a-time thresholds around config1 (configs 1-11)"),
    (range(12, 14), "normalization control: legacy double-normalised .raw input (configs 12-13)"),
    (range(14, 16), "topic-number control: 15 and 30 topics (configs 14-15)"),
    (range(16, 29), "Latin-hypercube sweep of the motif-stage thresholds (configs 16-28)"),
    (range(29, 34), "GABA/LGE region-set removal test and size-matched control (configs 29-33)"),
    (range(34, 36), "revised cisTopic object: 40 topics, no song-pair DARs (configs 34-35)"),
]

def classify(p):
    """returns (group, category, rationale)"""
    parts = p.split("/")
    d, base = parts[1], parts[-1]
    if d == "cistarget":
        return "cistarget build scripts", "PUBLICATION", "README-listed custom cisTarget database build scripts"
    if d in ("ra-arco-hvc-nc_glut", "ra-arco-hvc-nc_gaba", "ra-arco-hvc-nc_astro-oligo", "glut-ra_glut-arco-1"):
        why = ("README-listed per-subset pipeline on the earlier (pre-hybrid) labels; which label set is published? Their SCENIC+ outputs are read by "
               "overexpression/chick_mafb-emx2_run1.qmd via the SCENICPLUS_*_CONFIG_DIR constants in config/paths.R")
        if p.endswith(".ipynb"):
            why += " Notebook carries embedded outputs."
        return f"older-label subset run: {d}", "DECIDE", why
    if d == "ra-arco-hvc-nc_glut_hybrid":
        return "glut hybrid reference run", "DEVELOPMENT", "glut-only hybrid-label run, topic-grid exploration and perturbation set-up; not the config37 result the MAFB analyses use"
    if d == "ra-arco-hvc-nc":
        sub = parts[2]
        if sub == "anndata_rna":
            if base == "make_adata_doublenorm.py":
                return "hybrid all-cell: inputs", "DEVELOPMENT", "builds the deliberately wrong double-normalised input (control for the old .raw bug)"
            if base == "split_adata_by_region.ipynb":
                return "hybrid all-cell: inputs", "DEVELOPMENT", "legacy region split (identical copies in 5 directories)"
            return "hybrid all-cell: inputs", "PUBLICATION", "hybrid-label AnnData export / preparation for SCENIC+"
        if sub == "pycisTopic":
            if base == "pycisTopic.ipynb":
                return "hybrid all-cell: pycisTopic", "PUBLICATION", "topic modelling and region sets for the SCENIC+ input (6 MB: strip outputs)"
            if base in ("make_song_pair_dars_archr.R", "pycistopic_qc_commands.txt"):
                return "hybrid all-cell: pycisTopic", "PUBLICATION", "ArchR song-pair DARs used by config37 / pycisTopic QC commands"
            if base == "make_interneuron_dars_archr.R":
                return "hybrid all-cell: pycisTopic", "PUBLICATION", "ArchR interneuron DARs (PVALB-1/-2 vs LAMP5): called by scenicplus_mafb_up_both_dar_overlap.py, scenicplus_regulon_scatter_plots.R and the DAR heatmap"
            if base in ("song_pair_dars_archr_summary.csv", "interneuron_dars_archr_summary.csv"):
                return "hybrid all-cell: pycisTopic", "DECIDE", "record of the ArchR DAR counts written by make_song_pair_dars_archr.R / make_interneuron_dars_archr.R (small)"
            return "hybrid all-cell: pycisTopic", "DEVELOPMENT", "topic-number / region-set-variant / DAR-validation exploration"
        if sub == "cisTarget":
            return "hybrid all-cell: cisTarget", "PUBLICATION", "build of the hybrid cisTarget database that config37 uses"
        if sub == "scenicplus":
            m = re.match(r"grn/ra-arco-hvc-nc/scenicplus/config(\d+)/", p)
            if m:
                n = int(m.group(1))
                if n == 37:
                    return "hybrid all-cell: config37", "PUBLICATION", "config37 = the main comparison config (commit 'Make config37 the main comparison config'); hybrid cisTarget DB + ArchR song pairs"
                if n == 34:
                    return "hybrid all-cell: other configs", "DECIDE", "config34 builds the shared_k40 upstream (ACC_GEX, region-to-gene) that config36-39 reuse, so config37 cannot be rebuilt without it"
                if n == 35:
                    return "hybrid all-cell: other configs", "DECIDE", "config35 (40-topic object, loose thresholds) appears as a comparison row in the config37 MAFB shared-region tables"
                if n in (36, 38, 39):
                    return "hybrid all-cell: other configs", "DECIDE", "sibling of config37: 36 = no song-pair DARs, 38/39 = loose thresholds; controls for the 37 comparison"
                why = next(w for r, w in CONFIG_PURPOSE if n in r)
                return "hybrid all-cell: other configs", "DEVELOPMENT", why
            if base in ("genome_annotation.tsv", "chromsizes.tsv"):
                return "hybrid all-cell: scenicplus", "PUBLICATION", "reference inputs of the SCENIC+ run (identical copy in 6 directories)"
            if p.endswith("workflow/Snakefile"):
                return "hybrid all-cell: scenicplus", "PUBLICATION", "the Snakemake workflow"
            if base == "make_configs.py":
                return "hybrid all-cell: scenicplus", "DECIDE", "generates all 39 configs (and holds the sweep definitions); config37's definition is what a reader needs"
            if base == "config_parameters.tsv":
                return "hybrid all-cell: scenicplus", "DEVELOPMENT", "parameter table of the 39-config sweep"
            return "hybrid all-cell: scenicplus", "DECIDE", "cluster operations helper (prepare_shared_upstream.sh is needed to run config37)"
    if d == "analysis":
        stem = base
        if re.search(r"(^all\.qmd|^gaba\.qmd|^glut\.qmd|^glut-ra\.qmd|^astro-oligo\.qmd|_template\.qmd$)", stem):
            return "analysis: report templates", "DECIDE", "per-subset SCENIC+ result report (the older template serves the README-listed subsets; the hybrid ones are newer)"
        if stem.startswith("assemble_slim"):
            return "analysis: local helpers", "DECIDE", "builds the slim scplusmdata the report templates read"
        if re.search(r"(sweep_summary|normalization_control|gaba_removal|homeodomain|focal_tf)", stem):
            return "analysis: sweep / control diagnostics", "DEVELOPMENT", "cross-config diagnostic of the sweep, controls or GABA/LGE-removal test"
        if "coherence" in stem:
            return "analysis: coherence tiers", "DECIDE", "descriptive concordance tiers for every eRegulon in configs 1-13 (1.1 MB table)"
        if re.search(r"^(scenicplus_mafb_targets\.py|scenicplus_mafb_target_membership\.csv|scenicplus_mafb_shared_program(_vs_PVALB1|_vs_SST1)?\.csv)$", stem):
            return "analysis: cross-config MAFB analysis", "DECIDE", "compares MAFB across configs 1-35 (needs the sweep results and config_parameters.tsv); the config37-suffixed files are the single-config versions"
        if re.search(r"(mafb|regulon_|shared_program|shared_regions|dar_overlap|accessibility_concordance|region_accessibility|pv_fs_literature|fs_gene_aliases|loc_gene_labels)", stem):
            return "analysis: MAFB / regulon figure code", "PUBLICATION", "config37 MAFB shared-program analysis (manuscript figure code; the unmerged branch worktree-grn-k40-nosongpairs continues it: 'Panel A', 'Fig. 4D')"
        return "analysis: other", "DECIDE", "not matched by any rule"
    return "other", "DECIDE", "not matched by any rule"

out = []
for r in rows:
    g, c, why = classify(r["path"])
    out.append(dict(path=r["path"], bytes=r["bytes"], category=c, group=g, rationale=why, first_added=r["first_date"], last_modified=r["last_date"], last_commit_subject=r["last_subject"][:90]))
with open(T + "inventory.csv", "w", newline="") as fh:   # copied into grn/cleanup/ by the next step
    w = csv.DictWriter(fh, fieldnames=list(out[0])); w.writeheader(); w.writerows(out)
tot = collections.defaultdict(lambda: [0, 0])
grp = collections.defaultdict(lambda: collections.defaultdict(lambda: [0, 0]))
for o in out:
    tot[o["category"]][0] += 1; tot[o["category"]][1] += int(o["bytes"])
    grp[o["group"]][o["category"]][0] += 1; grp[o["group"]][o["category"]][1] += int(o["bytes"])
print("TOTAL", len(out), "files,", round(sum(int(o['bytes']) for o in out) / 1e6, 1), "MB")
for c, (n, b) in sorted(tot.items()): print(f"  {c:12} {n:4d} files {b/1e6:6.1f} MB")
print()
for g in sorted(grp, key=lambda g: (-sum(v[1] for v in grp[g].values()))):
    print(f"{g:55} " + "  ".join(f"{c[:4]} {n}f/{b/1e6:.2f}MB" for c, (n, b) in sorted(grp[g].items())))

# ---- 3. dependency scan
import csv, os, re, subprocess, collections
T = OUT
inv = list(csv.DictReader(open(T + "inventory.csv")))
cat = {r["path"]: r["category"] for r in inv}
text = {}
for p in cat:
    if re.search(r"\.(R|py|qmd|sh|yaml|md|txt|tsv|csv|ipynb)$", p) and int([r for r in inv if r["path"] == p][0]["bytes"]) < 600000:
        text[p] = subprocess.run(["git", "show", f"{REF}:{p}"], capture_output=True, text=True).stdout
trivial = {"Snakefile", "config.yaml", "run_snakemake.sbatch", "README.md"}
dev = [p for p, c in cat.items() if c == "DEVELOPMENT" and os.path.basename(p) not in trivial and "/scenicplus/config" not in p]
dev_cfg = {int(m.group(1)) for p, c in cat.items() if c == "DEVELOPMENT" for m in [re.search(r"/scenicplus/config(\d+)/", p)] if m}
pub = [p for p, c in cat.items() if c == "PUBLICATION"]
hits = collections.defaultdict(list)
for q in pub:
    t = text.get(q, "")
    for p in dev:
        b = os.path.basename(p)
        if len(b) > 8 and b in t:
            hits[q].append(b)
    for n in sorted(dev_cfg):
        if re.search(rf"\bconfig{n}\b", t) and n != 37:
            hits[q].append(f"config{n}")
print("PUBLICATION files that mention a DEVELOPMENT file / config:")
for q, hs in sorted(hits.items()):
    cfgs = sorted({h for h in hs if h.startswith("config")}, key=lambda s: int(s[6:])); files = sorted({h for h in hs if not h.startswith("config")})
    print(f"  {q.replace('grn/', '')}: " + (", ".join(files) + ("; " if files and cfgs else "")) + (("configs " + ",".join(c[6:] for c in cfgs)) if cfgs else ""))
print(len(hits), "of", len(pub), "publication files have such a reference")
# which publication files read the sweep summary / parameter table specifically
for needle in ("scenicplus_sweep_summary.csv", "config_parameters.tsv", "scenicplus_coherence_tf_tiers.csv", "scenicplus_coherence_eregulons.csv"):
    r = [q.replace("grn/", "") for q, t in text.items() if needle in t and q != "grn/analysis/" + needle]
    print(f"  uses {needle}: {r[:6]}")

os.remove(OUT + 'files_raw.csv')
