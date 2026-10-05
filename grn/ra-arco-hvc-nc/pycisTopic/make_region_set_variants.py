#!/usr/bin/env python3
"""Region-set folders with interneuron sets removed, to test how much they drive the homeodomain eRegulons.

Motif enrichment runs once per region set (DAR sets, topic sets). For homeodomain TFs the enriched cistromes come from
GABA sets (the motif is genuinely active in interneuron chromatin), and every TF annotated to that shared motif
inherits those LGE-open regions even where its own targets are expressed elsewhere (ALX4 in HVCx, EMX2 in astro).
These variants take the GABA signal out of the motif-enrichment input to see whether the same TFs then get cistromes
from the lineages where they are expressed. Same production cisTopic object, same expression; only the folder differs.

    region_sets_noGABA   all 9 GABA DAR sets, and the GABA topics (>= 50% of the topic's mass in GABA clusters: topics
                         1, 15, 20) from both topic sets. Topics 19 and 6 (41% and 40% GABA, diffuse; 6 is a Glut-NB /
                         GABA-Im progenitor mix) are borderline and stay in.
    region_sets_noLGE    only the two LGE DAR sets (GABA-LGE-1 and GABA-LGE-2) and Topic20, the topic carried by exactly
                         those two clusters (31% LGE-1, 30% LGE-2). Other interneuron sets stay, so this isolates LGE.
    region_sets_ctrl     size-matched control for noGABA: 9 randomly chosen non-GABA, non-focal DAR sets and 3 randomly
                         chosen non-GABA (< 30% GABA mass), non-focal topics. Removes as much as noGABA but nothing from the lineages of
                         interest (HVCx, HVCra, RA, astro/NSC), so a difference from noGABA is not just "fewer sets".

DARs_song-pairs is always kept. Each folder gets a REMOVED_SETS.tsv. Topic composition comes from the production
20-topic model (topic_cluster_composition_k20.tsv, kept in this directory).

    python make_region_set_variants.py
"""
import random
import shutil
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
P = Path("/hdd/jupyter/brad/scenicplus/motor-pathway_multiome/"
         "motor-pathway_multiome_seurat_cellbender.0.05_preprocess_cr/ra-arco-hvc-nc_hybrid/pycisTopic")
SRC = P / "region_sets"
SEED = 2026
SUBDIRS = ("DARs_all", "DARs_song-pairs", "Topics_otsu", "Topics_top_3k")
FOCAL_DARS = {"Glut-DACH2-HVCx", "Glut-DACH2-HVCra", "Glut-CACNA1H-RA", "Astro", "Epen", "Glut-NSC"}

comp = pd.read_csv(HERE / "topic_cluster_composition_k20.tsv", sep="\t", index_col=0)  # topics x clusters, share of topic mass
gaba_cols = [c for c in comp.columns if c.startswith("GABA")]
gaba_share = comp[gaba_cols].sum(1)
top = comp.idxmax(1)
gaba_topics = sorted([t for t in comp.index if gaba_share[t] >= 0.5], key=lambda t: int(t[5:]))
lge_topics = sorted([t for t in comp.index if top[t] in ("GABA-LGE-1", "GABA-LGE-2")], key=lambda t: int(t[5:]))
focal_topics = sorted([t for t in comp.index if top[t] in FOCAL_DARS or top[t] in ("Glut-DACH2-HVCx", "Glut-DACH2-HVCra", "Glut-CACNA1H-RA", "Astro")],
                      key=lambda t: int(t[5:]))
dar_all = sorted(p.stem for p in (SRC / "DARs_all").glob("*.bed"))
gaba_dars = [d for d in dar_all if d.startswith("GABA")]

rng = random.Random(SEED)
ctrl_dars = sorted(rng.sample([d for d in dar_all if d not in gaba_dars and d not in FOCAL_DARS], len(gaba_dars)))
ctrl_topics = sorted(rng.sample([t for t in comp.index if gaba_share[t] < 0.3 and t not in focal_topics], len(gaba_topics)), key=lambda t: int(t[5:]))  # clearly non-GABA only

VARIANTS = {
    "region_sets_noGABA": {"DARs_all": gaba_dars, "Topics_otsu": gaba_topics, "Topics_top_3k": gaba_topics},
    "region_sets_noLGE": {"DARs_all": ["GABA-LGE-1", "GABA-LGE-2"], "Topics_otsu": lge_topics, "Topics_top_3k": lge_topics},
    "region_sets_ctrl": {"DARs_all": ctrl_dars, "Topics_otsu": ctrl_topics, "Topics_top_3k": ctrl_topics},
}
print("GABA topics:", gaba_topics, "| LGE topics:", lge_topics, "| focal topics (protected in ctrl):", focal_topics)
print("control removes DAR sets:", ctrl_dars, "and topics:", ctrl_topics)

for name, removed in VARIANTS.items():
    dst = P / name
    shutil.rmtree(dst, ignore_errors=True)
    rows = []
    for sub in SUBDIRS:
        (dst / sub).mkdir(parents=True)
        drop = set(removed.get(sub, []))
        for bed in sorted((SRC / sub).glob("*.bed")):
            if bed.stem in drop:
                rows.append((sub, bed.stem, name))
            else:
                shutil.copy2(bed, dst / sub / bed.name)
        assert {r[1] for r in rows if r[0] == sub} == drop, (name, sub, drop)
    pd.DataFrame(rows, columns=["subdir", "set", "removed_by"]).to_csv(dst / "REMOVED_SETS.tsv", sep="\t", index=False)
    counts = {s: len(list((dst / s).glob("*.bed"))) for s in SUBDIRS}
    print(f"{name}: kept {counts}; removed {len(rows)} sets")
