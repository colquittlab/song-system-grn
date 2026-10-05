#!/usr/bin/env python3
"""Did the earlier run's double-transformed `.raw` drive the difference between the old and new SCENIC+ networks?

Controls: config12 / config13 have the same parameters as config1 / config2 and the same new labels, cells and
region sets, but read adata_doublenorm.h5ad (.raw = log1p(normalize_total(SCT data)), the earlier run's matrix).
Only the expression input differs, so, on positive (+/+) eRegulons:

    12 vs 1  (13 vs 2)    effect of the expression input alone, everything else equal
    12 vs old config3     does legacy input + new labels reproduce the earlier network? (old config3 used config1's
                          thresholds; the earlier run's other configs are shown as a range)
    1 vs old config3      the gap to explain (edge Jaccard 0.11 in the first comparison)

Reading it: if config12 lands near old config3 (edge Jaccard well above 0.11) and far from config1, the input
explains most of the gap; if config12 stays near config1, labels / cells / region sets do. Jaccard on TF sets and on
TF-gene edges; MAFB and EMX2 status alongside. Run after config12/13 results are in the results directory and
all.qmd has exported their tables (all.qmd discovers configN automatically).

    python scenicplus_normalization_control.py [new_out_dir] [old_out_dir]
"""
import sys
from pathlib import Path

import pandas as pd

BASE = Path("~/ssd/rstudio/multiome/motor-pathway/scenicplus").expanduser()
NEW = Path(sys.argv[1]).expanduser() if len(sys.argv) > 1 else BASE / "motor-pathway_scenicplus_v2_hybrid_all"
OLD = Path(sys.argv[2]).expanduser() if len(sys.argv) > 2 else BASE / "motor-pathway_scenicplus_v2_cluster-snrna-cr_all"
OLD_REF = "config3"  # earlier run's config with config1's thresholds


def load(d, cfg):
    f = d / cfg / "scenicplus_eRegulons.txt"
    if not f.exists():
        return None
    t = pd.read_csv(f, sep="\t", usecols=["TF", "Gene", "eRegulon_name"])
    t = t[t.eRegulon_name.str.contains(r"\+/\+")]
    return set(t.TF), set(map(tuple, t[["TF", "Gene"]].drop_duplicates().values))


def jac(a, b):
    return len(a & b) / len(a | b)


nets = {f"new {c}": load(NEW, c) for c in ("config1", "config2", "config12", "config13")}
nets[f"old {OLD_REF}"] = load(OLD, OLD_REF)
old_all = {p.name: load(OLD, p.name) for p in OLD.glob("config[0-9]*") if (p / "scenicplus_eRegulons.txt").exists()}

print("positive-eRegulon networks:")
for k, v in nets.items():
    print(f"  {k:14s}", "not available yet" if v is None else f"{len(v[0])} TFs, {len(v[1])} edges, MAFB {'MAFB' in v[0]}, EMX2 {'EMX2' in v[0]}")
print(f"  earlier run, all {len(old_all)} configs: EMX2 in {sum('EMX2' in v[0] for v in old_all.values())}, "
      f"edges {min(len(v[1]) for v in old_all.values())}-{max(len(v[1]) for v in old_all.values())}")

pairs = [("new config12", "new config1", "input alone (cfg1 params)"),
         ("new config13", "new config2", "input alone (cfg2 params)"),
         ("new config12", f"old {OLD_REF}", "legacy input + new labels vs earlier run"),
         ("new config1", f"old {OLD_REF}", "the gap to explain")]
rows = []
for a, b, what in pairs:
    if nets[a] is None or nets[b] is None:
        rows.append({"comparison": f"{a} vs {b}", "meaning": what, "TF_jaccard": "n/a", "edge_jaccard": "n/a"})
        continue
    rows.append({"comparison": f"{a} vs {b}", "meaning": what,
                 "TF_jaccard": round(jac(nets[a][0], nets[b][0]), 2), "edge_jaccard": round(jac(nets[a][1], nets[b][1]), 2)})
# config13 against every earlier config: the earlier run's closest match, since its thresholds are not config2's
if nets["new config13"] is not None:
    best = max(old_all.items(), key=lambda kv: jac(nets["new config13"][1], kv[1][1]))
    rows.append({"comparison": f"new config13 vs best earlier ({best[0]})", "meaning": "legacy input + new labels vs earlier run",
                 "TF_jaccard": round(jac(nets["new config13"][0], best[1][0]), 2),
                 "edge_jaccard": round(jac(nets["new config13"][1], best[1][1]), 2)})
out = pd.DataFrame(rows)
pd.set_option("display.width", 200, "display.max_colwidth", 60)
print("\n" + out.to_string(index=False))
