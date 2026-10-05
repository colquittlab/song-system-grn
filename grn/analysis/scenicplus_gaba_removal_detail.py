#!/usr/bin/env python3
"""Homeodomain regulons under GABA / LGE region-set removal (configs 29-33 vs config1 / config11): where region-activity peaks go
and which homeodomain TFs lose their eRegulon. See scenicplus_homeodomain_region_check.py for the per-config summary."""
import re
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

R = Path("/hdd/jupyter/brad/scenicplus/motor-pathway_multiome/motor-pathway_multiome_seurat_cellbender.0.05_preprocess_cr/ra-arco-hvc-nc_hybrid/results")
A = Path("/ssd/brad/rstudio/multiome/song-system-grn/.claude/worktrees/grn-hybrid-rerun/")
s = lambda x: x.decode() if isinstance(x, bytes) else x
dbd = dict(zip(*pd.read_csv(A / "data/lambert_tfs/tfs.csv")[["Name", "DBD"]].values.T))
is_hd = lambda t: "homeodomain" in str(dbd.get(t, "")).lower()
ARMS = {"config1": "baseline strict", "config33": "ctrl strict", "config31": "noLGE strict", "config29": "noGABA strict",
        "config11": "baseline loose", "config32": "noLGE loose", "config30": "noGABA loose"}


def region_peaks(cfg):
    outs = R / cfg / "outs"
    full = outs / "scplusmdata.h5mu"
    f = h5py.File(full if full.exists() else outs / "scplusmdata_slim.h5mu", "r")
    o = f["mod"]["scRNA_counts"]["obs"]
    c2c = dict(zip((s(x) for x in o["_index"][:]), np.array([s(c) for c in o["cluster"]["categories"][:]])[o["cluster"]["codes"][:]]))
    rows = []
    for kind in ("direct", "extended"):
        m = f["mod"][f"{kind}_region_based_AUC"]
        names = [s(x) for x in m["var"]["_index"][:]]
        cl = np.array([c2c[s(x)] for x in m["obs"]["Cell"][:]])
        X = np.asarray(m["X"][:], dtype=np.float32)
        clusters = sorted(set(cl))
        M = np.vstack([X[cl == c].mean(0) for c in clusters])
        for j, n in enumerate(names):
            rows.append({"TF": n.split("_")[0], "reg": n.rsplit("_(", 1)[0], "peak": clusters[int(np.argmax(M[:, j]))]})
    f.close()
    return pd.DataFrame(rows)


def grp(c):
    return ("LGE" if c in ("GABA-LGE-1", "GABA-LGE-2") else "other GABA" if c.startswith("GABA") else "glutamatergic" if c.startswith("Glut") else "glia/other")


res = {c: region_peaks(c) for c in ARMS}
print("region-activity peak of HOMEODOMAIN regulons (share by group):")
print(f"{'arm':16s} {'n':>4s}  " + "  ".join(f"{g:>13s}" for g in ("LGE", "other GABA", "glutamatergic", "glia/other")) + "   top clusters")
for c, lab in ARMS.items():
    d = res[c]
    d = d[d.TF.map(is_hd)]
    g = d.peak.map(grp).value_counts(normalize=True)
    top = ", ".join(f"{k.replace('Glut-', '').replace('GABA-', 'G-')} {v:.0%}" for k, v in d.peak.value_counts(normalize=True).head(4).items())
    print(f"{lab:16s} {len(d):4d}  " + "  ".join(f"{g.get(k, 0):13.0%}" for k in ("LGE", "other GABA", "glutamatergic", "glia/other")) + f"   {top}")

print("\nhomeodomain TFs with any eRegulon:")
tfs = {c: {t for t in res[c].TF if is_hd(t)} for c in ARMS}
for c, lab in ARMS.items():
    print(f"  {lab:16s} {len(tfs[c])}")
for a, b, lab in (("config1", "config29", "strict: baseline -> noGABA"), ("config1", "config31", "strict: baseline -> noLGE"), ("config11", "config30", "loose: baseline -> noGABA"), ("config11", "config32", "loose: baseline -> noLGE")):
    lost, kept = sorted(tfs[a] - tfs[b]), sorted(tfs[a] & tfs[b])
    print(f"\n{lab}: {len(kept)} kept, {len(lost)} lost, {len(tfs[b] - tfs[a])} gained")
    print("   lost:", ", ".join(lost[:45]))
