"""Check the song-pair DAR sets in region_sets_k40_archr against the raw fragments of the cisTopic object: share of DARs genuinely
more accessible in the foreground (>=5% of cells and >=1.5x), per-cluster detection, and DAR/random-region enrichment restricted to
10-40k-fragment cells and within a single library (so depth and library cannot explain it). Run: python validate_song_pair_dars.py"""
import pickle
import numpy as np
import pandas as pd

S = "/hdd/jupyter/brad/scenicplus/motor-pathway_multiome/motor-pathway_multiome_seurat_cellbender.0.05_preprocess_cr/ra-arco-hvc-nc_hybrid/pycisTopic/"
o = pickle.load(open(S + "cistopic_obj_glut_k40.pkl", "rb"))
cd = o.cell_data.loc[o.cell_names]
cl = cd["cluster"].astype(str).values
lib = cd["sample_short"].astype(str).values
fm = o.fragment_matrix.tocsr()
det = (fm > 0).astype(np.float32).tocsr()
depth = np.asarray(fm.sum(0)).ravel()
rn = pd.Series(np.arange(len(o.region_names)), index=o.region_names)
base = np.random.default_rng(1).choice(len(o.region_names), 20000, replace=False)
rnd = np.asarray(det[base].mean(0)).ravel()
focal = ["Glut-CACNA1H-RA", "Glut-CACNA1H-1", "Glut-CACNA1H-2", "Glut-DACH2-HVCra", "Glut-DACH2-HVCx", "Glut-DACH2-1", "Glut-DACH2-4"]
dm = (depth >= 10000) & (depth <= 40000)

for fg, bg in (("Glut-CACNA1H-RA", "Glut-CACNA1H-1"), ("Glut-DACH2-HVCra", "Glut-DACH2-1"), ("Glut-DACH2-HVCx", "Glut-DACH2-4")):
    name = f"{fg}_VS_{bg}"
    b = pd.read_csv(S + f"region_sets_k40_archr/DARs_song-pairs/{name}.bed", sep="\t", header=None, names=["c", "s", "e"])
    idx = rn.reindex(b.c + ":" + b.s.astype(str) + "-" + b.e.astype(str)).dropna().astype(int).values
    fgi, bgi = np.where(cl == fg)[0], np.where(cl == bg)[0]
    pf = np.asarray(det[idx][:, fgi].mean(1)).ravel()
    pb = np.asarray(det[idx][:, bgi].mean(1)).ravel()
    up = (pf >= 0.05) & (pf >= 1.5 * pb)
    dar = np.asarray(det[idx].mean(0)).ravel()
    enr = dar / np.maximum(rnd, 1e-3)
    print(f"\n===== {name}: {len(idx)} DARs | genuinely more accessible in {fg} by raw detection (>=5% of cells, >=1.5x): {up.mean():.1%}")
    prof = pd.Series(dar, index=o.cell_names).groupby(cl).mean().sort_values(ascending=False)
    print("   raw detection of the set, top clusters:", dict(prof.head(5).round(3)), "| fg %.3f bg %.3f" % (prof[fg], prof[bg]))
    t = pd.DataFrame({"cl": cl, "lib": lib, "enr": enr, "dm": dm})
    rows = {c: round(t[(t.cl == c) & t.dm].enr.median(), 2) for c in focal if ((t.cl == c) & t.dm).sum() >= 15}
    print("   depth-restricted (10-40k frags) median DAR/random enrichment:", rows)
    for lb in sorted(set(lib[fgi])):
        a = t[(t.lib == lb) & (t.cl == fg) & t.dm].enr
        c1 = t[(t.lib == lb) & (t.cl == bg) & t.dm].enr
        if len(a) >= 15 and len(c1) >= 10:
            print(f"   same library {lb}: {fg} {a.median():.2f} (n={len(a)}) vs {bg} {c1.median():.2f} (n={len(c1)})")
