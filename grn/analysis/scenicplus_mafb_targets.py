#!/usr/bin/env python3
"""MAFB target genes across configs: are fast-spiking (FS) genes such as PVALB and KCNC1 included, and under what settings?

For every config with results (all sign classes, direct + extended eRegulons, read from the merged/slim file), records
which genes are MAFB targets. Then:
  * membership of a curated FS gene list (Kv3/Kv1 channels, Nav, HCN, Na/K-ATPase, release machinery, ...) across the 20
    correct-input configs (1-11 and the looser-regime sweep), and in the legacy-input controls (12, 13)
  * for PVALB, KCNC1 and every FS gene: relation of inclusion to the six swept parameters AND to the size of MAFB's
    eRegulon (a bigger regulon includes more of everything)
  * a data-driven shortlist: MAFB targets that are high in both RA and MGE-PVALB clusters

Writes scenicplus_mafb_target_membership.csv (one row per gene, small) and prints a report.
The FS list is curated from the literature (see FS_GENES), not derived from these data.

    python scenicplus_mafb_targets.py
"""
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix

HERE = Path(__file__).resolve().parent
B = ("/hdd/jupyter/brad/scenicplus/motor-pathway_multiome/"
     "motor-pathway_multiome_seurat_cellbender.0.05_preprocess_cr/ra-arco-hvc-nc_hybrid/")
RES = Path(B + "results")
s = lambda x: x.decode() if isinstance(x, bytes) else x

# Curated: genes tied to the fast-spiking phenotype of PV interneurons (and of RA projection neurons where known)
FS_GENES = {
    "Kv3 channels": ["KCNC1", "KCNC2", "KCNC3", "KCNC4"],
    "Kv1 / other K+": ["KCNA1", "KCNA2", "KCNAB1", "KCNS3", "KCNMA1"],
    "Na+ channels": ["SCN1A", "SCN8A", "SCN1B", "SCN2B"],
    "HCN": ["HCN1", "HCN2"],
    "Na+/K+ pump": ["ATP1A3", "ATP1B1"],
    "Ca2+ buffer / channel": ["PVALB", "CACNA1A", "CACNB4"],
    "release machinery": ["SYT2", "CPLX1", "CPLX2", "SYT1", "STX1B", "SNAP25"],
    "fast AMPA/NMDA": ["GRIA1", "GRIA4", "GRIN2A"],
    "interneuron identity": ["ERBB4", "LHX6", "SOX6", "NKX2-1", "GAD1", "GAD2"],
}
# Several classic FS genes carry no symbol in lonStrDom2 (NCBI names them LOCxxxx); see fs_gene_aliases_lonStrDom2.tsv
# (mapped by sequence similarity and synteny against the newer ASM5065582v1 assembly).
ALIASES = dict(pd.read_csv(HERE / "fs_gene_aliases_lonStrDom2.tsv", sep="\t")[["symbol", "lonStrDom2_gene"]].values)
FS_FLAT = [g for v in FS_GENES.values() for g in v] + list(ALIASES.values())
PARAMS = ["ctx_nes_threshold", "ctx_auc_threshold", "ctx_rank_threshold", "dem_adj_pval_thr", "dem_log2fc_thr", "motif_similarity_fdr"]

P = pd.read_csv(HERE.parent / "ra-arco-hvc-nc/scenicplus/config_parameters.tsv", sep="\t").set_index("config")
cfgs = sorted((p.name for p in RES.glob("config*") if (p / "outs" / "scplusmdata.h5mu").exists() or (p / "outs" / "scplusmdata_slim.h5mu").exists()),
              key=lambda c: int(c[6:]))
correct = [c for c in cfgs if int(c[6:]) in list(range(1, 12)) + list(range(16, 29))]
assoc = [c for c in correct if c not in ("config9", "config10")]  # these two change non-motif parameters

# 1. MAFB rows per config
rows = []
for cfg in cfgs:
    outs = RES / cfg / "outs"
    full = outs / "scplusmdata.h5mu"
    f = h5py.File(full if full.exists() else outs / "scplusmdata_slim.h5mu", "r")
    for kind in ("direct", "extended"):
        g = f["uns"][f"{kind}_e_regulon_metadata"]
        tf = np.array([s(x) for x in g["TF"][:]])
        m = np.where(tf == "MAFB")[0]
        if not len(m):
            continue
        d = pd.DataFrame({k: [s(x) for x in g[k][:][m]] for k in ("Gene", "eRegulon_name", "Region")})
        for k in ("importance_TF2G", "rho_TF2G", "importance_R2G", "rho_R2G"):
            d[k] = g[k][:][m]
        d["config"], d["type"] = cfg, kind
        d["sign"] = d.eRegulon_name.str.extract(r"_(\+/\+|-/-|\+/-|-/\+)$")[0]
        rows.append(d)
    f.close()
M = pd.concat(rows, ignore_index=True)
size = M.groupby("config").Gene.nunique()                       # MAFB regulon size per config (union over classes)
print(f"MAFB eRegulon rows: {len(M)} over {M.config.nunique()} configs; target genes per config: {int(size.min())}-{int(size.max())}")

# per gene x config membership
M["pp"], M["pm"], M["is_direct"] = M.sign == "+/+", M.sign == "+/-", M.type == "direct"
mem = M.groupby(["Gene", "config"]).agg(pp=("pp", "any"), pm=("pm", "any"), direct=("is_direct", "any")).reset_index()
gc = mem[mem.config.isin(correct)].groupby("Gene").agg(n_configs=("config", "nunique"), n_pp=("pp", "sum"),
                                                      n_pm=("pm", "sum"), n_direct=("direct", "sum"))
genes = gc.join(mem[mem.config.isin(["config12", "config13"])].groupby("Gene").config.nunique().rename("in_legacy"), how="outer").fillna(0).astype(int).reset_index()

# 2. expression across clusters (new run's adata, .raw counts -> per-cell 1e4 -> log1p -> cluster mean) and detection
a = ad.read_h5ad(B + "anndata_rna/adata.h5ad")
raw = csr_matrix(a.raw.X)
tot = np.asarray(raw.sum(1)).ravel()
nrm = raw.multiply(1e4 / tot[:, None]).tocsr()
nrm.data = np.log1p(nrm.data)
clusters = sorted(a.obs.cluster.unique())
oh = csr_matrix((np.ones(len(a)), (np.arange(len(a)), [clusters.index(c) for c in a.obs.cluster])), shape=(len(a), len(clusters)))
n_cl = np.asarray(oh.sum(0)).ravel()
EX = pd.DataFrame((oh.T @ nrm).toarray() / n_cl[:, None], index=clusters, columns=list(a.raw.var_names))
Zx = (EX - EX.mean()) / EX.std()
PV, RA = ["GABA-MGE-PVALB-1", "GABA-MGE-PVALB-2"], "Glut-CACNA1H-RA"
genes["z_RA"] = genes.Gene.map(Zx.loc[RA]).round(2)
genes["z_PVALB"] = genes.Gene.map(Zx.loc[PV].mean()).round(2)
genes["expr_RA"] = genes.Gene.map(EX.loc[RA]).round(2)
genes["expr_PVALB"] = genes.Gene.map(EX.loc[PV].mean()).round(2)
genes["FS_curated"] = genes.Gene.isin(FS_FLAT)
genes = genes.sort_values(["n_configs", "n_pp"], ascending=False)
genes.to_csv(HERE / "scenicplus_mafb_target_membership.csv", index=False)

# 3. report: curated genes
N = len(correct)
print(f"\nMAFB target membership across {N} correct-input configs (and the 2 legacy-input controls):")
print(f"{'gene':8s} {'any':>4s} {'+/+':>4s} {'+/-':>4s} {'direct':>6s} {'legacy':>6s}   {'z RA':>5s} {'z PV':>5s}  group")
for grp, gl in FS_GENES.items():
    for gene in gl:
        label = gene
        if gene not in EX.columns and gene in ALIASES:
            gene = ALIASES[gene]                      # symbol absent in lonStrDom2: use its LOC id
        r = genes[genes.Gene == gene]
        if gene not in EX.columns:
            print(f"{label:8s}  (not in the finch gene set / expression matrix)  {grp}")
            continue
        r = r.iloc[0] if len(r) else None
        n = (int(r.n_configs), int(r.n_pp), int(r.n_pm), int(r.n_direct), int(r.in_legacy)) if r is not None else (0, 0, 0, 0, 0)
        print(f"{label:8s} {n[0]:4d} {n[1]:4d} {n[2]:4d} {n[3]:6d} {n[4]:6d}   {Zx.loc[RA, gene]:5.1f} {Zx.loc[PV, gene].mean():5.1f}  {grp}" + (f"  [= {gene}]" if gene != label else ""))

# 4. what explains inclusion across configs? parameters and MAFB regulon size
X = P.loc[assoc, PARAMS].astype(float)
sz = size.reindex(assoc).astype(float)


def auc(v, pres):
    a_, b_ = v[pres], v[~pres]
    return float(np.mean([(x > y) + 0.5 * (x == y) for x in a_ for y in b_])) if len(a_) and len(b_) else np.nan


print(f"\nWhat tracks inclusion (AUC; >0.5 = present configs have the higher value), {len(assoc)} configs:")
print(f"{'gene':7s} {'n':>3s} {'size':>5s} " + " ".join(f"{p.replace('ctx_', '').replace('_threshold', '').replace('dem_', '').replace('motif_similarity_', 'sim_')[:8]:>8s}" for p in PARAMS))
for gene in ("PVALB", "KCNC1"):
    present = pd.Series([bool(len(mem[(mem.Gene == gene) & (mem.config == c)])) for c in assoc], index=assoc)
    if present.sum() in (0, len(assoc)):
        print(gene, "present in", int(present.sum()), "of", len(assoc), "- no contrast")
        continue
    print(f"{gene:7s} {int(present.sum()):3d} {auc(sz.values, present.values):5.2f} " + " ".join(f"{auc(X[p].values, present.values):8.2f}" for p in PARAMS))
print(f"(MAFB regulon size spans {int(sz.min())}-{int(sz.max())} genes across these configs)")

# 5. data-driven: MAFB targets high in BOTH RA and MGE-PVALB
sh = genes[(genes.z_RA > 1) & (genes.z_PVALB > 1) & (genes.n_configs >= 10)].copy()
sh["shared"] = np.minimum(sh.z_RA, sh.z_PVALB)
print(f"\nMAFB targets (in >= 10 of {N} configs) high in both RA and MGE-PVALB (z > 1 in each): {len(sh)}")
print(sh.sort_values("shared", ascending=False).head(30)[["Gene", "n_configs", "n_pp", "z_RA", "z_PVALB", "expr_RA", "expr_PVALB", "FS_curated"]].to_string(index=False))


# 6. expression and detection of MAFB and the FS genes in RA, the two PVALB clusters, SST and the rest
nrm = nrm.tocsc()
names = list(a.raw.var_names)
cl = a.obs.cluster.values
groups = {"RA": ["Glut-CACNA1H-RA"], "PVALB-1": ["GABA-MGE-PVALB-1"], "PVALB-2": ["GABA-MGE-PVALB-2"], "SST": ["GABA-MGE-SST-1"],
          "other glut": [c for c in sorted(set(cl)) if c.startswith("Glut") and c != "Glut-CACNA1H-RA"],
          "other GABA": [c for c in sorted(set(cl)) if c.startswith("GABA") and "PVALB" not in c]}
print(f"{'gene':7s} " + "  ".join(f"{g:>16s}" for g in groups) + "   (mean log-norm expression, % cells detected)")
for gene in ("MAFB", "PVALB", "KCNC1", "KCNC2", "KCNS3", "LHX6", "CPLX1", "HCN1", "SCN1B", "SYT2"):
    col = nrm[:, names.index(gene)].toarray().ravel()
    out = []
    for g, cs in groups.items():
        m = np.isin(cl, cs)
        out.append(f"{col[m].mean():5.2f} ({100 * (col[m] > 0).mean():3.0f}%)")
    print(f"{gene:7s} " + "  ".join(f"{o:>16s}" for o in out))

# 2. threshold-free: rank of PVALB / KCNC1 among MAFB's TF-to-gene links
print("\nMAFB TF-to-gene links (all genes the model scored for MAFB): rank of the gene by importance, and its rho")
print(f"{'config':9s} {'links':>6s}  {'PVALB rank (pct)':>20s} {'rho':>6s}   {'KCNC1 rank (pct)':>20s} {'rho':>6s}")
for n in (1, 2, 16, 17, 18, 19, 21, 23, 24, 25, 26):
    p = RES / f"config{n}" / "outs" / "tf_to_gene_adj.tsv"
    if not p.exists():
        continue
    t = pd.read_csv(p, sep="\t")
    t = t[t.TF == "MAFB"].sort_values("importance", ascending=False).reset_index(drop=True)
    row = [f"config{n:<3d}", f"{len(t):6d}"]
    for gene in ("PVALB", "KCNC1"):
        r = t.index[t.target == gene]
        row.append(f"{int(r[0]) + 1:>9d} ({100 * (int(r[0]) + 1) / len(t):4.1f}%)  {t.loc[r[0], 'rho']:6.2f}" if len(r) else "      not scored")
    print("  ".join(row))

# 7. what decides inclusion? PVALB and KCNC1 versus their importance rank in MAFB's TF-to-gene links. Their correlation
# with MAFB is the same in every config, but the TF-to-gene model is refit per config with that config's candidate-TF list
# (tfs.txt, 482-669 TFs), so the importance rank moves; inclusion is a rank cutoff inside that ranking.
print("\nInclusion vs TF-to-gene importance percentile (lower = more strongly predicted by MAFB), configs with the table:")
rows7 = []
for n in (1, 2) + tuple(range(16, 29)):
    cfg = f"config{n}"
    adj = pd.read_csv(RES / cfg / "outs" / "tf_to_gene_adj.tsv", sep="\t")
    adj = adj[adj.TF == "MAFB"].sort_values("importance", ascending=False).reset_index(drop=True)
    inc = set(mem[mem.config == cfg].Gene)
    r = {"config": cfg, "candidate_TFs": sum(1 for line in open(RES / cfg / "outs" / "tfs.txt") if line.strip())}
    for gene in ("PVALB", "KCNC1"):
        hit = adj.index[adj.target == gene]
        r[f"{gene}_pct"] = round(100 * (hit[0] + 1) / len(adj), 1) if len(hit) else np.nan   # nan: gene not scored for MAFB
        r[f"{gene}_included"] = gene in inc
    rows7.append(r)
d7 = pd.DataFrame(rows7).dropna(subset=["PVALB_pct", "KCNC1_pct"]).sort_values("KCNC1_pct")
print(d7.to_string(index=False))
for gene in ("PVALB", "KCNC1"):
    i, o = d7[d7[f"{gene}_included"]][f"{gene}_pct"], d7[~d7[f"{gene}_included"]][f"{gene}_pct"]
    print(f"  {gene}: percentile when included {i.min():.1f}-{i.max():.1f} (n={len(i)}); when not included "
          + (f"{o.min():.1f}-{o.max():.1f} (n={len(o)})" if len(o) else "never"))
print("  Spearman, candidate-TF count vs KCNC1 percentile: %.2f; vs PVALB percentile: %.2f" % (
    d7.candidate_TFs.corr(d7.KCNC1_pct, method="spearman"), d7.candidate_TFs.corr(d7.PVALB_pct, method="spearman")))

# 8. topic-model configs, outside the aggregates above (different cisTopic object): same readout, each next to its reference.
# 14/15 = 15/30 topics; 34/35 = 40 topics with the song-pair DAR sets removed (34 at config1's thresholds, 35 at config11's).
print("\nTopic-model configs (not in the aggregates above):")
rows8 = []
for n in (1, 14, 15, 34, 11, 35):
    cfg = f"config{n}"
    if not (RES / cfg / "outs" / "tf_to_gene_adj.tsv").exists():
        print(f"  {cfg}: tf_to_gene_adj.tsv not transferred, skipped")
        continue
    adj = pd.read_csv(RES / cfg / "outs" / "tf_to_gene_adj.tsv", sep="\t")
    n_tfs = adj.TF.nunique()   # the candidate-TF list (tfs.txt is not always transferred)
    adj = adj[adj.TF == "MAFB"].sort_values("importance", ascending=False).reset_index(drop=True)
    inc = set(mem[mem.config == cfg].Gene)
    r = {"config": cfg, "n_topics": P.loc[cfg, "n_topics"], "MAFB_targets": int(size.get(cfg, 0)), "candidate_TFs": n_tfs}
    for gene in ("PVALB", "KCNC1"):
        hit = adj.index[adj.target == gene]
        r[f"{gene}_pct"] = round(100 * (hit[0] + 1) / len(adj), 1) if len(hit) else np.nan
        r[f"{gene}_included"] = gene in inc
    rows8.append(r)
print(pd.DataFrame(rows8).to_string(index=False))
