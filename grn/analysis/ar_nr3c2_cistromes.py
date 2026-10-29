#!/usr/bin/env python3
"""AR vs NR3C2 in the RA vs C1H-1 runs: how similar are their cistromes (the regions with motif hits that SCENIC+ assigns to each TF)?

For configs 1 and 2 of glut-ra_glut-c1h-1_hybrid (run3 excluded) and glut-ra_glut-c1h-1_allruns_hybrid, read cistromes_direct/extended.h5ad (regions x TFs),
compare AR's cistrome with every other TF's (Jaccard) and with NR3C2's, and check whether the regions linked to MAFB are in AR's / NR3C2's cistromes.
Finding: AR, NR3C2, VDR, NR1I2 and ZNF626 have IDENTICAL direct cistromes (they share a motif annotation), PGR is ~0.98-0.99, NR3C1 0.37-0.44; so
the binding-site evidence cannot tell AR from NR3C2, and the two eRegulons differ only through the TF-to-gene (expression) step.

    python ar_nr3c2_cistromes.py        # writes ar_nr3c2_cistromes.csv next to this script
Note: the cistrome matrices must be cast to float before the product (a boolean sparse product collapses every overlap to 1).
"""
warnings.filterwarnings("ignore")
HB = "/hdd/jupyter/brad/scenicplus/motor-pathway_multiome/motor-pathway_multiome_seurat_cellbender.0.05_preprocess_cr/"
D = {"RA_run2": "glut-ra_glut-c1h-1_hybrid", "RA_all": "glut-ra_glut-c1h-1_allruns_hybrid"}
cols = ["TF", "Gene", "Region"]
out = []
for ds, d in D.items():
    for n in (1, 2):
        er = pd.concat([pd.read_table(f"{HB}{d}/results/config{n}/outs/{f}", usecols=cols) for f in ("eRegulon_direct.tsv", "eRegulons_extended.tsv")])
        mafb_regions = set(er[er.Gene == "MAFB"].Region); nr3c2_mafb = set(er[(er.TF == "NR3C2") & (er.Gene == "MAFB")].Region)
        for kind in ("direct", "extended"):
            a = ad.read_h5ad(f"{HB}{d}/results/config{n}/outs/cistromes_{kind}.h5ad")
            X = sp.csc_matrix((a.X > 0).astype(np.float32)); tfs = list(a.var_names); regs = pd.Index(a.obs_names)
            if "AR" not in tfs or "NR3C2" not in tfs:
                print(f"{ds} c{n} {kind}: AR in cistromes {'AR' in tfs}, NR3C2 in cistromes {'NR3C2' in tfs}"); continue
            size = np.asarray(X.sum(0)).ravel()
            iAR, iNR = tfs.index("AR"), tfs.index("NR3C2")
            ov = np.asarray((X[:, iAR].T @ X).todense()).ravel()                 # |AR cistrome  intersect  each TF's cistrome|
            jac = ov / (size[iAR] + size - ov)
            frac_of_AR = ov / size[iAR]
            order = np.argsort(-jac); rank = {tfs[i]: int(np.where(order == i)[0][0]) for i in range(len(tfs)) if i != iAR}   # rank among others (0 = most similar)
            others = [i for i in range(len(tfs)) if i != iAR]
            ranks = {t: rank[t] for t in ("NR3C2", "NR3C1", "PGR") if t in tfs}
            top5 = [tfs[i] for i in order if i != iAR][:5]
            inAR = set(regs[np.asarray(X[:, iAR].todense()).ravel() > 0]); inNR = set(regs[np.asarray(X[:, iNR].todense()).ravel() > 0])
            row = dict(dataset=ds, config=n, kind=kind, n_regions=len(regs), n_TF=len(tfs), AR_cistrome=int(size[iAR]), NR3C2_cistrome=int(size[iNR]), shared=int(ov[iNR]),
                       pct_of_AR_in_NR3C2=round(100 * frac_of_AR[iNR], 1), jaccard=round(jac[iNR], 3), rank_of_NR3C2_by_jaccard=f"{ranks['NR3C2'] + 1}/{len(tfs) - 1}",
                       median_jaccard_AR_vs_all=round(float(np.median(jac[others])), 4), top5_similar_to_AR=",".join(top5),
                       MAFB_regions=len(mafb_regions), MAFB_regions_in_AR_cistrome=len(mafb_regions & inAR), MAFB_regions_in_NR3C2_cistrome=len(mafb_regions & inNR),
                       NR3C2toMAFB_regions=len(nr3c2_mafb), NR3C2toMAFB_in_AR_cistrome=len(nr3c2_mafb & inAR))
            out.append(row)
R = pd.DataFrame(out)
pd.set_option("display.width", 250); pd.set_option("display.max_columns", 30); pd.set_option("display.max_colwidth", 40)
print(R[["dataset", "config", "kind", "n_TF", "AR_cistrome", "NR3C2_cistrome", "shared", "pct_of_AR_in_NR3C2", "jaccard", "rank_of_NR3C2_by_jaccard", "median_jaccard_AR_vs_all"]].to_string(index=False))
print()
print(R[["dataset", "config", "kind", "top5_similar_to_AR"]].to_string(index=False))
print()
print(R[["dataset", "config", "kind", "MAFB_regions", "MAFB_regions_in_AR_cistrome", "MAFB_regions_in_NR3C2_cistrome", "NR3C2toMAFB_regions", "NR3C2toMAFB_in_AR_cistrome"]].to_string(index=False))
R.to_csv(os.path.join(os.path.dirname(os.path.abspath(__file__)), "ar_nr3c2_cistromes.csv"), index=False)
