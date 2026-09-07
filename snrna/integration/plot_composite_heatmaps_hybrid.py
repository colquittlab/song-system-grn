"""Re-render every composite / per-method heatmap for the two full-suite hybrid-label
composites (finch x chicken: gg_adult_hybrid; finch x mouse: yao_adult_hybrid) with the
vendored toolkit/plot_rank_heatmap.py.

Self-contained within this repo: all paths are relative to this file. Inputs are the
matrices written by assemble_gg_adult_hybrid.py / assemble_yao_adult_hybrid.py into
composite_scoring/results/<tag>/ (tracked CSVs), the Colquitt-2021-method GSI variants in
composite_scoring/results/, and the reference-label annotation tables in
composite_scoring/annotations/. Python env: envs/integration_plots.yaml. Outputs (PDF+PNG, gitignored)
overwrite the files of the same name in each results dir. Titles/flags reproduce the original 2025-08-19 renders (recovered
from the PDF title text and page geometry; the --scale on the pt5/pt6 miniatures is an
estimate from page size, ~0.6).

Usage:  python plot_composite_heatmaps_hybrid.py [gg_adult_hybrid|yao_adult_hybrid] [name ...]
        (no args = everything)
"""
import subprocess, sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
TOOL = HERE / "toolkit" / "plot_rank_heatmap.py"
RES = HERE / "composite_scoring" / "results"
ANN_DIR = HERE / "composite_scoring" / "annotations"
ANNOT = {
    "gg_adult_hybrid": ANN_DIR / "gg_adult_label_annotation.csv",
    "yao_adult_hybrid": ANN_DIR / "yao_label_annotation.csv",
}
METHOD_TITLES = {
    "gsi": ("GSI (correlation)", ["--signed"]),
    "samap": ("SAMap alignment score", []),
    "cca": ("Seurat CCA transfer score", []),
    "saturn": ("SATURN transfer score", []),
}
GSI_C2021 = {
    "all": "GSI, Colquitt 2021 method (DEG markers, all genes)",
    "nontf": "GSI, Colquitt 2021 method (variable genes, non-TF)",
    "tf": "GSI, Colquitt 2021 method (variable genes, TF-only)",
}


def jobs(tag):
    d = RES / tag
    conf = ["--matrix", d / "composite_confidence_matrix.csv", "--cbar_label", "mapping confidence",
            "--dots_csv", d / "composite_agreement_count_matrix.csv"]
    J = {
        "rank_score_clustermap": ["--matrix", d / "composite_rank_score_matrix.csv",
                                  "--title", "Composite rank-aggregate score"],
    }
    if tag == "gg_adult_hybrid":
        J["confidence_clustermap"] = conf + ["--title", "Composite mapping confidence"]
        for k, t in GSI_C2021.items():
            J[f"method_gsi_colquitt2021_{k}_heatmap"] = [
                "--matrix", RES / f"gsi_corr_gg_adult_hybrid_colquitt2021_{k}.csv", "--title", t, "--signed"]
    else:
        t = ["--title", "Composite mapping confidence (finch x Yao mouse)"]
        J["confidence_clustermap"] = conf + t
        J["confidence_clustermap_transposed"] = conf + t + ["--transpose", "--scale", "0.8"]
        J["confidence_clustermap_topk1_anchored"] = conf + t + ["--top_k", "1"]
        J["confidence_clustermap_topk1_anchored_transposed"] = conf + t + ["--top_k", "1", "--transpose"]
        for pt in (5, 6):
            J[f"confidence_clustermap_topk1_anchored_transposed_pt{pt}"] = conf + t + [
                "--top_k", "1", "--transpose", "--scale", "0.6", "--label_pt_override", str(pt)]
    for m, (title, extra) in METHOD_TITLES.items():
        J[f"method_{m}_heatmap"] = ["--matrix", d / f"method_{m}_matrix.csv", "--title", title] + extra
    return J


def main(argv):
    tags = [a for a in argv if a in ANNOT] or list(ANNOT)
    only = [a for a in argv if a not in ANNOT]
    for tag in tags:
        for name, args in jobs(tag).items():
            if only and name not in only:
                continue
            cmd = [sys.executable, TOOL, "--out_prefix", RES / tag / name, "--annot_csv", ANNOT[tag]] + args
            print(f"\n>>> {tag}/{name}", flush=True)
            subprocess.run([str(c) for c in cmd], check=True)


if __name__ == "__main__":
    main(sys.argv[1:])
