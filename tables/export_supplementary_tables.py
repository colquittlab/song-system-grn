#!/usr/bin/env python
"""Export the paper's supplementary tables as CSV (default) or one .xlsx workbook per table.

The table list follows the Google Sheet "adult-multiome_science_tables". CSV output is a
folder per table with one CSV per sheet and a README.txt. Each
workbook opens with a README sheet giving the source files, the method and
thresholds, and the git commit of this script, followed by one or more data
sheets.

Inputs are read read-only. Most of them are gitignored or untracked analysis
outputs that exist only in the main working tree, so paths resolve against
DATA_ROOT (the main checkout by default), not the checkout this script runs from.

    python tables/export_supplementary_tables.py [--format csv|xlsx|both] [--out DIR] [--only 01,07]

Environment overrides: SSG_DATA_ROOT, SSG_SCENICPLUS_ROOT, SSG_TABLE_OUT.
Requires pandas and xlsxwriter (the `cpdb` conda env has both).
"""

import argparse
import datetime
import os
import subprocess
from pathlib import Path

import pandas as pd
import xlsxwriter

SCRIPT_DIR = Path(__file__).resolve().parent
DATA_ROOT = Path(os.environ.get(
    "SSG_DATA_ROOT", "/ssd/brad/rstudio/multiome/song-system-grn"))
SCPLUS_ROOT = Path(os.environ.get(
    "SSG_SCENICPLUS_ROOT",
    "/hdd/jupyter/brad/scenicplus/motor-pathway_multiome/"
    "motor-pathway_multiome_seurat_cellbender.0.05_preprocess_cr"))
DEFAULT_OUT = Path(os.environ.get("SSG_TABLE_OUT", DATA_ROOT / "tables" / "output"))

INTEG = DATA_ROOT / "snrna/integration"
COMPOSITE = INTEG / "composite_scoring/results"
GLUT_COMPOSITE = INTEG / "gg_glutonly_hybrid/composite/results/composite_plusSATB2_noPre"
LABEL_LOOKUP = DATA_ROOT / "snrna/naming/hybrid_division_naming/hybrid_label_lookup.csv"
DEG_DIR = DATA_ROOT / "snrna/deg/song-surround_deg_glut_hybrid"
DAR_BEDS = DATA_ROOT / "multiome/archr/peaks_glut_hybrid/beds"

# Song type -> non-song counterpart, as tested in song-surround_deg_glut_hybrid.qmd
SONG_PAIRS = {
    "Glut-CACNA1H-RA": "Glut-CACNA1H-1",
    "Glut-DACH2-HVCra": "Glut-DACH2-1",
    "Glut-DACH2-HVCx": "Glut-DACH2-4",
    "Glut-DACH2-LMANco": "Glut-DACH2-8",
    "Glut-DACH2-LMANsh": "Glut-DACH2-5",
}
SHORT = {"Glut-CACNA1H-RA": "RA", "Glut-DACH2-HVCra": "HVCra", "Glut-DACH2-HVCx": "HVCx",
         "Glut-DACH2-LMANco": "LMANco", "Glut-DACH2-LMANsh": "LMANsh"}

EREG_SIGNS = ("+/+", "-/-")


# --------------------------------------------------------------------------- helpers

def git_commit():
    try:
        sha = subprocess.run(["git", "-C", str(SCRIPT_DIR), "rev-parse", "--short", "HEAD"],
                             capture_output=True, text=True, check=True).stdout.strip()
        dirty = subprocess.run(["git", "-C", str(SCRIPT_DIR), "status", "--porcelain", "--",
                                str(Path(__file__).resolve())],
                               capture_output=True, text=True).stdout.strip()
        return sha + (" (script modified)" if dirty else "")
    except Exception:
        return "unknown"


def relabel(names):
    """Old cluster IDs -> hybrid labels, as relabel_finch_clusters() does for the figures."""
    lut = pd.read_csv(LABEL_LOOKUP).dropna(subset=["new_label"])
    m = dict(zip(lut["cluster"], lut["new_label"]))
    out = [m.get(n, n) for n in names]
    dup = pd.Series(out)[pd.Series(out).duplicated()].tolist()
    if dup:
        raise ValueError(f"relabeling produced duplicate finch clusters: {dup}")
    return out


def need(*paths):
    for p in paths:
        if not Path(p).exists():
            raise FileNotFoundError(p)
    return paths


class XlsxBook:
    """One workbook: README first, then data sheets written row by row."""

    def __init__(self, path, title, description, sources, notes=()):
        self.path = path
        self.wb = xlsxwriter.Workbook(str(path), {"constant_memory": True,
                                                  "strings_to_numbers": False,
                                                  "strings_to_urls": False})
        self.wb.formats[0].set_font_name("Arial")
        self.wb.formats[0].set_font_size(10)
        self.bold = self.wb.add_format({"bold": True, "font_name": "Arial", "font_size": 10})
        self.title_fmt = self.wb.add_format({"bold": True, "font_name": "Arial", "font_size": 12})
        self.wrap = self.wb.add_format({"font_name": "Arial", "font_size": 10, "text_wrap": True,
                                        "valign": "top"})
        self.readme = self.wb.add_worksheet("README")
        self.readme.set_column(0, 0, 18)
        self.readme.set_column(1, 1, 110)
        self.readme.write(0, 0, title, self.title_fmt)
        self.readme.write(2, 0, "Description", self.bold)
        self.readme.write(2, 1, description, self.wrap)
        self._row = 4
        self._block("Sources", [str(s) for s in sources])
        self._block("Notes", list(notes))
        self._sheets = []
        self.sources = sources

    def _block(self, label, lines):
        if not lines:
            return
        self.readme.write(self._row, 0, label, self.bold)
        for line in lines:
            self.readme.write(self._row, 1, line, self.wrap)
            self._row += 1
        self._row += 1

    def sheet(self, name, df, description=""):
        assert len(name) <= 31, name
        ws = self.wb.add_worksheet(name)
        cols = list(df.columns)
        ws.write_row(0, 0, cols, self.bold)
        for i, c in enumerate(cols):
            ws.set_column(i, i, min(max(len(str(c)) + 2, 10), 40))
        ws.freeze_panes(1, 1 if cols and str(cols[0]) in ("finch_cluster", "query", "gene") else 0)
        ws.autofilter(0, 0, max(len(df), 1), len(cols) - 1)
        for r, row in enumerate(df.itertuples(index=False, name=None), start=1):
            ws.write_row(r, 0, [None if (isinstance(v, float) and v != v) else v for v in row])
        self._sheets.append((name, len(df), len(cols), description))

    def close(self):
        lines = [f"{n}: {r:,} rows x {c} columns" + (f". {d}" if d else "")
                 for n, r, c, d in self._sheets]
        self._block("Sheets", lines)
        self._block("Provenance", [
            f"Generated {datetime.date.today().isoformat()} by tables/export_supplementary_tables.py "
            f"(song-system-grn, commit {git_commit()}).",
            "Missing values are left blank."])
        self.wb.close()
        print(f"wrote {self.path} ({self.path.stat().st_size / 1e6:.1f} MB): "
              + ", ".join(f"{n}={r}" for n, r, _, _ in self._sheets))


class CsvBook:
    """Same interface as XlsxBook: a folder per table, one UTF-8 CSV per sheet, README.txt."""

    def __init__(self, path, title, description, sources, notes=()):
        self.dir = path.with_suffix("")
        self.dir.mkdir(parents=True, exist_ok=True)
        for old in self.dir.glob("*.csv"):
            old.unlink()
        self.lines = [title, "=" * len(title), "", description, ""]
        self._block("Sources", [str(s) for s in sources])
        self._block("Notes", list(notes))
        self._sheets = []

    def _block(self, label, lines):
        if lines:
            self.lines += [label + ":"] + [f"  - {l}" for l in lines] + [""]

    def sheet(self, name, df, description=""):
        f = self.dir / f"{name}.csv"
        df.to_csv(f, index=False, na_rep="", encoding="utf-8")
        self._sheets.append((f.name, len(df), len(df.columns), description))

    def close(self):
        self._block("Files", [f"{n}: {r:,} rows x {c} columns" + (f". {d}" if d else "")
                              for n, r, c, d in self._sheets])
        self._block("Provenance", [
            f"Generated {datetime.date.today().isoformat()} by tables/export_supplementary_tables.py "
            f"(song-system-grn, commit {git_commit()}).",
            "Missing values are empty fields."])
        (self.dir / "README.txt").write_text("\n".join(self.lines), encoding="utf-8")
        size = sum(p.stat().st_size for p in self.dir.iterdir()) / 1e6
        print(f"wrote {self.dir}/ ({size:.1f} MB): "
              + ", ".join(f"{n}={r}" for n, r, _, _ in self._sheets))


BOOK_CLASS = CsvBook


def Book(*args, **kwargs):
    return BOOK_CLASS(*args, **kwargs)


# --------------------------------------------------------------------------- tables

CALLS_NOTE = ("composite_calls: one row per finch cluster. top_rank_agg is the composite best match "
              "(rank aggregation across methods); confidence/confidence_tier score that call; "
              "conf_<method>/top_<method>/nrecip_<method> give each method's own confidence, best "
              "match and reciprocal support.")


def integration_book(out, key, species, ref, figure):
    d = COMPOSITE / f"{key}_hybrid"
    (calls,) = need(d / "composite_calls.csv")
    b = Book(out, f"Finch-{species} integration labels",
             f"Cross-species correspondence of adult zebra finch snRNA-seq cell types ({figure}) with "
             f"{ref}. Calls combine four methods (GSI, SAMap, Seurat CCA, SATURN at weight 0.5) into a "
             f"composite score per finch cluster x {species} cluster.",
             [calls, LABEL_LOOKUP],
             [CALLS_NOTE,
              "Finch GABAergic clusters were renamed to the hybrid labels with hybrid_label_lookup.csv, "
              "as the figures do (relabel_finch_clusters()); the source CSVs keep the raw cluster IDs."])
    df = pd.read_csv(calls)
    df["query"] = relabel(df["query"])
    b.sheet("composite_calls", df.rename(columns={"query": "finch_cluster"}))
    b.close()


def glut_integration_book(out):
    (calls,) = need(GLUT_COMPOSITE / "composite_calls_summary.csv")
    b = Book(out, "Glutamatergic finch-chicken label transfers",
             "Correspondence of finch glutamatergic neuron types with chicken (Zaremba et al. 2025) "
             "excitatory types (fig2-s3). Composite of GSI, SAMap (SCT) and Seurat RPCA label transfer "
             "on the plusSATB2_noPre glutamatergic subset (Glut-DACH2-HVCra-Int excluded).",
             [calls],
             [CALLS_NOTE + " group marks song vs non-song finch types."])
    b.sheet("composite_calls", pd.read_csv(calls).rename(columns={"query": "finch_cluster"}))
    b.close()


def deg_book(out):
    (src,) = need(DEG_DIR / "song_surround_paired_deg.csv")
    df = pd.read_csv(src)
    assert set(df["compare"]) == set(SONG_PAIRS), set(df["compare"])
    b = Book(out, "Glutamatergic song vs non-song differentially expressed genes",
             "Differential expression between each song glutamatergic neuron type and its non-song "
             "counterpart (fig2-s1). DESeq2 on raw-count pseudobulks (cell type x bird), design "
             "~ individual + group, paired within birds (>= 20 cells per side); log2 fold changes "
             "are apeglm-shrunken, positive = higher in the song type.",
             [src, DEG_DIR.parent / "song-surround_deg_glut_hybrid.qmd"],
             ["Significance threshold used in the paper: p_val_adj < 0.1 (heatmaps additionally "
              "require |avg_log2FC| >= 1).",
              "All tested genes are included, in a single table. contrast names the song type "
              "and the non-song type it was tested against (song vs non-song); rows are grouped "
              "by contrast and sorted by p_val_adj within each. The five contrasts are:"]
             + [f"{s} vs {n}" for s, n in SONG_PAIRS.items()])
    order = {s: i for i, s in enumerate(SONG_PAIRS)}
    df = (df.assign(contrast=df["compare"].map(lambda s: f"{s} vs {SONG_PAIRS[s]}"),
                    _o=df["compare"].map(order))
          .sort_values(["_o", "p_val_adj"], na_position="last", kind="stable"))
    df = df[["contrast", "gene", "avg_log2FC", "lfcSE", "baseMean", "p_val", "p_val_adj"]]
    b.sheet("glut_deg", df)
    b.close()


def core_book(out):
    core4, core5 = need(DEG_DIR / "core_no_lmansh_gene_sets.csv", DEG_DIR / "core_gene_sets.csv")
    b = Book(out, "Song system core genes",
             "Genes concordantly differentially expressed across song glutamatergic neuron types "
             "(fig2-s1), derived from the paired song vs non-song DESeq2 fits.",
             [core4, core5],
             ["core_no_LMANsh (shown in the figure): tested in the four RA, HVCra, HVCx and LMANco "
              "contrasts, same log2FC sign in all four, p_val_adj < 0.1 in at least three. lmansh "
              "describes the LMANsh contrast for that gene; in_5contrast_core marks genes that are "
              "also core across all five contrasts.",
              "five_contrast_sets: core = concordant sign across all five contrasts and p_val_adj < 0.1 "
              "in at least three; HVC-/LMAN-common = significant in both HVC (or both LMAN) pairs; "
              "-exclusive = additionally not significant outside that set. One row per gene per set.",
              "log2FC_<type>/padj_<type>: values from that type's song vs non-song contrast."])
    b.sheet("core_no_LMANsh", pd.read_csv(core4).sort_values(["direction", "gene"]))
    b.sheet("five_contrast_sets", pd.read_csv(core5).sort_values(["set", "direction", "gene"]))
    b.close()


def dar_book(out):
    pairs = [("Glut-CACNA1H-RA", "Glut-CACNA1H-1"), ("Glut-DACH2-HVCra", "Glut-DACH2-1"),
             ("Glut-DACH2-HVCx", "Glut-DACH2-4")]
    files = []
    for s, n in pairs:
        for d in ("1", "-1"):
            files.append((s, n, d, DAR_BEDS / f"{s}_vs_{n}_{d}.bed"))
    inter = DAR_BEDS / "Glut-CACNA1H-RA_Glut-DACH2-HVCra_Glut-DACH2-HVCx_intersect.bed"
    need(inter, *[f for *_, f in files])
    b = Book(out, "Song vs non-song differentially accessible regions",
             "Differentially accessible regions (DARs) between song glutamatergic neuron types and "
             "their non-song counterparts in the multiome data (fig3-s1). ArchR getMarkerFeatures on "
             "the PeakMatrix, Wilcoxon test with TSS enrichment and log10(nFrags) bias matching; "
             "DARs are FDR < 0.01 and |Log2FC| >= 0.6. Coordinates are 0-based, half-open (BED) on "
             "the zebra finch genome used for the multiome analysis.",
             [*[f for *_, f in files], inter, DATA_ROOT / "multiome/archr/peaks_glut.qmd"],
             ["One table (CSV file or workbook sheet) per pair and direction: _up = more accessible in the song type, _down = less.",
              "Log2FC is song type / non-song type.",
              "intersect: regions shared across the three song types' DAR sets (both directions).",
              "The repository's *_union.bed is a plain concatenation of the six pair files and is "
              "not repeated here.",
              "There is no LMAN pair in the multiome data."])
    for s, n, d, f in files:
        df = pd.read_csv(f, sep="\t", header=None,
                         names=["chrom", "start", "end", "name", "Log2FC", "strand"])
        df.insert(4, "song_type", s)
        df.insert(5, "nonsong_type", n)
        b.sheet(f"{SHORT[s]}_{'up' if d == '1' else 'down'}", df)
    df = pd.read_csv(inter, sep="\t", header=None, names=["chrom", "start", "end", "name"])
    b.sheet("intersect", df)
    b.close()


EREG_RUNS = {
    "full": ("ra-arco-hvc-nc_hybrid", "config37",
             "all cell types of RA, Arco, HVC and NC (hybrid labels)"),
    "ra_c1h1": ("glut-ra_glut-c1h-1_allruns_hybrid", "config2",
                "Glut-CACNA1H-RA and Glut-CACNA1H-1 neurons, all sequencing runs"),
    "gaba": ("ra-arco-hvc-nc_gaba_hybrid", "config2",
             "GABAergic neurons of RA, Arco, HVC and NC (8 clusters, GABA-Im excluded)"),
}


def ereg_book(out, which, title, figure):
    run, cfg, cells = EREG_RUNS[which]
    o = SCPLUS_ROOT / run / "results" / cfg / "outs"
    direct, ext = need(o / "eRegulon_direct.tsv", o / "eRegulons_extended.tsv")
    grn_dir = "ra-arco-hvc-nc" if which == "full" else run
    yaml = f"grn/{grn_dir}/scenicplus/{cfg}/Snakemake/config/config.yaml (song-system-grn repository)"
    frames = []
    for f in (direct, ext):
        df = pd.read_csv(f, sep="\t")
        sign = df["eRegulon_name"].str.rsplit("_", n=1).str[-1]
        frames.append(df[sign.isin(EREG_SIGNS)])
    df = pd.concat(frames, ignore_index=True)
    df = df.sort_values(["TF", "is_extended", "eRegulon_name", "triplet_rank"])
    n_ereg = df["eRegulon_name"].nunique()
    b = Book(out, title,
             f"SCENIC+ enhancer-driven regulons (eRegulons) for {cells} ({figure}). SCENIC+ run "
             f"{run}, {cfg}. One row per TF - region - target gene triplet.",
             [direct, ext, f"Run configuration: {yaml}"],
             [f"eRegulon signs are TF-to-gene / region-to-gene correlation. Restricted to +/+ "
              f"(activating TF, accessibility rises with expression) and -/- (repressing TF, "
              f"accessibility falls with expression) eRegulons ({n_ereg} eRegulons); +/- and -/+ "
              f"eRegulons are omitted. Direct and extended "
              f"eRegulons are pooled; is_extended = TRUE marks those built from orthology-extended "
              f"motif annotations.",
              "importance_R2G/rho_R2G: region-to-gene importance and correlation; importance_TF2G/"
              "rho_TF2G: TF-to-gene importance and correlation; regulation: 1 = activating, "
              "-1 = repressing; triplet_rank: SCENIC+ ranking of the TF-region-gene triplet.",
              "Gene_signature_name/Region_signature_name give the eRegulon with its gene (g) and "
              "region (r) counts."])
    b.sheet("eRegulons", df)
    b.close()


TABLES = {
    "01": ("01_finch-mouse_integration_labels.xlsx",
           lambda p: integration_book(p, "yao_adult", "mouse",
                                      "the adult mouse whole-brain atlas (Yao et al. 2023)",
                                      "fig1-s3")),
    "02": ("02_finch-chicken_integration_labels.xlsx",
           lambda p: integration_book(p, "gg_adult", "chicken",
                                      "the adult chicken pallium (Zaremba et al. 2025)", "fig1-s2")),
    "03": ("03_glut_differentially_expressed_genes.xlsx", deg_book),
    "04": ("04_song_system_core_genes.xlsx", core_book),
    "05": ("05_glut_finch-chicken_label_transfers.xlsx", glut_integration_book),
    "06": ("06_song-nonsong_DARs.xlsx", dar_book),
    "07": ("07_eRegulons_full_dataset.xlsx",
           lambda p: ereg_book(p, "full", "eRegulons for the full dataset", "fig4-s1")),
    "08": ("08_eRegulons_RA_C1H-1.xlsx",
           lambda p: ereg_book(p, "ra_c1h1", "eRegulons for Glut-CACNA1H-RA / Glut-CACNA1H-1",
                               "fig6")),
    "09": ("09_eRegulons_GABA.xlsx",
           lambda p: ereg_book(p, "gaba", "eRegulons for GABAergic neurons", "fig6")),
}


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--only", default="", help="comma-separated table numbers, e.g. 01,07")
    ap.add_argument("--format", choices=["csv", "xlsx", "both"], default="csv",
                    help="csv: <out>/csv/<table>/<sheet>.csv + README.txt (default); "
                         "xlsx: <out>/<table>.xlsx")
    a = ap.parse_args()
    keys = a.only.split(",") if a.only else list(TABLES)
    global BOOK_CLASS
    for fmt in (["csv", "xlsx"] if a.format == "both" else [a.format]):
        BOOK_CLASS = CsvBook if fmt == "csv" else XlsxBook
        out = a.out / "csv" if fmt == "csv" else a.out
        out.mkdir(parents=True, exist_ok=True)
        for k in keys:
            fname, fn = TABLES[k]
            fn(out / fname)


if __name__ == "__main__":
    main()
