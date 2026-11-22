# Figure 6 — standalone panel scripts

One script per scripted panel of Figure 6 (`fig6_mafb_locus`, "Song neuron-specific regulatory landscape of *MAFB*"). Each regenerates its panel from the
analysis tables, without running any other panel. The analysis code is in `comparative_genomics/phyloP_matched/loci_focus/` (and, for panel B,
`grn/analysis/`); the panel scripts run it and copy what it draws to `figures/figure6/output/` (not tracked; set `FIG6_OUT` to write elsewhere).
Paths and interpreters are in `common.sh`, which every panel script sources; each can be overridden with the environment variable named there.

| Panel | Content | Script | Output (`output/`) |
|---|---|---|---|
| A | Model of *cis*-regulatory co-option | — (drawn by hand) | — |
| B | SCENIC+ inputs to MAFB in C1H-RA (vs C1H-1) and GABAergic neurons | `panel_B_mafb_inputs_network.sh` | `panel_B_mafb_inputs_network.{svg,pdf,png}` |
| C | MAFB domain: peak-to-gene links, accessibility (RA, C1H-1, PVALB-1, LGE-1), genes, zoom-ins 1–4 | `panel_C_mafb_domain_overview.sh` | `panel_C_mafb_domain_overview.{pdf,svg}` |
| D | Avian and oscine phyloP of the DARs, by phyloP class | `panel_D_element_constraint_heatmap.sh` | `panel_D_element_constraint_heatmap.{pdf,png}` |
| E | phyloP class composition of RA, PVALB-1 and shared DARs | `panel_E_constraint_class_bars.sh` | `panel_E_constraint_class_bars.{pdf,png}` |
| F | Zoom-ins 1–4: accessibility, regulator motif hits, phyloP | `panel_F_element_zooms.sh` | `panel_F_element_zooms.{pdf,svg,png}` |
| G | ARE-full site of element 3 across the 363-bird alignment | `panel_G_are_site_alignment.sh` | `panel_G_are_site_alignment.{pdf,svg,png}` |

## Running

```bash
bash figures/figure6/make_inputs.sh          # once: the analysis tables (see below); ~5 min, plus ~1 h for the alignment extract if absent
bash figures/figure6/make_all.sh             # panels C-G (~1 min); or name panels: make_all.sh D E
bash figures/figure6/panel_B_mafb_inputs_network.sh   # on lark (SCENIC+ tables under /hdd/jupyter)
```

Python 3 with numpy, pandas, scipy, pysam, pyBigWig and matplotlib (`FIG6_PY`); for panel G and the browser clean-up also lxml (`FIG6_PY_ALN`);
R with data.table, ggplot2, ComplexHeatmap, circlize and cowplot (`FIG6_RSCRIPT`). Panels C and F read the accessibility and phyloP tracks from the
lab track hub (`colquitt-track-hub.pbsci.ucsc.edu`) and need network access; panel C is rendered by the UCSC Genome Browser and needs an API key in
`~/.ucsc_api_key` and Inkscape on the PATH of a login shell. `make_inputs.sh` needs `bedtools` and, to build the alignment extract, Cactus' `hal2maf`
and the 363-avian HAL file.

## Inputs (`make_inputs.sh`, in order)

| Step | Script (`loci_focus/`) | Writes | Used by |
|---|---|---|---|
| Elements: SnapATAC2 peaks of C1H-RA and PVALB-1 overlapping their ArchR DARs (vs C1H-1, LGE-1; FDR < 0.1, log2FC > 0) in the MAFB domain, curated, split into lobes; accessibility-weighted phyloP; alignment depth | `element_profiles.py` | `results/mafb_domain/element_summaries.tsv`, `element_profiles_*.tsv.gz` | D, E (via 2–3) |
| Per-element summaries for the class step | `element_rule_classes.py` | `element_rule_classes.tsv` | D |
| phyloP classes (avian 0.52, oscine 0.27, strong avian 3.0) | `element_threshold_classes.py` | `element_threshold_classes.tsv` | D, E, F |
| AR half-site (FIMO) | `element_ar_sites.py` | `element_ar_sites.tsv` | F |
| avian363 alignment with ancestors over every element ± 100 bp | `hal2maf` | `maf_switch/data/lobes_gap_pvalb1.maf` | F (via 6), G |
| Regulator motif sites, presence across species, oscine gain and placement | `maf_switch/lobe_site_gain.py` | `element_motif_gain.tsv` | F |
| (optional, `FIG6_NULLS=1`) shuffled and background-motif expectations | `maf_switch/lobe_gain_null.py`, `_report.py` | `lobe_gain_null_*` | — (text) |

Upstream of these: the ArchR DAR tables (`tracks_lark/panelE_inputs/dar_{ra,pvalb1_lge1}.tsv`; `multiome/archr/dar_pvalb1_lge1.R` for PVALB-1),
the SnapATAC2 MACS3 peaks per cluster (`tracks_lark/snapatac_peaks/`; `multiome/snapatac/call_peaks.ipynb`), the fragment bigWigs
(`tracks_lark/fragment/`), the per-base phyloP tracks (`phyloP_matched/perbase/`: `phyloP --mode CONACC --method LRT`, all branches and
`--subtree birdAnc222`), the per-base alignment depth tracks (`domain_depth.py`), and the 1-Mb peak-to-gene links
(`multiome/archr/peaks2geneslinks_mafb_domain.R`, on the track hub).

Curation is in files next to the analysis code, read by `element_profiles.py`: `excluded_peaks.bed`, `manual_lobes.bed`, `lobe_count.bed`,
`excluded_lobes.bed`, `class_override.bed`, and the PVALB-1 review (`class_override_pvalb1.bed`, `excluded_peaks_pvalb1.bed`, which take precedence).
`GABA_COMP=PVALB2` rebuilds the earlier PVALB-2 element set exactly (kept in `results/mafb_domain_pvalb2/`).

## Notes on individual panels

- **B.** KLF15, weakly expressed in C1H-RA, is left out (`MAFB_NET_EXCLUDE=RA_all:KLF15`). AR -> MAFB is a hand-added dashed edge (AR is not an
  eRegulon input to MAFB in the all-runs tables); see the docstring of `mafb_inputs_network_circle.py`.
- **C.** `ucsc_dar_batch.py overview zooms`: the browser render (`*_browser.pdf`) is converted to SVG with editable text, track bars merged
  (`scripts/browser_pdf_to_svg.py`), labels right-aligned, the ruler replaced by a 100-kb scale bar, and the zoom windows of panel F
  (`zoom_panels.windows()`) drawn as numbered grey bars.
- **D, E.** Elements aligned in ≥ 80 % of oscine and non-oscine genomes over every 50-bp window (66 elements, 61 peaks). E: two-sided Fisher's exact
  tests on the full class composition, Holm-adjusted over the three pairs (RA vs PVALB-1 p = 0.00049, PVALB-1 vs shared 0.038, RA vs shared 0.21).
- **F.** `zoom_panels.py`: zoom-ins and the boxed site are set in `ZOOMS` and `BOX` at the top of the script; flanks 400 bp (`FLANK`).
- **G.** `plot_site_alignment.py --paper 143 --row-mm 2.18 --flank 8 --compact --show Ptilonorhynchus_violaceus=Climacteris_rufus`.

## Check against the figure files

Run 2026-11-22; outputs compared with the files placed in the figure: F and G pixel-identical (PNG, 600 dpi); C, D, E are the same files the
analysis scripts write for the figure.
