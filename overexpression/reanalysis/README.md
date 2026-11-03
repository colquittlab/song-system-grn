# MAFB / EMX2 overexpression in chicken embryo pallium — reanalysis of run 1

Reanalysis of the PIPseq run `260312_chicken_mafb_emx2`, originally analyzed in `overexpression/chick_mafb-emx2_run1.qmd`
(left untouched). The question is the same: does dox-induced MAFB (or EMX2) push chicken glutamatergic neurons toward a
song-neuron expression profile? What changes is how it is tested:

| | original (`chick_mafb-emx2_run1.qmd`) | this reanalysis |
|---|---|---|
| unit of replication | cell (Wilcoxon `FindMarkers`) | embryo (souporcell, 2–4 per library) |
| tg+ | MAFB-tg ≥ 1, **pooled over the −dox and +dox libraries** | within one library; any plasmid read, or any inducible read (TF-tg or mScarlet) |
| control | tg− cells | tg− cells of the **same embryo**, and the same contrast in the −dox library (difference-in-differences) |
| DE | per-cell Wilcoxon | pseudobulk DESeq2 (project convention, `figures/figure4/common.R`) |
| song comparison | RRHO / KS vs SCENIC+ eRegulon genes and old song-vs-surround DE | fgsea + UCell/AUCell scoring of gated adult signatures against an expression-matched null |
| song axes | RA, HVC (old labels) | RA vs CACNA1H-1, HVCra vs DACH2-1, HVCx vs DACH2-4 (hybrid labels) |

## Scripts (run in order)

```bash
Rscript overexpression/reanalysis/01_reference_signatures.R   # adult song-axis signatures (finch snRNA, hybrid labels)
Rscript overexpression/reanalysis/02_process_oe.R             # chicken cells, embryos, classes, plasmid status
Rscript overexpression/reanalysis/03_pseudobulk_de.R          # detection model, pseudobulk DE, fgsea
Rscript overexpression/reanalysis/04_signature_scoring.R      # UCell/AUCell vs expression-matched null
OE_TG_DEF=inducible Rscript overexpression/reanalysis/03_pseudobulk_de.R      # tg+ = MAFB-tg/EMX2-tg/mScarlet read
OE_TG_DEF=inducible Rscript overexpression/reanalysis/04_signature_scoring.R
```

Shared helpers are in `oe_helpers.R`. The scoring helpers (`fast_auroc`, `bin_genes_by_expression`,
`expression_matched_sets` and the seeded AUCell rankings) are copied unchanged from the songbird_dev repo
(`lonStrDom3/seurat/dev-full1_cellbender.0.1_multiome_seurat_preprocess/nido/`, PR #5). That is the method from the session
"glutamatergic neuron clustering seurat". UCell and AUCell are not in this project's renv, so `oe_helpers.R` appends
dev_full1's R 4.6 library (`DEV_FULL1_RLIB`). Small tables go to `tables/` (tracked). Objects, full DE tables and figures go
to `output/` (gitignored).

### Design choices

- **Signatures** (`01`): pseudobulk DESeq2 on the adult reference
  (`snrna/naming/hybrid_division_naming/obj_hybrid_labels.qs2`). Each replicate is a library × bird × cell type with
  ≥ 10 nuclei, and the design is `~ lib + grp`. Genes need padj < 0.05 and |log2FC| ≥ 0.25. They then pass the
  songbird_dev specificity gate: a song-up gene must separate the song cluster from all 19 adult glutamatergic clusters
  (AUROC ≥ 0.70), and a sister-up gene must sit low in it (≤ 0.30). Genes are mapped to chicken by OMA 1:1 orthologs, the
  `snrna/integration/` rule. MAFB and EMX2 are excluded from every set. Gated, orthologous set sizes are below; the
  gate keeps 5–37% of each DE list.

  | axis | song-up | sister-up |
  |---|---|---|
  | RA vs CACNA1H-1 | 161 | 74 |
  | HVCra vs DACH2-1 | 71 | 76 |
  | HVCx vs DACH2-4 | 59 | 34 |

- **Scoring** (`04`): UCell up-only (maxRank 1500), with AUCell as a second scorer (seeded rankings, aucMaxRank = 5%).
  - Each cell's score is depth-adjusted, using the residual of `score ~ library + log10 UMI`. This matters because plasmid
    detection selects deeper cells.
  - Per embryo, the shift is the mean score of tg+ cells minus the mean of tg− cells.
  - Each shift is compared with the same shift for 1,000 expression-matched random gene sets. The adult-anchor calibration
    is **not** used: adult finch 10x anchors do not transfer to embryonic chicken PIPseq nuclei, so the claim is relative
    (tg+ vs tg− in the same embryo).
- **Minimum cells**: only one embryo has ≥ 10 plasmid+ glutamatergic cells. tg+ pseudobulk samples therefore use ≥ 3 cells,
  a documented exception to the project's ≥ 10.

## Results

### 1. The data available to test the question are very thin

Plasmid+ cells after QC (souporcell singlets, mito < 2%, 500–20,000 UMI, ≥ 300 genes; 3,387 cells, 861 glutamatergic):

| library | dox | embryos (Glut cells) | Glut plasmid+ | Glut inducible+ | Glut TF-tg+ |
|---|---|---|---|---|---|
| MAFBn | − | 118, 105 | 15 (8, 7) | 13 | 9 |
| MAFBp | + | 179, 80, 5 | 25 (18, 7, 0) | 25 | 22 |
| EMX2p | + | 116, 100, 86, 72 | 20 (8, 6, 4, 2) | 19 | 14 |

The three libraries share no embryos. Depth-filtered souporcell genotypes are equally discordant (0.32–0.45) for every
cross-library cluster pair.

### 2. There is no detectable dox induction of the inducible transcripts (`tables/oe_plasmid_by_library.csv`)

| | MAFB −dox | MAFB +dox |
|---|---|---|
| MAFB-tg UMI per 10k UMI | 0.151 | 0.117 |
| mScarlet UMI per 10k UMI | 0.042 | 0.021 |
| % cells MAFB-tg+ | 3.4 | 4.5 |
| mNG (constitutive) % cells | 1.4 | 0.3 |

The −dox library carries as much inducible transcript as +dox. MAFB-tg has no BLAST hit to endogenous chicken MAFB, so
this is not mis-mapping. It is leaky expression and/or ambient RNA. Whether +dox raised MAFB *protein* cannot be seen here.
Until induction is confirmed independently, a "+dox vs −dox" contrast in this run is not an OE-vs-control contrast.

### 3. The reviewer's point: tg− cells that are really positive (`tables/oe_plasmid_detection.csv`)

A depth-offset zero-truncated Poisson was fitted to the plasmid UMIs in positive cells (1.3–1.8 UMI per positive cell). A
plasmid-carrying glutamatergic cell at typical depth shows zero plasmid reads 35–57% of the time. The implied
contamination of the tg− group is **4–15%** (MAFB −dox 4%, MAFB +dox 15%, EMX2 7%). Hidden positives of this size
*attenuate* a tg+ vs tg− difference; they cannot create one.

The larger error runs the other way. Souporcell estimates 20% ambient RNA in MAFBn, 8% in MAFBp and 27% in EMX2p. At
those levels, an estimated **36% (MAFBn) and 44% (EMX2p)** of plasmid+ calls could come from ambient RNA alone, versus 10%
in MAFBp. Single-UMI tg+ calls in those libraries are unreliable.

### 4. Pseudobulk DE: essentially nothing within embryos (`tables/oe_de_summary*.csv`, `oe_de_top25*.csv`)

| contrast (glutamatergic, embryo-paired) | padj < 0.05 (plasmid / inducible tg+) |
|---|---|
| MAFB +dox, tg+ vs tg− | 1 / 1 |
| MAFB −dox, tg+ vs tg− | 0 / 0 |
| EMX2 +dox, tg+ vs tg− | 2 / 2 |
| MAFB difference-in-differences (+dox vs −dox tg effect) | 0 / 0 |
| all Glut, +dox vs −dox embryos (2 vs 2) | 14 / 16 |

The +dox vs −dox library contrast picks up OTP, RORB, BNC2, CDH13 and ABI3BP. These look like a difference in which
tissue each dissection sampled, not a MAFB response.

### 5. Song axes: no song-like shift (`tables/oe_fgsea_song_axes*.csv`, `oe_signature_scores*.csv`)

- **fgsea** on the within-embryo rankings finds song-up *and* sister-up sets lower in tg+ cells, in every library
  including MAFB −dox (NES −1.2 to −2.0). That is a general difference between electroporated and neighboring cells, not
  a song axis.
- **Matched-null scoring** removes it: once depth is adjusted and random expression-matched genes are subtracted, no
  song-up set shifts in any condition.

  | UCell z (plasmid tg+) | RA song-up | HVCra song-up | HVCx song-up |
  |---|---|---|---|
  | MAFB +dox | +0.17 | −0.27 | +0.68 |
  | MAFB −dox | −0.64 | −1.20 | −0.84 |
  | EMX2 +dox | −1.15 | −0.76 | −0.97 |
  | MAFB DiD | +0.56 | +0.69 | +1.02 |

  AUCell agrees, and the inducible tg+ definition agrees (RA song-up DiD z = +0.28).
- **The only nominal signals point away from song identity.** In MAFB +dox, tg+ cells are higher on *sister-up* genes:
  HVCx z = 2.31 (p = 0.024), with DiD z = 2.1–2.4 for HVCra and HVCx sisters. None survives correction across the 24
  set × statistic tests, and the signs are inconsistent between embryos (`output/signature_shift_per_embryo*.png`).
- **A TF-transgene-only tg+** (MAFB-tg / EMX2-tg reads only) gave the same picture in step 3 before it was replaced by the
  inducible definition: 0–1 DE genes within embryos, and an RA song-up DiD NES of −0.97.

**Bottom line.** With embryo as the replicate, this run does not support MAFB- or EMX2-induced RA- or HVC-like
expression. Read plainly, the earlier RA-like result came from per-cell testing, and its "MAFB+" group pooled −dox
(leaky/ambient) cells with +dox cells. The design can only detect large effects: 2 embryos per arm, 7–18 tg+ cells
each, ~1 UMI per transgene call and no measurable induction. So this is "not detectable here", not evidence of no effect.

## For the repeat experiment

1. **Confirm induction before sequencing.** Check mScarlet fluorescence ± dox in the same batch, or run qPCR on dissected
   tissue. In run 1, −dox and +dox carried the same amount of inducible transcript.
2. **Make the inducible transcript detectable.** PIPseq is 3′-biased, and the TF CDS sits far from the poly(A). A short
   unique barcode in the 3′ UTR just upstream of the poly(A) site of both the inducible and the constitutive transcript
   would raise detection well above 1–2 UMI per cell. That reduces both hidden positives and ambient false positives.
3. **Enrich electroporated cells.** Sort mNG+ nuclei, or at least dissect only the electroporated region. Run 1 had 5–9%
   plasmid+ glutamatergic cells per library.
4. **Replicates where the test needs them.** Plan ≥ 3–4 embryos per arm with ≥ 10–20 tg+ glutamatergic cells each, with
   the −dox arm processed in parallel. Souporcell demultiplexing worked well, so pooling embryos per library is fine.
   Ideally, split one embryo cohort across ± dox libraries so that the DiD is not confounded with dissection.
5. **Remove ambient RNA** (CellBender, as for the finch data) before calling tg+. Single-UMI calls in libraries with
   20–27% ambient are unreliable.

The scripts take the library table (`LIBS` in `02_process_oe.R`) and souporcell k as inputs. Adding the new libraries
there should rerun the whole analysis on the repeat.
