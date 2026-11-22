suppressMessages(library(SummarizedExperiment))
G <- "/private/home/bcolquit/group/song-system-grn/multiome/archr"
O <- "/private/warm-archive/colquittlab/assembly/tracks_lark/panelE_inputs"
for (f in c(glut = "peaks_glut_hybrid/peaks_nonsong_glut_se.rds", gaba = "peaks_gaba_hybrid/peaks_nonsong_se.rds")) {
  x <- readRDS(file.path(G, f))
  cat(f, " class:", class(x), " contrasts:", paste(names(x), collapse = ", "), "\n")
  for (nm in names(x)) {
    se <- x[[nm]]; rd <- rowData(se)
    df <- data.frame(chrom = rd$seqnames, start = rd$start, end = rd$end,
                     Log2FC = assay(se, "Log2FC")[, 1], FDR = assay(se, "FDR")[, 1])
    write.table(df, file.path(O, paste0("dar_", nm, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
    cat(sprintf("  %-14s %d peaks; width %s; FDR<0.1&LFC>1: %d\n", nm, nrow(df),
                paste(range(df$end - df$start), collapse = "-"), sum(df$FDR < 0.1 & df$Log2FC > 1, na.rm = TRUE)))
  }
}
