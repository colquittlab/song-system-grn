
# Neurogenesis -------------------------------------------------------------------------
## https://cshperspectives.cshlp.org/content/8/5/a018820.full

genes_b1_cells = c("SOX2", "PRRX1", "NR2E1", "ID1", "SRRT", "BMI1", "EZH2", "HDAC1")
genes_c_cells = c("NEUROG2","OLIG2", "ASCL1", "TBR1", "TBR2", "NEUROD1")

## Summarized https://dev.biologists.org/content/146/4/dev156059.long
genes_qNSC = c("BMP1", "BMP6", "SOX9", "SLC1A3", "NOTCH2","FABP7", "HES5", "NR2E1")
genes_aNSC = c("ASCL1", "DLL1", "DLL3", "EGFR1", "FOS", "CCND1", "E2F1", "FOXM1","MYC", "MYCN")


# Neuropeptides
nps = c("ADCYAP1", "CCK", "CRH", "CRHBP", "NPY", "PDYN", "PENK", "PNOC", "SST", "TAC1", "TAC3", "VIP" )
nprs = c("ADCYAP1R1", "CCKBR", "CRHR1", "CRHR2", "NPY1R", "NPY2R", "NPY5R", "OPKR1", "OPRD1", "OPRM1", "OPRL1", 
         "SSTR1", "SSTR2", "SSTR3", "SSTR4","SSTR5",
         "TACR1", "TACR2", "TACR3", "VIPR1", "VIPR2"
         )

# Cell type markers
ct_markers = list(
  Glut=c("SLC17A6", "SYT1"), 
  GABA=c("GAD1", "GAD2", "DLX1"), 
  Astro=c("SLC1A3", "GFAP"), 
  Oligo=c("MBP", "PLP1"), 
  OPC=c("PDGFRA"), 
  Micro=c("CSF1R"), 
  Endo=c("LUM", "RGS5", "FLI1"), 
  Epen=c("SPEF2", "FOXJ1"),
  Mito=c("TOP2A", "NUF2", "MCM5", "PCNA", "TYMS"),
  NP=c("SOX2", "PAX6", "NES", "CCND2", "CDCA7"), 
  GABA_NP=c("VAX1", "NKX2-1", "OLIG2", "ISL1"), 
  IP=c("EOMES")
)

## Arcopallium 
arco_marker_genes = list(
  ARCO=c("C1QL3", "NR4A2"),
  RA=c("SIX2", "NDNF", "PVALB", "SYF2", "AR", "GABRE", "KCNS1"),
  AMD=c("ADRA1D", "CYP19A1"),
  AV=c("NECAB2", "AQP1"),
  AP=c("CDH4", "ZEB2", "SV2B", "PLD1"),
  AId=c("CABP1", "ADAM23", "CAMTA1", "SYNGR3", "HTR2A", "RCAN2"),
  AD=c("CBLN2",   "CCK", "CRHR1", "CBLN2", "ATP2A3", "QSOX1"),
  AMV=c("ESR2",  "VIP", "FAM163B", "SLC6A6", "ID2",  "ZBTB20", "DAPK1", "FABP7", "KCND2", "MAOA", "MGAT4C"),
  AA=c("KCTD20", "ID2", "CRHR2", "GLRA4", "KCTD12", "SCUBE1")
)
marker_genes_AV = c("NECAB2", "AQP1", "SV2B")
marker_genes_AP = c("NECAB2", "C1QL3",
                    "SV2B" # moderate
)
marker_genes_AId = c("HTR2A", "KCNS1", "PVALB", "CNTN4")
marker_genes_AMV = c("ESR2", "ZBTB20", "KCNQ5")
marker_genes_AMD = c("CYP19A1",
                     "ZBTB20", #negative,
                     "AQP1", "NECAB2")