#!/bin/env Rscript

# Motif analysis with Signac https://satijalab.org/signac/articles/motif_vignette.html
# Load necessary packages
library(Signac)   #v1.5.0
packageVersion('Signac')

library(Seurat)   #v4.0.5
packageVersion('Seurat')

library(JASPAR2020)  #BiocManager::install("JASPAR2020")
library(TFBSTools)  
library(BSgenome.Hsapiens.UCSC.hg38)  #BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
library(ggplot2)
library(svglite)

# Load named, harmony integrated object
load("output/FC_ATAC_harmony_int_named.Rdata")
FC.atac #38474 nuclei

# Get a list of motif positions from the JASPAR database
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)

# Add motif information to the FC.atac object
FC.atac <- AddMotifs(
  object = FC.atac,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  pfm = pfm
)

# Confirm that we see enrichment of MEF2C in PVALB+ InN (positive control)
da_peaks_InN <- FindMarkers(
  object = FC.atac,
  ident.1 = 'InN.7',  #PVALB+
  ident.2 = 'InN.8',  #SST+
  only.pos = TRUE,
  test.use = 'LR',
  min.pct = 0.1,
  latent.vars = 'nCount_peaks'
)

## Sort differentially accessible peaks
top_da_peaks_InN <- rownames(da_peaks_InN[da_peaks_InN$p_val < 0.005, ])

## Find overrepresented motifs
enriched_motifs_InN <- FindMotifs(
  object = FC.atac,
  features = top_da_peaks_InN
)

## Sort motifs
enriched_motifs_InN <- enriched_motifs_InN[order(enriched_motifs_InN[, 7], -enriched_motifs_InN[, 6]), ]
head(enriched_motifs_InN, n = 15)  # We do see enrichment of MEF2C as well as MEF2A, MEF2D, MEF2B in the PVALB+ InN


# Find motifs enriched in MGL vs all other cell types
da_peaks_MGL <- FindMarkers(
  object = FC.atac,
  ident.1 = 'MGL.4',
  only.pos = TRUE,
  test.use = 'LR',
  min.pct = 0.1,
  latent.vars = 'nCount_peaks'
)

top_da_peaks_MGL <- rownames(da_peaks_MGL[da_peaks_MGL$p_val < 0.005, ])

enriched_motifs_MGL <- FindMotifs(
  object = FC.atac,
  features = top_da_peaks_MGL
)

enriched_motifs_MGL <- enriched_motifs_MGL[order(enriched_motifs_MGL[, 7], -enriched_motifs_MGL[, 6]), ]
write.csv(file = "output/Set12_motifs_enriched_MGL4_vs_allothers.csv", enriched_motifs_MGL, row.names = F, quote = F)
head(enriched_motifs_MGL, n = 30)

p1 <- MotifPlot(
  object = FC.atac,
  motifs = head(rownames(enriched_motifs_MGL), n = 12)
)
ggsave(file = "plots/Set12_FCatac_MGL4_top_motifs.png", plot = p1, width = 10, height = 8)

# Make plot showing TF most likely to bind overlying rs1491942 (in LRRK2 TSS peak), using JASPAR 2020 track in ucsc genome browser
## https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr12%3A40226978%2D40227033&hgsid=1234207631_rBXLDbPG8FlD9cyjbAksid2PHE7w
## These four TF also on list of enriched motifs in MGL4
motifs_rs1942 <- c("MA0154.4", "MA0116.1", "MA0774.1", "MA0775.1")
p2 <- MotifPlot(
  object = FC.atac,
  motifs = motifs_rs1942,
  ncol = 4
)

ggsave("plots/Set12_FCatac_enriched_MGL4_over_rs1491942.png", plot = p2, width = 10, height = 2.67)


# Make plot showing TF most likely to bind overlying rs6581593, using JASPAR 2022 track in ucsc genome browser
## This TF also on list of enriched motifs in MGL4, and differentially expression in FCtx snRNAseq MGL13
p3 <- MotifPlot(
  object = FC.atac,
  motifs = "MA0473.3"
)

ggsave("plots/Set12_FCatac_enriched_MGL4_over_rs6581593_ELF1.png", plot = p3, width = 5, height = 2.5)
ggsave("plots/Set12_FCatac_enriched_MGL4_over_rs6581593_ELF1.svg", plot = p3, width = 5, height = 2.5)






