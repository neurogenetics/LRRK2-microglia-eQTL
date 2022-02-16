#!/bin/env Rscript

# Load necessary packages
library(Signac)  #v1.5.0
packageVersion('Signac')

library(Seurat)  #v4.0.5
packageVersion('Seurat')

library(tidyverse)
library(svglite)

# Load named integrated object
load("output/SN2960_SN2972_SN2942_Signac_int_predicted_ids.Rdata")
SN.int  #8001 nuclei

# Find celltype markers in GEX dataset (using SCT normalized values)
SN.int.markers <- FindAllMarkers(object = SN.int, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(SN.int.markers, file = "output/SN3_multiGEX_celltype_markers.csv")


# Find motifs enriched in SN microglia
library(JASPAR2020)  #BiocManager::install("JASPAR2020")
library(TFBSTools)  
library(BSgenome.Hsapiens.UCSC.hg38)  #BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")

# Switch default assay to peaks
DefaultAssay(SN.int) <- "peaks"

# Exclude scaffolds
main.chroms <- standardChromosomes(BSgenome.Hsapiens.UCSC.hg38)
keep.peaks <- which(as.character(seqnames(granges(SN.int))) %in% main.chroms)
SN.int[["peaks"]] <- subset(SN.int[["peaks"]], features = rownames(SN.int[["peaks"]])[keep.peaks])

# Get a list of motif positions from the JASPAR database
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)

# Add motif information to the FC.atac object
SN.int <- AddMotifs(
  object = SN.int,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  pfm = pfm
)

# Find motifs enriched in MGL vs all other cell types
da_peaks_MGL <- FindMarkers(
  object = SN.int,
  ident.1 = 'MGL.6',
  only.pos = TRUE,
  test.use = 'LR',
  min.pct = 0.1,
  latent.vars = 'nCount_peaks'
)

top_da_peaks_MGL <- rownames(da_peaks_MGL[da_peaks_MGL$p_val < 0.005, ])

enriched_motifs_MGL <- FindMotifs(
  object = SN.int,
  features = top_da_peaks_MGL
)

enriched_motifs_MGL <- enriched_motifs_MGL[order(enriched_motifs_MGL[, 7], -enriched_motifs_MGL[, 6]), ]
write.csv(file = "output/SN3_motifs_enriched_MGL6_vs_allothers.csv", enriched_motifs_MGL, row.names = F, quote = F)
head(enriched_motifs_MGL, n = 30)

p1 <- MotifPlot(
  object = SN.int,
  motifs = head(rownames(enriched_motifs_MGL), n = 12)
)
ggsave(file = "plots/SN3_multiATAC_MGL6_top_motifs.png", plot = p1, width = 10, height = 8)
ggsave(file = "plots/SN3_multiATAC_MGL6_top_motifs.svg", plot = p1, width = 10, height = 8)
