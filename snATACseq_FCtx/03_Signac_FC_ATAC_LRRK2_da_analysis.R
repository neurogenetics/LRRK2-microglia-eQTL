#!/bin/env Rscript

# Load necessary packages
library(Signac)  #confirm v1.1.1
packageVersion('Signac')

library(Seurat)  #confirm v4.0.0
packageVersion('Seurat')

library(patchwork)
library(tidyverse)
library(svglite)

# Load harmony integrated object
load("output/FC_ATAC_clustered_harmony_int.Rdata")
FC.atac  #38474 nuclei

DefaultAssay(FC.atac) <- "peaks"

# Name cell types based on marker gene expression
FC.atac$num.ident <- Idents(FC.atac)

FC.atac <- RenameIdents(
  object = FC.atac,
  '0' = 'ODC.0',
  '1' = 'ODC.1',
  '2' = 'ExN.2',
  '3' = 'AST.3',
  '4' = 'MGL.4',
  '5' = 'OPC.5',
  '6' = 'ExN.6',
  '7' = 'InN.7',
  '8' = 'InN.8',
  '9' = 'ExN.9',
  '10' = 'ExN.10',
  '11' = 'InN.11',
  '12' = 'ExN.12',
  '13' = 'InN.13',
  '14' = 'AST.14',
  '15' = 'VC.15',
  '16' = 'ExN.16',
  '17' = 'ExN.17',
  '18' = 'ExN.18',
  '19' = 'ExN.19',
  '20' = 'InN.20',
  '21' = 'ExN.21'
)

FC.atac$celltype <- Idents(FC.atac)

p1 <- DimPlot(FC.atac, label = TRUE) + NoLegend()
ggsave(file = "plots/FC_ATAC_UMAP_harmony_int_named.png", plot = p1, width = 5.1, height = 5, dpi = 300)
ggsave(file = "plots/FC_ATAC_UMAP_harmony_int_named.svg", plot = p1, width = 5.1, height = 5, dpi = 300)

save(FC.atac, file = "output/FC_ATAC_harmony_int_named.Rdata")


# Write tables describing contribution of nuclei to each cluster
write.csv(table(FC.atac$dataset), file = "output/FC_ATAC_hm_nuclei_per_sample.csv")
write.csv(table(Idents(FC.atac), FC.atac$dataset), file = "output/FC_ATAC_hm_nuclei_per_sample_per_celltype.csv")
write.csv(prop.table(table(Idents(FC.atac), FC.atac$dataset), margin = 2), file = "output/FC_ATAC_hm_ATAC_prop_nuclei_per_sample_per_celltype.csv")

# Make a marker gene heatmap
## Set plotting order
levels(FC.atac) <- c("MGL.4", "AST.3", "AST.14", "OPC.5", "ODC.0", "ODC.1", "VC.15", "ExN.2", "ExN.6", 
                        "ExN.9", "ExN.10", "ExN.12", "ExN.16", "ExN.17", "ExN.18", "ExN.19", "ExN.21", 
                        "InN.7", "InN.8", "InN.11", "InN.13", "InN.20")

DefaultAssay(FC.atac) <- "SCT"
marker.genes <- c("CSF1R", "P2RY12", "ITGAM", "CX3CR1", "SLC1A2", "AQP4", "GFAP", "PDGFRA", "MYT1", "CSPG4", "PLP1", "MBP", "OPALIN", "CLDN5", "ICAM2", "COLEC12", "GGT5", "SATB2", "SLC17A7", "SLC17A6", "THEMIS", "FEZF2", "VAT1L", "GAD1", "GAD2", "PVALB", "SST", "VIP", "ADARB2")
cluster.averages <- AverageExpression(FC.atac, features = marker.genes, assay = "SCT", return.seurat = TRUE)


p2 <- DoHeatmap(cluster.averages, features = marker.genes, size = 3, draw.lines = F) +
  scale_fill_gradientn(colors = c("#0000FF", "#CCFFFF", "#FF3300"))
ggsave(file = "plots/NEW_FC_ATAC_marker_heatmap.png", plot = p2, height = 7.5, width = 7, dpi = 300)
ggsave(file = "plots/NEW_FC_ATAC_marker_heatmap.svg", plot = p2, height = 7.5, width = 7)


# Find differentially accessible peaks between MGL and all other clusters
DefaultAssay(FC.atac) <- "peaks"

da_peaks.MGL4 <- FindMarkers(
  object = FC.atac,
  ident.1 = "MGL.4", 
  #ident.2 = ,
  min.pct = 0.3,
  test.use = 'LR',
  latent.vars = 'peak_region_fragments'
)
write.csv(da_peaks.MGL4, file = "output/FC_ATAC_hm_MGL4_da_peaks.csv")

# Export da peaks in MGL on chr12
da_peaks.df <- da_peaks.MGL4 %>% rownames_to_column(var = "peak_coord")
da_peaks.df <- da_peaks.df %>% separate(peak_coord, c("chr", "start", "end"), sep = "-")
da_peaks.chr12 <- da_peaks.df %>% filter(chr == "chr12")

#write bed file to view in IGV
write.table(da_peaks.chr12[ , 1:3], file = "output/FC_ATAC_MGL4_da_peaks_chr12.txt", quote = F, row.names = F, col.names = F, sep = "\t")
write.table(da_peaks.chr12[ , 1:3], file = "output/FC_ATAC_MGL4_da_peaks_chr12.bed", quote = F, row.names = F, col.names = F)

# Find da_peaks in MGL near LRRK2
da_peaks.df <- da_peaks.MGL4 %>% rownames_to_column(var = "peak_coord")
da_peaks.MGL.LRRK2 <- filter(da_peaks.df, grepl("chr12-402", peak_coord))
da_peaks.MGL.LRRK2

# Determine which genes are closest to these da peaks
#open_MGL <- rownames(da_peaks.MGL4[da_peaks.MGL4$avg_log2FC > 0.25, ])
#closest_MGL <- ClosestFeature(FC.atac.hm, open_MGL)

# View LRRK2 promoter region
p3 <- CoveragePlot(
  object = FC.atac,
  region = "chr12-40206700-40228700",
  extend.upstream = 0,
  extend.downstream = 0,
  ncol = 1
)
ggsave(file = "plots/FC_ATAC_LRRK2_region_pileups_all_celltypes.png", plot = p3, width = 7, height = 18, dpi = 300)
ggsave(file = "plots/FC_ATAC_LRRK2_region_pileups_all_celltypes.svg", plot = p3, width = 7, height = 18)


# Visualize da peak at LRRK2 TSS
p4 <- VlnPlot(
  object = FC.atac,
  features = da_peaks.MGL.LRRK2$peak_coord[1],
  pt.size = 0
) + NoLegend()

p5 <- FeaturePlot(
  object = FC.atac,
  features = da_peaks.MGL.LRRK2$peak_coord[1],
  pt.size = 0.1,
  cols = c("#CCFFFF", "lightgrey", "#FF0000"),
  label = TRUE,
  order = T,
  label.size = 4
)

p6 <- p4 | p5
ggsave("plots/FC_ATAC_vis_MGL_da_peak_near_LRRK2.png", plot = p6, width = 10, height = 5, dpi = 300)


p5 <- FeaturePlot(
  object = FC.atac,
  features = da_peaks.MGL.LRRK2$peak_coord[1],
  pt.size = 0.1,
  cols = c("#CCFFFF", "lightgrey", "#FF0000"),
  label = TRUE,
  order = T,
  label.size = 3
)
ggsave(file = "plots/NEW_LRRK2_DApeak_MGL4_vs_AllOthers_Feature.png", plot = p5, height = 4, width = 4, dpi = 300)
ggsave(file = "plots/NEW_LRRK2_DApeak_MGL4_vs_AllOthers_Feature.svg", plot = p5, height = 4, width = 4)


# Any difference between "CC" and "TT" MGL.4?
DefaultAssay(FC.atac) <- 'SCT'

p7 <- VlnPlot(FC.atac, "LRRK2", pt.size = 0, split.by = "rs76904798", cols = c("yellow", "orchid1", "turquoise3"))
ggsave("plots/FC_ATAC_LRRK2_Vln_hm.png", plot = p7, height = 3, width = 10, dpi = 300)
ggsave("plots/FC_ATAC_LRRK2_Vln_hm.svg", plot = p7, height = 3, width = 10)

# Write table describing genotype info about nuclei contributing to each cluster
FC.atac$celltype.rs76904798 <- paste(Idents(FC.atac), FC.atac$rs76904798, sep = "_")
Idents(FC.atac) <- "celltype.rs76904798"
write.csv(table(Idents(FC.atac)), file = "output/FC_ATAC_nuclei_per_celltype_per_genotype_new.csv")

# Check if LRRK2 is DE gene between CC and TT cells of any cell type
rs.effectODC.0_CCvsTT <- FindMarkers(FC.atac, ident.1 = "ODC.0_CC", ident.2 = "ODC.0_TT", verbose = FALSE)
write.csv(rs.effectODC.0_CCvsTT, file = "output/rs_effect_ODC0_CCvsTT_markers.csv")

rs.effectODC.1_CCvsTT <- FindMarkers(FC.atac, ident.1 = "ODC.1_CC", ident.2 = "ODC.1_TT", verbose = FALSE)
write.csv(rs.effectODC.1_CCvsTT, file = "output/rs_effect_ODC1_CCvsTT_markers.csv")

rs.effectExN.2_CCvsTT <- FindMarkers(FC.atac, ident.1 = "ExN.2_CC", ident.2 = "ExN.2_TT", verbose = FALSE)
write.csv(rs.effectExN.2_CCvsTT, file = "output/rs_effect_ExN2_CCvsTT_markers.csv")

rs.effectExN.6_CCvsTT <- FindMarkers(FC.atac, ident.1 = "ExN.6_CC", ident.2 = "ExN.6_TT", verbose = FALSE)
write.csv(rs.effectExN.6_CCvsTT, file = "output/rs_effect_ExN6_CCvsTT_markers.csv")

rs.effectExN.9_CCvsTT <- FindMarkers(FC.atac, ident.1 = "ExN.9_CC", ident.2 = "ExN.9_TT", verbose = FALSE)
write.csv(rs.effectExN.9_CCvsTT, file = "output/rs_effect_ExN9_CCvsTT_markers.csv")

rs.effectExN.10_CCvsTT <- FindMarkers(FC.atac, ident.1 = "ExN.10_CC", ident.2 = "ExN.10_TT", verbose = FALSE)
write.csv(rs.effectExN.10_CCvsTT, file = "output/rs_effect_ExN10_CCvsTT_markers.csv")

rs.effectExN.12_CCvsTT <- FindMarkers(FC.atac, ident.1 = "ExN.12_CC", ident.2 = "ExN.12_TT", verbose = FALSE)
write.csv(rs.effectExN.12_CCvsTT, file = "output/rs_effect_ExN12_CCvsTT_markers.csv")

rs.effectExN.16_CCvsTT <- FindMarkers(FC.atac, ident.1 = "ExN.16_CC", ident.2 = "ExN.16_TT", verbose = FALSE)
write.csv(rs.effectExN.16_CCvsTT, file = "output/rs_effect_ExN16_CCvsTT_markers.csv")

rs.effectExN.17_CCvsTT <- FindMarkers(FC.atac, ident.1 = "ExN.17_CC", ident.2 = "ExN.17_TT", verbose = FALSE)
write.csv(rs.effectExN.17_CCvsTT, file = "output/rs_effect_ExN17_CCvsTT_markers.csv")

rs.effectExN.18_CCvsTT <- FindMarkers(FC.atac, ident.1 = "ExN.18_CC", ident.2 = "ExN.18_TT", verbose = FALSE)
write.csv(rs.effectExN.18_CCvsTT, file = "output/rs_effect_ExN18_CCvsTT_markers.csv")

rs.effectExN.19_CCvsTT <- FindMarkers(FC.atac, ident.1 = "ExN.19_CC", ident.2 = "ExN.19_TT", verbose = FALSE)
write.csv(rs.effectExN.19_CCvsTT, file = "output/rs_effect_ExN19_CCvsTT_markers.csv")

rs.effectExN.21_CCvsTT <- FindMarkers(FC.atac, ident.1 = "ExN.21_CC", ident.2 = "ExN.21_TT", verbose = FALSE)
write.csv(rs.effectExN.21_CCvsTT, file = "output/rs_effect_ExN21_CCvsTT_markers.csv")

rs.effectAST.3_CCvsTT <- FindMarkers(FC.atac, ident.1 = "AST.3_CC", ident.2 = "AST.3_TT", verbose = FALSE)
write.csv(rs.effectAST.3_CCvsTT, file = "output/rs_effect_AST3_CCvsTT_markers.csv")

rs.effectAST.14_CCvsTT <- FindMarkers(FC.atac, ident.1 = "AST.14_CC", ident.2 = "AST.14_TT", verbose = FALSE)
write.csv(rs.effectAST.14_CCvsTT, file = "output/rs_effect_AST14_CCvsTT_markers.csv")

rs.effectMGL.4_CCvsTT <- FindMarkers(FC.atac, ident.1 = "MGL.4_CC", ident.2 = "MGL.4_TT", verbose = FALSE)
write.csv(rs.effectMGL.4_CCvsTT, file = "output/rs_effect_MGL4_CCvsTT_markers.csv")

rs.effectOPC.5_CCvsTT <- FindMarkers(FC.atac, ident.1 = "OPC.5_CC", ident.2 = "OPC.5_TT", verbose = FALSE)
write.csv(rs.effectOPC.5_CCvsTT, file = "output/rs_effect_OPC5_CCvsTT_markers.csv")

rs.effectInN.7_CCvsTT <- FindMarkers(FC.atac, ident.1 = "InN.7_CC", ident.2 = "InN.7_TT", verbose = FALSE)
write.csv(rs.effectInN.7_CCvsTT, file = "output/rs_effect_InN7_CCvsTT_markers.csv")

rs.effectInN.8_CCvsTT <- FindMarkers(FC.atac, ident.1 = "InN.8_CC", ident.2 = "InN.8_TT", verbose = FALSE)
write.csv(rs.effectInN.8_CCvsTT, file = "output/rs_effect_InN8_CCvsTT_markers.csv")

rs.effectInN.11_CCvsTT <- FindMarkers(FC.atac, ident.1 = "InN.11_CC", ident.2 = "InN.11_TT", verbose = FALSE)
write.csv(rs.effectInN.11_CCvsTT, file = "output/rs_effect_InN11_CCvsTT_markers.csv")

rs.effectInN.13_CCvsTT <- FindMarkers(FC.atac, ident.1 = "InN.13_CC", ident.2 = "InN.13_TT", verbose = FALSE)
write.csv(rs.effectInN.13_CCvsTT, file = "output/rs_effect_InN13_CCvsTT_markers.csv")

rs.effectInN.20_CCvsTT <- FindMarkers(FC.atac, ident.1 = "InN.20_CC", ident.2 = "InN.20_TT", verbose = FALSE)
write.csv(rs.effectInN.20_CCvsTT, file = "output/rs_effect_InN20_CCvsTT_markers.csv")

rs.effectVC.15_CCvsTT <- FindMarkers(FC.atac, ident.1 = "VC.15_CC", ident.2 = "VC.15_TT", verbose = FALSE)
write.csv(rs.effectVC.15_CCvsTT, file = "output/rs_effect_VC15_CCvsTT_markers.csv")


# Find differentially accessible peaks between "TT" and "CC" microglia
## Switch back to working with peaks instead of gene activities
DefaultAssay(FC.atac) <- 'peaks'

da_peaks_CCvsTT <- FindMarkers(
  object = FC.atac,
  ident.1 = "MGL.4_TT",
  ident.2 = "MGL.4_CC",
  min.pct = 0.2,
  test.use = 'LR',
  latent.vars = 'peak_region_fragments'
)

da_peaks.df <- da_peaks_CCvsTT %>% rownames_to_column(var = "peak_coord")
da_peaks.CCvsTT.LRRK2 <- filter(da_peaks.df, grepl("chr12-402", peak_coord))
da_peaks.CCvsTT.LRRK2 #none at LRRK2


# LRRK2 region CC vs TT MGL
p8 <- CoveragePlot(
  object = FC.atac,
  region = "chr12-40206700-40400000",
  idents = c("MGL.4_TT", "MGL.4_CC"),
  extend.upstream = 0,
  extend.downstream = 0,
  ncol = 1
)
ggsave(file = "plots/FC_ATAC_LRRK2_pileups_CCvsTT_MGL4.png", plot = p8, width = 8, height = 3, dpi = 300)


# Pull average expression of LRRK2 per donor in each cluster
## Note this time used SCT normalization, previously used lognormalize method for RNA assay
Idents(FC.atac) <- "celltype"
LRRK2_av.expression <- AverageExpression(FC.atac, assay = "SCT", features = "LRRK2", add.ident = "dataset", return.seurat = F)
LRRK2_av.df <- as.data.frame(LRRK2_av.expression$SCT)
LRRK2_av.df <- tibble::rownames_to_column(LRRK2_av.df, var = "Gene")
write.csv(file = "output/Set12_av_LRRK2_accessibility_per_donor.csv", quote = F, row.names = F, LRRK2_av.df)





