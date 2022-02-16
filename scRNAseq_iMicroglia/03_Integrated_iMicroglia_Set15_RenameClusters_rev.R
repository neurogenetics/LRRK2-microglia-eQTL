#!/bin/env Rscript

# sinteractive --mem=100g --gres=lscratch:50

# Load necessary packages
library(Seurat)  #v4.0.1
library(ggplot2)
library(svglite)

# Load object produced by script iMicroglia_vs_Set15FCsnRNAseq.R
load("output/FCSet15_iMGL_integrated_clustered.Rdata")

# Switch to SCT assay
DefaultAssay(FC_iMGL.integrated) <- "SCT"
FC_iMGL.integrated # 127,662 cells

# Save numerical cluster IDs before renaming based on marker gene expression
FC_iMGL.integrated[["num.ident"]] <- Idents(object = FC_iMGL.integrated)

# Exclude clusters with unclear/multiple cell type markers - 25, 41, 44, 46, 47
FC_iMGL.integrated <- subset(FC_iMGL.integrated, idents = c(25, 41, 44, 46, 47), invert = TRUE)
FC_iMGL.integrated # 125,036 cells

# Name OPC and ODC (PDGFRA, OLIG2 and PLP1, MOBP)
FC_iMGL.integrated <- RenameIdents(FC_iMGL.integrated, `6` = "OPC.6", `0` = "ODC.0", `5` = "ODC.5", `9` = "ODC.9")

# Name MGL and AST
FC_iMGL.integrated <- RenameIdents(FC_iMGL.integrated, `3` = "AST.3", `24` = "AST.24", `1` = "iMGL.1", `10` = "MGL.10", `39` = "MGL.39", `42` = "MGL.42")

# Name VC
FC_iMGL.integrated <- RenameIdents(FC_iMGL.integrated, `28` = "VC.28", `29` = "VC.29", `37` = "VC.37", `43` = "VC.43")

# Name ExN
FC_iMGL.integrated <- RenameIdents(FC_iMGL.integrated, `2` = "ExN.2", `4` = "ExN.4", `7` = "ExN.7", `8` = "ExN.8", `11` = "ExN.11", `14` = "ExN.14", `15` = "ExN.15", `16` = "ExN.16", `17` = "ExN.17", `19` = "ExN.19",`20` = "ExN.20", `21` = "ExN.21", `30` = "ExN.30", `31` = "ExN.31", `32` = "ExN.32", `38` = "ExN.38", `40` = "ExN.40", `45` = "ExN.45")

# Name InN
FC_iMGL.integrated <- RenameIdents(FC_iMGL.integrated, `12` = "InN.12", `13` = "InN.13", `18` = "InN.18", `22` = "InN.22", `23` = "InN.23", `26` = "InN.26", `27` = "InN.27", `33` = "InN.33", `34` = "InN.34", `35` = "InN.35", `36` = "InN.36")


p8 <- DimPlot(FC_iMGL.integrated, reduction = 'umap', label = TRUE) + NoLegend()
ggsave("plots/FCSet15_iMGL_labeled_UMAP_rev.png", plot = p8, width = 7, height = 5)

FC_iMGL.integrated$celltype <- Idents(FC_iMGL.integrated)

write.csv(table(Idents(FC_iMGL.integrated), FC_iMGL.integrated$orig.ident), file = "output/FCSet15_iMGL_cells_per_sample_per_celltype_rev.csv")
write.csv(table(Idents(FC_iMGL.integrated), FC_iMGL.integrated$type), file = "output/FCSet15_iMGL_cells_per_type_per_celltype_rev.csv")


# Name clusters by broad cell type
Idents(FC_iMGL.integrated) <- "num.ident"

# Name OPC and ODC (PDGFRA, OLIG2 and PLP1, MOBP)
FC_iMGL.integrated <- RenameIdents(FC_iMGL.integrated, `6` = "OPC", `0` = "ODC", `5` = "ODC", `9` = "ODC")

# Name MGL and AST
FC_iMGL.integrated <- RenameIdents(FC_iMGL.integrated, `3` = "AST", `24` = "AST", `1` = "iMGL", `10` = "MGL", `39` = "MGL", `42` = "MGL")

# Name VC
FC_iMGL.integrated <- RenameIdents(FC_iMGL.integrated, `28` = "VC", `29` = "VC", `37` = "VC", `43` = "VC")

# Name ExN
FC_iMGL.integrated <- RenameIdents(FC_iMGL.integrated, `2` = "ExN", `4` = "ExN", `7` = "ExN", `8` = "ExN", `11` = "ExN", `14` = "ExN", `15` = "ExN", `16` = "ExN", `17` = "ExN", `19` = "ExN",`20` = "ExN", `21` = "ExN", `30` = "ExN", `31` = "ExN", `32` = "ExN", `38` = "ExN", `40` = "ExN", `45` = "ExN")

# Name InN
FC_iMGL.integrated <- RenameIdents(FC_iMGL.integrated, `12` = "InN", `13` = "InN", `18` = "InN", `22` = "InN", `23` = "InN", `26` = "InN", `27` = "InN", `33` = "InN", `34` = "InN", `35` = "InN", `36` = "InN")

FC_iMGL.integrated$broad_celltype <- Idents(FC_iMGL.integrated)

#save(file = "output/FCSet15_iMGL_clustered_named.Rdata", object = FC_iMGL.integrated)

p1 <- DimPlot(FC_iMGL.integrated, reduction = 'umap', group.by = "broad_celltype", pt.size = 0.08, label = TRUE, label.size = 5, raster = F) + NoLegend()
ggsave("plots/FCSet15_iMGL_UMAP_by_broadcelltype_rev.png", plot = p1, width = 6, height = 6.2)
ggsave("plots/FCSet15_iMGL_UMAP_by_broadcelltype_rev.svg", plot = p1, width = 6, height = 6.2)

p2 <- DimPlot(FC_iMGL.integrated, reduction = 'umap', group.by = "type", cols = c("aquamarine2", "#FF00FF"), pt.size = 0.08, label = TRUE, label.size = 5)
ggsave("plots/FCSet15_iMGL_UMAP_by_sourcetype_rev.png", plot = p2, width = 6.5, height = 6)
ggsave("plots/FCSet15_iMGL_UMAP_by_sourcetype_rev.svg", plot = p2, width = 6.5, height = 6)

FC_iMGL.integrated$broad_celltype <- Idents(FC_iMGL.integrated)

write.csv(table(Idents(FC_iMGL.integrated), FC_iMGL.integrated$orig.ident), file = "output/FCSet15_iMGL_cells_per_sample_per_broadcelltype_rev.csv")
write.csv(table(Idents(FC_iMGL.integrated), FC_iMGL.integrated$type), file = "output/FCSet15_iMGL_cells_per_type_per_broadcelltype_rev.csv")

# Make heatmap showing Pearson correlation between broad celltypes
av.exp <- AverageExpression(FC_iMGL.integrated, "SCT")$SCT
write.csv(av.exp, file = "output/FC_iMGL_integrated_SCT_AvExp_rev.csv") #Used in FCSet15_iMGL_integration_plots_rev.Rmd



