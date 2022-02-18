#!/bin/env Rscript

# Load necessary packages
library(Signac)  #confirm v1.1.1
packageVersion('Signac')

library(Seurat)  #confirm v4.0.3
packageVersion('Seurat')

library(tidyverse)
library(svglite)

# Load integrated object
load("output/SN2960_SN2972_SN2942_Signac_int.Rdata")
SN.int  #14902 nuclei

# Subset to exclude cells with very low quality GEX data
SN.int <- subset(
  x = SN.int,
  subset = nCount_RNA < 25000 &
    nCount_RNA > 200
)

SN.int  #8001 nuclei

# Make QC plot
p1 <- VlnPlot(
  object = SN.int,
  features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),
  ncol = 4,
  pt.size = 0
)
ggsave(file = "plots/SN2960_SN2972_SN2942_Signac_int_QC_Vln_filt.png", plot = p1, width = 10, height = 3)
ggsave(file = "plots/SN2960_SN2972_SN2942_Signac_int_QC_Vln_filt.svg", plot = p1, width = 10, height = 3)


# Normalize GEX data
DefaultAssay(SN.int) <- "RNA"
SN.int <- PercentageFeatureSet(SN.int, pattern = "^MT-", col.name = "percent.mt")
SN.int <- SCTransform(SN.int, vars.to.regress = "percent.mt", verbose = FALSE)
SN.int <- RunPCA(SN.int, verbose = FALSE)


DefaultAssay(SN.int) <- "SCT"


# Annotate cell types based on publicly available datasets
## Agarwal 2020, Welch 2019 - object made with v4 Seurat
load("/data/langstonrg/snRNAseq/Outside_SN/output/Agarwal_Broad_SN_Seuratv4_named.Rdata")
DefaultAssay(SN.integrated) <- "integrated"
SN.integrated

# Transfer cell type labels from reference to query
transfer_anchors <- FindTransferAnchors(
  reference = SN.integrated,
  query = SN.int,
  normalization.method = "SCT",
  reference.reduction = "pca",
  recompute.residuals = FALSE,
  dims = 1:50
)

predictions <- TransferData(
  anchorset = transfer_anchors, 
  refdata = SN.integrated$celltype,
  weight.reduction = SN.int[['pca']],
  dims = 1:50
)

SN.int <- AddMetaData(
  object = SN.int,
  metadata = predictions
)

# Set the cell identities to the cell type predictions, save object, and visualize
Idents(SN.int) <- "predicted.id"

save(SN.int, file = "output/SN2960_SN2972_SN2942_Signac_int_predicted_ids.Rdata")

p2 <- DimPlot(SN.int, reduction = 'umap', label = TRUE) + NoLegend()
ggsave(file = "plots/SN2960_SN2972_SN2942_UMAP_predicted_ids_filt.png", plot = p2, width = 5.1, height = 5, dpi = 300)
ggsave(file = "plots/SN2960_SN2972_SN2942_UMAP_predicted_ids_filt.svg", plot = p2, width = 5.1, height = 5)

# Confirm cell type assignments
levels(SN.int) <- c("DaN.7", "InN.10", "ExN.13", "ODC.0", "ODC.1", "ODC.3", "ODC.5", "ODC.9", "OPC.2", "MGL.6", "AST.4", "VC.8", "VC.12")

p3 <- VlnPlot(SN.int, c("TH", "GAD1", "SLC17A7", "MOBP", "PDGFRA", "CSF1R", "GFAP", "CLDN5", "CYP1B1"), pt.size = 0, ncol = 1)
ggsave("plots/SN2960_SN2972_SN2942_predicted_id_marker_gene_Vln.png", plot = p3, height = 18, width = 8, dpi = 300)
ggsave("plots/SN2960_SN2972_SN2942_predicted_id_marker_gene_Vln.svg", plot = p3, height = 18, width = 8)

# Write tables describing contribution of nuclei to each cluster
write.csv(table(SN.int$orig.ident), file = "output/SN2960_SN2972_SN2942_nuclei_per_sample.csv")
write.csv(table(Idents(SN.int), SN.int$orig.ident), file = "output/SN2960_SN2972_SN2942_nuclei_per_sample_per_celltype.csv")
write.csv(prop.table(table(Idents(SN.int), SN.int$orig.ident), margin = 2), file = "output/SN2960_SN2972_SN2942_prop_nuclei_per_sample_per_celltype.csv")


# Look at LRRK2
p4 <- FeaturePlot(SN.int, reduction = "umap", "LRRK2", cols = c("#CCFFFF", "lightgrey", "#FF0000"), label = TRUE, order = T, pt.size = 0.1)
ggsave("plots/SN2960_SN2972_SN2942_Feature_predicted_id_LRRK2.png", plot = p4, height = 5, width = 5.5)

# Write table of average LRRK2 RNA expression per cluster
cl.av.LRRK2 <- AverageExpression(SN.int, assays = "SCT", features = "LRRK2")
write.csv(file = "output/SN2960_SN2972_SN2942_av_LRRK2_per_cluster.csv", cl.av.LRRK2$SCT)

# Make coverage plot over LRRK2 promoter region
DefaultAssay(SN.int) <- "peaks"

levels(SN.int) <- c("MGL.6", "AST.4", "OPC.2", "ODC.0", "ODC.1", "ODC.3", "ODC.5", "ODC.9", "VC.8", "VC.12", "DaN.7", "InN.10", "ExN.13")

p5 <- CoveragePlot(
  object = SN.int,
  region = "chr12-40206700-40228700",
  features = "LRRK2",
  expression.assay = "SCT"
)

ggsave(file = "plots/SN2960_SN2972_SN2942_CoveragePlot_LRRK2_promoter.png", plot = p5, width = 8, height = 11, dpi = 300)
ggsave(file = "plots/SN2960_SN2972_SN2942_CoveragePlot_LRRK2_promoter.svg", plot = p5, width = 8, height = 11)

# Find differentially accessible peaks between MGL.6 and all other clusters
da_peaks.MGL <- FindMarkers(
  object = SN.int,
  ident.1 = "MGL.6", 
  #ident.2 = ,
  min.pct = 0.3,
  test.use = 'LR',
  latent.vars = 'atac_peak_region_fragments'
)
write.csv(da_peaks.MGL, file = "output/SN2960_SN2972_SN2942_MGL6_da_peaks.csv")

# Visualize da peak at LRRK2 TSS
da_peaks.MGL.LRRK2 <- da_peaks.MGL[grepl( "chr12-402", rownames(da_peaks.MGL) ), ]
da_peaks.MGL.LRRK2[ , c(2, 5)]

p6 <- FeaturePlot(
  object = SN.int,
  features = rownames(da_peaks.MGL.LRRK2)[1],
  #pt.size = 0.1,
  cols = c("#CCFFFF", "lightgrey", "#FF0000"),
  label = TRUE,
  order = T,
  label.size = 4
)
ggsave("plots/SN2960_SN2972_SN2942_Feature_MGL6_da_peak_near_LRRK2.png", plot = p6, width = 5.5, height = 5, dpi = 300)
ggsave("plots/SN2960_SN2972_SN2942_Feature_MGL6_da_peak_near_LRRK2.svg", plot = p6, width = 5.5, height = 5)





# sbatch --mem=200g --gres=lscratch:100 --time=8:00:00


