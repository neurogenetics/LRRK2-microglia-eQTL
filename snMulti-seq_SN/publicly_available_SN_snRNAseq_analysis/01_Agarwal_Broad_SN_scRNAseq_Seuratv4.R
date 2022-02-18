#!/bin/env Rscript

# Load necessary packages
library(Seurat)
packageVersion('Seurat')

library(sctransform)
library(dplyr)
library(ggplot2)
library(cowplot)

# Set up future for parallelization
library(future)
library(future.apply)
plan("multiprocess", workers = 10)
options(future.globals.maxSize = 40000 * 1024^2)

# Read in data
## Retrieved from: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE140231  PMID: 32826893
SN3.data <- Read10X(data.dir = "Agarwal_SN/data/GSM4157068_Sample_6_N3/")
SN4.data <- Read10X(data.dir = "Agarwal_SN/data/GSM4157069_Sample_7_N4/")
SN5.data <- Read10X(data.dir = "Agarwal_SN/data/GSM4157070_Sample_8_N5/")
SN1B.data <- Read10X(data.dir = "Agarwal_SN/data/GSM4157072_Sample_10_N1B/")
SN2B.data <- Read10X(data.dir = "Agarwal_SN/data/GSM4157074_Sample_12_N2B/")
SN4B.data <- Read10X(data.dir = "Agarwal_SN/data/GSM4157076_Sample_14_N4B/")
SN5B.data <- Read10X(data.dir = "Agarwal_SN/data/GSM4157078_Sample_16_N5B/")

## Retrieved from: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE126836  PMID: 31178122
S5534.data <- Read10X(data.dir = "Broad_SN/data/SN_MD5534_bc_gene_matrix/")
S5828.data <- Read10X(data.dir = "Broad_SN/data/SN_MD5828_bc_gene_matrix/")
S5840.data <- Read10X(data.dir = "Broad_SN/data/SN_MD5840_bc_gene_matrix/")
S5862.data <- Read10X(data.dir = "Broad_SN/data/SN_MD5862_bc_gene_matrix/")
S5893.data <- Read10X(data.dir = "Broad_SN/data/SN_MD5893_bc_gene_matrix/")
S6060.data <- Read10X(data.dir = "Broad_SN/data/SN_MD6060_bc_gene_matrix/")
S6063.data <- Read10X(data.dir = "Broad_SN/data/SN_MD6063_bc_gene_matrix/")


# Make Seurat objects
SN3 <- CreateSeuratObject(counts = SN3.data, min.cells = 3, min.features = 500, project = "SN3")
SN4 <- CreateSeuratObject(counts = SN4.data, min.cells = 3, min.features = 500, project = "SN4")
SN5 <- CreateSeuratObject(counts = SN5.data, min.cells = 3, min.features = 500, project = "SN5")
SN1B <- CreateSeuratObject(counts = SN1B.data, min.cells = 3, min.features = 500, project = "SN1B")
SN2B <- CreateSeuratObject(counts = SN2B.data, min.cells = 3, min.features = 500, project = "SN2B")
SN4B <- CreateSeuratObject(counts = SN4B.data, min.cells = 3, min.features = 500, project = "SN4B")
SN5B <- CreateSeuratObject(counts = SN5B.data, min.cells = 3, min.features = 500, project = "SN5B")

S5534 <- CreateSeuratObject(counts = S5534.data, min.cells = 3, min.features = 500, project = "S5534")
S5828 <- CreateSeuratObject(counts = S5828.data, min.cells = 3, min.features = 500, project = "S5828")
S5840 <- CreateSeuratObject(counts = S5840.data, min.cells = 3, min.features = 500, project = "S5840")
S5862 <- CreateSeuratObject(counts = S5862.data, min.cells = 3, min.features = 500, project = "S5862")
S5893 <- CreateSeuratObject(counts = S5893.data, min.cells = 3, min.features = 500, project = "S5893")
S6060 <- CreateSeuratObject(counts = S6060.data, min.cells = 3, min.features = 500, project = "S6060")
S6063 <- CreateSeuratObject(counts = S6063.data, min.cells = 3, min.features = 500, project = "S6063")


# Normalize (SCTransform) and Integrate
SN.list <- c(SN3, SN4, SN5, SN1B, SN2B, SN4B, SN5B, S5534, S5828, S5840, S5862, S5893, S6060, S6063)
SN.list <- future_lapply(X = SN.list, FUN = function(x) {
  x <- PercentageFeatureSet(x, pattern = "^MT-", col.name = "percent.mt")
  x <- SCTransform(x, vars.to.regress = "percent.mt", verbose = FALSE)
})

SN.features <- SelectIntegrationFeatures(object.list = SN.list)
SN.list <- future_lapply(X = SN.list, FUN = function(x) {
  x <- RunPCA(x, features = SN.features, verbose = FALSE)
})

SN.list <- PrepSCTIntegration(object.list = SN.list, anchor.features = SN.features)
SN.anchors <- FindIntegrationAnchors(object.list = SN.list, normalization.method = "SCT", anchor.features = SN.features, reduction = "rpca", dims = 1:30)
SN.integrated <- IntegrateData(anchorset = SN.anchors, normalization.method = "SCT", dims = 1:30)
SN.integrated


#Check out integrated object
p1 <- VlnPlot(object = SN.integrated, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave("plots/SN_integrated_characterization1_v4.png", plot = p1, width = 12, height = 5, units = "in", dpi = 300)

p2 <- FeatureScatter(object = SN.integrated, feature1 = "nCount_RNA", feature2 = "percent.mt") + NoLegend()
p3 <- FeatureScatter(object = SN.integrated, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
ggsave("plots/SN_integrated_characterization2_v4.png", plot = (plot_grid(p2, p3)), width = 10, height = 5, units = "in", dpi = 300)


#Cluster
SN.integrated <- RunPCA(SN.integrated, verbose = FALSE)
SN.integrated <- RunUMAP(SN.integrated, dims = 1:30, verbose = FALSE)
SN.integrated <- FindNeighbors(SN.integrated, dims = 1:30, verbose = FALSE)
SN.integrated <- FindClusters(SN.integrated, resolution = 0.3, verbose = FALSE)
save(file = "output/Agarwal_Broad_SN_Seuratv4_clustered.Rdata", object = SN.integrated)

p4 <- DimPlot(SN.integrated, reduction = "umap", group.by = "orig.ident", pt.size = 0.1)
p5 <- DimPlot(SN.integrated, reduction = "umap", label = TRUE, pt.size = 0.1) + NoLegend()
ggsave("plots/Agarwal_Broad_SN_UMAP_Clusters_v4.png", plot = (plot_grid(p4, p5)), width = 11, height = 5, units = "in", dpi = 300)


# sbatch --mem=200g --cpus-per-task=10 --gres=lscratch:100 --time=4:00:00

