#!/bin/env Rscript

# Load necessary packages
library(Seurat)
library(tidyverse)
library(sctransform)
library(cowplot)

# Set up future for parallelization
library(future)
library(future.apply)
plan("multiprocess", workers = 4)
options(future.globals.maxSize = 20000 * 1024^2)

# Read in FC snRNAseq data (Set15)
S1027.data <- Read10X(data.dir = "/data/langstonrg/snRNAseq/data/FC_Set15_rs76904798/S1027/filtered_feature_bc_matrix/")
S1672.data <- Read10X(data.dir = "/data/langstonrg/snRNAseq/data/FC_Set15_rs76904798/S1672/filtered_feature_bc_matrix/")
S1363.data <- Read10X(data.dir = "/data/langstonrg/snRNAseq/data/FC_Set15_rs76904798/S1363/filtered_feature_bc_matrix/")
S5123.data <- Read10X(data.dir = "/data/langstonrg/snRNAseq/data/FC_Set15_rs76904798/S5123/filtered_feature_bc_matrix/")
S630.data <- Read10X(data.dir = "/data/langstonrg/snRNAseq/data/FC_Set15_rs76904798/S630/filtered_feature_bc_matrix/")
S1584.data <- Read10X(data.dir = "/data/langstonrg/snRNAseq/data/FC_Set15_rs76904798/S1584/filtered_feature_bc_matrix/")
S4022.data <- Read10X(data.dir = "/data/langstonrg/snRNAseq/data/FC_Set15_rs76904798/S4022/filtered_feature_bc_matrix/")
S5079.data <- Read10X(data.dir = "/data/langstonrg/snRNAseq/data/FC_Set15_rs76904798/S5079/filtered_feature_bc_matrix/")
S4924.data <- Read10X(data.dir = "/data/langstonrg/snRNAseq/data/FC_Set15_rs76904798/S4924/filtered_feature_bc_matrix/")
S1135.data <- Read10X(data.dir = "/data/langstonrg/snRNAseq/data/FC_Set15_rs76904798/S1135/filtered_feature_bc_matrix/")
S4724.data <- Read10X(data.dir = "/data/langstonrg/snRNAseq/data/FC_Set15_rs76904798/S4724/filtered_feature_bc_matrix/")
S794.data <- Read10X(data.dir = "/data/langstonrg/snRNAseq/data/FC_Set15_rs76904798/S794/filtered_feature_bc_matrix/")
S1209.data <- Read10X(data.dir = "/data/langstonrg/snRNAseq/data/FC_Set15_rs76904798/S1209/filtered_feature_bc_matrix/")
S1230.data <- Read10X(data.dir = "/data/langstonrg/snRNAseq/data/FC_Set15_rs76904798/S1230/filtered_feature_bc_matrix/")
S1274.data <- Read10X(data.dir = "/data/langstonrg/snRNAseq/data/FC_Set15_rs76904798/S1274/filtered_feature_bc_matrix/")


# Create Seurat objects, add rs76904798 genotype group and batch
S1027 <- CreateSeuratObject(counts = S1027.data, min.cells = 3, min.features = 500, project = "FC_S1027")
S1027$rs76904798 <- "TT"
S1027$batch <- 2
S1027

S1672 <- CreateSeuratObject(counts = S1672.data, min.cells = 3, min.features = 500, project = "FC_S1672")
S1672$rs76904798 <- "CT"
S1672$batch <- 2
S1672

S1363 <- CreateSeuratObject(counts = S1363.data, min.cells = 3, min.features = 500, project = "FC_S1363")
S1363$rs76904798 <- "CC"
S1363$batch <- 2
S1363

S5123 <- CreateSeuratObject(counts = S5123.data, min.cells = 3, min.features = 500, project = "FC_S5123")
S5123$rs76904798 <- "CT"
S5123$batch <- 2
S5123

S630 <- CreateSeuratObject(counts = S630.data, min.cells = 3, min.features = 500, project = "FC_S630")
S630$rs76904798 <- "CC"
S630$batch <- 3
S630

S1584 <- CreateSeuratObject(counts = S1584.data, min.cells = 3, min.features = 500, project = "FC_S1584")
S1584$rs76904798 <- "CC"
S1584$batch <- 3
S1584

S4022 <- CreateSeuratObject(counts = S4022.data, min.cells = 3, min.features = 500, project = "FC_S4022")
S4022$rs76904798 <- "TT"
S4022$batch <- 3
S4022

S5079 <- CreateSeuratObject(counts = S5079.data, min.cells = 3, min.features = 500, project = "FC_S5079")
S5079$rs76904798 <- "TT"
S5079$batch <- 3
S5079

S4924 <- CreateSeuratObject(counts = S4924.data, min.cells = 3, min.features = 500, project = "FC_S4924")
S4924$rs76904798 <- "TT"
S4924$batch <- 4
S4924

S1135 <- CreateSeuratObject(counts = S1135.data, min.cells = 3, min.features = 500, project = "FC_S1135")
S1135$rs76904798 <- "CT"
S1135$batch <- 4
S1135

S4724 <- CreateSeuratObject(counts = S4724.data, min.cells = 3, min.features = 500, project = "FC_S4724")
S4724$rs76904798 <- "CT"
S4724$batch <- 4
S4724

S794 <- CreateSeuratObject(counts = S794.data, min.cells = 3, min.features = 500, project = "FC_S794")
S794$rs76904798 <- "CC"
S794$batch <- 4
S794

S1209 <- CreateSeuratObject(counts = S1209.data, min.cells = 3, min.features = 500, project = "FC_S1209")
S1209$rs76904798 <- "CT"
S1209$batch <- 1
S1209

S1230 <- CreateSeuratObject(counts = S1230.data, min.cells = 3, min.features = 500, project = "FC_S1230")
S1230$rs76904798 <- "TT"
S1230$batch <- 1
S1230

S1274 <- CreateSeuratObject(counts = S1274.data, min.cells = 3, min.features = 500, project = "FC_S1274")
S1274$rs76904798 <- "CC"
S1274$batch <- 1
S1274

# Read in iMicroglia scRNAseq data and create Seurat Objects
PPMI3453.data <- Read10X(data.dir = "data/PPMI3453_iMicroglia/filtered_feature_bc_matrix/")
PPMI4101.data <- Read10X(data.dir = "data/PPMI4101_iMicroglia/filtered_feature_bc_matrix/")

PPMI3453 <- CreateSeuratObject(counts = PPMI3453.data, min.cells = 3, min.features = 500, project = "iMGL_3453")
PPMI3453$rs76904798 <- "CC"
PPMI3453$batch <- 5
PPMI3453  #5,951 cells

PPMI4101 <- CreateSeuratObject(counts = PPMI4101.data, min.cells = 3, min.features = 500, project = "iMGL_4101")
PPMI4101$rs76904798 <- "CT"
PPMI4101$batch <- 5
PPMI4101  #4,079 cells


# Add source type
PPMI3453$type <- "iPSC_derived"
PPMI4101$type <- "iPSC_derived"

S1027$type <- "brain_FC"
S1672$type <- "brain_FC"
S1363$type <- "brain_FC"
S5123$type <- "brain_FC"
S630$type <- "brain_FC"
S1584$type <- "brain_FC"
S4022$type <- "brain_FC"
S5079$type <- "brain_FC"
S4924$type <- "brain_FC"
S1135$type <- "brain_FC"
S4724$type <- "brain_FC"
S794$type <- "brain_FC"
S1209$type <- "brain_FC"
S1230$type <- "brain_FC"
S1274$type <- "brain_FC"


# Normalize (SCTransform) and Integrate
FC_iMGL.list <- c(S1027, S1672, S1363, S5123, S630, S1584, S4022, S5079, S4924, S1135, S4724, S794, S1209, S1230, S1274, PPMI3453, PPMI4101)
FC_iMGL.list <- future_lapply(X = FC_iMGL.list, FUN = function(x) {
  x <- PercentageFeatureSet(x, pattern = "^MT-", col.name = "percent.mt")
  x <- SCTransform(x, vars.to.regress = "percent.mt", verbose = FALSE)
})

features <- SelectIntegrationFeatures(object.list = FC_iMGL.list)
FC_iMGL.list <- future_lapply(X = FC_iMGL.list, FUN = function(x) {
  x <- RunPCA(x, features = features, verbose = FALSE)
})

FC_iMGL.list <- PrepSCTIntegration(object.list = FC_iMGL.list, anchor.features = features)
anchors <- FindIntegrationAnchors(object.list = FC_iMGL.list, normalization.method = "SCT", anchor.features = features, reduction = "rpca", dims = 1:30)
FC_iMGL.integrated <- IntegrateData(anchorset = anchors, normalization.method = "SCT", dims = 1:30)
FC_iMGL.integrated

save(file = "output/FCSet15_iMGL_integrated.Rdata", object = FC_iMGL.integrated)

#Check out integrated object
p1 <- VlnPlot(object = FC_iMGL.integrated, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
ggsave("plots/FCSet15_iMGL_integration_QC1.png", plot = p1, width = 12, height = 5, units = "in")

p2 <- FeatureScatter(object = FC_iMGL.integrated, feature1 = "nCount_RNA", feature2 = "percent.mt") + NoLegend()
p3 <- FeatureScatter(object = FC_iMGL.integrated, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
ggsave("plots/FCSet15_iMGL_integration_QC2.png", plot = (plot_grid(p2, p3)), width = 10, height = 5, units = "in")

#Cluster
FC_iMGL.integrated <- RunPCA(FC_iMGL.integrated, verbose = FALSE)
FC_iMGL.integrated <- RunUMAP(FC_iMGL.integrated, dims = 1:30, verbose = FALSE)
FC_iMGL.integrated <- FindNeighbors(FC_iMGL.integrated, dims = 1:30, verbose = FALSE)
FC_iMGL.integrated <- FindClusters(FC_iMGL.integrated, resolution = 1, verbose = FALSE)
save(file = "output/FCSet15_iMGL_integrated_clustered.Rdata", object = FC_iMGL.integrated)

p3 <- DimPlot(FC_iMGL.integrated, reduction = "umap", label = TRUE) + NoLegend()
ggsave("plots/FCSet15_iMGL_Clusters_1.png", plot = p3, width = 7, height = 5, units = "in")

p4 <- DimPlot(FC_iMGL.integrated, reduction = "umap", group.by = "orig.ident")
ggsave("plots/FCSet15_iMGL_Clusters_2.png", plot = p4, width = 9, height = 5, units = "in")

p5 <- DimPlot(FC_iMGL.integrated, reduction = "umap", group.by = "type")
ggsave("plots/FCSet15_iMGL_Clusters_3.png", plot = p5, width = 8, height = 5, units = "in")

FC_iMGL.integrated.markers <- FindAllMarkers(object = FC_iMGL.integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(FC_iMGL.integrated.markers, file = "output/FCSet15_iMGL_markers.csv")
top10.FC_iMGL.integrated <- FC_iMGL.integrated.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
write.csv(top10.FC_iMGL.integrated, file = "output/top10_FCSet15_iMGL_markers.csv")


