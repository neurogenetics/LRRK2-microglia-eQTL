#!/bin/env Rscript

# Load necessary packages
library(Seurat)
packageVersion('Seurat')

library(dplyr)
library(ggplot2)
library(cowplot)

# Set up future for parallelization
library(future)
library(future.apply)
plan("multiprocess", workers = 10)
options(future.globals.maxSize = 40000 * 1024^2)

# Read in data
S1363.data <- Read10X(data.dir = "/data/langstonrg/snRNAseq/data/FC_Set15_rs76904798/S1363/filtered_feature_bc_matrix/")
S5123.data <- Read10X(data.dir = "/data/langstonrg/snRNAseq/data/FC_Set15_rs76904798/S5123/filtered_feature_bc_matrix/")
S630.data <- Read10X(data.dir = "/data/langstonrg/snRNAseq/data/FC_Set15_rs76904798/S630/filtered_feature_bc_matrix/")
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
S1363 <- CreateSeuratObject(counts = S1363.data, min.cells = 3, min.features = 500, project = "S1363")
S1363$rs76904798 <- "CC"
S1363$batch <- 2

S5123 <- CreateSeuratObject(counts = S5123.data, min.cells = 3, min.features = 500, project = "S5123")
S5123$rs76904798 <- "CT"
S5123$batch <- 2

S630 <- CreateSeuratObject(counts = S630.data, min.cells = 3, min.features = 500, project = "S630")
S630$rs76904798 <- "CC"
S630$batch <- 3

S4022 <- CreateSeuratObject(counts = S4022.data, min.cells = 3, min.features = 500, project = "S4022")
S4022$rs76904798 <- "TT"
S4022$batch <- 3

S5079 <- CreateSeuratObject(counts = S5079.data, min.cells = 3, min.features = 500, project = "S5079")
S5079$rs76904798 <- "TT"
S5079$batch <- 3

S4924 <- CreateSeuratObject(counts = S4924.data, min.cells = 3, min.features = 500, project = "S4924")
S4924$rs76904798 <- "TT"
S4924$batch <- 4

S1135 <- CreateSeuratObject(counts = S1135.data, min.cells = 3, min.features = 500, project = "S1135")
S1135$rs76904798 <- "CT"
S1135$batch <- 4

S4724 <- CreateSeuratObject(counts = S4724.data, min.cells = 3, min.features = 500, project = "S4724")
S4724$rs76904798 <- "CT"
S4724$batch <- 4

S794 <- CreateSeuratObject(counts = S794.data, min.cells = 3, min.features = 500, project = "S794")
S794$rs76904798 <- "CC"
S794$batch <- 4

S1209 <- CreateSeuratObject(counts = S1209.data, min.cells = 3, min.features = 500, project = "S1209")
S1209$rs76904798 <- "CT"
S1209$batch <- 1

S1230 <- CreateSeuratObject(counts = S1230.data, min.cells = 3, min.features = 500, project = "S1230")
S1230$rs76904798 <- "TT"
S1230$batch <- 1

S1274 <- CreateSeuratObject(counts = S1274.data, min.cells = 3, min.features = 500, project = "S1274")
S1274$rs76904798 <- "CC"
S1274$batch <- 1


# Normalize (SCTransform) and Integrate
FC_Set12.list <- c(S1363, S5123, S630, S4022, S5079, S4924, S1135, S4724, S794, S1209, S1230, S1274)
FC_Set12.list <- future_lapply(X = FC_Set12.list, FUN = function(x) {
  x <- PercentageFeatureSet(x, pattern = "^MT-", col.name = "percent.mt")
  x <- SCTransform(x, vars.to.regress = "percent.mt", verbose = FALSE)
})

frontal.features <- SelectIntegrationFeatures(object.list = FC_Set12.list)
FC_Set12.list <- future_lapply(X = FC_Set12.list, FUN = function(x) {
  x <- RunPCA(x, features = frontal.features, verbose = FALSE)
})

FC_Set12.list <- PrepSCTIntegration(object.list = FC_Set12.list, anchor.features = frontal.features)
frontal.anchors <- FindIntegrationAnchors(object.list = FC_Set12.list, normalization.method = "SCT", anchor.features = frontal.features, reduction = "rpca", dims = 1:30)
frontal.integrated <- IntegrateData(anchorset = frontal.anchors, normalization.method = "SCT", dims = 1:30)
frontal.integrated

#Check out integrated object
p1 <- VlnPlot(object = frontal.integrated, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave("plots/FC_Set12_characterization1.png", plot = p1, width = 12, height = 5, units = "in", dpi = 300)

p2 <- FeatureScatter(object = frontal.integrated, feature1 = "nCount_RNA", feature2 = "percent.mt") + NoLegend()
p3 <- FeatureScatter(object = frontal.integrated, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
ggsave("plots/FC_Set12_characterization2.png", plot = (plot_grid(p2, p3)), width = 10, height = 5, units = "in", dpi = 300)


#Cluster
frontal.integrated <- RunPCA(frontal.integrated, verbose = FALSE)
frontal.integrated <- RunUMAP(frontal.integrated, dims = 1:30, verbose = FALSE)
frontal.integrated <- FindNeighbors(frontal.integrated, dims = 1:30, verbose = FALSE)
frontal.integrated <- FindClusters(frontal.integrated, resolution = 0.5)
save(file = "output/FC_Set12_SeuratAnalysis_clustered.Rdata", object = frontal.integrated)

p4 <- DimPlot(frontal.integrated, reduction = "umap", group.by = "orig.ident")
p5 <- DimPlot(frontal.integrated, reduction = "umap", label = TRUE) + NoLegend()
ggsave("plots/FC_Set12_UMAP_Clusters_1.png", plot = (plot_grid(p4, p5)), width = 10, height = 5, units = "in", dpi = 300)

# Switch to SCT assay
DefaultAssay(frontal.integrated) <- "SCT"
frontal.integrated 

# Take a look at clusters
p6 <- VlnPlot(frontal.integrated, c("SLC17A7", "RORB", "THEMIS", "FEZF2" ), pt.size = 0, ncol = 2) + NoLegend()
ggsave("plots/markers/SET12_Vln_ExN_SLC17A7_RORB_THEMIS_FEZF2.png", plot = p6, width = 12, height = 5)

p7 <- VlnPlot(frontal.integrated, c("GAD1", "SST", "VIP", "PVALB"), pt.size = 0, ncol = 2) + NoLegend()
ggsave("plots/markers/SET12_Vln_InN_GAD1_SST_VIP_PVALB.png", plot = p7, width = 12, height = 5)

p8 <- VlnPlot(frontal.integrated, c("PDGFRA", "OLIG2", "PLP1", "MOBP"), pt.size = 0, ncol = 2) + NoLegend()
ggsave("plots/markers/SET12_Vln_OPC_ODC_PDGFRA_MOBP.png", plot = p8, width = 12, height = 5)

p9 <- VlnPlot(frontal.integrated, c("SLC1A2", "SLC1A3", "GFAP", "AQP4"), pt.size = 0, ncol = 2) + NoLegend()
ggsave("plots/markers/SET12_Vln_AST_SLC1A2_SLC1A3_GFAP_AQP4.png", plot = p9, width = 12, height = 5)

p10 <- VlnPlot(frontal.integrated, c("P2RY12", "ITGAM", "CSF1R", "CX3CR1"), pt.size = 0, ncol = 2) + NoLegend()  #ITGAM is CD11b
ggsave("plots/markers/SET12_Vln_P2RY12_CSF1R_CX3CR1_ITGAM.png", plot = p10, width = 12, height = 5)

p11 <- VlnPlot(frontal.integrated, c("CLDN5", "COLEC12"), pt.size = 0, ncol = 1) + NoLegend()
ggsave("plots/markers/SET12_Vln_CLDN5_COLEC12.png", plot = p11, width = 8, height = 6)

p12 <- VlnPlot(frontal.integrated, c("CUX2", "SNCA"), pt.size = 0, ncol = 1) + NoLegend()
ggsave("plots/markers/SET12_Vln_CUX2_SNCA.png", plot = p12, width = 8, height = 6)

frontal.integrated.markers <- FindAllMarkers(object = frontal.integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5)
write.csv(frontal.integrated.markers, file = "output/FC_Set12_snRNAseq_markers.csv")
top10.frontal.integrated <- frontal.integrated.markers %>% group_by(cluster) %>% top_n(10, avg_log2FC)
write.csv(top10.frontal.integrated, file = "output/top10_FC_Set12_snRNAseq_markers.csv")



# sbatch --partition=largemem --mem=375g --cpus-per-task=10 --gres=lscratch:100 --time=16:00:00

