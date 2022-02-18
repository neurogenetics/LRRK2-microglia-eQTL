#!/bin/env Rscript

# Load necessary packages
library(Signac)  #confirm v1.1.1
packageVersion('Signac')

library(Seurat)  #confirm v4.0.0
packageVersion('Seurat')

library(GenomicRanges)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)

library(patchwork)
library(tidyverse)

library(future)
library(future.apply)
plan("multiprocess", workers = 4)
options(future.globals.maxSize = 40000 * 1024^2)
set.seed(1234)

# Load the harmony integrated Set12 FC snATACseq object
load("output/FC_ATAC_clustered_harmony_int.Rdata")
FC.atac

# Load the Seuratv4 integrated Set12 FC snRNAseq object
load("/data/langstonrg/snRNAseq/FrontalCortex/FC_Set15_Integrated/output/FC_Set12_SeuratAnalysis_clustered.Rdata")
frontal.integrated #91810 nuclei

## Name clusters based on marker gene expression
frontal.integrated$num.ident <- Idents(frontal.integrated)

## Remove cluster 17 (markers of multiple glial celltypes AST, MGL, ODC)
frontal.integrated <- subset(frontal.integrated, idents = 17, invert = T)
frontal.integrated #90035 nuclei

frontal.integrated <- RenameIdents(
  object = frontal.integrated,
  '0' = 'ODC.0',
  '1' = 'ExN.1',
  '2' = 'ExN.2',
  '3' = 'ExN.3',
  '4' = 'ExN.4',
  '5' = 'AST.5',
  '6' = 'ODC.6',
  '7' = 'OPC.7',
  '8' = 'ExN.8',
  '9' = 'InN.9',
  '10' = 'InN.10',
  '11' = 'InN.11',
  '12' = 'MGL.12',
  '13' = 'ExN.13',
  '14' = 'ExN.14',
  '15' = 'ExN.15',
  '16' = 'AST.16',
  '18' = 'VC.18',   
  '19' = 'InN.19',
  '20' = 'ExN.20',
  '21' = 'InN.21',
  '22' = 'InN.22',
  '23' = 'VC.23',  
  '24' = 'ExN.24',
  '25' = 'ExN.25',
  '26' = 'ODC.26',
  '27' = 'ExN.27',
  '28' = 'InN.28',
  '29' = 'InN.29',
  '30' = 'ExN.30',
  '31' = 'ExN.31',
  '32' = 'VC.32',
  '33' = 'InN.33',
  '34' = 'ExN.34'
)

frontal.integrated$celltype <- Idents(frontal.integrated)

p1 <- DimPlot(frontal.integrated, label = TRUE) + NoLegend()
ggsave(file = "plots/FC_RNA_SET12_UMAP_named.png", plot = p1, width = 6, height = 6, dpi = 300)

# Annotate cell types in ATAC dataset based on ids in RNA dataset
#DefaultAssay(FC.atac) <- "SCT"
DefaultAssay(frontal.integrated) <- "integrated"

# transfer cell type labels from reference to query
transfer_anchors <- FindTransferAnchors(
  reference = frontal.integrated,
  query = FC.atac,
  normalization.method = "SCT",
  reference.reduction = "pca",
  recompute.residuals = FALSE,
  dims = 1:50
)

predictions <- TransferData(
  anchorset = transfer_anchors, 
  refdata = frontal.integrated$celltype,
  weight.reduction = FC.atac[['pca']],
  dims = 1:50
)

FC.atac <- AddMetaData(
  object = FC.atac,
  metadata = predictions
)

# set the cell identities to the cell type predictions
Idents(FC.atac) <- "predicted.id"
DefaultAssay(FC.atac) <- "peaks"

p2 <- DimPlot(FC.atac, reduction = 'umap', label = TRUE) + NoLegend()
ggsave(file = "plots/FC_atac_hm_UMAP_predicted_ids.png", plot = p2, width = 6, height = 6, dpi = 300)

p3 <- DimPlot(
  object = frontal.integrated,
  group.by = 'celltype',
  label = TRUE,
  repel = F) + NoLegend() + ggtitle('snRNA-seq')

p4 <- DimPlot(
  object = FC.atac,
  group.by = 'predicted.id',
  label = TRUE,
  repel = T) + NoLegend() + ggtitle('snATAC-seq')

p5 <- p3 | p4
ggsave(file = "plots/NEW_FC_snATACseq_snRNAseq_UMAP_combined.png", plot = p5, width = 12, height = 6, dpi = 300)
ggsave(file = "plots/NEW_FC_snATACseq_snRNAseq_UMAP_combined.svg", plot = p5, width = 12, height = 6)

# sbatch --mem=200g --cpus-per-task=4 --gres=lscratch:100 --time=8:00:00
# or sinteractive









