#!/bin/env Rscript

# Integration of three substantia nigra snMultiome datasets based on ATAC data

# Load necessary packages
library(Signac)  #confirm v1.1.1
packageVersion('Signac')

library(Seurat)  #confirm v4.0.3
packageVersion('Seurat')

library(GenomicRanges)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(patchwork)
library(future)
library(tidyverse)

plan("multiprocess", workers = 8)
options(future.globals.maxSize = 50000 * 1024^2)
set.seed(1234)

# load the RNA and ATAC data
SN2960.counts <- Read10X_h5("./SN_arc/SN2960_outs/filtered_feature_bc_matrix.h5")
SN2960.fragpath <- "./SN_arc/SN2960_outs/atac_fragments.tsv.gz"
SN2960.md <- read.table("./SN_arc/SN2960_outs/per_barcode_metrics.csv", stringsAsFactors = FALSE, sep = ",",
                        header = TRUE, row.names = 1)

SN2972.counts <- Read10X_h5("./SN_arc/SN2972_outs/filtered_feature_bc_matrix.h5")
SN2972.fragpath <- "./SN_arc/SN2972_outs/atac_fragments.tsv.gz"
SN2972.md <- read.table("./SN_arc/SN2972_outs/per_barcode_metrics.csv", stringsAsFactors = FALSE, sep = ",",
                        header = TRUE, row.names = 1)

SN2942.counts <- Read10X_h5("./SN_arc/SN2942_outs/filtered_feature_bc_matrix.h5")
SN2942.fragpath <- "./SN_arc/SN2942_outs/atac_fragments.tsv.gz"
SN2942.md <- read.table("./SN_arc/SN2942_outs/per_barcode_metrics.csv", stringsAsFactors = FALSE, sep = ",",
                        header = TRUE, row.names = 1)

# get gene annotations for hg38
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotation) <- "UCSC"
genome(annotation) <- "hg38"

# Create Seurat objects
SN2960 <- CreateSeuratObject(
  counts = SN2960.counts$`Gene Expression`,
  assay = "RNA",
  meta.data = SN2960.md,
  project = "SN2960"
)

SN2960[["ATAC"]] <- CreateChromatinAssay(
  counts = SN2960.counts$Peaks,
  sep = c(":", "-"),
  fragments = SN2960.fragpath,
  annotation = annotation
)
SN2960


SN2972 <- CreateSeuratObject(
  counts = SN2972.counts$`Gene Expression`,
  assay = "RNA",
  meta.data = SN2972.md,
  project = "SN2972"
)

SN2972[["ATAC"]] <- CreateChromatinAssay(
  counts = SN2972.counts$Peaks,
  sep = c(":", "-"),
  fragments = SN2972.fragpath,
  annotation = annotation
)
SN2972


SN2942 <- CreateSeuratObject(
  counts = SN2942.counts$`Gene Expression`,
  assay = "RNA",
  meta.data = SN2942.md,
  project = "SN2942"
)

SN2942[["ATAC"]] <- CreateChromatinAssay(
  counts = SN2942.counts$Peaks,
  sep = c(":", "-"),
  fragments = SN2942.fragpath,
  annotation = annotation
)
SN2942

# Perform basic QC filtering based on ATAC data
DefaultAssay(SN2960) <- "ATAC"
DefaultAssay(SN2972) <- "ATAC"
DefaultAssay(SN2942) <- "ATAC"

## SN2960
SN2960 <- NucleosomeSignal(SN2960)
SN2960 <- TSSEnrichment(SN2960)

SN2960 <- subset(
  x = SN2960,
  subset = nCount_ATAC < 100000 &
    nCount_ATAC > 1000 &
    nucleosome_signal < 2 &
    TSS.enrichment > 1
)
SN2960

## SN2972
SN2972 <- NucleosomeSignal(SN2972)
SN2972 <- TSSEnrichment(SN2972)

SN2972 <- subset(
  x = SN2972,
  subset = nCount_ATAC < 100000 &
    nCount_ATAC > 1000 &
    nucleosome_signal < 2 &
    TSS.enrichment > 1
)
SN2972

## SN2942
SN2942 <- NucleosomeSignal(SN2942)
SN2942 <- TSSEnrichment(SN2942)

SN2942 <- subset(
  x = SN2942,
  subset = nCount_ATAC < 100000 &
    nCount_ATAC > 1000 &
    nucleosome_signal < 2 &
    TSS.enrichment > 1
)
SN2942

# Call peaks with macs2
## SN2960
macs2.path <- "/usr/local/apps/macs/2.2.7.1/bin/macs2"

peaks <- CallPeaks(SN2960, macs2.path = macs2.path)
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg38_unified, invert = TRUE)

macs2_counts <- FeatureMatrix(
  fragments = Fragments(SN2960),
  features = peaks,
  cells = colnames(SN2960)
)

SN2960[["peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  fragments = SN2960.fragpath,
  annotation = annotation
)

## SN2972
peaks <- CallPeaks(SN2972, macs2.path = macs2.path)
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg38_unified, invert = TRUE)

macs2_counts <- FeatureMatrix(
  fragments = Fragments(SN2972),
  features = peaks,
  cells = colnames(SN2972)
)

SN2972[["peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  fragments = SN2972.fragpath,
  annotation = annotation
)

## SN2942
peaks <- CallPeaks(SN2942, macs2.path = macs2.path)
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg38_unified, invert = TRUE)

macs2_counts <- FeatureMatrix(
  fragments = Fragments(SN2942),
  features = peaks,
  cells = colnames(SN2942)
)

SN2942[["peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  fragments = SN2942.fragpath,
  annotation = annotation
)

# Create a consensus peak set
SN.list <- list(SN2960, SN2972, SN2942)
combined.peaks <- UnifyPeaks(object.list = SN.list, mode = "reduce")

## Filter
peakwidths <- width(combined.peaks)
combined.peaks <- combined.peaks[peakwidths < 10000 & peakwidths > 20]
combined.peaks

# Count fragments per cell overlapping the set of peaks in each dataset and add to object
## SN2960
SN2960.cp <- FeatureMatrix(
  fragments = Fragments(SN2960),
  features = combined.peaks,
  cells = colnames(SN2960)
)

SN2960[['peaks']] <- CreateChromatinAssay(
  counts = SN2960.cp,
  fragments = SN2960.fragpath,
  ranges = combined.peaks,
  annotation = annotation
)

## SN2972
SN2972.cp <- FeatureMatrix(
  fragments = Fragments(SN2972),
  features = combined.peaks,
  cells = colnames(SN2972)
)

SN2972[['peaks']] <- CreateChromatinAssay(
  counts = SN2972.cp,
  fragments = SN2972.fragpath,
  ranges = combined.peaks,
  annotation = annotation
)

## SN2942
SN2942.cp <- FeatureMatrix(
  fragments = Fragments(SN2942),
  features = combined.peaks,
  cells = colnames(SN2942)
)

SN2942[['peaks']] <- CreateChromatinAssay(
  counts = SN2942.cp,
  fragments = SN2942.fragpath,
  ranges = combined.peaks,
  annotation = annotation
)

# Perform dimensional reduction on each object
DefaultAssay(SN2960) <- "peaks"
DefaultAssay(SN2972) <- "peaks"
DefaultAssay(SN2942) <- "peaks"

SN2960 <- FindTopFeatures(SN2960, min.cutoff = 10)
SN2960 <- RunTFIDF(SN2960)
SN2960 <- RunSVD(SN2960)

SN2972 <- FindTopFeatures(SN2972, min.cutoff = 10)
SN2972 <- RunTFIDF(SN2972)
SN2972 <- RunSVD(SN2972)

SN2942 <- FindTopFeatures(SN2942, min.cutoff = 10)
SN2942 <- RunTFIDF(SN2942)
SN2942 <- RunSVD(SN2942)


# add dataset-identifying metadata
SN2960$dataset <- "SN2960"
SN2972$dataset <- "SN2972"
SN2942$dataset <- "SN2942"

# merge
SN.combined <- merge(SN2960, list(SN2972, SN2942))

# process the combined dataset
SN.combined <- FindTopFeatures(SN.combined, min.cutoff = 10)
SN.combined <- RunTFIDF(SN.combined)
SN.combined <- RunSVD(SN.combined)
SN.combined <- RunUMAP(SN.combined, reduction = "lsi", dims = 2:30)
p1 <- DimPlot(SN.combined, group.by = "dataset", pt.size = 0.1)


# find integration anchors
integration.anchors <- FindIntegrationAnchors(
  object.list = list(SN2960, SN2972, SN2942),
  anchor.features = rownames(SN2960),
  reduction = "rlsi",
  dims = 2:30
)

# integrate LSI embeddings
SN.int <- IntegrateEmbeddings(
  anchorset = integration.anchors,
  reductions = SN.combined[["lsi"]],
  new.reduction.name = "integrated_lsi",
  dims.to.integrate = 1:30
)

# create a new UMAP using the integrated embeddings
SN.int <- RunUMAP(SN.int, reduction = "integrated_lsi", dims = 2:30)
p2 <- DimPlot(SN.int, group.by = "dataset", pt.size = 0.1)

# Visualize the difference between merged and Integrated
p3 <- (p1 + ggtitle("Merged")) | (p2 + ggtitle("Integrated"))
ggsave(file = "plots/SN2960_SN2972_SN2942_Signac_merged_vs_int.png", plot = p3, width = 12.5, height = 6)

# Save integrated object
save(SN.int, file = "output/SN2960_SN2972_SN2942_Signac_int.Rdata")



# sbatch --mem=200g --cpus-per-task=8 --gres=lscratch:100 --time=8:00:00


