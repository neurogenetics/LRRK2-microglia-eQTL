#!/bin/env Rscript

# Load necessary packages
library(Signac)  #confirm v1.1.1
packageVersion('Signac')

library(Seurat)  #confirm v4.0.0
packageVersion('Seurat')

library(harmony)
packageVersion('harmony')

library(GenomicRanges)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)

library(patchwork)
library(future)
library(tidyverse)
library(svglite)

library(future)
library(future.apply)
plan("multiprocess", workers = 10)
options(future.globals.maxSize = 50000 * 1024^2)
set.seed(1234)


# Get gene annotations for hg38
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotation) <- "UCSC"
genome(annotation) <- "hg38"

# Define a function to load each dataset and create Seurat objects
create_obj <- function(dir) {
  count.path <- list.files(path = dir, pattern = "*filtered_peak_bc_matrix.h5", full.names = TRUE)
  frag.path <- list.files(path = dir, pattern = "*fragments.tsv.gz", full.names = TRUE)[1]
  counts <- Read10X_h5(count.path)
  md.path <- list.files(path = dir, pattern = "*singlecell.csv", full.names = TRUE)
  md <- read.table(file = md.path, stringsAsFactors = FALSE, sep = ",", header = TRUE, row.names = 1)
  chrom.assay <- CreateChromatinAssay(counts = counts, sep = c(":", "-"), fragments = frag.path, annotation = annotation)
  obj <- CreateSeuratObject(counts = chrom.assay, assay = "peaks", meta.data = md)
  return(obj)
}


# Create Seurat objects, add additional metadata
S794 <- create_obj("data/UMARY794_hg38/outs/")
S794$dataset <- 'S794'
S794$rs76904798 <- 'CC'
S794$age <- 51
S794$sex <- 'Female'
S794$batch <- 1
S794

S1230 <- create_obj("data/UMARY1230_hg38/outs/")
S1230$dataset <- 'S1230'
S1230$rs76904798 <- 'TT'
S1230$age <- 16
S1230$sex <- 'Female'
S1230$batch <- 1
S1230

S1274 <- create_obj("data/UMARY1274_hg38/outs/")
S1274$dataset <- 'S1274'
S1274$rs76904798 <- 'CC'
S1274$age <- 49
S1274$sex <- 'Male'
S1274$batch <- 1
S1274

S5079 <- create_obj("data/UMARY5079_hg38/outs/")
S5079$dataset <- 'S5079'
S5079$rs76904798 <- 'TT'
S5079$age <- 33
S5079$sex <- 'Male'
S5079$batch <- 1
S5079

S630 <- create_obj("data/FC_ATAC_630_hg38/outs/")
S630$dataset <- 'S630'
S630$rs76904798 <- 'CC'
S630$age <- 19
S630$sex <- 'Male'
S630$batch <- 2
S630

S1135 <- create_obj("data/FC_ATAC_1135_hg38/outs/")
S1135$dataset <- 'S1135'
S1135$rs76904798 <- 'CT'
S1135$age <- 42
S1135$sex <- 'Male'
S1135$batch <- 2
S1135

S1209 <- create_obj("data/FC_ATAC_1209_hg38/outs/")
S1209$dataset <- 'S1209'
S1209$rs76904798 <- 'CT'
S1209$age <- 39
S1209$sex <- 'Female'
S1209$batch <- 2
S1209

S4022 <- create_obj("data/FC_ATAC_4022_hg38/outs/")
S4022$dataset <- 'S4022'
S4022$rs76904798 <- 'TT'
S4022$age <- 57
S4022$sex <- 'Female'
S4022$batch <- 2
S4022

S4924 <- create_obj("data/FC_ATAC_4924_hg38/outs/")
S4924$dataset <- 'S4924'
S4924$rs76904798 <- 'TT'
S4924$age <- 48
S4924$sex <- 'Male'
S4924$batch <- 3
S4924

S4724 <- create_obj("data/FC_ATAC_4724_hg38/outs/")
S4724$dataset <- 'S4724'
S4724$rs76904798 <- 'CT'
S4724$age <- 16
S4724$sex <- 'Female'
S4724$batch <- 3
S4724

S5123 <- create_obj("data/FC_ATAC_5123_hg38/outs/")
S5123$dataset <- 'S5123'
S5123$rs76904798 <- 'CT'
S5123$age <- 61
S5123$sex <- 'Male'
S5123$batch <- 3
S5123

S1363 <- create_obj("data/FC_ATAC_1363_hg38/outs/")
S1363$dataset <- 'S1363'
S1363$rs76904798 <- 'CC'
S1363$age <- 40
S1363$sex <- 'Female'
S1363$batch <- 3
S1363


# Perform basic QC filtering on each object
S794 <- NucleosomeSignal(S794)
S794 <- TSSEnrichment(S794)
S794$pct_reads_in_peaks <- S794$peak_region_fragments / S794$passed_filters * 100
S794$blacklist_ratio <- S794$blacklist_region_fragments / S794$peak_region_fragments

S1274 <- NucleosomeSignal(S1274)
S1274 <- TSSEnrichment(S1274)
S1274$pct_reads_in_peaks <- S1274$peak_region_fragments / S1274$passed_filters * 100
S1274$blacklist_ratio <- S1274$blacklist_region_fragments / S1274$peak_region_fragments

S1230 <- NucleosomeSignal(S1230)
S1230 <- TSSEnrichment(S1230)
S1230$pct_reads_in_peaks <- S1230$peak_region_fragments / S1230$passed_filters * 100
S1230$blacklist_ratio <- S1230$blacklist_region_fragments / S1230$peak_region_fragments

S5079 <- NucleosomeSignal(S5079)
S5079 <- TSSEnrichment(S5079)
S5079$pct_reads_in_peaks <- S5079$peak_region_fragments / S5079$passed_filters * 100
S5079$blacklist_ratio <- S5079$blacklist_region_fragments / S5079$peak_region_fragments

S4022 <- NucleosomeSignal(S4022)
S4022 <- TSSEnrichment(S4022)
S4022$pct_reads_in_peaks <- S4022$peak_region_fragments / S4022$passed_filters * 100
S4022$blacklist_ratio <- S4022$blacklist_region_fragments / S4022$peak_region_fragments

S4924 <- NucleosomeSignal(S4924)
S4924 <- TSSEnrichment(S4924)
S4924$pct_reads_in_peaks <- S4924$peak_region_fragments / S4924$passed_filters * 100
S4924$blacklist_ratio <- S4924$blacklist_region_fragments / S4924$peak_region_fragments

S1135 <- NucleosomeSignal(S1135)
S1135 <- TSSEnrichment(S1135)
S1135$pct_reads_in_peaks <- S1135$peak_region_fragments / S1135$passed_filters * 100
S1135$blacklist_ratio <- S1135$blacklist_region_fragments / S1135$peak_region_fragments

S1209 <- NucleosomeSignal(S1209)
S1209 <- TSSEnrichment(S1209)
S1209$pct_reads_in_peaks <- S1209$peak_region_fragments / S1209$passed_filters * 100
S1209$blacklist_ratio <- S1209$blacklist_region_fragments / S1209$peak_region_fragments

S4724 <- NucleosomeSignal(S4724)
S4724 <- TSSEnrichment(S4724)
S4724$pct_reads_in_peaks <- S4724$peak_region_fragments / S4724$passed_filters * 100
S4724$blacklist_ratio <- S4724$blacklist_region_fragments / S4724$peak_region_fragments

S5123 <- NucleosomeSignal(S5123)
S5123 <- TSSEnrichment(S5123)
S5123$pct_reads_in_peaks <- S5123$peak_region_fragments / S5123$passed_filters * 100
S5123$blacklist_ratio <- S5123$blacklist_region_fragments / S5123$peak_region_fragments

S630 <- NucleosomeSignal(S630)
S630 <- TSSEnrichment(S630)
S630$pct_reads_in_peaks <- S630$peak_region_fragments / S630$passed_filters * 100
S630$blacklist_ratio <- S630$blacklist_region_fragments / S630$peak_region_fragments

S1363 <- NucleosomeSignal(S1363)
S1363 <- TSSEnrichment(S1363)
S1363$pct_reads_in_peaks <- S1363$peak_region_fragments / S1363$passed_filters * 100
S1363$blacklist_ratio <- S1363$blacklist_region_fragments / S1363$peak_region_fragments


# Filter out low quality nuclei
S794 <- subset(
  x = S794,
  subset = peak_region_fragments > 2500 &
    peak_region_fragments < 25000 &
    pct_reads_in_peaks > 15 &
    blacklist_ratio < 0.05 &
    nucleosome_signal < 3 &
    TSS.enrichment > 2
)
S794

S1274 <- subset(
  x = S1274,
  subset = peak_region_fragments > 2500 &
    peak_region_fragments < 25000 &
    pct_reads_in_peaks > 15 &
    blacklist_ratio < 0.05 &
    nucleosome_signal < 3 &
    TSS.enrichment > 2
)
S1274

S630 <- subset(
  x = S630,
  subset = peak_region_fragments > 2500 &
    peak_region_fragments < 25000 &
    pct_reads_in_peaks > 15 &
    blacklist_ratio < 0.05 &
    nucleosome_signal < 3 &
    TSS.enrichment > 2
)
S630

S1363 <- subset(
  x = S1363,
  subset = peak_region_fragments > 2500 &
    peak_region_fragments < 25000 &
    pct_reads_in_peaks > 15 &
    blacklist_ratio < 0.05 &
    nucleosome_signal < 3 &
    TSS.enrichment > 2
)
S1363

S4022 <- subset(
  x = S4022,
  subset = peak_region_fragments > 2500 &
    peak_region_fragments < 25000 &
    pct_reads_in_peaks > 15 &
    blacklist_ratio < 0.05 &
    nucleosome_signal < 3 &
    TSS.enrichment > 2
)
S4022

S1230 <- subset(
  x = S1230,
  subset = peak_region_fragments > 2500 &
    peak_region_fragments < 25000 &
    pct_reads_in_peaks > 15 &
    blacklist_ratio < 0.05 &
    nucleosome_signal < 3 &
    TSS.enrichment > 2
)
S1230

S4924 <- subset(
  x = S4924,
  subset = peak_region_fragments > 2500 &
    peak_region_fragments < 25000 &
    pct_reads_in_peaks > 15 &
    blacklist_ratio < 0.05 &
    nucleosome_signal < 3 &
    TSS.enrichment > 2
)
S4924

S5079 <- subset(
  x = S5079,
  subset = peak_region_fragments > 2500 &
    peak_region_fragments < 25000 &
    pct_reads_in_peaks > 15 &
    blacklist_ratio < 0.05 &
    nucleosome_signal < 3 &
    TSS.enrichment > 2
)
S5079

S1135 <- subset(
  x = S1135,
  subset = peak_region_fragments > 2500 &
    peak_region_fragments < 25000 &
    pct_reads_in_peaks > 15 &
    blacklist_ratio < 0.05 &
    nucleosome_signal < 3 &
    TSS.enrichment > 2
)
S1135

S1209 <- subset(
  x = S1209,
  subset = peak_region_fragments > 2500 &
    peak_region_fragments < 25000 &
    pct_reads_in_peaks > 15 &
    blacklist_ratio < 0.05 &
    nucleosome_signal < 3 &
    TSS.enrichment > 2
)
S1209

S4724 <- subset(
  x = S4724,
  subset = peak_region_fragments > 2500 &
    peak_region_fragments < 25000 &
    pct_reads_in_peaks > 15 &
    blacklist_ratio < 0.05 &
    nucleosome_signal < 3 &
    TSS.enrichment > 2
)
S4724

S5123 <- subset(
  x = S5123,
  subset = peak_region_fragments > 2500 &
    peak_region_fragments < 25000 &
    pct_reads_in_peaks > 15 &
    blacklist_ratio < 0.05 &
    nucleosome_signal < 3 &
    TSS.enrichment > 2
)
S5123


# Call peaks with macs2
macs2.path <- "/usr/local/apps/macs/2.2.7.1/bin/macs2"

## S794
peaks <- CallPeaks(S794, macs2.path = macs2.path)
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg38_unified, invert = TRUE)

macs2_counts <- FeatureMatrix(
  fragments = Fragments(S794),
  features = peaks,
  cells = colnames(S794)
)

S794[["peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  fragments = "data/UMARY794_hg38/outs/fragments.tsv.gz",
  annotation = annotation
)

## S1274
peaks <- CallPeaks(S1274, macs2.path = macs2.path)
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg38_unified, invert = TRUE)

macs2_counts <- FeatureMatrix(
  fragments = Fragments(S1274),
  features = peaks,
  cells = colnames(S1274)
)

S1274[["peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  fragments = "data/UMARY1274_hg38/outs/fragments.tsv.gz",
  annotation = annotation
)

## S1230
peaks <- CallPeaks(S1230, macs2.path = macs2.path)
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg38_unified, invert = TRUE)

macs2_counts <- FeatureMatrix(
  fragments = Fragments(S1230),
  features = peaks,
  cells = colnames(S1230)
)

S1230[["peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  fragments = "data/UMARY1230_hg38/outs/fragments.tsv.gz",
  annotation = annotation
)

## S5079
peaks <- CallPeaks(S5079, macs2.path = macs2.path)
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg38_unified, invert = TRUE)

macs2_counts <- FeatureMatrix(
  fragments = Fragments(S5079),
  features = peaks,
  cells = colnames(S5079)
)

S5079[["peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  fragments = "data/UMARY5079_hg38/outs/fragments.tsv.gz",
  annotation = annotation
)

## S4022
peaks <- CallPeaks(S4022, macs2.path = macs2.path)
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg38_unified, invert = TRUE)

macs2_counts <- FeatureMatrix(
  fragments = Fragments(S4022),
  features = peaks,
  cells = colnames(S4022)
)

S4022[["peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  fragments = "data/FC_ATAC_4022_hg38/outs/fragments.tsv.gz",
  annotation = annotation
)

## S4924
peaks <- CallPeaks(S4924, macs2.path = macs2.path)
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg38_unified, invert = TRUE)

macs2_counts <- FeatureMatrix(
  fragments = Fragments(S4924),
  features = peaks,
  cells = colnames(S4924)
)

S4924[["peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  fragments = "data/FC_ATAC_4924_hg38/outs/fragments.tsv.gz",
  annotation = annotation
)

## S1135
peaks <- CallPeaks(S1135, macs2.path = macs2.path)
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg38_unified, invert = TRUE)

macs2_counts <- FeatureMatrix(
  fragments = Fragments(S1135),
  features = peaks,
  cells = colnames(S1135)
)

S1135[["peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  fragments = "data/FC_ATAC_1135_hg38/outs/fragments.tsv.gz",
  annotation = annotation
)

## S1209
peaks <- CallPeaks(S1209, macs2.path = macs2.path)
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg38_unified, invert = TRUE)

macs2_counts <- FeatureMatrix(
  fragments = Fragments(S1209),
  features = peaks,
  cells = colnames(S1209)
)

S1209[["peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  fragments = "data/FC_ATAC_1209_hg38/outs/fragments.tsv.gz",
  annotation = annotation
)

## S4724
peaks <- CallPeaks(S4724, macs2.path = macs2.path)
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg38_unified, invert = TRUE)

macs2_counts <- FeatureMatrix(
  fragments = Fragments(S4724),
  features = peaks,
  cells = colnames(S4724)
)

S4724[["peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  fragments = "data/FC_ATAC_4724_hg38/outs/fragments.tsv.gz",
  annotation = annotation
)

## S5123
peaks <- CallPeaks(S5123, macs2.path = macs2.path)
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg38_unified, invert = TRUE)

macs2_counts <- FeatureMatrix(
  fragments = Fragments(S5123),
  features = peaks,
  cells = colnames(S5123)
)

S5123[["peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  fragments = "data/FC_ATAC_5123_hg38/outs/fragments.tsv.gz",
  annotation = annotation
)

## S630
peaks <- CallPeaks(S630, macs2.path = macs2.path)
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg38_unified, invert = TRUE)

macs2_counts <- FeatureMatrix(
  fragments = Fragments(S630),
  features = peaks,
  cells = colnames(S630)
)

S630[["peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  fragments = "data/FC_ATAC_630_hg38/outs/fragments.tsv.gz",
  annotation = annotation
)

## S1363
peaks <- CallPeaks(S1363, macs2.path = macs2.path)
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg38_unified, invert = TRUE)

macs2_counts <- FeatureMatrix(
  fragments = Fragments(S1363),
  features = peaks,
  cells = colnames(S1363)
)

S1363[["peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  fragments = "data/FC_ATAC_1363_hg38/outs/fragments.tsv.gz",
  annotation = annotation
)


# Perform dimensional reduction on each object
DefaultAssay(S794) <- "peaks"
DefaultAssay(S1274) <- "peaks"
DefaultAssay(S1230) <- "peaks"
DefaultAssay(S5079) <- "peaks"
DefaultAssay(S4022) <- "peaks"
DefaultAssay(S4924) <- "peaks"
DefaultAssay(S1135) <- "peaks"
DefaultAssay(S1209) <- "peaks"
DefaultAssay(S4724) <- "peaks"
DefaultAssay(S5123) <- "peaks"
DefaultAssay(S630) <- "peaks"
DefaultAssay(S1363) <- "peaks"

S794 <- FindTopFeatures(S794, min.cutoff = 'q0')
S794 <- RunTFIDF(S794)
S794 <- RunSVD(S794)
S794 <- RunUMAP(S794, reduction = 'lsi', dims = 2:30)

S1274 <- FindTopFeatures(S1274, min.cutoff = 'q0')
S1274 <- RunTFIDF(S1274)
S1274 <- RunSVD(S1274)
S1274 <- RunUMAP(S1274, reduction = 'lsi', dims = 2:30)

S1230 <- FindTopFeatures(S1230, min.cutoff = 'q0')
S1230 <- RunTFIDF(S1230)
S1230 <- RunSVD(S1230)
S1230 <- RunUMAP(S1230, reduction = 'lsi', dims = 2:30)

S5079 <- FindTopFeatures(S5079, min.cutoff = 'q0')
S5079 <- RunTFIDF(S5079)
S5079 <- RunSVD(S5079)
S5079 <- RunUMAP(S5079, reduction = 'lsi', dims = 2:30)

S4022 <- FindTopFeatures(S4022, min.cutoff = 'q0')
S4022 <- RunTFIDF(S4022)
S4022 <- RunSVD(S4022)
S4022 <- RunUMAP(S4022, reduction = 'lsi', dims = 2:30)

S4924 <- FindTopFeatures(S4924, min.cutoff = 'q0')
S4924 <- RunTFIDF(S4924)
S4924 <- RunSVD(S4924)
S4924 <- RunUMAP(S4924, reduction = 'lsi', dims = 2:30)

S1135 <- FindTopFeatures(S1135, min.cutoff = 'q0')
S1135 <- RunTFIDF(S1135)
S1135 <- RunSVD(S1135)
S1135 <- RunUMAP(S1135, reduction = 'lsi', dims = 2:30)

S1209 <- FindTopFeatures(S1209, min.cutoff = 'q0')
S1209 <- RunTFIDF(S1209)
S1209 <- RunSVD(S1209)
S1209 <- RunUMAP(S1209, reduction = 'lsi', dims = 2:30)

S4724 <- FindTopFeatures(S4724, min.cutoff = 'q0')
S4724 <- RunTFIDF(S4724)
S4724 <- RunSVD(S4724)
S4724 <- RunUMAP(S4724, reduction = 'lsi', dims = 2:30)

S5123 <- FindTopFeatures(S5123, min.cutoff = 'q0')
S5123 <- RunTFIDF(S5123)
S5123 <- RunSVD(S5123)
S5123 <- RunUMAP(S5123, reduction = 'lsi', dims = 2:30)

S630 <- FindTopFeatures(S630, min.cutoff = 'q0')
S630 <- RunTFIDF(S630)
S630 <- RunSVD(S630)
S630 <- RunUMAP(S630, reduction = 'lsi', dims = 2:30)

S1363 <- FindTopFeatures(S1363, min.cutoff = 'q0')
S1363 <- RunTFIDF(S1363)
S1363 <- RunSVD(S1363)
S1363 <- RunUMAP(S1363, reduction = 'lsi', dims = 2:30)


# Create a consensus peakset
FC.atac.list <- list(S794, S1230, S1274, S5079, S4022, S4924, S1135, S1209, S4724, S5123, S630, S1363)
combined.peaks <- UnifyPeaks(object.list = FC.atac.list, mode = "reduce")

## Filter
peakwidths <- width(combined.peaks)
combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20]
combined.peaks

# Count fragments per cell overlapping the set of peaks in each dataset and add to object
## S794
## S794
S794.cp <- FeatureMatrix(
  fragments = Fragments(S794),
  features = combined.peaks,
  cells = colnames(S794)
)

S794[['peaks']] <- CreateChromatinAssay(
  counts = S794.cp,
  min.cells = 1,
  ranges = combined.peaks,
  annotation = annotation
)

## S1230
S1230.cp <- FeatureMatrix(
  fragments = Fragments(S1230),
  features = combined.peaks,
  cells = colnames(S1230)
)

S1230[['peaks']] <- CreateChromatinAssay(
  counts = S1230.cp,
  min.cells = 1,
  ranges = combined.peaks,
  annotation = annotation
)

## S1274
S1274.cp <- FeatureMatrix(
  fragments = Fragments(S1274),
  features = combined.peaks,
  cells = colnames(S1274)
)

S1274[['peaks']] <- CreateChromatinAssay(
  counts = S1274.cp,
  min.cells = 1,
  ranges = combined.peaks,
  annotation = annotation
)

## S5079
S5079.cp <- FeatureMatrix(
  fragments = Fragments(S5079),
  features = combined.peaks,
  cells = colnames(S5079)
)

S5079[['peaks']] <- CreateChromatinAssay(
  counts = S5079.cp,
  min.cells = 1,
  ranges = combined.peaks,
  annotation = annotation
)


## S4022
S4022.cp <- FeatureMatrix(
  fragments = Fragments(S4022),
  features = combined.peaks,
  cells = colnames(S4022)
)

S4022[['peaks']] <- CreateChromatinAssay(
  counts = S4022.cp,
  min.cells = 1,
  ranges = combined.peaks,
  annotation = annotation
)


## S4924
S4924.cp <- FeatureMatrix(
  fragments = Fragments(S4924),
  features = combined.peaks,
  cells = colnames(S4924)
)

S4924[['peaks']] <- CreateChromatinAssay(
  counts = S4924.cp,
  min.cells = 1,
  ranges = combined.peaks,
  annotation = annotation
)


## S1135
S1135.cp <- FeatureMatrix(
  fragments = Fragments(S1135),
  features = combined.peaks,
  cells = colnames(S1135)
)

S1135[['peaks']] <- CreateChromatinAssay(
  counts = S1135.cp,
  min.cells = 1,
  ranges = combined.peaks,
  annotation = annotation
)


## S1209
S1209.cp <- FeatureMatrix(
  fragments = Fragments(S1209),
  features = combined.peaks,
  cells = colnames(S1209)
)

S1209[['peaks']] <- CreateChromatinAssay(
  counts = S1209.cp,
  min.cells = 1,
  ranges = combined.peaks,
  annotation = annotation
)


## S4724
S4724.cp <- FeatureMatrix(
  fragments = Fragments(S4724),
  features = combined.peaks,
  cells = colnames(S4724)
)

S4724[['peaks']] <- CreateChromatinAssay(
  counts = S4724.cp,
  min.cells = 1,
  ranges = combined.peaks,
  annotation = annotation
)

## S5123
S5123.cp <- FeatureMatrix(
  fragments = Fragments(S5123),
  features = combined.peaks,
  cells = colnames(S5123)
)

S5123[['peaks']] <- CreateChromatinAssay(
  counts = S5123.cp,
  min.cells = 1,
  ranges = combined.peaks,
  annotation = annotation
)


## S630
S630.cp <- FeatureMatrix(
  fragments = Fragments(S630),
  features = combined.peaks,
  cells = colnames(S630)
)

S630[['peaks']] <- CreateChromatinAssay(
  counts = S630.cp,
  min.cells = 1,
  ranges = combined.peaks,
  annotation = annotation
)


## S1363
S1363.cp <- FeatureMatrix(
  fragments = Fragments(S1363),
  features = combined.peaks,
  cells = colnames(S1363)
)

S1363[['peaks']] <- CreateChromatinAssay(
  counts = S1363.cp,
  min.cells = 1,
  ranges = combined.peaks,
  annotation = annotation
)

# Merge
FC.atac.list <- list(S1230, S1274, S5079, S4022, S4924, S1135, S1209, S4724, S5123, S630, S1363)
FC.atac <- merge(
  x = S794, 
  y = FC.atac.list, 
  add.cell.ids = c("S794", "S1230", "S1274", "S5079", "S4022", "S4924", "S1135", "S1209", "S4724", "S5123", "S630", "S1363")
  )

FC.atac

FC.atac[["peaks"]]

# Add merged fragment file to the chromatin assay
fragments <- CreateFragmentObject(
  path = "data/Set12_FC_fragments.tsv.gz",
  cells = colnames(FC.atac), 
  validate.fragments = TRUE
)
fragments

Fragments(FC.atac) <- fragments

# Normalize merged object
DefaultAssay(FC.atac) <- "peaks" 

FC.atac <- RunTFIDF(FC.atac)
FC.atac <- FindTopFeatures(FC.atac, min.cutoff = 20)
FC.atac <- RunSVD(FC.atac)
FC.atac <- RunUMAP(FC.atac, dims = 2:50, reduction = 'lsi')

# Run harmony
FC.atac <- RunHarmony(
  object = FC.atac,
  group.by.vars = 'dataset',
  reduction = 'lsi',
  assay.use = 'peaks',
  project.dim = FALSE
)

# re-compute the UMAP using corrected LSI embeddings
FC.atac <- RunUMAP(FC.atac, dims = 2:30, reduction = 'harmony')

# Check out harmony integrated object
p1 <- VlnPlot(
  object = FC.atac,
  features = c('pct_reads_in_peaks', 'peak_region_fragments', 'nucleosome_signal', 'blacklist_ratio'),
  group.by = 'dataset',
  pt.size = 0,
  ncol = 4) + NoLegend()
ggsave(file = "plots/FC_ATAC_QC_Vln_harmony_int.png", plot = p1, width = 12, height = 4, dpi = 300)
ggsave(file = "plots/FC_ATAC_QC_Vln_harmony_int.svg", plot = p1, width = 12, height = 4)

# Cluster
FC.atac <- FindNeighbors(FC.atac, reduction = 'harmony', dims = 2:30)
FC.atac <- FindClusters(FC.atac, verbose = FALSE, algorithm = 3, resolution = 0.4)

# Visualize
p2 <- DimPlot(FC.atac, label = TRUE) + NoLegend()
p3 <- DimPlot(FC.atac, group.by = 'dataset', pt.size = 0.1)
p4 <- p2 + p3
ggsave(file = "plots/FC_ATAC_UMAP_harmony_int.png", plot = p4, width = 12.5, height = 6, dpi = 300)

# Create gene activity matrix and normalize
gene.activities <- GeneActivity(FC.atac)
FC.atac[['RNA']] <- CreateAssayObject(counts = gene.activities)

DefaultAssay(FC.atac) <- "RNA"
FC.atac <- SCTransform(FC.atac)
FC.atac <- RunPCA(FC.atac)

save(FC.atac, file = "output/FC_ATAC_clustered_harmony_int.Rdata")

# Look at cell type cluster markers
DefaultAssay(FC.atac) <- "SCT"

FC.atac.markers <- FindAllMarkers(object = FC.atac, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(FC.atac.markers, file = "output/FC_ATAC_markers_harmony_int.csv")

p5 <- FeaturePlot(FC.atac, reduction = "umap", c("PDGFRA", "GFAP", "CSF1R", "MOBP", "CLDN5", "COLEC12"), cols = c("#CCFFFF", "lightgrey", "#FF0000"), label = TRUE, order = T, ncol = 3, label.size = 4)
ggsave("plots/FC_ATAC_markers_glia_Feature_harmony_int.png", plot = p5, height = 8, width = 12)

p6 <- VlnPlot(FC.atac, c("PDGFRA", "GFAP", "CSF1R", "MOBP", "CLDN5", "COLEC12"), pt.size = 0, ncol = 2)
ggsave("plots/FC_ATAC_markers_glia_Vln_harmony_int.png", plot = p6, height = 8, width = 12)

p7 <- FeaturePlot(FC.atac, reduction = "umap", c("SLC17A7", "SLC17A6", "THEMIS", "FEZF2", "LINC00507", "RORB"), cols = c("#CCFFFF", "lightgrey", "#FF0000"), label = TRUE, order = T, ncol = 3, label.size = 4)
ggsave("plots/FC_ATAC_markers_ExN_Feature_harmony_int.png", plot = p7, height = 8, width = 12)

p8 <- VlnPlot(FC.atac, c("SLC17A7", "SLC17A6", "THEMIS", "FEZF2", "LINC00507", "RORB"), pt.size = 0, ncol = 2)
ggsave("plots/FC_ATAC_markers_ExN_Vln_harmony_int.png", plot = p8, height = 8, width = 12)

p9 <- FeaturePlot(FC.atac, reduction = "umap", c("GAD1", "GAD2", "VIP", "SST", "PVALB", "ADARB2"), cols = c("#CCFFFF", "lightgrey", "#FF0000"), label = TRUE, order = T, ncol = 3, label.size = 4)
ggsave("plots/FC_ATAC_markers_InN_Feature_harmony_int.png", plot = p9, height = 8, width = 12)

p10 <- VlnPlot(FC.atac, c("GAD1", "GAD2", "VIP", "SST", "PVALB", "ADARB2"), pt.size = 0, ncol = 2)
ggsave("plots/FC_ATAC_markers_InN_Vln_harmony_int.png", plot = p10, height = 8, width = 12)



# sbatch --partition=largemem --mem=375g --cpus-per-task=10 --gres=lscratch:200 --time=36:00:00

