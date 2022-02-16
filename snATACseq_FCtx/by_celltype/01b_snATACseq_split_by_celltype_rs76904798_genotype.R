#!/bin/env Rscript

# Load necessary packages
library(Signac)  #confirm v1.1.1
packageVersion('Signac')

library(Seurat)  #confirm v4.0.0
packageVersion('Seurat')

library(dplyr)
library(tidyr)
library(stringr)

# Load integrated, named FC.atac object
load("output/FC_ATAC_harmony_int_named.Rdata")
FC.atac

# Split celltypes by genotype at rs76904798 and create new ident
FC.atac$celltype.rs76904798 <- paste(Idents(FC.atac), FC.atac$rs76904798, sep = "_")
Idents(FC.atac) <- "celltype.rs76904798"

# Fetch cell barcodes + celltype classification for cell type of interest, and split by sample ID
named_cell_bcs <- FetchData(object = FC.atac, vars = "celltype.rs76904798")
named_cell_bcs_all <- tibble::rownames_to_column(named_cell_bcs, var = "barcode")
named_cell_bcs_all <- separate(named_cell_bcs_all, col = barcode, into = c("ID", "barcode"), sep = "_")

## For population MGL.4
named_cell_bcs_MGL <- named_cell_bcs_all %>% filter(str_detect(celltype.rs76904798, "MGL.4"))

named_cell_bcs_S794 <- named_cell_bcs_MGL %>% filter(ID == "S794")
named_cell_bcs_S1230 <- named_cell_bcs_MGL %>% filter(ID == "S1230")
named_cell_bcs_S1274 <- named_cell_bcs_MGL %>% filter(ID == "S1274")
named_cell_bcs_S5079 <- named_cell_bcs_MGL %>% filter(ID == "S5079")
named_cell_bcs_S4022 <- named_cell_bcs_MGL %>% filter(ID == "S4022")
named_cell_bcs_S4724 <- named_cell_bcs_MGL %>% filter(ID == "S4724")
named_cell_bcs_S1135 <- named_cell_bcs_MGL %>% filter(ID == "S1135")
named_cell_bcs_S5123 <- named_cell_bcs_MGL %>% filter(ID == "S5123")
named_cell_bcs_S1209 <- named_cell_bcs_MGL %>% filter(ID == "S1209")
named_cell_bcs_S4924 <- named_cell_bcs_MGL %>% filter(ID == "S4924")
named_cell_bcs_S630 <- named_cell_bcs_MGL %>% filter(ID == "S630")
named_cell_bcs_S1363 <- named_cell_bcs_MGL %>% filter(ID == "S1363")


# Remove sample IDs and write to tab-delimited files
write.table(named_cell_bcs_S794[2:3], "output/UMARY794_bcs_celltype_MGL4_CCvsCTvsTT.tsv", sep="\t", col.names=F, row.names=F, quote=F)
write.table(named_cell_bcs_S1230[2:3], "output/UMARY1230_bcs_celltype_MGL4_CCvsCTvsTT.tsv", sep="\t", col.names=F, row.names=F, quote=F)
write.table(named_cell_bcs_S1274[2:3], "output/UMARY1274_bcs_celltype_MGL4_CCvsCTvsTT.tsv", sep="\t", col.names=F, row.names=F, quote=F)
write.table(named_cell_bcs_S5079[2:3], "output/UMARY5079_bcs_celltype_MGL4_CCvsCTvsTT.tsv", sep="\t", col.names=F, row.names=F, quote=F)
write.table(named_cell_bcs_S4022[2:3], "output/FC_ATAC_4022_bcs_celltype_MGL4_CCvsCTvsTT.tsv", sep="\t", col.names=F, row.names=F, quote=F)
write.table(named_cell_bcs_S4724[2:3], "output/FC_ATAC_4724_bcs_celltype_MGL4_CCvsCTvsTT.tsv", sep="\t", col.names=F, row.names=F, quote=F)
write.table(named_cell_bcs_S1135[2:3], "output/FC_ATAC_1135_bcs_celltype_MGL4_CCvsCTvsTT.tsv", sep="\t", col.names=F, row.names=F, quote=F)
write.table(named_cell_bcs_S5123[2:3], "output/FC_ATAC_5123_bcs_celltype_MGL4_CCvsCTvsTT.tsv", sep="\t", col.names=F, row.names=F, quote=F)
write.table(named_cell_bcs_S1209[2:3], "output/FC_ATAC_1209_bcs_celltype_MGL4_CCvsCTvsTT.tsv", sep="\t", col.names=F, row.names=F, quote=F)
write.table(named_cell_bcs_S4924[2:3], "output/FC_ATAC_4924_bcs_celltype_MGL4_CCvsCTvsTT.tsv", sep="\t", col.names=F, row.names=F, quote=F)
write.table(named_cell_bcs_S630[2:3], "output/FC_ATAC_630_bcs_celltype_MGL4_CCvsCTvsTT.tsv", sep="\t", col.names=F, row.names=F, quote=F)
write.table(named_cell_bcs_S1363[2:3], "output/FC_ATAC_1363_bcs_celltype_MGL4_CCvsCTvsTT.tsv", sep="\t", col.names=F, row.names=F, quote=F)

# sbatch --mem=50g --gres=lscratch:100 --time=6:00:00




