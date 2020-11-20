#!/bin/env Rscript

library(DiffBind)

dba <- dba(sampleSheet = "./data/sampleSheet_DiffBind_with_iMGL_iFbN.csv")

# Export raw read counts
dba <- dba.count(dba, score=DBA_SCORE_READS, bUseSummarizeOverlaps=TRUE)     #requires indexed bam
allReads <- dba.peakset(dba, bRetrieve=TRUE, DataType = DBA_DATA_FRAME)
write.csv(allReads, "./output/snATACseq_bulkATAC_celltype_counts_summ.csv", row.names=F, quote=F)

