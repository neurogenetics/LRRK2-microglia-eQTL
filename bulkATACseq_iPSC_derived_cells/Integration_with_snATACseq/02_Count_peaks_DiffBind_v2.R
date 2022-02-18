#!/bin/env Rscript

# Load necessary package
library(DiffBind)
packageVersion("DiffBind")

dba <- dba(sampleSheet = "sampleSheet_DiffBind_Set12_iMGL_iFbN.csv")

# Export raw read counts
dba <- dba.count(dba, score=DBA_SCORE_READS, bUseSummarizeOverlaps=TRUE)   #requires .bam.bai index
allReads <- dba.peakset(dba, bRetrieve=TRUE, DataType = DBA_DATA_FRAME)
write.csv(allReads, "output/Set12_snATACseq_bulkATAC_celltype_counts_summ.csv", row.names=F, quote=F)

#sbatch --mem=250g --time=6:00:00 --gres=lscratch:100 BatchSubmit_CountPeaks.sh


