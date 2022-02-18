#!/bin/bash

# sbatch --cpus-per-task=20 --mem=50g --time=1:00:00 --mail-type=BEGIN,END Sinto_filterbarcodes_p2_v2.sh CELLTYPE

CELLTYPE=$1
echo "$CELLTYPE"

# Change into snATACseq data folder
cd /data/langstonrg/snATACseq/data

# Merge bam files for each celltype
module load samtools

samtools merge -@ 20 ./celltype_bams/Set12_FC_"$CELLTYPE".bam \
	./UMARY794_hg38/outs/split_frags/"$CELLTYPE".bam \
	./UMARY1230_hg38/outs/split_frags/"$CELLTYPE".bam \
	./UMARY1274_hg38/outs/split_frags/"$CELLTYPE".bam \
	./UMARY5079_hg38/outs/split_frags/"$CELLTYPE".bam \
	./FC_ATAC_630_hg38/outs/split_frags/"$CELLTYPE".bam \
	./FC_ATAC_1135_hg38/outs/split_frags/"$CELLTYPE".bam \
	./FC_ATAC_1209_hg38/outs/split_frags/"$CELLTYPE".bam \
	./FC_ATAC_1363_hg38/outs/split_frags/"$CELLTYPE".bam \
	./FC_ATAC_4022_hg38/outs/split_frags/"$CELLTYPE".bam \
	./FC_ATAC_4724_hg38/outs/split_frags/"$CELLTYPE".bam \
	./FC_ATAC_4924_hg38/outs/split_frags/"$CELLTYPE".bam \
	./FC_ATAC_5123_hg38/outs/split_frags/"$CELLTYPE".bam

# Index merged bam
samtools index ./celltype_bams/Set12_FC_"$CELLTYPE".bam -@ 20



	