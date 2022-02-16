#!/bin/bash

# On Biowulf cluster http://hpc.nih.gov
# sbatch --cpus-per-task=20 --mem=100g --mail-type=BEGIN,END --time=2:00:00 script.sh CELLTYPE
# Using bam files produced from by_celltype scripts in snATACseq folder

CELLTYPE=$1
echo "$CELLTYPE"

# Change into directory containing bam files split by cell type
cd /data/langstonrg/snATACseq/data/celltype_bams

# Remove mitochondrial reads with removeChrom.py from https://github.com/harvardinformatics/ATAC-seq
module load samtools
module load python/3.8

samtools view -h  ./Set12_FC_"$CELLTYPE".bam  \
	|  python /data/langstonrg/ATACseq/removeChrom.py - - chrM  |  \
	samtools view -b -  >  ./Set12_FC_"$CELLTYPE"_noChrM.bam

# Remove PCR duplicates
module load picard

java -XX:ParallelGCThreads=20 -jar $PICARDJARPATH/picard.jar MarkDuplicates \
	I=./Set12_FC_"$CELLTYPE"_noChrM.bam \
	O=./Set12_FC_"$CELLTYPE"_noChrM_noDup.bam \
	M=./Set12_FC_"$CELLTYPE"_dups.txt \
	REMOVE_DUPLICATES=true
	
# Remove non-unique alignments
samtools view -b  -q 10 -@ 20  ./Set12_FC_"$CELLTYPE"_noChrM_noDup.bam \
	>  ./Set12_FC_"$CELLTYPE"_filtered.bam
	
# Remove intermediate files
rm ./Set12_FC_"$CELLTYPE"_noChrM.bam ./Set12_FC_"$CELLTYPE"_noChrM_noDup.bam

# Index
samtools index ./Set12_FC_"$CELLTYPE"_filtered.bam -@ 20