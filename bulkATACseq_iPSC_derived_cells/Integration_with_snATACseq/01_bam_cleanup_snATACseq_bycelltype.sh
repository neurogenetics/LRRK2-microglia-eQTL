#!/bin/bash

# On Biowulf cluster http://hpc.nih.gov
# sbatch --cpus-per-task=21 --mem=100g --mail-type=BEGIN,END --time=2:00:00 script.sh CELLTYPE
# Using bam files produced from by_celltype scripts in snATACseq folder

CELLTYPE=$1
echo "$CELLTYPE"

# Remove mitochondrial reads with removeChrom.py from https://github.com/harvardinformatics/ATAC-seq
module load samtools
module load python/3.7

samtools view -h  ./"$CELLTYPE".bam  \
	|  python /data/langstonrg/ATACseq/removeChrom.py - - chrM  |  \
	samtools view -b -  >  ./"$CELLTYPE"_noChrM.bam

# Remove PCR duplicates
module load picard

java -XX:ParallelGCThreads=20 -jar $PICARDJARPATH/picard.jar MarkDuplicates \
	I=./"$CELLTYPE"_noChrM.bam \
	O=./"$CELLTYPE"_noChrM_noDup.bam \
	M=./"$CELLTYPE"_dups.txt \
	REMOVE_DUPLICATES=true
	
# Remove non-unique alignments
samtools view -b  -q 10 -@ 20  ./"$CELLTYPE"_noChrM_noDup.bam \
	>  ./"$CELLTYPE"_filtered.bam
	
# Remove intermediate files
rm ./"$CELLTYPE"_noChrM.bam ./"$CELLTYPE"_noChrM_noDup.bam

# Index
samtools index ./"$CELLTYPE"_filtered.bam -@ 20
