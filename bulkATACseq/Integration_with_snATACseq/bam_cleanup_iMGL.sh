#!/bin/bash

# On Biowulf cluster http://hpc.nih.gov
# sbatch --cpus-per-task=21 --mem=100g --mail-type=BEGIN,END script.sh SAMPLE

SAMPLE=$1
echo "$SAMPLE"

# Remove mitochondrial reads with removeChrom.py from https://github.com/harvardinformatics/ATAC-seq
module load samtools
module load python/3.7

samtools view -h -@ 20  ./mapped_reads/PPMI_"$SAMPLE"_iMGL.bam  \
	|  python ../removeChrom.py - - chrM  |  \
	samtools view -b -@ 20 -  >  ./mapped_reads/PPMI_"$SAMPLE"_iMGL_noChrM.bam
  
# Remove PCR duplicates
module load picard

java -XX:ParallelGCThreads=20 -jar $PICARDJARPATH/picard.jar MarkDuplicates \
	I=./mapped_reads/PPMI_"$SAMPLE"_iMGL_noChrM.bam \
	O=./mapped_reads/PPMI_"$SAMPLE"_iMGL_noChrM_noDup.bam \
	M=./mapped_reads/"$SAMPLE"_dups.txt \
	REMOVE_DUPLICATES=true

# Remove non-unique alignments
samtools view -b  -q 10 -@ 20 ./mapped_reads/PPMI_"$SAMPLE"_iMGL_noChrM_noDup.bam \
	>  ./mapped_reads/PPMI_"$SAMPLE"_iMGL_filtered.bam
	
# Remove intermediate files
rm ./mapped_reads/PPMI_"$SAMPLE"_iMGL_noChrM.bam
rm ./mapped_reads/PPMI_"$SAMPLE"_iMGL_noChrM_noDup.bam
	
# Index
samtools index ./mapped_reads/PPMI_"$SAMPLE"_iMGL_filtered.bam -@ 20
