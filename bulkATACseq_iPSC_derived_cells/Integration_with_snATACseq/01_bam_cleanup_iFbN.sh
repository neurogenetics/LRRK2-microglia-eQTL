#!/bin/bash

# On Biowulf cluster http://hpc.nih.gov
# sbatch --cpus-per-task=21 --mem=100g --mail-type=BEGIN,END script.sh SAMPLE

SAMPLE=$1
echo "$SAMPLE"

# Remove mitochondrial reads with removeChrom.py from https://github.com/harvardinformatics/ATAC-seq
module load samtools
module load python/3.7

samtools view -h  ./FBn_comparison_XR/mapped_reads/33i_FBn_"$SAMPLE".bam  \
	|  python removeChrom.py - - chrM  |  \
	samtools view -b -  >  ./FBn_comparison_XR/mapped_reads/33i_FBn_"$SAMPLE"_noChrM.bam

# Remove PCR duplicates
module load picard

java -XX:ParallelGCThreads=20 -jar $PICARDJARPATH/picard.jar MarkDuplicates \
	I=./FBn_comparison_XR/mapped_reads/33i_FBn_"$SAMPLE"_noChrM.bam \
	O=./FBn_comparison_XR/mapped_reads/33i_FBn_"$SAMPLE"_noChrM_noDup.bam \
	M=./FBn_comparison_XR/mapped_reads/"$SAMPLE"_dups.txt \
	REMOVE_DUPLICATES=true
	
# Remove non-unique alignments
samtools view -b  -q 10 -@ 20  ./FBn_comparison_XR/mapped_reads/33i_FBn_"$SAMPLE"_noChrM_noDup.bam \
	>  ./FBn_comparison_XR/mapped_reads/33i_FBn_"$SAMPLE"_filtered.bam

# Remove intermediate files
rm ./FBn_comparison_XR/mapped_reads/33i_FBn_"$SAMPLE"_noChrM.bam \
	./FBn_comparison_XR/mapped_reads/33i_FBn_"$SAMPLE"_noChrM_noDup.bam
	
# Index
samtools index ./FBn_comparison_XR/mapped_reads/33i_FBn_"$SAMPLE"_filtered.bam -@ 20

  

