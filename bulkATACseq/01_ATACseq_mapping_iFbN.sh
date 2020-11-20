#!/bin/bash

# On Biowulf cluster http://hpc.nih.gov
# sbatch --cpus-per-task=20 --mem=100g --mail-type=BEGIN,END script.sh

SAMPLE=$1
echo "$SAMPLE"

# Trim adapter sequences with NGmerge https://github.com/jsh58/NGmerge, then align
module load bowtie/2
module load samtools

export BOWTIE2_INDEXES=/fdb/igenomes/Homo_sapiens/NCBI/GRCh38/Sequence/Bowtie2Index/

./NGmerge-master/NGmerge  -a  \
	-1 ./FBn_comparison_XR/data/33i_FBn_"$SAMPLE"_R1.fastq.gz  \
	-2 ./FBn_comparison_XR/data/33i_FBn_"$SAMPLE"_R2.fastq.gz  \
	-o ./FBn_comparison_XR/data/33i_FBn_"$SAMPLE"_noadapters

bowtie2 -p $SLURM_CPUS_PER_TASK --very-sensitive  -x genome \
	-1 ./FBn_comparison_XR/data/33i_FBn_"$SAMPLE"_noadapters_1.fastq.gz \
	-2 ./FBn_comparison_XR/data/33i_FBn_"$SAMPLE"_noadapters_2.fastq.gz \
  |  samtools view -u -  \
  |  samtools sort -  >  ./FBn_comparison_XR/mapped_reads/33i_FBn_"$SAMPLE".bam


  