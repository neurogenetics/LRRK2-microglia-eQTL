#!/bin/bash

# On Biowulf cluster http://hpc.nih.gov
# sbatch --cpus-per-task=20 --mem=100g --mail-type=BEGIN,END

# Generate qc reports
module load fastqc

fastqc ./data/PPMI_iMGL_3411_1x_R1.fastq.gz
fastqc ./data/PPMI_iMGL_3411_1x_R2.fastq.gz

fastqc ./data/PPMI_iMGL_3446_1x_R1.fastq.gz
fastqc ./data/PPMI_iMGL_3446_1x_R2.fastq.gz

fastqc ./data/PPMI_iMGL_3453_1x_R1.fastq.gz
fastqc ./data/PPMI_iMGL_3453_1x_R2.fastq.gz

fastqc ./data/PPMI_iMGL_4101_1x_R1.fastq.gz
fastqc ./data/PPMI_iMGL_4101_1x_R2.fastq.gz


# Trim adapter sequences with NGmerge https://github.com/jsh58/NGmerge
./NGmerge-master/NGmerge  -a  -1 ./data/PPMI_iMGL_3411_1x_R1.fastq.gz  \
	-2 ./data/PPMI_iMGL_3411_1x_R2.fastq.gz  -o ./data/PPMI_iMGL_3411_1x_noadapters
	
./NGmerge-master/NGmerge  -a  -1 ./data/PPMI_iMGL_3446_1x_R1.fastq.gz  \
	-2 ./data/PPMI_iMGL_3446_1x_R2.fastq.gz  -o ./data/PPMI_iMGL_3446_1x_noadapters
	
./NGmerge-master/NGmerge  -a  -1 ./data/PPMI_iMGL_3453_1x_R1.fastq.gz  \
	-2 ./data/PPMI_iMGL_3453_1x_R2.fastq.gz  -o ./data/PPMI_iMGL_3453_1x_noadapters
	
./NGmerge-master/NGmerge  -a  -1 ./data/PPMI_iMGL_4101_1x_R1.fastq.gz  \
	-2 ./data/PPMI_iMGL_4101_1x_R2.fastq.gz  -o ./data/PPMI_iMGL_4101_1x_noadapters


# Align
module load bowtie/2
module load samtools

export BOWTIE2_INDEXES=/fdb/igenomes/Homo_sapiens/NCBI/GRCh38/Sequence/Bowtie2Index/

bowtie2 -p 10 --very-sensitive  -x genome \
	-1 ./data/PPMI_iMGL_3411_1x_noadapters_1.fastq.gz \
	-2 ./data/PPMI_iMGL_3411_1x_noadapters_2.fastq.gz \
  |  samtools view -u -  \
  |  samtools sort -  >  ./mapped_reads/PPMI_iMGL_3411_1x.bam
  
  
bowtie2 -p 10 --very-sensitive  -x genome \
	-1 ./data/PPMI_iMGL_3446_1x_noadapters_1.fastq.gz \
	-2 ./data/PPMI_iMGL_3446_1x_noadapters_2.fastq.gz \
  |  samtools view -u -  \
  |  samtools sort -  >  ./mapped_reads/PPMI_iMGL_3446_1x.bam
  
bowtie2 -p 10 --very-sensitive  -x genome \
	-1 ./data/PPMI_iMGL_3453_1x_noadapters_1.fastq.gz \
	-2 ./data/PPMI_iMGL_3453_1x_noadapters_2.fastq.gz \
  |  samtools view -u -  \
  |  samtools sort -  >  ./mapped_reads/PPMI_iMGL_3453_1x.bam

bowtie2 -p 10 --very-sensitive  -x genome \
	-1 ./data/PPMI_iMGL_4101_1x_noadapters_1.fastq.gz \
	-2 ./data/PPMI_iMGL_4101_1x_noadapters_2.fastq.gz \
  |  samtools view -u -  \
  |  samtools sort -  >  ./mapped_reads/PPMI_iMGL_4101_1x.bam
  
  
