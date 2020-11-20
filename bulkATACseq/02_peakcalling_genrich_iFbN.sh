#!/bin/bash

# On Biowulf cluster http://hpc.nih.gov
# sbatch --mem=100g --cpus-per-task=20 --mail-type=BEGIN,END --time=8:00:00 script.sh SAMPLE

SAMPLE=$1
echo "$SAMPLE"

# Align trimmed reads
module load bowtie/2
module load samtools

export BOWTIE2_INDEXES=/fdb/igenomes/Homo_sapiens/NCBI/GRCh38/Sequence/Bowtie2Index/


bowtie2 -p $SLURM_CPUS_PER_TASK --very-sensitive  -x genome \
	-1 ./FBn_comparison_XR/data/"$SAMPLE"_noadapters_1.fastq.gz \
	-2 ./FBn_comparison_XR/data/"$SAMPLE"_noadapters_2.fastq.gz \
  |  samtools view -u -  \
  |  samtools sort -  >  ./FBn_comparison_XR/mapped_reads/"$SAMPLE".bam

# Sort bam file by queryname
samtools sort -o ./FBn_comparison_XR/mapped_reads/"$SAMPLE"_sorted.bam  -n  -@ 20  ./FBn_comparison_XR/mapped_reads/"$SAMPLE".bam

# Call peaks with genrich  (https://github.com/jsh58/Genrich)
module load genrich

Genrich  -t ./FBn_comparison_XR/mapped_reads/"$SAMPLE"_sorted.bam  \
	-o ./FBn_comparison_XR/peaks_genrich/"$SAMPLE".narrowPeak -j  -y  -r  -v  -e chrM  \
	-E ./ENCODEhg38blacklist.v2.bed.gz  \
	-k ./FBn_comparison_XR/peaks_genrich/"$SAMPLE"_pileups.bed  \
	-b ./FBn_comparison_XR/peaks_genrich/"$SAMPLE"_fragments.bed
	
# Process pileups for IGV viewing
awk 'NR>2 { print $1, $2, $3, $4 }' ./FBn_comparison_XR/peaks_genrich/"$SAMPLE"_pileups.bed >  \
	./FBn_comparison_XR/peaks_genrich/"$SAMPLE"_pileups.bdg

module load IGV IGVTools

igvtools toTDF ./FBn_comparison_XR/peaks_genrich/"$SAMPLE"_pileups.bdg ./FBn_comparison_XR/peaks_genrich/"$SAMPLE"_pileups.bdg.tdf hg38

