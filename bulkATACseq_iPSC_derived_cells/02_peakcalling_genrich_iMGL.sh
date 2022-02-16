#!/bin/bash

# On Biowulf cluster http://hpc.nih.gov
# sbatch --mem=100g --cpus-per-task=20 --mail-type=BEGIN,END --time=3:00:00 script.sh SAMPLE

SAMPLE=$1
echo "$SAMPLE"

# Sort bam file by queryname
module load samtools

samtools sort -o ./mapped_reads/PPMI_iMGL_"$SAMPLE"_sorted.bam  -n  -@ 20  ./mapped_reads/PPMI_iMGL_"$SAMPLE"_1x.bam

# Call peaks with genrich (https://github.com/jsh58/Genrich)
module load genrich

Genrich  -t ./mapped_reads/PPMI_iMGL_"$SAMPLE"_sorted.bam  \
	-o ./peaks_genrich/PPMI_iMGL_"$SAMPLE".narrowPeak -j  -y  -r  -v  -e chrM  \
	-E ./ENCODEhg38blacklist.v2.bed.gz  \
	-k ./peaks_genrich/PPMI_iMGL_"$SAMPLE"_pileups.bed  \
	-b ./peaks_genrich/PPMI_iMGL_"$SAMPLE"_fragments.bed
	
# Process pileups for IGV viewing
awk 'NR>2 { print $1, $2, $3, $4 }' ./peaks_genrich/PPMI_iMGL_"$SAMPLE"_pileups.bed >  \
	./peaks_genrich/PPMI_iMGL_"$SAMPLE"_pileups.bdg

module load IGV IGVTools

igvtools toTDF ./peaks_genrich/PPMI_iMGL_"$SAMPLE"_pileups.bdg ./peaks_genrich/PPMI_iMGL_"$SAMPLE"_pileups.bdg.tdf hg38
