#!/bin/bash

# On Biowulf cluster http://hpc.nih.gov
# sbatch --mem=100g --cpus-per-task=20 --mail-type=BEGIN,END --time=1:00:00 scripts/Celltype_peakcalling_genrich.sh CELLTYPE

CELLTYPE=$1

echo "$CELLTYPE"

# Change into data directory
cd ./data

# Sort bam file by queryname
module load samtools

samtools sort -o ./FC_"$CELLTYPE"_sorted.bam  -n  -@ 20  ./FC_"$CELLTYPE".bam

# Call peaks with genrich https://github.com/jsh58/Genrich
module load genrich

Genrich  -t ./FC_"$CELLTYPE"_sorted.bam  -o ./FC_peaks_genrich/FC_"$CELLTYPE".narrowPeak -j  -y  -r  -v  \
	-e chrM  \
	-E ./ENCODEhg38blacklist.v2.bed.gz  \
	-k ./FC_peaks_genrich/FC_"$CELLTYPE"_pileups.bed  \
	-b ./FC_peaks_genrich/FC_"$CELLTYPE"_fragments.bed
	
# Process pileups for IGV viewing
awk 'NR>2 { print $1, $2, $3, $4 }' ./FC_peaks_genrich/FC_"$CELLTYPE"_pileups.bed > ./FC_peaks_genrich/FC_"$CELLTYPE"_pileups.bdg

module load IGV IGVTools

igvtools toTDF ./FC_peaks_genrich/FC_"$CELLTYPE"_pileups.bdg ./FC_peaks_genrich/FC_"$CELLTYPE"_pileups.bdg.tdf hg38

