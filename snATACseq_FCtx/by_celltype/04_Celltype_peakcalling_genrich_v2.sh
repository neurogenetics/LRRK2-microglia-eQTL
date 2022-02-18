#!/bin/bash

# sbatch --mem=100g --cpus-per-task=20 --mail-type=BEGIN,END --time=1:00:00 Celltype_peakcalling_genrich_v2.sh CELLTYPE

CELLTYPE=$1
echo "$CELLTYPE"

# Change into snATACseq data folder
cd /data/langstonrg/snATACseq

# sort bam file by queryname
module load samtools #v1.11

samtools sort -o ./data/celltype_bams/Set12_FC_"$CELLTYPE"_sorted.bam  -n  -@ 20  ./data/celltype_bams/Set12_FC_"$CELLTYPE".bam

# call peaks with genrich
## https://github.com/harvardinformatics/ATAC-seq
## https://github.com/jsh58/Genrich#pileup
module load genrich  #v0.6

Genrich  -t ./data/celltype_bams/Set12_FC_"$CELLTYPE"_sorted.bam \
	-o ./Set12_FC_peaks_genrich/Set12_FC_"$CELLTYPE".narrowPeak \
	-j  -y  -r  -v  \
	-e chrM  \
	-E ./ENCODEhg38blacklist.v2.bed.gz  \
	-k ./Set12_FC_peaks_genrich/Set12_FC_"$CELLTYPE"_pileups.bed  \
	-b ./Set12_FC_peaks_genrich/Set12_FC_"$CELLTYPE"_fragments.bed
	
# Process pileups for IGV viewing
awk 'NR>2 { print $1, $2, $3, $4 }' ./Set12_FC_peaks_genrich/Set12_FC_"$CELLTYPE"_pileups.bed > ./Set12_FC_peaks_genrich/Set12_FC_"$CELLTYPE"_pileups.bdg

module load IGV IGVTools

igvtools toTDF ./Set12_FC_peaks_genrich/Set12_FC_"$CELLTYPE"_pileups.bdg ./Set12_FC_peaks_genrich/Set12_FC_"$CELLTYPE"_pileups.bdg.tdf hg38
