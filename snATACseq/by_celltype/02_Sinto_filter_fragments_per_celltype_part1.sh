#!/bin/bash

# On Biowulf cluster http://hpc.nih.gov
# sbatch --cpus-per-task=20 --mem=100g --time=1:00:00 --mail-type=BEGIN,END scripts/Sinto_filterbarcodes.sh SAMPLE

SAMPLE=$1

echo "$SAMPLE"

# Change directories into data folder
cd ./data/"$SAMPLE"_hg38/outs

module load python/3.7

# pip install --user sinto  #if needed

# Filter cell barcodes by broadcelltype for input into bulk ATACseq pipeline
sinto filterbarcodes -b ./possorted_bam.bam -c ./"$SAMPLE"_bcs_celltype.tsv -p 20

