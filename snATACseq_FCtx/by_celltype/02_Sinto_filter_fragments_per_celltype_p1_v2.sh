#!/bin/bash

# sbatch --cpus-per-task=20 --mem=100g --gres=lscratch:100 --time=1:00:00 --mail-type=BEGIN,END Sinto_filterbarcodes_p1.sh UMARYID

SAMPLE=$1
echo "$SAMPLE"

module load python/3.8
module load samtools  #v1.11

# pip install --user sinto  #required once only
## Successfully installed sinto-0.7.2.2

# Change into outs/split_frags folder
cd /data/langstonrg/snATACseq/data/"$SAMPLE"_hg38/outs/split_frags

# Filter cell barcodes by broadcelltype for input into bulk ATACseq pipeline
sinto filterbarcodes \
	-b ../possorted_bam.bam \
	-c /data/langstonrg/snATACseq/output/"$SAMPLE"_bcs_celltype.tsv \
	-p 20


	