#!/bin/bash

# sbatch --cpus-per-task=8 --mem=50g --time=4:00:00 --mail-type=BEGIN,END

cd ./data
module load samtools

# decompress files and add the same cell prefix as was added to the Seurat object
gzip -dc UMARY794_hg38/fragments.tsv.gz | awk 'BEGIN {FS=OFS="\t"} {print $1,$2,$3,"S794_"$4,$5}' - > S794_fragments.tsv
gzip -dc UMARY1230_hg38/fragments.tsv.gz | awk 'BEGIN {FS=OFS="\t"} {print $1,$2,$3,"S1230_"$4,$5}' - > S1230_fragments.tsv
gzip -dc UMARY1274_hg38/fragments.tsv.gz | awk 'BEGIN {FS=OFS="\t"} {print $1,$2,$3,"S1274_"$4,$5}' - > S1274_fragments.tsv
gzip -dc UMARY5079_hg38/fragments.tsv.gz | awk 'BEGIN {FS=OFS="\t"} {print $1,$2,$3,"S5079_"$4,$5}' - > S5079_fragments.tsv

## merge files (avoids having to re-sort)
sort -m -k 1,1V -k2,2n S794_fragments.tsv S1230_fragments.tsv S1274_fragments.tsv S5079_fragments.tsv > FC_fragments.tsv

## block gzip compress the merged file
bgzip -@ 8 FC_fragments.tsv # -@ 8 uses 8 threads

## index the bgzipped file
tabix -p bed FC_fragments.tsv.gz

## remove intermediate files
rm S794_fragments.tsv S1230_fragments.tsv S1274_fragments.tsv S5079_fragments.tsv

