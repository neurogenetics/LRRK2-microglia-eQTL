#!/bin/bash

# sbatch --cpus-per-task=20 --mem=100g --time=8:00:00 --mail-type=BEGIN,END

cd ./data
module load samtools

# decompress files and add the same cell prefix as was added to the Seurat object
gzip -dc UMARY794_hg38/outs/fragments.tsv.gz | awk 'BEGIN {FS=OFS="\t"} {print $1,$2,$3,"S794_"$4,$5}' - > S794_fragments.tsv
gzip -dc UMARY1230_hg38/outs/fragments.tsv.gz | awk 'BEGIN {FS=OFS="\t"} {print $1,$2,$3,"S1230_"$4,$5}' - > S1230_fragments.tsv
gzip -dc UMARY1274_hg38/outs/fragments.tsv.gz | awk 'BEGIN {FS=OFS="\t"} {print $1,$2,$3,"S1274_"$4,$5}' - > S1274_fragments.tsv
gzip -dc UMARY5079_hg38/outs/fragments.tsv.gz | awk 'BEGIN {FS=OFS="\t"} {print $1,$2,$3,"S5079_"$4,$5}' - > S5079_fragments.tsv

gzip -dc FC_ATAC_4022_hg38/outs/fragments.tsv.gz | awk 'BEGIN {FS=OFS="\t"} {print $1,$2,$3,"S4022_"$4,$5}' - > S4022_fragments.tsv
gzip -dc FC_ATAC_4924_hg38/outs/fragments.tsv.gz | awk 'BEGIN {FS=OFS="\t"} {print $1,$2,$3,"S4924_"$4,$5}' - > S4924_fragments.tsv
gzip -dc FC_ATAC_4724_hg38/outs/fragments.tsv.gz | awk 'BEGIN {FS=OFS="\t"} {print $1,$2,$3,"S4724_"$4,$5}' - > S4724_fragments.tsv
gzip -dc FC_ATAC_1209_hg38/outs/fragments.tsv.gz | awk 'BEGIN {FS=OFS="\t"} {print $1,$2,$3,"S1209_"$4,$5}' - > S1209_fragments.tsv

gzip -dc FC_ATAC_5123_hg38/outs/fragments.tsv.gz | awk 'BEGIN {FS=OFS="\t"} {print $1,$2,$3,"S5123_"$4,$5}' - > S5123_fragments.tsv
gzip -dc FC_ATAC_1363_hg38/outs/fragments.tsv.gz | awk 'BEGIN {FS=OFS="\t"} {print $1,$2,$3,"S1363_"$4,$5}' - > S1363_fragments.tsv
gzip -dc FC_ATAC_630_hg38/outs/fragments.tsv.gz | awk 'BEGIN {FS=OFS="\t"} {print $1,$2,$3,"S630_"$4,$5}' - > S630_fragments.tsv
gzip -dc FC_ATAC_1135_hg38/outs/fragments.tsv.gz | awk 'BEGIN {FS=OFS="\t"} {print $1,$2,$3,"S1135_"$4,$5}' - > S1135_fragments.tsv


## merge files (avoids having to re-sort)
sort -m -k 1,1V -k2,2n S794_fragments.tsv S1230_fragments.tsv S1274_fragments.tsv S5079_fragments.tsv \
	S4022_fragments.tsv S4924_fragments.tsv S4724_fragments.tsv S1209_fragments.tsv \
	S5123_fragments.tsv S1363_fragments.tsv S630_fragments.tsv S1135_fragments.tsv > Set12_FC_fragments.tsv

## block gzip compress the merged file
bgzip -@ 20 Set12_FC_fragments.tsv # -@ 8 uses 8 threads

## index the bgzipped file
tabix -p bed Set12_FC_fragments.tsv.gz

## remove intermediate files
rm S794_fragments.tsv S1230_fragments.tsv S1274_fragments.tsv S5079_fragments.tsv \
	S4022_fragments.tsv S4924_fragments.tsv S4724_fragments.tsv S1209_fragments.tsv \
	S5123_fragments.tsv S1363_fragments.tsv S630_fragments.tsv S1135_fragments.tsv

