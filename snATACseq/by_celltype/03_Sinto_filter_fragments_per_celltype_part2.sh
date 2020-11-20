#!/bin/bash

# On Biowulf cluster http://hpc.nih.gov
# sbatch --cpus-per-task=20 --mem=150g --time=4:00:00 --mail-type=BEGIN,END scripts/Sinto_filterbarcodes_part2.sh


# Change into snATACseq data directory
cd ./data

# Merge bam files for each celltype and index
module load samtools

samtools merge -@ 20 ./FC_ODC.0.bam \
	./UMARY794_hg38/outs/ODC.0.bam \
	./UMARY1230_hg38/outs/ODC.0.bam \
	./UMARY1274_hg38/outs/ODC.0.bam \
	./UMARY5079_hg38/outs/ODC.0.bam
	
samtools merge -@ 20 ./FC_ODC.2.bam \
	./UMARY794_hg38/outs/ODC.2.bam \
	./UMARY1230_hg38/outs/ODC.2.bam \
	./UMARY1274_hg38/outs/ODC.2.bam \
	./UMARY5079_hg38/outs/ODC.2.bam
	
samtools merge -@ 20 ./FC_ODC.5.bam \
	./UMARY794_hg38/outs/ODC.5.bam \
	./UMARY1274_hg38/outs/ODC.5.bam \
	./UMARY5079_hg38/outs/ODC.5.bam
	
samtools merge -@ 20 ./FC_OPC.6.bam \
	./UMARY794_hg38/outs/OPC.6.bam \
	./UMARY1230_hg38/outs/OPC.6.bam \
	./UMARY1274_hg38/outs/OPC.6.bam \
	./UMARY5079_hg38/outs/OPC.6.bam
	
samtools merge -@ 20 ./FC_EC.13.bam \
	./UMARY794_hg38/outs/EC.13.bam \
	./UMARY1230_hg38/outs/EC.13.bam \
	./UMARY1274_hg38/outs/EC.13.bam \
	./UMARY5079_hg38/outs/EC.13.bam
	
samtools merge -@ 20 ./FC_AST.3.bam \
	./UMARY794_hg38/outs/AST.3.bam \
	./UMARY1230_hg38/outs/AST.3.bam \
	./UMARY1274_hg38/outs/AST.3.bam \
	./UMARY5079_hg38/outs/AST.3.bam
	
samtools merge -@ 20 ./FC_MGL.4.bam \
	./UMARY794_hg38/outs/MGL.4.bam \
	./UMARY1230_hg38/outs/MGL.4.bam \
	./UMARY1274_hg38/outs/MGL.4.bam \
	./UMARY5079_hg38/outs/MGL.4.bam
	
samtools merge -@ 20 ./FC_ExN.1.bam \
	./UMARY794_hg38/outs/ExN.1.bam \
	./UMARY1230_hg38/outs/ExN.1.bam \
	./UMARY1274_hg38/outs/ExN.1.bam \
	./UMARY5079_hg38/outs/ExN.1.bam
	
samtools merge -@ 20 ./FC_ExN.7.bam \
	./UMARY794_hg38/outs/ExN.7.bam \
	./UMARY1230_hg38/outs/ExN.7.bam \
	./UMARY1274_hg38/outs/ExN.7.bam \
	./UMARY5079_hg38/outs/ExN.7.bam
	
samtools merge -@ 20 ./FC_ExN.12.bam \
	./UMARY794_hg38/outs/ExN.12.bam \
	./UMARY1230_hg38/outs/ExN.12.bam \
	./UMARY1274_hg38/outs/ExN.12.bam \
	./UMARY5079_hg38/outs/ExN.12.bam
	
samtools merge -@ 20 ./FC_InN.8.bam \
	./UMARY794_hg38/outs/InN.8.bam \
	./UMARY1230_hg38/outs/InN.8.bam \
	./UMARY1274_hg38/outs/InN.8.bam \
	./UMARY5079_hg38/outs/InN.8.bam
	
samtools merge -@ 20 ./FC_InN.9.bam \
	./UMARY794_hg38/outs/InN.9.bam \
	./UMARY1230_hg38/outs/InN.9.bam \
	./UMARY1274_hg38/outs/InN.9.bam \
	./UMARY5079_hg38/outs/InN.9.bam
	
samtools merge -@ 20 ./FC_InN.10.bam \
	./UMARY794_hg38/outs/InN.10.bam \
	./UMARY1230_hg38/outs/InN.10.bam \
	./UMARY1274_hg38/outs/InN.10.bam \
	./UMARY5079_hg38/outs/InN.10.bam
	
samtools merge -@ 20 ./FC_InN.11.bam \
	./UMARY794_hg38/outs/InN.11.bam \
	./UMARY1230_hg38/outs/InN.11.bam \
	./UMARY1274_hg38/outs/InN.11.bam \
	./UMARY5079_hg38/outs/InN.11.bam


samtools index ./FC_ODC.0.bam -@ 20
samtools index ./FC_ODC.2.bam -@ 20
samtools index ./FC_ODC.5.bam -@ 20
samtools index ./FC_OPC.6.bam -@ 20
samtools index ./FC_EC.13.bam -@ 20
samtools index ./FC_AST.3.bam -@ 20
samtools index ./FC_MGL.4.bam -@ 20
samtools index ./FC_ExN.1.bam -@ 20
samtools index ./FC_ExN.7.bam -@ 20
samtools index ./FC_ExN.12.bam -@ 20
samtools index ./FC_InN.8.bam -@ 20
samtools index ./FC_InN.9.bam -@ 20
samtools index ./FC_InN.10.bam -@ 20
samtools index ./FC_InN.11.bam -@ 20

