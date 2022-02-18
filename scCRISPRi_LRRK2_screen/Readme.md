## Analysis for CRISPR-dcas9 experiment towards understanding the cis-regultator elements of LRRK2 in microglia differentiated from human iPSCs (iMicroglia)

- current analysis suggests that common PD risk variant captures variation that leads to increased expression of LRRK2
- typically we can see LRRK2 eQTL in blood, but more diffficult in human brain expression cohorts
    - we previously seen suggestions of LRRK2 eQTL in NABEC DLPFC (need to double check if can see in Metabrain)
- a snRNAseq experiment was done in selected subjects from NABEC DLPFC samples
    - where selection was based on balanced number of genotypes for PD risk index variant near LRRK2, rs76904798
    - results showed that PD risk variant effect was an increase of LRRK2 in microglia
    - snATAC was performed to look for LRRK2 CRE locations but results were inconclusive
- a CRISPR-dcas9 experiment was done in iMicroglia to try and determine CREs for LRRK2 in microglia
    - note that the cell line used here is reference homozygous for the PD risk index variant, so this is more about LRRK2 CREs in microglia but may not be conclusive in regards to PD risk; as data interagated thus far is the dosage of alternate allele leads to increase of LRRK2 expression in microglia

### guides were designed against elements near LRRK2
| #	| Gene Symbol/Target Name	| # of sgRNA
| - | --------------------------| ----------
| 1	| Non_Targeting_Human_CRi	| 50
| 2	| peak_38433	| 46
| 3	| peak_38437	| 34
| 4	| peak_38439	| 11
| 5	| rs7294619	| 3
| 6	| peak_38440	| 22
| 7	| peak_38441	| 367

### 10X Genomics did processing of the three experiments using a modified version of their pipeline

### we did analysis of the processed single-cell data
- combine_experiments.ipynb, combines the experiments and adds in the guides info
- label_cells_maca.ipynb, uses MACA and couple different cell marker DBs to label the cell types
- scanpy_processing.ipynb, adds MACA generated cell labels and runs the more typical QC and clustering of the data
- reg_disease_gene_cell.ipynb, looks at effects of guides on disease relevant gene expression in specified tissue
    - looks at effects of designed guides on the expressoin of LRRK2 in iMicroglia

