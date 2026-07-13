Overview
This repository contains the scripts used for RNA-seq analysis in [ ]
The pipeline includes quality control, read preprocessing, alignment, gene quantification, differential expression analysis, PCA, and WGCNA.

Requirements
- FastQC v0.12.1 
- fastp v0.23.4 
- HISAT2 v2.2.1 
- samtools v1.21 
- StringTie v2.2.3 
- DESeq2 v1.44.0 
- WGCNA v1.73 

Input
- Raw FASTQ files (50-bp single-end) 
- UCSC hg38 reference genome 
- Corresponding gene annotation (GTF)

Usage
Run the scripts in the following order. 
1. ｀01_fastqc.sh 
2. ｀02_fastp.sh 
3. ｀03_hisat2.sh 
4. ｀04_stringtie.sh 
5. ｀05_merge_counts.R 
6. ｀06_deseq2.R 
7. ｀07_pca.R 
8. ｀08_wgcna.R


