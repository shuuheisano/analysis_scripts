# Analysis scripts

This repository contains the analysis scripts used in the following study:

"Helios-Dependent Chromatin Remodeling Drives IFN-α–Responsive
Plasmablast Differentiation in NMOSD Naïve B Cells"

## Repository structure

### RNAseq

The RNA-seq analysis pipeline includes:

1. Quality control (FastQC, MultiQC)
2. Read trimming (fastp)
3. Read alignment (HISAT2)
4. Gene quantification (StringTie and prepDE.py3)
5. Differential expression analysis (DESeq2)
6. Principal component analysis (PCA)
7. Gene Ontology (GO) enrichment analysis (Metascape)
8. Gene set enrichment analysis (fgsea)
9. Weighted gene co-expression network analysis (WGCNA)

### ATACseq

The ATAC-seq analysis pipeline includes:

1. Quality control (FastQC)
2. Adapter trimming (Trimmomatic)
3. Read alignment (BWA-MEM)
4. Peak calling (MACS2)
5. Differential accessibility analysis (DiffBind)
6. Peak annotation (ChIPseeker)
7. Motif analysis (chromVAR)
8. BigWig generation and IGV visualization


