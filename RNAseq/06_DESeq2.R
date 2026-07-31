###################################################
### Differential expression analysis using DESeq2
###################################################

###################################################
### Sample group definitions ----
#
# Code names used in this script:
#
# NMO_plus   : NMO CD25hi NB
# NMO_minus  : NMO CD25− NB
# HC_plus    : HC CD25hi NB
# HC_minus   : HC CD25− NB
###################################################

library(DESeq2)

###################################################
### Input and output paths ----

# Input gene count matrix
gene_count_file <- "path/to/gene_count_matrix.csv"

# Output directory
output_dir <- "path/to/output_directory"

# Thresholds used to define differentially expressed genes
padj_cutoff <- 0.05
fold_change_cutoff <- 1.5
log2fc_cutoff <- log2(fold_change_cutoff)


###################################################
### Import raw gene count table and metadata ----

# (A) Read and check raw gene count data
gene_count_table <- read.csv(
  gene_count_file,
  header = TRUE,
  check.names = FALSE
)

# Basic checks
dim(gene_count_table)
colnames(gene_count_table)
class(gene_count_table)

# Check whether gene symbols are unique
gene_count_table |>
  dplyr::count(Symbol, sort = TRUE)

# Summarize the count table
skimr::skim(gene_count_table)


# (B) Create sample metadata
#
# The order of row names in coldata must match the order of
# sample columns in the gene count table.

coldata <- data.frame(
  sample = c(
    "NMO3_minus",
    "NMO2_minus",
    "NMO4_minus",
    "NMO8_minus",
    "NMO7_minus",
    "NMO5_plus",
    "NMO2_plus",
    "NMO3_plus",
    "NMO6_plus",
    "NMO8_plus",
    "NMO7_plus",
    "NMO5_minus",
    "NMO4_plus",
    "NMO6_minus",
    "HC3_plus",
    "HC4_minus",
    "HC2_plus",
    "HC5_minus",
    "HC4_plus",
    "HC3_minus",
    "HC2_minus",
    "HC5_plus"
  ),
  condition = c(
    "NMO_minus",
    "NMO_minus",
    "NMO_minus",
    "NMO_minus",
    "NMO_minus",
    "NMO_plus",
    "NMO_plus",
    "NMO_plus",
    "NMO_plus",
    "NMO_plus",
    "NMO_plus",
    "NMO_minus",
    "NMO_plus",
    "NMO_minus",
    "HC_plus",
    "HC_minus",
    "HC_plus",
    "HC_minus",
    "HC_plus",
    "HC_minus",
    "HC_minus",
    "HC_plus"
  ),
  row.names = "sample"
)

# Basic checks
dim(coldata)
colnames(coldata)
class(coldata)

# Convert condition to a factor
coldata$condition <- factor(coldata$condition)

# Confirm the number of samples in each condition
table(coldata$condition)


###################################################
### Prepare count matrix ----

# Move gene symbols from a column to row names
data_DE <- gene_count_table |>
  tibble::column_to_rownames("Symbol")


###################################################
### Create DESeq2 object ----

DESeq2Object <- DESeq2::DESeqDataSetFromMatrix(
  countData = data_DE,
  colData = coldata,
  design = ~ condition
)


###################################################
### Filter count data ----

dim(DESeq2Object)

# Keep genes expressed in at least one sample
keep <- rowSums(
  DESeq2::counts(DESeq2Object)
) > 0

DESeq2Object <- DESeq2Object[keep, ]

dim(DESeq2Object)


###################################################
### Variance-stabilizing transformation ----

vst <- DESeq2::vst(DESeq2Object, blind = TRUE)
hist(assay(vst))

###################################################
### Run DESeq2 ----

DESeq2Object <- DESeq2::DESeq(DESeq2Object)

# Display fitted DESeq2 object
DESeq2Object


###################################################
### NMO_plus vs HC_plus ----

my_contrast_1 <- c(
  "condition",
  "NMO_plus",
  "HC_plus"
)

deseq2Results_1 <- DESeq2::results(
  DESeq2Object,
  contrast = my_contrast_1,
  alpha = padj_cutoff
)

DESeq2::summary(deseq2Results_1)


###################################################
### Export all DEGs: NMO_plus vs HC_plus ----

# Convert DESeq2 results to a data frame
DEG_NMO_plus_vs_HC_plus <- as.data.frame(
  deseq2Results_1
)

# Move gene symbols from row names to a column
DEG_NMO_plus_vs_HC_plus$Symbol <- rownames(
  DEG_NMO_plus_vs_HC_plus
)

# Exclude genes with NA adjusted P-values and retain genes with:
# adjusted P-value < 0.05 and absolute fold change >= 1.5
DEG_NMO_plus_vs_HC_plus <- DEG_NMO_plus_vs_HC_plus[
  !is.na(DEG_NMO_plus_vs_HC_plus$padj) &
    DEG_NMO_plus_vs_HC_plus$padj < padj_cutoff &
    abs(DEG_NMO_plus_vs_HC_plus$log2FoldChange) >= log2fc_cutoff,
]

# Add direction of differential expression
DEG_NMO_plus_vs_HC_plus$regulation <- ifelse(
  DEG_NMO_plus_vs_HC_plus$log2FoldChange >= log2fc_cutoff,
  "Upregulated in NMO_plus",
  "Downregulated in NMO_plus"
)

# Add fold change while preserving its direction
DEG_NMO_plus_vs_HC_plus$foldChange <- 2^
  DEG_NMO_plus_vs_HC_plus$log2FoldChange

# Reorder columns
DEG_NMO_plus_vs_HC_plus <- DEG_NMO_plus_vs_HC_plus[
  ,
  c(
    "Symbol",
    "baseMean",
    "log2FoldChange",
    "foldChange",
    "lfcSE",
    "stat",
    "pvalue",
    "padj",
    "regulation"
  )
]

# Sort by adjusted P-value
DEG_NMO_plus_vs_HC_plus <- DEG_NMO_plus_vs_HC_plus[
  order(DEG_NMO_plus_vs_HC_plus$padj),
]

# Reset row names
rownames(DEG_NMO_plus_vs_HC_plus) <- NULL

# Check output
dim(DEG_NMO_plus_vs_HC_plus)
head(DEG_NMO_plus_vs_HC_plus)
table(DEG_NMO_plus_vs_HC_plus$regulation)

# Export DEGs
write.csv(
  DEG_NMO_plus_vs_HC_plus,
  file = file.path(
    output_dir,
    "DEG_NMO_plus_vs_HC_plus_FC1.5_padj0.05.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)


###################################################
### NMO_minus vs HC_minus ----

my_contrast_2 <- c(
  "condition",
  "NMO_minus",
  "HC_minus"
)

deseq2Results_2 <- DESeq2::results(
  DESeq2Object,
  contrast = my_contrast_2,
  alpha = padj_cutoff
)

DESeq2::summary(deseq2Results_2)


###################################################
### Export all DEGs: NMO_minus vs HC_minus ----

DEG_NMO_minus_vs_HC_minus <- as.data.frame(
  deseq2Results_2
)

DEG_NMO_minus_vs_HC_minus$Symbol <- rownames(
  DEG_NMO_minus_vs_HC_minus
)

DEG_NMO_minus_vs_HC_minus <- DEG_NMO_minus_vs_HC_minus[
  !is.na(DEG_NMO_minus_vs_HC_minus$padj) &
    DEG_NMO_minus_vs_HC_minus$padj < padj_cutoff &
    abs(DEG_NMO_minus_vs_HC_minus$log2FoldChange) >= log2fc_cutoff,
]

DEG_NMO_minus_vs_HC_minus$regulation <- ifelse(
  DEG_NMO_minus_vs_HC_minus$log2FoldChange >= log2fc_cutoff,
  "Upregulated in NMO_minus",
  "Downregulated in NMO_minus"
)

DEG_NMO_minus_vs_HC_minus$foldChange <- 2^
  DEG_NMO_minus_vs_HC_minus$log2FoldChange

DEG_NMO_minus_vs_HC_minus <- DEG_NMO_minus_vs_HC_minus[
  ,
  c(
    "Symbol",
    "baseMean",
    "log2FoldChange",
    "foldChange",
    "lfcSE",
    "stat",
    "pvalue",
    "padj",
    "regulation"
  )
]

DEG_NMO_minus_vs_HC_minus <- DEG_NMO_minus_vs_HC_minus[
  order(DEG_NMO_minus_vs_HC_minus$padj),
]

rownames(DEG_NMO_minus_vs_HC_minus) <- NULL

# Check output
dim(DEG_NMO_minus_vs_HC_minus)
head(DEG_NMO_minus_vs_HC_minus)
table(DEG_NMO_minus_vs_HC_minus$regulation)

# Export DEGs
write.csv(
  DEG_NMO_minus_vs_HC_minus,
  file = file.path(
    output_dir,
    "DEG_NMO_minus_vs_HC_minus_FC1.5_padj0.05.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)
