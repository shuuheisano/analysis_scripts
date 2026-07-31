###################################################
### DEG overlap and transcription factor analysis
###################################################

library(dplyr)
library(ggplot2)
library(ggVennDiagram)


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


###################################################
### Input and output paths ----

# DESeq2 DEG tables
deg_minus_file <- "path/to/DEG_NMO_minus_vs_HC_minus_FC1.5_padj0.05.csv"

deg_plus_file <- "path/to/DEG_NMO_plus_vs_HC_plus_FC1.5_padj0.05.csv"

# Human transcription factor list obtained from AnimalTFDB
tf_file <- "path/to/Homo_sapiens_TF.txt"

# Output directory
output_dir <- "path/to/output_directory"

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


###################################################
### Analysis parameters ----

padj_cutoff <- 0.05
fold_change_cutoff <- 1.5
log2fc_cutoff <- log2(fold_change_cutoff)


###################################################
### Import DEG tables ----

deg_minus <- read.csv(
  deg_minus_file,
  header = TRUE,
  check.names = FALSE
)

deg_plus <- read.csv(
  deg_plus_file,
  header = TRUE,
  check.names = FALSE
)


###################################################
### Extract significant DEGs ----

# NMO_minus vs HC_minus
genes_deg_minus <- deg_minus |>
  dplyr::filter(
    !is.na(padj),
    padj < padj_cutoff,
    abs(log2FoldChange) >= log2fc_cutoff,
    !is.na(Symbol),
    Symbol != ""
  ) |>
  dplyr::pull(Symbol) |>
  unique()

# NMO_plus vs HC_plus
genes_deg_plus <- deg_plus |>
  dplyr::filter(
    !is.na(padj),
    padj < padj_cutoff,
    abs(log2FoldChange) >= log2fc_cutoff,
    !is.na(Symbol),
    Symbol != ""
  ) |>
  dplyr::pull(Symbol) |>
  unique()

# Check the number of DEGs
length(genes_deg_minus)
length(genes_deg_plus)


###################################################
### Venn diagram ----

venn_list <- list(
  "NMO_minus vs HC_minus" = genes_deg_minus,
  "NMO_plus vs HC_plus" = genes_deg_plus
)

p_venn <- ggVennDiagram::ggVennDiagram(
  venn_list,
  label_alpha = 0,
  category.names = c(
    "NMO CD25− NB\nvs HC CD25− NB",
    "NMO CD25hi NB\nvs HC CD25hi NB"
  )
) +
  ggplot2::scale_fill_gradient(
    low = "white",
    high = "red"
  ) +
  ggplot2::theme(
    legend.position = "none",
    text = ggplot2::element_text(size = 9)
  )

print(p_venn)

ggplot2::ggsave(
  filename = file.path(
    output_dir,
    "Venn_DEGs_NMO_minus_vs_HC_minus_and_NMO_plus_vs_HC_plus.pdf"
  ),
  plot = p_venn,
  width = 6,
  height = 6,
  device = "pdf"
)


###################################################
### NMO_plus-specific DEGs ----

# Genes detected as DEGs only in NMO_plus vs HC_plus
specific_genes_NMO_plus <- setdiff(
  genes_deg_plus,
  genes_deg_minus
)

# Check the number of NMO_plus-specific genes
length(specific_genes_NMO_plus)

specific_genes_NMO_plus_df <- data.frame(
  Symbol = specific_genes_NMO_plus
)

# Export NMO_plus-specific genes
write.csv(
  specific_genes_NMO_plus_df,
  file = file.path(
    output_dir,
    "NMO_plus_vs_HC_plus_specific_genes.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)


###################################################
### Extract transcription factors ----

# Human transcription factor list obtained from AnimalTFDB
#
# This script assumes that the TF table contains a column named "Symbol".
tfdb <- read.delim(
  tf_file,
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# Identify transcription factors among NMO_plus-specific DEGs
tf_genes_NMO_plus <- specific_genes_NMO_plus_df |>
  dplyr::inner_join(
    tfdb,
    by = "Symbol"
  )

# Check extracted transcription factors
dim(tf_genes_NMO_plus)
print(tf_genes_NMO_plus)

# Export transcription factors
write.csv(
  tf_genes_NMO_plus,
  file = file.path(
    output_dir,
    "NMO_plus_vs_HC_plus_specific_transcription_factors.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)