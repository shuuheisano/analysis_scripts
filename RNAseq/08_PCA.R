###################################################
### PCA analysis using selected transcription factors
###################################################

# ---- 1) Required packages ----

library(dplyr)
library(tibble)
library(gprofiler2)
library(ggplot2)
library(DESeq2)


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

# Input gene count matrix
gene_count_file <- "path/to/gene_count_matrix.csv"

# Mouse transcription factor list
mouse_tf_csv <- "path/to/Fig7a_mouse_TF_list.csv"

# Output directory
output_dir <- "path/to/output_directory"

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


###################################################
### Import raw gene count table and metadata ----

# (A) Read and check raw gene count data
gene_count_table <- read.csv(
  gene_count_file,
  header = TRUE,
  check.names = FALSE
)

# checks
dim(gene_count_table)
colnames(gene_count_table)
class(gene_count_table)

gene_count_table |>
  dplyr::count(Symbol, sort = TRUE)


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

dim(coldata)
colnames(coldata)
class(coldata)

coldata$condition <- factor(coldata$condition)


###################################################
### Create DESeq2 object ----

data_DE <- gene_count_table |>
  tibble::column_to_rownames("Symbol")

DESeq2Objct <- DESeq2::DESeqDataSetFromMatrix(
  countData = data_DE,
  colData = coldata,
  design = ~ condition
)


###################################################
### Filter and transform count data ----

dim(DESeq2Objct)

# Keep genes expressed in at least one sample
keep <- rowSums(
  DESeq2::counts(DESeq2Objct)
) > 0

DESeq2Objct <- DESeq2Objct[keep, ]

dim(DESeq2Objct)

# Variance-stabilizing transformation
vst <- DESeq2::vst(
  DESeq2Objct,
  blind = TRUE
)

# Extract VST matrix
expr_mat <- SummarizedExperiment::assay(vst)

head(rownames(expr_mat))
head(colnames(expr_mat))


###################################################
### Import mouse transcription factor list ----

# The CSV file must contain a column named "gene".
mouse_tf <- read.csv(
  mouse_tf_csv,
  header = TRUE,
  check.names = FALSE
)

head(mouse_tf)
colnames(mouse_tf)


###################################################
### Convert mouse genes to human orthologs ----

mapped <- gprofiler2::gorth(
  mouse_tf$gene,
  source_organism = "mmusculus",
  target_organism = "hsapiens"
)

head(mapped)
colnames(mapped)

# Extract required columns
mapped_df <- mapped[
  ,
  c(
    "input",
    "ortholog_name"
  )
]

# Remove duplicated mouse gene entries
mapped_df <- mapped_df[
  !duplicated(mapped_df$input),
]

head(mapped_df)
colnames(mapped_df)

# Join mouse and human gene symbols
joined <- dplyr::left_join(
  mouse_tf,
  mapped_df,
  by = c(
    "gene" = "input"
  )
)

mouse_tf_hs <- joined |>
  dplyr::rename(
    human_symbol = ortholog_name
  ) |>
  dplyr::filter(
    !is.na(human_symbol),
    human_symbol != ""
  )

# Human TF gene symbols
tf_human <- unique(
  mouse_tf_hs$human_symbol
)

# Export ortholog mapping
write.csv(
  mouse_tf_hs,
  file = file.path(
    output_dir,
    "Fig7a_mouse_to_human_TF_list.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)


###################################################
### Row-wise Z-score normalization function ----

scale_rows <- function(x) {

  mu <- rowMeans(
    x,
    na.rm = TRUE
  )

  sdv <- apply(
    x,
    1,
    sd,
    na.rm = TRUE
  )

  sdv[
    is.na(sdv) |
      sdv == 0
  ] <- 1

  (x - mu) / sdv
}


###################################################
### PCA: NMO_plus vs HC_plus ----

# Select target conditions
target_conditions_plus <- c(
  "NMO_plus",
  "HC_plus"
)

# Select target samples
target_samples_plus <- rownames(coldata)[
  coldata$condition %in% target_conditions_plus
]

# Subset VST matrix
expr_mat_subset_plus <- expr_mat[
  ,
  target_samples_plus,
  drop = FALSE
]

# Extract selected TF genes
tf_human_in_data_plus <- intersect(
  rownames(expr_mat_subset_plus),
  tf_human
)

vst_tf_plus <- expr_mat_subset_plus[
  tf_human_in_data_plus,
  ,
  drop = FALSE
]

# Z-score normalization
vst_tf_z_plus <- scale_rows(
  vst_tf_plus
)

# PCA
pca_plus <- prcomp(
  t(vst_tf_z_plus),
  scale. = FALSE
)

# Proportion of variance explained
pve_plus <- (
  pca_plus$sdev^2
) / sum(
  pca_plus$sdev^2
)

pc1_pct_plus <- round(
  pve_plus[1] * 100,
  1
)

pc2_pct_plus <- round(
  pve_plus[2] * 100,
  1
)

# Convert PCA results to data frame
pca_df_plus <- as.data.frame(
  pca_plus$x[, 1:2]
)

pca_df_plus$sample <- rownames(
  pca_df_plus
)

# Add condition information
pca_df_plus$condition <- coldata[
  pca_df_plus$sample,
  "condition"
]

# PCA plot
PCA_plus <- ggplot(
  pca_df_plus,
  aes(
    PC1,
    PC2,
    color = condition
  )
) +
  geom_point(
    shape = 17,
    size = 3.5,
    alpha = 0.9
  ) +
  scale_color_manual(
    values = c(
      "NMO_plus" = "red",
      "HC_plus" = "blue"
    ),
    name = "Condition"
  ) +
  theme_classic() +
  labs(
    title = "Principal component analysis in CD25hi NB cells",
    x = paste0(
      "PC1 (",
      pc1_pct_plus,
      "%)"
    ),
    y = paste0(
      "PC2 (",
      pc2_pct_plus,
      "%)"
    )
  )

print(PCA_plus)

ggsave(
  filename = file.path(
    output_dir,
    "PCA_NMO_plus_vs_HC_plus.pdf"
  ),
  plot = PCA_plus,
  width = 8,
  height = 6,
  device = "pdf"
)


###################################################
### PCA: NMO_minus vs HC_minus ----

# Select target conditions
target_conditions_minus <- c(
  "NMO_minus",
  "HC_minus"
)

# Select target samples
target_samples_minus <- rownames(coldata)[
  coldata$condition %in% target_conditions_minus
]

# Subset VST matrix
expr_mat_subset_minus <- expr_mat[
  ,
  target_samples_minus,
  drop = FALSE
]

# Extract selected TF genes
tf_human_in_data_minus <- intersect(
  rownames(expr_mat_subset_minus),
  tf_human
)

vst_tf_minus <- expr_mat_subset_minus[
  tf_human_in_data_minus,
  ,
  drop = FALSE
]

# Z-score normalization
vst_tf_z_minus <- scale_rows(
  vst_tf_minus
)

# PCA
pca_minus <- prcomp(
  t(vst_tf_z_minus),
  scale. = FALSE
)

# Proportion of variance explained
pve_minus <- (
  pca_minus$sdev^2
) / sum(
  pca_minus$sdev^2
)

pc1_pct_minus <- round(
  pve_minus[1] * 100,
  1
)

pc2_pct_minus <- round(
  pve_minus[2] * 100,
  1
)

# Convert PCA results to data frame
pca_df_minus <- as.data.frame(
  pca_minus$x[, 1:2]
)

pca_df_minus$sample <- rownames(
  pca_df_minus
)

# Add condition information
pca_df_minus$condition <- coldata[
  pca_df_minus$sample,
  "condition"
]

# PCA plot
PCA_minus <- ggplot(
  pca_df_minus,
  aes(
    PC1,
    PC2,
    color = condition
  )
) +
  geom_point(
    shape = 17,
    size = 3.5,
    alpha = 0.9
  ) +
  scale_color_manual(
    values = c(
      "NMO_minus" = "red",
      "HC_minus" = "blue"
    ),
    name = "Condition"
  ) +
  theme_classic() +
  labs(
    title = "Principal component analysis in CD25− NB cells",
    x = paste0(
      "PC1 (",
      pc1_pct_minus,
      "%)"
    ),
    y = paste0(
      "PC2 (",
      pc2_pct_minus,
      "%)"
    )
  )

print(PCA_minus)

ggsave(
  filename = file.path(
    output_dir,
    "PCA_NMO_minus_vs_HC_minus.pdf"
  ),
  plot = PCA_minus,
  width = 8,
  height = 6,
  device = "pdf"
)