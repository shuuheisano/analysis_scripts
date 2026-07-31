###################################################
### Ranked gene set enrichment analysis using fgsea
###################################################

library(fgsea)
library(msigdbr)
library(dplyr)
library(ggplot2)


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
### Output path ----

# Output directory
output_dir <- "path/to/output_directory"

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


###################################################
### Ranking: NMO_plus vs HC_plus ----

# Convert DESeq2 results to a data frame
rank_df_plus <- as.data.frame(deseq2Results_1)

# Add gene symbols as a column
rank_df_plus$Symbol <- rownames(rank_df_plus)

# Remove rows with missing statistics, missing gene symbols,
# or empty gene symbols
rank_df_plus <- rank_df_plus[
  !is.na(rank_df_plus$stat) &
    !is.na(rank_df_plus$Symbol) &
    rank_df_plus$Symbol != "",
  ,
  drop = FALSE
]

# If duplicated gene symbols are present, retain the row
# with the largest absolute Wald statistic
rank_df_plus <- rank_df_plus |>
  dplyr::group_by(Symbol) |>
  dplyr::slice_max(
    order_by = abs(stat),
    n = 1,
    with_ties = FALSE
  ) |>
  dplyr::ungroup()

# Create a named numeric vector for fgsea
ranks_plus <- rank_df_plus$stat
names(ranks_plus) <- rank_df_plus$Symbol

# Sort from the largest to the smallest statistic
ranks_plus <- sort(
  ranks_plus,
  decreasing = TRUE
)

# Check ranking vector
length(ranks_plus)
head(ranks_plus)
tail(ranks_plus)

# Export ranking data
rank_export_plus <- data.frame(
  Symbol = names(ranks_plus),
  rank_metric = as.numeric(ranks_plus)
)

write.csv(
  rank_export_plus,
  file = file.path(
    output_dir,
    "GSEA_rank_NMO_plus_vs_HC_plus.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)


###################################################
### Hallmark gene sets ----

msig_hallmark <- msigdbr::msigdbr(
  species = "Homo sapiens",
  collection = "H"
)

# Check the MSigDB version
unique(msig_hallmark$db_version)

# Convert Hallmark gene sets to a named list
hallmark_pathways <- split(
  x = msig_hallmark$gene_symbol,
  f = msig_hallmark$gs_name
)

# Remove duplicated genes within each pathway
hallmark_pathways <- lapply(
  hallmark_pathways,
  unique
)


###################################################
### Run fgsea: Hallmark gene sets ----

set.seed(123)

fgsea_plus_hallmark <- fgsea::fgsea(
  pathways = hallmark_pathways,
  stats = ranks_plus,
  minSize = 15,
  maxSize = 500,
  scoreType = "std"
)

fgsea_plus_hallmark <- as.data.frame(
  fgsea_plus_hallmark
)

# Sort by adjusted P-value and absolute NES
fgsea_plus_hallmark <- fgsea_plus_hallmark[
  order(
    fgsea_plus_hallmark$padj,
    -abs(fgsea_plus_hallmark$NES)
  ),
  ,
  drop = FALSE
]

# Convert the leading-edge list column to text
fgsea_plus_hallmark$leadingEdge <- vapply(
  fgsea_plus_hallmark$leadingEdge,
  paste,
  collapse = ";",
  FUN.VALUE = character(1)
)

# Check results
head(fgsea_plus_hallmark)

# Export all Hallmark GSEA results
write.csv(
  fgsea_plus_hallmark,
  file = file.path(
    output_dir,
    "GSEA_Hallmark_NMO_plus_vs_HC_plus.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)


###################################################
### Significant Hallmark gene sets ----

fgsea_plus_hallmark_sig <- fgsea_plus_hallmark[
  !is.na(fgsea_plus_hallmark$padj) &
    fgsea_plus_hallmark$padj < 0.05,
  ,
  drop = FALSE
]

# Check significant results
nrow(fgsea_plus_hallmark_sig)

fgsea_plus_hallmark_sig[
  ,
  c(
    "pathway",
    "NES",
    "pval",
    "padj",
    "size"
  )
]

# Export significant Hallmark GSEA results
write.csv(
  fgsea_plus_hallmark_sig,
  file = file.path(
    output_dir,
    "GSEA_Hallmark_significant_NMO_plus_vs_HC_plus.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)


###################################################
### Enrichment plot ----

plot_pathway <- "HALLMARK_INTERFERON_ALPHA_RESPONSE"

# Extract NES and adjusted P-value for the selected pathway
plot_result <- fgsea_plus_hallmark[
  fgsea_plus_hallmark$pathway == plot_pathway,
  ,
  drop = FALSE
]

plot_nes <- plot_result$NES
plot_padj <- plot_result$padj

# Generate enrichment plot
p_ifna_plus <- fgsea::plotEnrichment(
  pathway = hallmark_pathways[[plot_pathway]],
  stats = ranks_plus
) +
  ggplot2::labs(
    title = plot_pathway,
    subtitle = paste0(
      "NMO_plus vs HC_plus; NES = ",
      round(plot_nes, 2),
      ", adjusted P = ",
      signif(plot_padj, 3)
    )
  )

p_ifna_plus

# Save enrichment plot
ggplot2::ggsave(
  filename = file.path(
    output_dir,
    "GSEA_Hallmark_IFN_alpha_NMO_plus_vs_HC_plus.pdf"
  ),
  plot = p_ifna_plus,
  width = 6,
  height = 5
)

###################################################
### Ranking: NMO_minus vs HC_minus ----

# Convert DESeq2 results to a data frame
rank_df_minus <- as.data.frame(deseq2Results_2)

# Add gene symbols as a column
rank_df_minus$Symbol <- rownames(rank_df_minus)

# Remove rows with missing statistics, missing gene symbols,
# or empty gene symbols
rank_df_minus <- rank_df_minus[
  !is.na(rank_df_minus$stat) &
    !is.na(rank_df_minus$Symbol) &
    rank_df_minus$Symbol != "",
  ,
  drop = FALSE
]

# If duplicated gene symbols are present, retain the row
# with the largest absolute Wald statistic
rank_df_minus <- rank_df_minus |>
  dplyr::group_by(Symbol) |>
  dplyr::slice_max(
    order_by = abs(stat),
    n = 1,
    with_ties = FALSE
  ) |>
  dplyr::ungroup()

# Create a named numeric vector for fgsea
ranks_minus <- rank_df_minus$stat
names(ranks_minus) <- rank_df_minus$Symbol

# Sort from the largest to the smallest statistic
ranks_minus <- sort(
  ranks_minus,
  decreasing = TRUE
)

# Check ranking vector
length(ranks_minus)
head(ranks_minus)
tail(ranks_minus)

# Export ranking data
rank_export_minus <- data.frame(
  Symbol = names(ranks_minus),
  rank_metric = as.numeric(ranks_minus)
)

write.csv(
  rank_export_minus,
  file = file.path(
    output_dir,
    "GSEA_rank_NMO_minus_vs_HC_minus.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)


###################################################
### Run fgsea: Hallmark gene sets
### NMO_minus vs HC_minus ----

set.seed(123)

fgsea_minus_hallmark <- fgsea::fgsea(
  pathways = hallmark_pathways,
  stats = ranks_minus,
  minSize = 15,
  maxSize = 500,
  scoreType = "std"
)

fgsea_minus_hallmark <- as.data.frame(
  fgsea_minus_hallmark
)

# Sort by adjusted P-value and absolute NES
fgsea_minus_hallmark <- fgsea_minus_hallmark[
  order(
    fgsea_minus_hallmark$padj,
    -abs(fgsea_minus_hallmark$NES)
  ),
  ,
  drop = FALSE
]

# Convert the leading-edge list column to text
fgsea_minus_hallmark$leadingEdge <- vapply(
  fgsea_minus_hallmark$leadingEdge,
  paste,
  collapse = ";",
  FUN.VALUE = character(1)
)

# Check results
head(fgsea_minus_hallmark)

# Export all Hallmark GSEA results
write.csv(
  fgsea_minus_hallmark,
  file = file.path(
    output_dir,
    "GSEA_Hallmark_NMO_minus_vs_HC_minus.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)


###################################################
### Significant Hallmark gene sets ----

fgsea_minus_hallmark_sig <- fgsea_minus_hallmark[
  !is.na(fgsea_minus_hallmark$padj) &
    fgsea_minus_hallmark$padj < 0.05,
  ,
  drop = FALSE
]

# Check significant results
nrow(fgsea_minus_hallmark_sig)

fgsea_minus_hallmark_sig[
  ,
  c(
    "pathway",
    "NES",
    "pval",
    "padj",
    "size"
  )
]

# Export significant Hallmark GSEA results
write.csv(
  fgsea_minus_hallmark_sig,
  file = file.path(
    output_dir,
    "GSEA_Hallmark_significant_NMO_minus_vs_HC_minus.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)


###################################################
### IFN-alpha enrichment plot:
### NMO_minus vs HC_minus ----

plot_pathway_minus <- "HALLMARK_INTERFERON_ALPHA_RESPONSE"

# Extract GSEA results for the selected pathway
ifna_minus_result <- fgsea_minus_hallmark[
  fgsea_minus_hallmark$pathway == plot_pathway_minus,
  ,
  drop = FALSE
]

# Check selected pathway result
ifna_minus_result[
  ,
  c(
    "pathway",
    "NES",
    "pval",
    "padj",
    "size",
    "leadingEdge"
  )
]

# Generate enrichment plot
p_ifna_minus <- fgsea::plotEnrichment(
  pathway = hallmark_pathways[[plot_pathway_minus]],
  stats = ranks_minus
) +
  ggplot2::labs(
    title = plot_pathway_minus,
    subtitle = paste0(
      "NMO_minus vs HC_minus; NES = ",
      sprintf("%.2f", ifna_minus_result$NES),
      ", adjusted P = ",
      formatC(
        ifna_minus_result$padj,
        format = "g",
        digits = 3
      )
    )
  )

p_ifna_minus

# Save enrichment plot
ggplot2::ggsave(
  filename = file.path(
    output_dir,
    "GSEA_Hallmark_IFN_alpha_NMO_minus_vs_HC_minus.pdf"
  ),
  plot = p_ifna_minus,
  width = 6,
  height = 5
)