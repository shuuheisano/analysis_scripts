###################################################
### Volcano plots
###################################################

library(ggplot2)
library(ggrepel)
library(ggtext)
library(msigdbr)
library(dplyr)


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

# DESeq2 result tables
plus_result_file <- "path/to/NMO_plus_vs_HC_plus_all_DESeq2_results.csv"

minus_result_file <- "path/to/NMO_minus_vs_HC_minus_all_DESeq2_results.csv"

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
fc_cutoff <- log2(1.5)


###################################################
### Hallmark interferon-alpha response gene set ----

hallmark_ifna <- msigdbr::msigdbr(
  species = "Homo sapiens",
  collection = "H"
) |>
  dplyr::filter(
    gs_name == "HALLMARK_INTERFERON_ALPHA_RESPONSE"
  ) |>
  dplyr::pull(gene_symbol) |>
  unique()


###################################################
### Genes selected for labeling ----

label_genes <- unique(
  c(
    hallmark_ifna,
    "PRDM1",
    "IFIT1",
    "MX1",
    "OAS3",
    "RPS27A",
    "IFITM1",
    "SOCS1",
    "IFITM2",
    "FANCA",
    "FKBP5",
    "FLNB",
    "PLCG1",
    "TARBP2",
    "BATF",
    "FOXP1",
    "IL6ST",
    "RELB",
    "CBLB",
    "HCK",
    "SRC",
    "CDK6",
    "TRAF6"
  )
)


###################################################
### Volcano plot: NMO_plus vs HC_plus ----

# Import DESeq2 results
res_plus <- read.csv(
  plus_result_file,
  header = TRUE,
  check.names = FALSE
)

# Remove genes with missing or zero adjusted P-values
res_plus_plot <- res_plus |>
  dplyr::filter(
    !is.na(padj),
    padj > 0
  )

# Classify genes according to significance and fold change
res_plus_plot <- res_plus_plot |>
  dplyr::mutate(
    sig = dplyr::case_when(
      padj < padj_cutoff &
        log2FoldChange > fc_cutoff ~ "Up",

      padj < padj_cutoff &
        log2FoldChange < -fc_cutoff ~ "Down",

      TRUE ~ "NS"
    )
  )

# Prepare genes selected for labeling
label_data_plus <- res_plus_plot |>
  dplyr::filter(
    padj < padj_cutoff,
    abs(log2FoldChange) > fc_cutoff,
    Symbol %in% label_genes
  ) |>
  dplyr::mutate(
    nudge_x_value = ifelse(
      log2FoldChange > 0,
      1.2,
      -1.2
    )
  )

# Draw volcano plot
p_plus <- ggplot2::ggplot(
  res_plus_plot,
  ggplot2::aes(
    x = log2FoldChange,
    y = -log10(padj),
    color = sig
  )
) +
  ggplot2::geom_point(
    alpha = 0.7,
    size = 1.8
  ) +
  ggplot2::scale_color_manual(
    values = c(
      "Up" = "red3",
      "Down" = "blue3",
      "NS" = "grey70"
    )
  ) +
  ggplot2::geom_vline(
    xintercept = c(
      -fc_cutoff,
      fc_cutoff
    ),
    linetype = "dashed",
    color = "black"
  ) +
  ggplot2::geom_hline(
    yintercept = -log10(padj_cutoff),
    linetype = "dashed",
    color = "black"
  ) +
  ggrepel::geom_text_repel(
    data = label_data_plus,
    ggplot2::aes(
      label = Symbol
    ),
    nudge_x = label_data_plus$nudge_x_value,
    size = 3.5,
    box.padding = 0.5,
    point.padding = 0.3,
    force = 3,
    force_pull = 0.5,
    segment.color = "grey40",
    segment.size = 0.3,
    min.segment.length = 0,
    max.overlaps = 30,
    max.iter = 20000,
    seed = 123,
    show.legend = FALSE
  ) +
  ggplot2::theme_classic() +
  ggplot2::labs(
    title = paste0(
      "<span style='color:red;'>NMOSD CD25hi NB</span>",
      " vs ",
      "<span style='color:blue;'>HC CD25hi NB</span>"
    ),
    x = expression(log[2] * " fold change"),
    y = expression(-log[10] * " FDR"),
    color = NULL
  ) +
  ggplot2::theme(
    plot.title = ggtext::element_markdown(
      hjust = 0.5
    ),
    legend.position = "right"
  )

print(p_plus)

# Export volcano plot
ggplot2::ggsave(
  filename = file.path(
    output_dir,
    "Volcano_plot_NMO_plus_vs_HC_plus.pdf"
  ),
  plot = p_plus,
  width = 9,
  height = 7,
  device = "pdf"
)


###################################################
### Volcano plot: NMO_minus vs HC_minus ----

# Import DESeq2 results
res_minus <- read.csv(
  minus_result_file,
  header = TRUE,
  check.names = FALSE
)

# Remove genes with missing or zero adjusted P-values
res_minus_plot <- res_minus |>
  dplyr::filter(
    !is.na(padj),
    padj > 0
  )

# Classify genes according to significance and fold change
res_minus_plot <- res_minus_plot |>
  dplyr::mutate(
    sig = dplyr::case_when(
      padj < padj_cutoff &
        log2FoldChange > fc_cutoff ~ "Up",

      padj < padj_cutoff &
        log2FoldChange < -fc_cutoff ~ "Down",

      TRUE ~ "NS"
    )
  )

# Prepare genes selected for labeling
label_data_minus <- res_minus_plot |>
  dplyr::filter(
    padj < padj_cutoff,
    abs(log2FoldChange) > fc_cutoff,
    Symbol %in% label_genes
  ) |>
  dplyr::mutate(
    nudge_x_value = ifelse(
      log2FoldChange > 0,
      1.2,
      -1.2
    )
  )

# Draw volcano plot
p_minus <- ggplot2::ggplot(
  res_minus_plot,
  ggplot2::aes(
    x = log2FoldChange,
    y = -log10(padj),
    color = sig
  )
) +
  ggplot2::geom_point(
    alpha = 0.7,
    size = 1.8
  ) +
  ggplot2::scale_color_manual(
    values = c(
      "Up" = "red3",
      "Down" = "blue3",
      "NS" = "grey70"
    )
  ) +
  ggplot2::geom_vline(
    xintercept = c(
      -fc_cutoff,
      fc_cutoff
    ),
    linetype = "dashed",
    color = "black"
  ) +
  ggplot2::geom_hline(
    yintercept = -log10(padj_cutoff),
    linetype = "dashed",
    color = "black"
  ) +
  ggrepel::geom_text_repel(
    data = label_data_minus,
    ggplot2::aes(
      label = Symbol
    ),
    nudge_x = label_data_minus$nudge_x_value,
    size = 3.5,
    box.padding = 0.5,
    point.padding = 0.3,
    force = 3,
    force_pull = 0.5,
    segment.color = "grey40",
    segment.size = 0.3,
    min.segment.length = 0,
    max.overlaps = 30,
    max.iter = 20000,
    seed = 123,
    show.legend = FALSE
  ) +
  ggplot2::theme_classic() +
  ggplot2::labs(
    title = paste0(
      "<span style='color:red;'>NMOSD CD25− NB</span>",
      " vs ",
      "<span style='color:blue;'>HC CD25− NB</span>"
    ),
    x = expression(log[2] * " fold change"),
    y = expression(-log[10] * " FDR"),
    color = NULL
  ) +
  ggplot2::theme(
    plot.title = ggtext::element_markdown(
      hjust = 0.5
    ),
    legend.position = "right"
  )

print(p_minus)

# Export volcano plot
ggplot2::ggsave(
  filename = file.path(
    output_dir,
    "Volcano_plot_NMO_minus_vs_HC_minus.pdf"
  ),
  plot = p_minus,
  width = 9,
  height = 7,
  device = "pdf"
)