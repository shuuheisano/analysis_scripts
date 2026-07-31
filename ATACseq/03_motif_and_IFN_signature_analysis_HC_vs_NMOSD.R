############################################################
# 03_motif_and_IFN_signature_analysis_HC_vs_NMOSD.R
#
# chromVAR analysis of primary-cell ATAC-seq data from
# healthy controls (HC) and patients with NMOSD
#
# Workflow:
# 1. Import consensus peaks and BAM files
# 2. Construct a chromVAR count matrix
# 3. Remove overlapping peaks and correct for GC bias
# 4. Calculate motif deviation Z-scores
# 5. Perform differential motif analysis between HC and NMOSD
# 6. Calculate a core IFN-associated motif activity score
# 7. Generate the core IFN boxplot and motif heatmap
############################################################

rm(list = ls())

library(chromVAR)
library(motifmatchr)
library(SummarizedExperiment)
library(S4Vectors)
library(BSgenome.Hsapiens.UCSC.hg38)
library(limma)
library(dplyr)
library(ggplot2)
library(pheatmap)

PROJECT_DIR <- "/path/to/ATACseq_project"

PEAK_FILE <- file.path(
  PROJECT_DIR,
  "data",
  "processed",
  "macs2",
  "merged_peaks_woNMO5_500bp.bed"
)

BAM_DIR <- file.path(
  PROJECT_DIR,
  "data",
  "processed",
  "bwa"
)

OUTPUT_DIR <- file.path(
  PROJECT_DIR,
  "results",
  "motif_and_IFN_signature"
)

dir.create(
  OUTPUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

peaks <- getPeaks(
  PEAK_FILE,
  sort_peaks = TRUE
)

bamfiles <- file.path(
  BAM_DIR,
  c(
    "HC4.bam",
    "HC5.bam",
    "HC6.bam",
    "HC7.bam",
    "HC8.bam",
    "NMO3.bam",
    "NMO4.bam",
    "NMO7.bam"
  )
)

missing_files <- c(PEAK_FILE, bamfiles)[!file.exists(c(PEAK_FILE, bamfiles))]

if (length(missing_files) > 0) {
  stop(
    paste(
      "The following input files were not found:",
      paste(missing_files, collapse = ", ")
    )
  )
}

sample_id <- c(
  "HC4", "HC5", "HC6", "HC7", "HC8",
  "NMO3", "NMO4", "NMO7"
)

condition <- factor(
  c("HC", "HC", "HC", "HC", "HC", "NMO", "NMO", "NMO"),
  levels = c("HC", "NMO")
)

sample_metadata <- DataFrame(
  sample_id = sample_id,
  condition = condition
)

rownames(sample_metadata) <- basename(bamfiles)

counts <- getCounts(
  bamfiles,
  peaks,
  paired = TRUE,
  format = "bam",
  colData = sample_metadata
)

if (!identical(rownames(colData(counts)), colnames(counts))) {
  stop("The sample metadata do not match the count-matrix columns.")
}

counts_filtered <- filterPeaks(
  counts,
  non_overlapping = TRUE
)

counts_filtered <- addGCBias(
  counts_filtered,
  genome = BSgenome.Hsapiens.UCSC.hg38
)

motifs <- getJasparMotifs()

motif_matches <- matchMotifs(
  motifs,
  counts_filtered,
  genome = BSgenome.Hsapiens.UCSC.hg38
)

deviations <- computeDeviations(
  object = counts_filtered,
  annotations = motif_matches
)

deviation_z <- deviationScores(deviations)

write.table(
  data.frame(
    motif = rownames(deviation_z),
    deviation_z,
    check.names = FALSE
  ),
  file = file.path(
    OUTPUT_DIR,
    "HC_vs_NMOSD_chromVAR_motif_deviation_Z_scores.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

metadata_df <- as.data.frame(colData(deviations))

metadata_df$condition <- factor(
  metadata_df$condition,
  levels = c("HC", "NMO")
)

design <- model.matrix(
  ~ condition,
  data = metadata_df
)

fit <- lmFit(
  deviation_z,
  design
)

fit <- eBayes(fit)

motif_results <- topTable(
  fit,
  coef = "conditionNMO",
  number = Inf,
  sort.by = "P"
)

motif_results$motif <- rownames(motif_results)

motif_results <- motif_results %>%
  select(
    motif,
    logFC,
    AveExpr,
    t,
    P.Value,
    adj.P.Val,
    B
  )

write.table(
  motif_results,
  file = file.path(
    OUTPUT_DIR,
    "HC_vs_NMOSD_all_motif_differential_analysis.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

clean_motif_label <- function(x) {
  sub("^[^_]+_", "", x)
}

motif_labels <- rownames(deviation_z)
clean_labels <- clean_motif_label(motif_labels)

idx_stat1_stat2 <- grep(
  "STAT1::STAT2",
  clean_labels,
  ignore.case = TRUE
)

idx_stat1_only <- grep(
  "(^|[^A-Z0-9])STAT1([^A-Z0-9]|$)",
  clean_labels,
  ignore.case = TRUE
)

idx_stat1_only <- setdiff(
  idx_stat1_only,
  idx_stat1_stat2
)

idx_stat1_only <- idx_stat1_only[
  !grepl(
    "::",
    clean_labels[idx_stat1_only]
  )
]

idx_irf9 <- grep(
  "(^|[^A-Z0-9])IRF9([^A-Z0-9]|$)",
  clean_labels,
  ignore.case = TRUE
)

core_ifn_index <- sort(
  unique(
    c(
      idx_stat1_stat2,
      idx_irf9,
      idx_stat1_only
    )
  )
)

if (length(core_ifn_index) == 0) {
  stop(
    paste0(
      "No core IFN-associated motifs were identified. ",
      "Check the JASPAR motif labels."
    )
  )
}

core_ifn_motifs <- data.frame(
  motif = motif_labels[core_ifn_index],
  motif_label = clean_labels[core_ifn_index]
)

print(core_ifn_motifs)

write.table(
  core_ifn_motifs,
  file = file.path(
    OUTPUT_DIR,
    "HC_vs_NMOSD_core_IFN_motifs_used_for_signature.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

core_ifn_results <- motif_results %>%
  filter(motif %in% core_ifn_motifs$motif) %>%
  left_join(
    core_ifn_motifs,
    by = "motif"
  ) %>%
  select(
    motif,
    motif_label,
    everything()
  )

write.table(
  core_ifn_results,
  file = file.path(
    OUTPUT_DIR,
    "HC_vs_NMOSD_core_IFN_motif_differential_analysis.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

core_ifn_score <- colMeans(
  deviation_z[
    core_ifn_index,
    ,
    drop = FALSE
  ],
  na.rm = TRUE
)

ifn_score_df <- data.frame(
  sample = colnames(deviation_z),
  IFNcore = as.numeric(core_ifn_score),
  condition = metadata_df$condition
)

write.table(
  ifn_score_df,
  file = file.path(
    OUTPUT_DIR,
    "HC_vs_NMOSD_core_IFN_signature_scores.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

ifn_test <- wilcox.test(
  IFNcore ~ condition,
  data = ifn_score_df,
  paired = FALSE,
  exact = FALSE
)

ifn_test_results <- data.frame(
  comparison = "NMOSD vs HC",
  test = "Wilcoxon rank-sum test",
  n_HC = sum(ifn_score_df$condition == "HC"),
  n_NMOSD = sum(ifn_score_df$condition == "NMO"),
  statistic = unname(ifn_test$statistic),
  p_value = ifn_test$p.value
)

write.table(
  ifn_test_results,
  file = file.path(
    OUTPUT_DIR,
    "HC_vs_NMOSD_core_IFN_signature_test.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

ifn_p_label <- if (ifn_test$p.value < 0.001) {
  "P < 0.001"
} else {
  paste0(
    "P = ",
    signif(ifn_test$p.value, 3)
  )
}

ifn_y_range <- range(
  ifn_score_df$IFNcore,
  na.rm = TRUE
)

ifn_y_span <- diff(ifn_y_range)

if (!is.finite(ifn_y_span) || ifn_y_span == 0) {
  ifn_y_span <- 1
}

p_ifn <- ggplot(
  ifn_score_df,
  aes(
    x = condition,
    y = IFNcore,
    fill = condition
  )
) +
  geom_boxplot(
    width = 0.55,
    outlier.shape = NA,
    alpha = 0.85
  ) +
  geom_jitter(
    width = 0.10,
    size = 2.2,
    alpha = 0.8
  ) +
  annotate(
    "segment",
    x = 1,
    xend = 2,
    y = ifn_y_range[2] + ifn_y_span * 0.12,
    yend = ifn_y_range[2] + ifn_y_span * 0.12
  ) +
  annotate(
    "text",
    x = 1.5,
    y = ifn_y_range[2] + ifn_y_span * 0.18,
    label = ifn_p_label,
    size = 3.5
  ) +
  labs(
    title = "Core IFN-associated motif activity",
    subtitle = "STAT1::STAT2 / IRF9 / STAT1",
    x = NULL,
    y = "Mean chromVAR deviation Z-score"
  ) +
  scale_fill_manual(
    values = c(
      "HC" = "grey",
      "NMO" = "red"
    )
  ) +
  coord_cartesian(
    ylim = c(
      ifn_y_range[1] - ifn_y_span * 0.05,
      ifn_y_range[2] + ifn_y_span * 0.25
    )
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 11),
    plot.subtitle = element_text(size = 9)
  )

print(p_ifn)

ggsave(
  filename = file.path(
    OUTPUT_DIR,
    "HC_vs_NMOSD_core_IFN_signature_boxplot.pdf"
  ),
  plot = p_ifn,
  width = 6,
  height = 6
)

core_ifn_matrix <- as.matrix(
  deviation_z[
    core_ifn_index,
    ,
    drop = FALSE
  ]
)

heatmap_annotation <- data.frame(
  condition = metadata_df$condition
)

rownames(heatmap_annotation) <- colnames(core_ifn_matrix)

if (!identical(
  rownames(heatmap_annotation),
  colnames(core_ifn_matrix)
)) {
  stop(
    "The heatmap annotation does not match the motif deviation matrix."
  )
}

pdf(
  file.path(
    OUTPUT_DIR,
    "HC_vs_NMOSD_core_IFN_motif_heatmap.pdf"
  ),
  width = 8,
  height = 6,
  useDingbats = FALSE
)

pheatmap(
  core_ifn_matrix,
  scale = "none",
  annotation_col = heatmap_annotation,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  labels_row = clean_labels[core_ifn_index],
  main = "Core IFN-related motif deviations"
)

dev.off()

writeLines(
  capture.output(sessionInfo()),
  con = file.path(
    OUTPUT_DIR,
    "HC_vs_NMOSD_chromVAR_sessionInfo.txt"
  )
)
