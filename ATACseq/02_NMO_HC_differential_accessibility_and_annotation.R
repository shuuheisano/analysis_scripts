############################################################
# 02_NMO_HC_differential_accessibility_and_annotation_primary.R
#
# Differential accessibility analysis and peak annotation
# for primary-cell ATAC-seq data
#
# Workflow:
# 1. Register MACS2 peak files and BAM files in DiffBind
# 2. Remove ENCODE blacklist regions
# 3. Count reads within consensus peaks
# 4. Perform quality-control visualization
# 5. Normalize read counts
# 6. Perform differential accessibility analysis
# 7. Annotate significant differentially accessible regions
# 8. Optionally visualize accessibility at a selected PRDM1 region
############################################################


############################################################
# 1. Load required packages
############################################################

rm(list = ls())

library(DiffBind)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(GenomicRanges)
library(IRanges)
library(dplyr)
library(tidyr)
library(ggplot2)


############################################################
# 2. User-defined paths
############################################################

# Root directory of the ATAC-seq project
PROJECT_DIR <- "/path/to/ATACseq_project"

# Input directories
MACS2_DIR <- file.path(
  PROJECT_DIR,
  "data",
  "processed",
  "macs2"
)

BAM_DIR <- file.path(
  PROJECT_DIR,
  "data",
  "processed",
  "bwa"
)

# Output directory
OUTPUT_DIR <- file.path(
  PROJECT_DIR,
  "results",
  "differential_accessibility"
)

dir.create(
  OUTPUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)


############################################################
# 3. Register peak and BAM files in DiffBind
############################################################

HC_vs_NMO <- NULL


# HC4
HC_vs_NMO <- dba.peakset(
  HC_vs_NMO,
  peaks = file.path(
    MACS2_DIR,
    "HC4_peak_peaks.xls"
  ),
  peak.caller = "macs",
  sampID = "HC4",
  condition = "HC",
  replicate = 1,
  bamReads = file.path(
    BAM_DIR,
    "HC4.bam"
  )
)


# HC5
HC_vs_NMO <- dba.peakset(
  HC_vs_NMO,
  peaks = file.path(
    MACS2_DIR,
    "HC5_peak_peaks.xls"
  ),
  peak.caller = "macs",
  sampID = "HC5",
  condition = "HC",
  replicate = 2,
  bamReads = file.path(
    BAM_DIR,
    "HC5.bam"
  )
)


# HC6
HC_vs_NMO <- dba.peakset(
  HC_vs_NMO,
  peaks = file.path(
    MACS2_DIR,
    "HC6_peak_peaks.xls"
  ),
  peak.caller = "macs",
  sampID = "HC6",
  condition = "HC",
  replicate = 3,
  bamReads = file.path(
    BAM_DIR,
    "HC6.bam"
  )
)


# HC7
HC_vs_NMO <- dba.peakset(
  HC_vs_NMO,
  peaks = file.path(
    MACS2_DIR,
    "HC7_peak_peaks.xls"
  ),
  peak.caller = "macs",
  sampID = "HC7",
  condition = "HC",
  replicate = 4,
  bamReads = file.path(
    BAM_DIR,
    "HC7.bam"
  )
)


# HC8
HC_vs_NMO <- dba.peakset(
  HC_vs_NMO,
  peaks = file.path(
    MACS2_DIR,
    "HC8_peak_peaks.xls"
  ),
  peak.caller = "macs",
  sampID = "HC8",
  condition = "HC",
  replicate = 5,
  bamReads = file.path(
    BAM_DIR,
    "HC8.bam"
  )
)


# NMO3
HC_vs_NMO <- dba.peakset(
  HC_vs_NMO,
  peaks = file.path(
    MACS2_DIR,
    "NMO3_peak_peaks.xls"
  ),
  peak.caller = "macs",
  sampID = "NMO3",
  condition = "NMO",
  replicate = 1,
  bamReads = file.path(
    BAM_DIR,
    "NMO3.bam"
  )
)


# NMO4
HC_vs_NMO <- dba.peakset(
  HC_vs_NMO,
  peaks = file.path(
    MACS2_DIR,
    "NMO4_peak_peaks.xls"
  ),
  peak.caller = "macs",
  sampID = "NMO4",
  condition = "NMO",
  replicate = 2,
  bamReads = file.path(
    BAM_DIR,
    "NMO4.bam"
  )
)


# NMO7
HC_vs_NMO <- dba.peakset(
  HC_vs_NMO,
  peaks = file.path(
    MACS2_DIR,
    "NMO7_peak_peaks.xls"
  ),
  peak.caller = "macs",
  sampID = "NMO7",
  condition = "NMO",
  replicate = 3,
  bamReads = file.path(
    BAM_DIR,
    "NMO7.bam"
  )
)


# Display registered samples
dba.show(HC_vs_NMO)


############################################################
# 4. Apply the ENCODE blacklist
############################################################

HC_vs_NMO <- dba.blacklist(
  HC_vs_NMO,
  blacklist = TRUE,
  greylist = FALSE
)


############################################################
# 5. Count reads in consensus peak regions
############################################################

HC_vs_NMO_counts <- dba.count(
  HC_vs_NMO,
  summits = 200,
  filter = 0,
  bRemoveDuplicates = FALSE,
  mapQCth = 15,
  bParallel = TRUE,
  score = DBA_SCORE_READS
)


############################################################
# 6. Quality-control plots
############################################################

# Correlation heatmap
pdf(
  file.path(
    OUTPUT_DIR,
    "HC_vs_NMO_correlation_heatmap.pdf"
  ),
  width = 8,
  height = 6
)

plot(HC_vs_NMO_counts)

dev.off()


# Principal component analysis
pdf(
  file.path(
    OUTPUT_DIR,
    "HC_vs_NMO_PCA.pdf"
  ),
  width = 8,
  height = 6
)

dba.plotPCA(HC_vs_NMO_counts)

dev.off()


############################################################
# 7. Normalization
############################################################

# MA plot before normalization
pdf(
  file.path(
    OUTPUT_DIR,
    "HC_vs_NMO_MA_plot_before_normalization.pdf"
  ),
  width = 8,
  height = 6
)

dba.plotMA(
  HC_vs_NMO_counts,
  bNormalized = FALSE,
  sub = "Non-normalized",
  contrast = list(
    NMO = HC_vs_NMO_counts$masks$NMO,
    HC = HC_vs_NMO_counts$masks$HC
  )
)

dev.off()


# Normalize read counts
HC_vs_NMO_counts_normalized <- dba.normalize(
  HC_vs_NMO_counts
)


# MA plot after normalization
pdf(
  file.path(
    OUTPUT_DIR,
    "HC_vs_NMO_MA_plot_after_normalization.pdf"
  ),
  width = 8,
  height = 6
)

dba.plotMA(
  HC_vs_NMO_counts_normalized,
  sub = "Normalized",
  contrast = list(
    NMO = HC_vs_NMO_counts_normalized$masks$NMO,
    HC = HC_vs_NMO_counts_normalized$masks$HC
  )
)

dev.off()


############################################################
# 8. Differential accessibility analysis
############################################################

# Set HC as the reference condition
HC_vs_NMO_model <- dba.contrast(
  HC_vs_NMO_counts_normalized,
  reorderMeta = list(
    Condition = "HC"
  )
)


# Run differential accessibility analysis
HC_vs_NMO_model <- dba.analyze(
  HC_vs_NMO_model
)


# Display contrast information
dba.show(
  HC_vs_NMO_model,
  bContrasts = TRUE
)


# Obtain significant differentially accessible regions
# using an FDR threshold of 0.05
HC_vs_NMO_DAR <- dba.report(
  HC_vs_NMO_model,
  th = 0.05
)


# Display significant regions
HC_vs_NMO_DAR


############################################################
# 9. Export differential accessibility results
############################################################

HC_vs_NMO_DAR_df <- as.data.frame(
  HC_vs_NMO_DAR
)

write.table(
  HC_vs_NMO_DAR_df,
  file = file.path(
    OUTPUT_DIR,
    "HC_vs_NMO_differentially_accessible_regions.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


############################################################
# 10. Peak annotation with ChIPseeker
############################################################

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene


# Annotate significant differentially accessible regions
peakAnno <- annotatePeak(
  HC_vs_NMO_DAR,
  tssRegion = c(-3000, 3000),
  TxDb = txdb,
  annoDb = "org.Hs.eg.db"
)


# Convert annotation results to a data frame
peakAnno_df <- as.data.frame(
  peakAnno
)


# Preview annotation results
head(peakAnno_df)


# Export annotation results
write.table(
  peakAnno_df,
  file = file.path(
    OUTPUT_DIR,
    "HC_vs_NMO_peak_annotation.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


############################################################
# 11. visualization of a selected PRDM1 region
############################################################

# Obtain all tested peaks, including nonsignificant peaks,
# together with normalized read counts
all_df <- as.data.frame(
  dba.report(
    HC_vs_NMO_model,
    th = 1,
    bCounts = TRUE
  )
)


# Genomic region of interest near PRDM1
roi <- GRanges(
  seqnames = "chr6",
  ranges = IRanges(
    start = 106103667,
    end = 106104067
  )
)


# Convert all tested peaks to GRanges
all_ranges <- GRanges(
  seqnames = all_df$seqnames,
  ranges = IRanges(
    start = all_df$start,
    end = all_df$end
  )
)


# Identify peaks overlapping the selected region
hits <- findOverlaps(
  all_ranges,
  roi,
  ignore.strand = TRUE
)


target_df <- all_df[
  queryHits(hits)[1],
  ,
  drop = FALSE
]


############################################################
# 13. Convert normalized counts to long format
############################################################

sample_cols <- c(
  "NMO3",
  "NMO4",
  "NMO7",
  "HC4",
  "HC5",
  "HC6",
  "HC7",
  "HC8"
)


# Confirm that all expected sample columns are present
missing_sample_cols <- setdiff(
  sample_cols,
  colnames(target_df)
)

if (length(missing_sample_cols) > 0) {

  stop(
    paste(
      "The following sample columns were not found:",
      paste(missing_sample_cols, collapse = ", ")
    )
  )
}


t_long <- target_df %>%
  pivot_longer(
    cols = all_of(sample_cols),
    names_to = "sample",
    values_to = "norm_count"
  ) %>%
  mutate(
    condition = ifelse(
      grepl("^NMO", sample),
      "NMOSD",
      "HC"
    ),
    condition = factor(
      condition,
      levels = c("HC", "NMOSD")
    )
  )


############################################################
# 14. Create a DiffBind statistical label
############################################################

diffbind_label <- sprintf(
  "DiffBind: log2FC = %.2f, FDR = %.3g",
  target_df$Fold[1],
  target_df$FDR[1]
)


############################################################
# 15. Plot normalized accessibility at the selected region
############################################################

p_roi <- ggplot(
  t_long,
  aes(
    x = condition,
    y = norm_count,
    fill = condition
  )
) +
  geom_boxplot(
    width = 0.55,
    outlier.shape = NA,
    alpha = 0.85
  ) +
  geom_jitter(
    width = 0.12,
    size = 2,
    alpha = 0.8
  ) +
  annotate(
    "text",
    x = 1.5,
    y = max(
      t_long$norm_count,
      na.rm = TRUE
    ) * 1.08,
    label = diffbind_label,
    hjust = 0.5,
    size = 3.5
  ) +
  labs(
    title = sprintf(
      paste0(
        "Accessibility at chr6:106,103,667–106,104,067\n",
        "(selected peak: %s:%d–%d)"
      ),
      as.character(target_df$seqnames[1]),
      target_df$start[1],
      target_df$end[1]
    ),
    x = NULL,
    y = "Normalized read counts"
  ) +
  scale_fill_manual(
    values = c(
      "HC" = "grey",
      "NMOSD" = "red"
    )
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    plot.title = element_text(
      size = 10
    )
  )


print(p_roi)


ggsave(
  filename = file.path(
    OUTPUT_DIR,
    "HC_vs_NMO_PRDM1_region_accessibility.pdf"
  ),
  plot = p_roi,
  width = 6,
  height = 6
)
