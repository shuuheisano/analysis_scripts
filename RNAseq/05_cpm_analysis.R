#!/usr/bin/env Rscript

library(tidyverse)
library(skimr)

# Input files
COUNT_FILE <- "gene_count_matrix.csv"
METADATA_FILE <- "metadata.csv"

# Output directories
COUNT_OUTPUT_DIR <- "gene_count_matrix"
CPM_OUTPUT_DIR <- "cpm"

# Create output directories
dir.create(
  COUNT_OUTPUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  CPM_OUTPUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

# Verify the current working directory
base::getwd()


#### R script for CPM ####

# Read the gene count matrix
counts <- read.csv(
  COUNT_FILE,
  header = TRUE,
  check.names = FALSE
)

# Rename the first column from gene_id to Symbol
colnames(counts)[1] <- "Symbol"

# Remove the duplicated part of each gene label
counts$Symbol <- sub("\\|.*", "", counts$Symbol)

# Check the result
head(counts)


# Read metadata
metadata <- read.csv(
  METADATA_FILE,
  header = TRUE,
  check.names = FALSE
)

# Remove ".fastq.gz" from FASTQ names and append "_trim_counts.gtf"
adjusted_metadata <- paste0(
  sub("\\.fastq\\.gz$", "", metadata$fastq),
  "_trim_counts.gtf"
)

# Match count-matrix column names to sample names in metadata
sample_match <- match(
  colnames(counts)[-1],
  adjusted_metadata
)

# Stop if any sample columns cannot be matched
if (any(is.na(sample_match))) {
  stop(
    "The following count-matrix columns could not be matched to metadata: ",
    paste(
      colnames(counts)[-1][is.na(sample_match)],
      collapse = ", "
    )
  )
}

colnames(counts)[-1] <- metadata$sample[sample_match]


# Summarize the count matrix
skimr::skim(counts)

# Check the number of rows corresponding to each Symbol
counts |>
  dplyr::count(Symbol, sort = TRUE)


# Display the last six rows
tail(counts)

# Remove the "Total read count" row
counts <- counts %>%
  filter(Symbol != "Total read count")

# Confirm that the final row was removed
tail(counts)


# Confirm that count values contain no missing values
if (anyNA(counts[, -1])) {
  stop("NA values were found in the count matrix.")
}


# Save the cleaned gene count matrix
write.csv(
  counts,
  file.path(
    COUNT_OUTPUT_DIR,
    "gene_count_matrix.csv"
  ),
  row.names = FALSE
)


# Calculate the total read count for each sample
total_counts <- colSums(counts[, -1])

# Stop if any sample has a library size of zero
if (any(total_counts == 0)) {
  stop("At least one sample has a total count of zero.")
}

# Calculate CPM
cpm_values <- t(
  t(counts[, -1]) /
    total_counts *
    1e6
)

# Create the CPM matrix
cpm_matrix <- data.frame(
  Symbol = counts$Symbol,
  cpm_values,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# Confirm that each sample column sums to 1,000,000
colSums(cpm_matrix[, -1])

# Save the CPM matrix
write.csv(
  cpm_matrix,
  file.path(
    COUNT_OUTPUT_DIR,
    "gene_count_matrix_CPM.csv"
  ),
  row.names = FALSE,
  quote = TRUE
)


#### Comparison of IKZF2 CPM values ####

# Extract the IKZF2 row
IKZF2_cpm <- cpm_matrix[
  cpm_matrix$Symbol == "IKZF2",
]

if (nrow(IKZF2_cpm) == 0) {
  stop("IKZF2 was not found in the CPM matrix.")
}

if (nrow(IKZF2_cpm) > 1) {
  warning(
    "Multiple rows corresponding to IKZF2 were found in the CPM matrix."
  )
}


# Convert the IKZF2 data into long format
IKZF2_long <- IKZF2_cpm %>%
  pivot_longer(
    cols = -Symbol,
    names_to = "sample",
    values_to = "CPM"
  )

# Add condition information from metadata
IKZF2_long <- IKZF2_long %>%
  left_join(
    metadata[, c("sample", "condition")],
    by = "sample"
  )

# Retain the four groups used for comparison
IKZF2_subset <- IKZF2_long %>%
  filter(
    condition %in%
      c(
        "HC_minus",
        "HC_plus",
        "NMO_minus",
        "NMO_plus"
      )
  )

# Relabel the conditions
IKZF2_subset$condition <- factor(
  IKZF2_subset$condition,
  levels = c(
    "HC_minus",
    "HC_plus",
    "NMO_minus",
    "NMO_plus"
  ),
  labels = c(
    "HC_CD25-",
    "HC_CD25hi",
    "NMOSD_CD25-",
    "NMOSD_CD25hi"
  )
)


# Plot IKZF2 CPM values
IKZF2_cpm_plot <- ggplot(
  IKZF2_subset,
  aes(
    x = condition,
    y = CPM,
    fill = condition
  )
) +
  geom_boxplot(
    outlier.shape = NA,
    alpha = 0.5
  ) +
  geom_jitter(
    position = position_jitter(
      width = 0.2,
      seed = 123
    ),
    size = 2,
    alpha = 0.7
  ) +
  scale_fill_manual(
    values = c(
      "HC_CD25-" = "grey70",
      "HC_CD25hi" = "grey40",
      "NMOSD_CD25-" = "pink",
      "NMOSD_CD25hi" = "red3"
    )
  ) +
  theme_bw() +
  labs(
    title = "IKZF2",
    y = "CPM",
    x = ""
  )

IKZF2_cpm_plot


# Save the plot as a PDF
ggsave(
  filename = file.path(
    CPM_OUTPUT_DIR,
    "IKZF2_cpm_boxplot.pdf"
  ),
  plot = IKZF2_cpm_plot,
  width = 4,
  height = 4,
  device = "pdf"
)

#### Comparison of TBX21 CPM values ####

# Extract the TBX21 row
TBX21_cpm <- cpm_matrix[
  cpm_matrix$Symbol == "TBX21",
]

if (nrow(TBX21_cpm) == 0) {
  stop("TBX21 was not found in the CPM matrix.")
}

if (nrow(TBX21_cpm) > 1) {
  warning(
    "Multiple rows corresponding to TBX21 were found in the CPM matrix."
  )
}

# Convert the TBX21 data into long format
TBX21_long <- TBX21_cpm %>%
  pivot_longer(
    cols = -Symbol,
    names_to = "sample",
    values_to = "CPM"
  ) %>%
  left_join(
    metadata[, c("sample", "condition")],
    by = "sample"
  ) %>%
  filter(
    condition %in%
      c(
        "HC_minus",
        "HC_plus",
        "NMO_minus",
        "NMO_plus"
      )
  )

# Relabel the conditions
TBX21_long$condition <- factor(
  TBX21_long$condition,
  levels = c(
    "HC_minus",
    "HC_plus",
    "NMO_minus",
    "NMO_plus"
  ),
  labels = c(
    "HC_CD25-",
    "HC_CD25hi",
    "NMOSD_CD25-",
    "NMOSD_CD25hi"
  )
)

# Plot TBX21 CPM values
TBX21_cpm_plot <- ggplot(
  TBX21_long,
  aes(
    x = condition,
    y = CPM,
    fill = condition
  )
) +
  geom_boxplot(
    outlier.shape = NA,
    alpha = 0.5
  ) +
  geom_jitter(
    position = position_jitter(
      width = 0.2,
      seed = 123
    ),
    size = 2,
    alpha = 0.7
  ) +
  scale_fill_manual(
    values = c(
      "HC_CD25-" = "grey70",
      "HC_CD25hi" = "grey40",
      "NMOSD_CD25-" = "pink",
      "NMOSD_CD25hi" = "red3"
    )
  ) +
  theme_bw() +
  labs(
    title = "TBX21",
    y = "CPM",
    x = ""
  )

TBX21_cpm_plot

# Save the plot as a PDF
ggsave(
  filename = file.path(
    CPM_OUTPUT_DIR,
    "TBX21_cpm_boxplot.pdf"
  ),
  plot = TBX21_cpm_plot,
  width = 4,
  height = 4,
  device = "pdf"
)

#### Comparison of ZEB2 CPM values ####

# Extract the ZEB2 row
ZEB2_cpm <- cpm_matrix[
  cpm_matrix$Symbol == "ZEB2",
]

if (nrow(ZEB2_cpm) == 0) {
  stop("ZEB2 was not found in the CPM matrix.")
}

if (nrow(ZEB2_cpm) > 1) {
  warning(
    "Multiple rows corresponding to ZEB2 were found in the CPM matrix."
  )
}

# Convert the ZEB2 data into long format
ZEB2_long <- ZEB2_cpm %>%
  pivot_longer(
    cols = -Symbol,
    names_to = "sample",
    values_to = "CPM"
  ) %>%
  left_join(
    metadata[, c("sample", "condition")],
    by = "sample"
  ) %>%
  filter(
    condition %in%
      c(
        "HC_minus",
        "HC_plus",
        "NMO_minus",
        "NMO_plus"
      )
  )

# Relabel the conditions
ZEB2_long$condition <- factor(
  ZEB2_long$condition,
  levels = c(
    "HC_minus",
    "HC_plus",
    "NMO_minus",
    "NMO_plus"
  ),
  labels = c(
    "HC_CD25-",
    "HC_CD25hi",
    "NMOSD_CD25-",
    "NMOSD_CD25hi"
  )
)

# Plot ZEB2 CPM values
ZEB2_cpm_plot <- ggplot(
  ZEB2_long,
  aes(
    x = condition,
    y = CPM,
    fill = condition
  )
) +
  geom_boxplot(
    outlier.shape = NA,
    alpha = 0.5
  ) +
  geom_jitter(
    position = position_jitter(
      width = 0.2,
      seed = 123
    ),
    size = 2,
    alpha = 0.7
  ) +
  scale_fill_manual(
    values = c(
      "HC_CD25-" = "grey70",
      "HC_CD25hi" = "grey40",
      "NMOSD_CD25-" = "pink",
      "NMOSD_CD25hi" = "red3"
    )
  ) +
  theme_bw() +
  labs(
    title = "ZEB2",
    y = "CPM",
    x = ""
  )

ZEB2_cpm_plot

# Save the plot as a PDF
ggsave(
  filename = file.path(
    CPM_OUTPUT_DIR,
    "ZEB2_cpm_boxplot.pdf"
  ),
  plot = ZEB2_cpm_plot,
  width = 4,
  height = 4,
  device = "pdf"
)