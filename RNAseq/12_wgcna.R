###################################################
### R packages & Project directory ----

# install.packages("WGCNA")
# install.packages("tidyverse")
# install.packages("gridExtra")
# install.packages("devtools")
# devtools::install_github("kevinblighe/CorLevelPlot")
# install.packages("ggdendro")
# install.packages("ggrepel")
# install.packages("ggplotify")
# install.packages("tidygraph")
# install.packages("ggraph")
# install.packages("matrixStats")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
# BiocManager::install("DESeq2", force = TRUE)
# BiocManager::install("GEOquery", force = TRUE)
# BiocManager::install("impute", force = TRUE)

rm(list = ls())

# Set the working directory as appropriate for your environment.
# setwd("path/to/project_directory")

library(WGCNA)
library(DESeq2)
library(GEOquery)
library(tidyverse)
library(gridExtra)
library(impute)
library(CorLevelPlot)
library(ggplot2)
library(ggdendro)
library(ggrepel)
library(ggplotify)
library(tidygraph)
library(ggraph)
library(dplyr)
library(matrixStats)

options(stringsAsFactors = FALSE)
allowWGCNAThreads()

###################################################
### 0. SETTINGS  ----

# Use a signed network
NETWORK_TYPE <- "signed"
TOM_TYPE     <- "signed"

# Number of highly variable genes
topN <- 5000

# min module size 
minModuleSize <- 50
deepSplit <- 1
mergeCutHeight <- 0.25
pamStage <- FALSE

# soft power
soft_power <- 14

###################################################
### 1. Fetch Data ----

data <- read.csv("path/to/gene_count_matrix.csv", header = TRUE)
phenoData <- read.csv("path/to/phenoData.csv", header = TRUE)
phenoData <- phenoData |> tibble::column_to_rownames("sample")

###################################################
### 2. QC - outlier detection ----

data <- data |> tibble::column_to_rownames("Symbol")

gsg <- goodSamplesGenes(t(data))
summary(gsg)
gsg$allOK
table(gsg$goodGenes)
table(gsg$goodSamples)

data <- data[gsg$goodGenes == TRUE, ]

# sample clustering
htree <- hclust(dist(t(data)), method = "average")
plot(htree)

dend <- as.dendrogram(htree)
p1 <- ggdendrogram(dend)
ggsave(plot = p1,
       filename = "path/to/results/htree.png",
       width = 8, height = 6, dpi = 300)

# PCA
pca <- prcomp(t(data))
pca.dat <- as.data.frame(pca$x)
pca.var <- pca$sdev^2
pca.var.percent <- round(pca.var/sum(pca.var)*100, digits = 2)

p2 <- ggplot(pca.dat, aes(PC1, PC2)) +
  geom_point() +
  geom_text_repel(aes(label = rownames(pca.dat))) +
  labs(x = paste0("PC1: ", pca.var.percent[1], " %"),
       y = paste0("PC2: ", pca.var.percent[2], " %"))
p2

ggsave(plot = p2,
       filename = "path/to/results/pca.png",
       width = 8, height = 6, dpi = 300)

###################################################
### 3. Sample subset & colData alignment ----
data.subset <- data
colData <- phenoData

# Ensure identical sample order
colData <- colData[colnames(data.subset), ]
stopifnot(all(rownames(colData) == colnames(data.subset)))

###################################################
### 4. Normalization (DESeq2 vst) ----

dds <- DESeqDataSetFromMatrix(countData = data.subset,
                              colData = colData,
                              design = ~ 1)

# Retain genes expressed (count ≥1) in at least 50% of samples
n_samples <- ncol(dds)
threshold <- ceiling(n_samples * 0.5)

dds_filt <- dds[rowSums(counts(dds) >= 1) >= threshold, ]
nrow(dds_filt)

dds_norm <- vst(dds_filt)
norm.counts <- assay(dds_norm) %>% t()

###################################################
### 5. Network Construction (SIGNED) ----

# soft threshold
power <- c(c(1:10), seq(from = 12, to = 50, by = 2))

sft <- pickSoftThreshold(norm.counts,
                         powerVector = power,
                         networkType = NETWORK_TYPE,
                         verbose = 5)

sft.data <- sft$fitIndices

a1 <- ggplot(sft.data, aes(Power, SFT.R.sq, label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, color = "red") +
  labs(x = "Power", y = "Scale free topology model fit, signed R^2") +
  theme_classic()

a2 <- ggplot(sft.data, aes(Power, mean.k., label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  labs(x = "Power", y = "Mean Connectivity") +
  theme_classic()

p3 <- grid.arrange(a1, a2, nrow = 2)
ggsave(plot = p3,
       filename = "path/to/results/power_signed.png",
       width = 8, height = 6, dpi = 300)

# dat: genes x samples
dat <- t(norm.counts)

# Select the top topN most variable genes
keep <- order(rowVars(dat), decreasing = TRUE)[seq_len(min(nrow(dat), topN))]
dat_top <- dat[keep, ]


# WGCNA::cor 
temp_cor <- cor
cor <- WGCNA::cor

bwnet_top <- blockwiseModules(
  t(dat_top),                   # samples x genes
  power = soft_power,
  deepSplit = deepSplit,
  minModuleSize = minModuleSize,
  mergeCutHeight = mergeCutHeight,
  pamStage = pamStage,
  networkType = NETWORK_TYPE,
  TOMType = TOM_TYPE,
  verbose = 3
)

cor <- temp_cor
###################################################
### 6. Module Eigengenes ----

module_eigengenes <- bwnet_top$MEs
head(module_eigengenes)
table(bwnet_top$colors)

for (i in 1:length(bwnet_top$dendrograms)) {
  plotDendroAndColors(
    bwnet_top$dendrograms[[i]],
    cbind(bwnet_top$unmergedColors[[i]],
          bwnet_top$colors[bwnet_top$blockGenes[[i]]]),
    c("unmerged", "merged"),
    main = paste("Gene dendrogram and module colors in block", i),
    dendroLabels = FALSE,
    addGuide = TRUE,
    hang = 0.03,
    guideHang = 0.05
  )
}

###################################################
### 7. Relate modules to traits ----
# ---- Create trait matrix ----
colData$cell_state <- factor(colData$cell_state,
                             levels = c("HC_plus", "NMO_plus", "HC_minus", "NMO_minus"))

traits_B <- model.matrix(~ 0 + cell_state, data = colData)
colnames(traits_B) <- gsub("cell_state", "", colnames(traits_B))

traits <- traits_B  

# Calculate module–trait correlations
nSamples <- nrow(module_eigengenes)  # module eigengenes は samples x modules
module.trait.corr <- cor(module_eigengenes, traits, use = "p")
module.trait.corr.pvals <- corPvalueStudent(module.trait.corr, nSamples)

# heatmap
heatmap.data <- merge(module_eigengenes, traits, by = "row.names")
heatmap.data <- heatmap.data %>% column_to_rownames("Row.names")

# CorLevelPlot
trait_cols  <- (ncol(module_eigengenes) + 1):ncol(heatmap.data)
module_cols <- 1:ncol(module_eigengenes)

p4 <- as.ggplot(
  CorLevelPlot(
    heatmap.data,
    x = names(heatmap.data)[trait_cols],
    y = names(heatmap.data)[module_cols],
    col = c("blue1", "skyblue", "white", "pink", "red"),
    cexLabX = 0.7,
    cexLabY = 1.0,
    cexCorval = 0.9
  )
)

p4

ggsave("path/to/results/heatmap_module_trait_signed.pdf",
       p4, width = 6, height = 6, units = "in", dpi = 300)

###################################################
### 8. Hub gene identification ---

bestModule <- "red"

if (is.matrix(traits) && "NMO_plus" %in% colnames(traits)) {
  trait <- traits[, "NMO_plus"]
  names(trait) <- rownames(traits)
} else {
  trait <- traits$disease_state_bin
  names(trait) <- rownames(traits)
}

datExpr <- dat_top  # genes x samples
MEs <- module_eigengenes

nSamples <- ncol(datExpr)

geneModuleMembership <- as.data.frame(cor(t(datExpr), MEs, use = "p"))
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))

geneTraitSignificance <- as.data.frame(cor(t(datExpr), trait, use = "p"))
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))

moduleGeneNames <- names(bwnet_top$colors[bwnet_top$colors == bestModule])

hubTable <- data.frame(
  Gene = moduleGeneNames,
  ModuleMembership = geneModuleMembership[moduleGeneNames, paste0("ME", bestModule)],
  MM_pvalue = MMPvalue[moduleGeneNames, paste0("ME", bestModule)],
  GeneSignificance = geneTraitSignificance[moduleGeneNames, 1],
  GS_pvalue = GSPvalue[moduleGeneNames, 1]
) %>%
  arrange(desc(abs(ModuleMembership)))

top15 <- hubTable %>% slice_max(order_by = abs(ModuleMembership), n = 15)
cat("Top 15 hub genes in ", bestModule, " module:\n", sep = "")
print(top15$Gene)

###################################################
### 9. Network visualization via TOM (SIGNED) ----

datExpr_for_TOM <- t(dat_top)  # samples x genes
TOM <- TOMsimilarityFromExpr(datExpr_for_TOM,
                             power = soft_power,
                             networkType = NETWORK_TYPE,
                             TOMType = TOM_TYPE)

colnames(TOM) <- rownames(dat_top)
rownames(TOM) <- rownames(dat_top)

module_of_interest <- bestModule
genes_in_module <- names(bwnet_top$colors)[bwnet_top$colors == module_of_interest]

TOM_sub <- TOM[genes_in_module, genes_in_module, drop = FALSE]
TOM_sub[lower.tri(TOM_sub, diag = TRUE)] <- NA

edges <- as.data.frame(as.table(TOM_sub)) %>%
  filter(!is.na(Freq)) %>%
  rename(from = Var1, to = Var2, weight = Freq)

thr <- quantile(edges$weight, 0.99, na.rm = TRUE)
edges <- edges %>% filter(weight >= thr)

kME_df <- hubTable[, c("Gene", "ModuleMembership")] %>%
  rename(name = Gene, kME = ModuleMembership)

top_genes <- kME_df %>%
  arrange(desc(abs(kME))) %>%
  slice_head(n = 30) %>%
  pull(name)

edges <- edges %>% filter(from %in% top_genes & to %in% top_genes)

nodes <- kME_df %>%
  filter(name %in% top_genes) %>%
  mutate(hub = rank(-abs(kME), ties.method = "first") <= 15)

graph_tbl <- tbl_graph(nodes = nodes, edges = edges, directed = FALSE)

set.seed(123)
layout_df <- create_layout(graph_tbl, layout = "fr")

p_net <- ggraph(layout_df) +
  geom_edge_link(aes(width = weight), alpha = 0.35, colour = "grey55", show.legend = FALSE) +
  scale_edge_width(range = c(0.4, 1.2)) +
  geom_node_point(aes(size = abs(kME), colour = hub), show.legend = FALSE) +
  scale_color_manual(values = c("FALSE" = "grey75", "TRUE" = "tomato")) +
  geom_label_repel(
    data = layout_df %>% dplyr::filter(hub),
    aes(x = x, y = y, label = name),
    size = 6,
    color = "tomato4",
    fill = scales::alpha("white", 0.9),
    fontface = "bold",
    label.size = 0.2,
    segment.color = "tomato4",
    segment.size = 0.3,
    box.padding = 0.5,
    point.padding = 0.9,
    force = 3,
    max.overlaps = Inf,
    min.segment.length = 0,
    direction = "both",
    seed = 123
  ) +
  ggtitle(paste0("SIGNED network | Module: ", module_of_interest,
                 " | Top ", length(top_genes), " genes by kME")) +
  theme_void()

p_net

ggsave("path/to/results/network_signed_top15.pdf",
       p_net, width = 9, height = 7, units = "in", dpi = 300)


###################################################
### 10. Intramodular connectivity and hub ranking
###################################################

bestModule <- "red"

# trait vector
if (is.matrix(traits) && "NMO_plus" %in% colnames(traits)) {
  
  trait <- traits[, "NMO_plus"]
  names(trait) <- rownames(traits)
  
} else {
  
  trait <- traits$disease_state_bin
  names(trait) <- rownames(traits)
}

# dat_top: genes x samples
datExpr <- dat_top

# WGCNA input: samples x genes
datExpr_wgcna <- as.data.frame(
  t(datExpr)
)

MEs <- module_eigengenes
nSamples <- nrow(datExpr_wgcna)

# Alignment checks
stopifnot(
  identical(
    rownames(datExpr_wgcna),
    rownames(MEs)
  )
)

stopifnot(
  identical(
    rownames(datExpr_wgcna),
    names(trait)
  )
)

stopifnot(
  identical(
    colnames(datExpr_wgcna),
    names(bwnet_top$colors)
  )
)

###################################################
# Module membership
###################################################

geneModuleMembership <- as.data.frame(
  WGCNA::cor(
    datExpr_wgcna,
    MEs,
    use = "pairwise.complete.obs"
  )
)

MMPvalue <- as.data.frame(
  WGCNA::corPvalueStudent(
    as.matrix(geneModuleMembership),
    nSamples
  )
)

###################################################
# Gene significance
###################################################

geneTraitSignificance <- as.data.frame(
  WGCNA::cor(
    datExpr_wgcna,
    trait,
    use = "pairwise.complete.obs"
  )
)

GSPvalue <- as.data.frame(
  WGCNA::corPvalueStudent(
    as.matrix(geneTraitSignificance),
    nSamples
  )
)

###################################################
# Intramodular connectivity
###################################################
connectivity_all <- WGCNA::intramodularConnectivity.fromExpr(
  datExpr = datExpr_wgcna,
  colors = bwnet_top$colors,
  power = soft_power,
  networkType = NETWORK_TYPE
)

###################################################
# Ranking table for selected module
###################################################
rownames(connectivity_all) <- colnames(datExpr_wgcna)

head(rownames(connectivity_all))

identical(
  rownames(connectivity_all),
  colnames(datExpr_wgcna)
)


moduleGeneNames <- names(
  bwnet_top$colors[
    bwnet_top$colors == bestModule
  ]
)

moduleME_name <- paste0(
  "ME",
  bestModule
)

ranking_table <- data.frame(
  Gene = moduleGeneNames,
  
  Module = bestModule,
  
  kTotal = connectivity_all[
    moduleGeneNames,
    "kTotal"
  ],
  
  kWithin = connectivity_all[
    moduleGeneNames,
    "kWithin"
  ],
  
  kOut = connectivity_all[
    moduleGeneNames,
    "kOut"
  ],
  
  kDiff = connectivity_all[
    moduleGeneNames,
    "kDiff"
  ],
  
  ModuleMembership = geneModuleMembership[
    moduleGeneNames,
    moduleME_name
  ],
  
  MM_pvalue = MMPvalue[
    moduleGeneNames,
    moduleME_name
  ],
  
  GeneSignificance = geneTraitSignificance[
    moduleGeneNames,
    1
  ],
  
  GS_pvalue = GSPvalue[
    moduleGeneNames,
    1
  ],
  
  row.names = NULL,
  check.names = FALSE
)

ranking_table <- ranking_table %>%
  dplyr::arrange(
    dplyr::desc(kWithin)
  ) %>%
  dplyr::mutate(
    IntramodularConnectivityRank = dplyr::row_number(),
    
    kME_rank = rank(
      -abs(ModuleMembership),
      ties.method = "min"
    ),
    
    GS_rank = rank(
      -abs(GeneSignificance),
      ties.method = "min"
    )
  ) %>%
  dplyr::select(
    IntramodularConnectivityRank,
    Gene,
    Module,
    kWithin,
    kTotal,
    kOut,
    kDiff,
    ModuleMembership,
    MM_pvalue,
    kME_rank,
    GeneSignificance,
    GS_pvalue,
    GS_rank
  )

###################################################
# Export
###################################################

write.csv(
  ranking_table,
  file = paste0(
    "path/to/results/",
    bestModule,
    "_module_gene_ranking_by_intramodular_connectivity.csv"
  ),
  row.names = FALSE
)