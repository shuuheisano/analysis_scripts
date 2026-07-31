###################################################
### Export gene lists for Metascape ----

### Upregulated genes in NMO_plus

# Extract upregulated DEGs
NMO_plus_up <- DEG_NMO_plus_vs_HC_plus[
  DEG_NMO_plus_vs_HC_plus$log2FoldChange >= log2fc_cutoff,
  ,
  drop = FALSE
]

# Extract gene symbols
NMO_plus_up_gene_list <- NMO_plus_up[
  ,
  "Symbol",
  drop = FALSE
]

# Remove missing, empty, and duplicated gene symbols
NMO_plus_up_gene_list <- NMO_plus_up_gene_list[
  !is.na(NMO_plus_up_gene_list$Symbol) &
    NMO_plus_up_gene_list$Symbol != "",
  ,
  drop = FALSE
]

NMO_plus_up_gene_list <- unique(
  NMO_plus_up_gene_list
)

# Check output
nrow(NMO_plus_up_gene_list)
head(NMO_plus_up_gene_list)

# Export gene list
write.csv(
  NMO_plus_up_gene_list,
  file = file.path(
    output_dir,
    "NMO_plus_upregulated_genes_for_Metascape.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)


###################################################
### Upregulated genes in NMO_minus ----

# Extract upregulated DEGs
NMO_minus_up <- DEG_NMO_minus_vs_HC_minus[
  DEG_NMO_minus_vs_HC_minus$log2FoldChange >= log2fc_cutoff,
  ,
  drop = FALSE
]

# Extract gene symbols
NMO_minus_up_gene_list <- NMO_minus_up[
  ,
  "Symbol",
  drop = FALSE
]

# Remove missing, empty, and duplicated gene symbols
NMO_minus_up_gene_list <- NMO_minus_up_gene_list[
  !is.na(NMO_minus_up_gene_list$Symbol) &
    NMO_minus_up_gene_list$Symbol != "",
  ,
  drop = FALSE
]

NMO_minus_up_gene_list <- unique(
  NMO_minus_up_gene_list
)

# Check output
nrow(NMO_minus_up_gene_list)
head(NMO_minus_up_gene_list)

# Export gene list
write.csv(
  NMO_minus_up_gene_list,
  file = file.path(
    output_dir,
    "NMO_minus_upregulated_genes_for_Metascape.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)