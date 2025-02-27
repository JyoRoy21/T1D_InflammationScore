### 27 genes + Original Genes Heatmap time point 4
## Antonio H.
# 10/21/2024



# Pre-Requirements & T1 Dataset ---------------------------------------------------------

#
source("Code/Raw Code/Matrix/EnrichR/27 EnrichR Genes/27 Genes + OG/27 genes with OG heatmap.R")

t4_dataset <- gene_expected_count %>% select(entrezgene_id,external_gene_name,contains("_t4"))

# Matrix & Plotting ----------------------------------------------------------------

#Required Matrix & changes 
gene_expected_matrix <- t4_dataset %>%
  # GSVA can't the Ensembl IDs so we should drop this column as well as the means
  dplyr::select(#-Duplicate,
    -entrezgene_id) %>%
  # We need to store our gene identifiers as row names
  tibble::column_to_rownames("external_gene_name") %>%
  # Now we can convert our object into a matrix
  as.matrix()

# Create a data frame with the column names classified as "NOD" or "NS"
my_gene_col <- data.frame(group = ifelse(grepl("NOD_",
                                               colnames(gene_expected_matrix)), 
                                         yes = "NOD", no = "NS"))

# Set row names of the metadata to match the column names of the gene matrix
rownames(my_gene_col) <- colnames(gene_expected_matrix)

# Ensure the 'group' column is treated as a factor
my_gene_col$group <- as.factor(my_gene_col$group)

# Normalizing the data
gene_normalized <- t(scale(t(gene_expected_matrix)))


# Plot the heatmap with reordered columns
pheatmap::pheatmap(gene_normalized,annotation_col = my_gene_col, cluster_cols = FALSE,   # Disable column clustering
                   cluster_rows = FALSE    # Disable row clustering (optional)
)




