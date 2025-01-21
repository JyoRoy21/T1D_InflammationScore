### Data Driven GS Algorithm (Bootstrapping)
## Antonio Holmes
# 06/07/2024



# Creating Matrix from Degree of Centrality [TESTING] ---------------------

# Using the pathways tbl, create a tbl with the cols = pathways & rows = genes (this dataset was created in GS Algorithm BS)

### Add a source function for data_tidying then change gene_expected to gene_expected_count


# Lib Package -------------------------------------------------------------

library(tidyr)
library(ReactomePA)
library(tidyverse,dplyr)

# Creating Matrix --------------------------------------------------------------

# Pathway enrichment analysis
enriched_pathways <- enrichPathway(gene = genes_of_interest, organism = "mouse")

# Extract pathway and gene information
pathways <- as.data.frame(enriched_pathways) # 18 out of 27 genes

#Picking out the columns I am interested in
pathway2 <-pathways %>% rownames_to_column() %>% dplyr::select(Description,geneID)

# Separate the gene_id column into multiple rows
data_long<- pathway2 %>% separate_rows(geneID, sep = "/") %>% mutate(presence = 1) %>% dplyr::rename("pathway"=Description)

# Spread the data into a wide format with pathways as columns
data_wide<-data_long %>% pivot_wider(names_from = pathway, values_from = presence, values_fill = list(presence = 0))

# Calculate the degree centrality
network_matrix <- data_wide %>%
  rowwise() %>%
  mutate(degree_centrality = sum(c_across(-geneID))) %>%
  ungroup() 

# Adding the gene symbol and reordering the columns
# network_matrix<-genesets_hallmark %>% dplyr::select(entrez_gene,gene_symbol) %>%
#   right_join(network_matrix,by=join_by(entrez_gene==geneID)) %>%
#   dplyr::select(entrez_gene,gene_symbol,degree_centrality, everything()) %>% unique()

# Create the gene-gene adjacency matrix
genes <- network_matrix$geneID
pathways <- network_matrix %>% dplyr::select(-degree_centrality, -geneID#,-gene_symbol
                                             )

# Initialize an empty adjacency matrix
gene_gene_matrix <- matrix(0, nrow = length(genes), ncol = length(genes))
rownames(gene_gene_matrix) <- genes
colnames(gene_gene_matrix) <- genes

# Populate the adjacency matrix
for (i in 1:nrow(pathways)) {
  for (j in 1:nrow(pathways)) {
    if (i != j) {
      common_pathways <- sum(pathways[i, ] * pathways[j, ])
      if (common_pathways > 0) {
        gene_gene_matrix[i, j] <- 1
      }
    }
  }
}

# Convert the matrix to a dataframe for easier viewing
gene_gene_df <- as.data.frame(gene_gene_matrix)

gene_gene_df<-gene_gene_df %>% rownames_to_column() %>% rowwise() %>%
  dplyr::rename("geneID"=rowname) %>% 
  mutate(degree_centrality = sum(c_across(-geneID))) %>%
  ungroup() %>% 
  dplyr::select(geneID,degree_centrality, everything())


# Creating a Cutoff -------------------------------------------------------

# only want genes connected to at least 4 other
gene_gene_df<-gene_gene_df %>% filter(degree_centrality>=4)

# Missing Values ----------------------------------------------------------
# 
# ###
# ###
# 
# # Convert the data frame to a vector of unique values
# df_values <- unique(as.vector(as.matrix(data_wide$geneID)))
# 
# # Identify missing values
# missing_values <- !genes_of_interest %in% df_values
# 
# #
# gene_expected_count %>% mutate("missing"=missing_values) %>% dplyr::select(entrezgene_id,missing,everything()) %>% filter(missing==TRUE)%>% view()
# 
# ###
# ###
