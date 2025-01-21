





# Install & Load Packages -------------------------------------------------

library(dplyr)
library(ggplot2)
library(readr)
library(readxl)
library(tidyverse)
library(tidyr)
library(ReactomePA)

# Loading & Tidying GEO Data (UP) --------------------------------------------------

# Setting the dir I need for this script to run
setwd("project/")

#Load GEO hyper expressed gene data 
GEO_Dataset_up <- read_excel("Data/Public Data/GEO Dataset.xlsx", 
                             sheet = "Disease_Perturbations_from_GEO_")

# Selecting only the columns I am interest in
GEO_Dataset_up<-GEO_Dataset_up %>% filter(...2=="X")

# Dataset with all UP genes in row 1
row_as_column_df <- GEO_Dataset_up[1, 3:403] %>% 
  t() %>%
  as.data.frame() %>%
  rename("Genes"=V1) 
# V1="Genes" 
# "Genes"=V1

# Dataset with all UP genes in row 2
more_genes<-GEO_Dataset_up[2, 3:403] %>% 
  t() %>%
  as.data.frame() %>%
  rename("Genes"=V1)

# Adding row 2 genes to the dataset
row_as_column_df<-row_as_column_df %>% full_join(more_genes)

# Dataset with all UP genes in row 3
more_genes<-GEO_Dataset_up[3, 3:403] %>% 
  t() %>%
  as.data.frame() %>%
  rename("Genes"=V1)

# Adding row 3 genes to the dataset
row_as_column_df<-row_as_column_df %>% full_join(more_genes)

# Preventing the script from crashing
row_as_column_df<-row_as_column_df %>% distinct()

# Dataset with all UP genes in row 4
more_genes<-GEO_Dataset_up[4, 3:403] %>% 
  t() %>%
  as.data.frame() %>%
  rename("Genes"=V1)

# Adding row 4 genes to the dataset
row_as_column_df<-row_as_column_df %>% full_join(more_genes)

# Dataset with all UP genes in row 5
more_genes<-GEO_Dataset_up[5, 3:403] %>% 
  t() %>%
  as.data.frame() %>%
  rename("Genes"=V1)

# Adding row 5 genes to the dataset
row_as_column_df<-row_as_column_df %>% full_join(more_genes)

# Preventing the script from crashing
row_as_column_df<-row_as_column_df %>% distinct()

# Dataset with all UP genes in row 6
more_genes<-GEO_Dataset_up[6, 3:403] %>% 
  t() %>%
  as.data.frame() %>%
  rename("Genes"=V1)

# Adding row 6 genes to the dataset
up_row_as_column_df<-row_as_column_df %>%
  full_join(more_genes) %>% mutate("Type"="up") #1470

#
rm(row_as_column_df)

# Loading & Tidying GEO Data (DOWN) --------------------------------------------------

#Load GEO hyper expressed gene data 
GEO_Dataset_down <- read_excel("Data/Public Data/GEO Dataset.xlsx", 
                               sheet = "Sheet2")

# Selecting only the columns I am interest in
GEO_Dataset_down<-GEO_Dataset_down %>% filter(...2=="x")

# Dataset with all DOWN genes in row 1
row_as_column_df <- GEO_Dataset_down[1, 3:403] %>% 
  t() %>%
  as.data.frame() %>%
  rename("Genes"=V1)

# Dataset with all DOWN genes in row 2
more_genes<-GEO_Dataset_down[2, 3:403] %>% 
  t() %>%
  as.data.frame() %>%
  rename("Genes"=V1)

# Adding row 2 genes to the dataset
row_as_column_df<-row_as_column_df %>% full_join(more_genes)

# Dataset with all DOWN genes in row 3
more_genes<-GEO_Dataset_down[3, 3:403] %>% 
  t() %>%
  as.data.frame() %>%
  rename("Genes"=V1)

# Adding row 3 genes to the dataset
row_as_column_df<-row_as_column_df %>% full_join(more_genes)

# Preventing the script from crashing
row_as_column_df<-row_as_column_df %>% distinct()

# Dataset with all DOWN genes in row 4
more_genes<-GEO_Dataset_down[4, 3:403] %>% 
  t() %>%
  as.data.frame() %>%
  rename("Genes"=V1)

# Adding row 4 genes to the dataset
row_as_column_df<-row_as_column_df %>% full_join(more_genes)

# Dataset with all DOWN genes in row 5
more_genes<-GEO_Dataset_down[5, 3:403] %>% 
  t() %>%
  as.data.frame() %>%
  rename("Genes"=V1)

# Adding row 5 genes to the dataset
row_as_column_df<-row_as_column_df %>% full_join(more_genes)

# Preventing the script from crashing
row_as_column_df<-row_as_column_df %>% distinct()

# Dataset with all DOWN genes in row 6
more_genes<-GEO_Dataset_down[6, 3:403] %>% 
  t() %>%
  as.data.frame() %>%
  rename("Genes"=V1)

# Adding row 6 genes to the dataset
down_row_as_column_df<-row_as_column_df %>% 
  full_join(more_genes) %>% mutate("Type"="down") #1470

#combining UP and DOWN Datasets
full_row_as_column_df<- up_row_as_column_df %>% full_join(down_row_as_column_df)


# Loading & Tidying Scaffold Data --------------------------------------------------

#Load gene expression data 
gene_expected_count <- read_csv("Data/gene_expected_count.annot.csv")

# Filtering Raw Data using GEO Data
gene_dataset<- gene_expected_count %>% filter(str_detect(external_gene_name, paste(full_row_as_column_df$Genes, collapse = "|"))) %>% filter(entrezgene_id!=".") #96 w/o filter, 27 w/ filter

#
rm(down_row_as_column_df,full_row_as_column_df,gene_expected_count, GEO_Dataset_up, GEO_Dataset_down,more_genes,row_as_column_df,up_row_as_column_df)

stop()
## heatmap_27genes = gene_dataset

##################

# To understand the matrix better
name_and_id<-gene_dataset %>% dplyr::select(entrezgene_id,external_gene_name)

# No need for this 
rm(row_as_column_df)#,row_as_column_df_down,row_as_column_df_up)

# Creating Matrix --------------------------------------------------------------

# The genes I am interested in (HAS TO BE ENTREZ ID)
genes_of_interest<- gene_dataset$entrezgene_id

# Pathway enrichment analysis
enriched_pathways <- enrichPathway(gene = genes_of_interest, organism = "mouse") 

# Extract pathway and gene information
pathways <- as.data.frame(enriched_pathways)

#Picking out the columns I am interested in
pathway2 <-pathways %>% rownames_to_column() %>% dplyr::select(Description,geneID)

# Separate the gene_id column into multiple rows
data_long<- pathway2 %>% separate_rows(geneID, sep = "/") %>% mutate(presence = 1) %>% dplyr::rename("pathway"=Description)

# Spread the data into a wide format with pathways as columns
data_wide<-data_long %>% pivot_wider(names_from = pathway, values_from = presence, values_fill = list(presence = 0)) # 14 out 27 genes

# Calculate the degree centrality
network_matrix <- data_wide %>%
  rowwise() %>%
  mutate(degree_centrality = sum(c_across(-geneID))) %>%
  ungroup() 

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

#
gene_gene_df %>% left_join(name_and_id,by = join_by("geneID"=="entrezgene_id")) %>% dplyr::select(geneID,external_gene_name, everything())



