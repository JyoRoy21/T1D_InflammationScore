### 28 genes Heatmap time point 3
## Antonio H.
# 10/22/2024

# packages ---------------------------------------------------------------


#Installing GSVA for this script
install.packages("BiocManager")
BiocManager::install("GSVA", version="3.19", force=TRUE)
BiocManager::install("DESeq2")
install.packages("pheatmap")

#Loading the required packages
library(GSVA)
library(msigdbr)
library(pheatmap)
library(dplyr)
library(ggplot2)
library(readr)
library(readxl)
library(tidyverse)
library(DESeq2)
library(RColorBrewer)

# Creating NOD & NS Datasets ---------------------------------------------------------

# RUN LINE 1 TO LINE 174 ON COMBINE (TEST) SCRIPT
#source("project/Code/Raw Code/Matrix/GEO/Combined (test).R")

# Setting the dir I need for this script to run
setwd("project/")

# Load gene expression data 
gene_expected_count <- read_csv("Data/Private Data (from Lab)/gene_expected_count.annot.csv") # all in one dataset

#Filtering only the genes of interest
gene_expected<-gene_expected_count %>% 
  filter(external_gene_name%in%c( "Foxp3","Ptpn2","Tigit","Tgfb1", "Il1b", "Il17a", "Ifng","Serpine1", "Ptgs2","Nos1", "Nos2", "Nos3","Sod1","Sod2","Sod3","Cat", "Gpx1","Gpx2","Gclc", "Gclm","Gsr","Gstm1", "Gstt1","Gstp1", "Gzma","Prf1","Klrg1","Kctd16"))

# NOD Sample 3
nod_sample_3 <-gene_expected %>% dplyr::select(entrezgene_id,external_gene_name,`6872-JK-3`,`6872-JK-8`,`6872-JK-13`,`6872-JK-18`,`6872-JK-23`) %>% 
  dplyr::rename("NOD_Sample3_t1"="6872-JK-3","NOD_Sample3_t2"="6872-JK-8","NOD_Sample3_t3"="6872-JK-13","NOD_Sample3_t4"="6872-JK-18","NOD_Sample3_t5"="6872-JK-23")

# NS Sample 3
ns_sample_3 <-gene_expected %>% dplyr::select(entrezgene_id,external_gene_name,`6872-JK-53`,`6872-JK-57`,`6872-JK-61`,`6872-JK-65`,`6872-JK-69`) %>% 
  dplyr::rename("NS_Sample3_t1"="6872-JK-53","NS_Sample3_t2"="6872-JK-57","NS_Sample3_t3"="6872-JK-61","NS_Sample3_t4"="6872-JK-65","NS_Sample3_t5"="6872-JK-69")

t3_dataset<- nod_sample_3 %>% full_join(ns_sample_3)

rm(gene_expected,nod_sample_3,ns_sample_3)

# Creating Matrix ---------------------------------------------------------

#Required Matrix & changes 
gene_expected_matrix <- t3_dataset %>%
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


# Creating Plot -----------------------------------------------------------


# Plot the heatmap with reordered columns
pheatmap::pheatmap(gene_normalized,annotation_col = my_gene_col, cluster_cols = FALSE,   # Disable column clustering
                   cluster_rows = FALSE    # Disable row clustering (optional)
)