### 28 genes Heatmap time point 2
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

# NOD Sample 2
nod_sample_2 <-gene_expected %>% dplyr::select(entrezgene_id,external_gene_name,`6872-JK-2`,`6872-JK-7`,`6872-JK-12`,`6872-JK-17`,`6872-JK-22`) %>% 
  dplyr::rename("NOD_Sample2_t1"="6872-JK-2","NOD_Sample2_t2"="6872-JK-7","NOD_Sample2_t3"="6872-JK-12","NOD_Sample2_t4"="6872-JK-17","NOD_Sample2_t5"="6872-JK-22")

# NS Sample 2 == Seq 52 missing or the real column is 52_REPREP
ns_sample_2 <-gene_expected %>% dplyr::select(entrezgene_id,external_gene_name,`6872-JK-52-REPREP`,`6872-JK-56`,`6872-JK-60`,`6872-JK-64`,`6872-JK-68`) %>% 
  dplyr::rename("NS_Sample2_t1"="6872-JK-52-REPREP","NS_Sample2_t2"="6872-JK-56","NS_Sample2_t3"="6872-JK-60","NS_Sample2_t4"="6872-JK-64","NS_Sample2_t5"="6872-JK-68")

t2_dataset<- nod_sample_2 %>% full_join(ns_sample_2)

rm(gene_expected,nod_sample_2,ns_sample_2)

# Creating Matrix ---------------------------------------------------------

#Required Matrix & changes 
gene_expected_matrix <- t2_dataset %>%
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