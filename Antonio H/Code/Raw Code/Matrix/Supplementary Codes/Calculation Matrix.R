### Calculation Matrix for Mutation Dataset
## Antonio Holmes
# 06/07/2024


## WHAT NEEDS TO HAPPEN? ---------------------------------------------------

# I NEED A DATASET WITH JUST GENE SYMBOL, LOG FC(optional), AND P-VALUE
# then
# INSERT THIS TIABLE INTO run_pathfindR(example_pathfindR_input, p_val_threshold = 0.01)


#PATHFINR WEBSITE:  https://cran.r-project.org/web/packages/pathfindR/vignettes/intro_vignette.html

# Install & Load Packages -------------------------------------------------

#install.packages("pathfindR")
#BiocManager::install("org.Hs.eg.db")

#library(pathfindR)
library(limma)
library(dplyr)
library(ggplot2)
library(readr)
library(readxl)
library(tidyverse)

# Creating Req Datasets ---------------------------------------------------

#Running Data Tidying Script
#source("/Users/antonioholmes/project/Code/Raw Code/data_tidying.R")

# Load gene expression data 
gene_expected<-read_csv("Data/Private Data (from Lab)/gene_expected_count.annot.csv")

###
### NOD Samples
###

# NOD sample 1
nod_sample_1 <-gene_expected %>%dplyr::select(entrezgene_id, external_gene_name,`6872-JK-1`,`6872-JK-6`,`6872-JK-11`,`6872-JK-16`,`6872-JK-21`) %>% 
  dplyr::rename("NOD_Sample1_t1"="6872-JK-1","NOD_Sample1_t2"="6872-JK-6","NOD_Sample1_t3"="6872-JK-11","NOD_Sample1_t4"="6872-JK-16","NOD_Sample1_t5"="6872-JK-21")

# NOD Sample 2
nod_sample_2 <-gene_expected %>% dplyr::select(entrezgene_id,external_gene_name,`6872-JK-2`,`6872-JK-7`,`6872-JK-12`,`6872-JK-17`,`6872-JK-22`) %>% 
  dplyr::rename("NOD_Sample2_t1"="6872-JK-2","NOD_Sample2_t2"="6872-JK-7","NOD_Sample2_t3"="6872-JK-12","NOD_Sample2_t4"="6872-JK-17","NOD_Sample2_t5"="6872-JK-22")

# NOD Sample 3
nod_sample_3 <-gene_expected %>% dplyr::select(entrezgene_id,external_gene_name,`6872-JK-3`,`6872-JK-8`,`6872-JK-13`,`6872-JK-18`,`6872-JK-23`) %>% 
  dplyr::rename("NOD_Sample3_t1"="6872-JK-3","NOD_Sample3_t2"="6872-JK-8","NOD_Sample3_t3"="6872-JK-13","NOD_Sample3_t4"="6872-JK-18","NOD_Sample3_t5"="6872-JK-23")

# NOD Sample 4
nod_sample_4 <-gene_expected %>% dplyr::select(entrezgene_id,external_gene_name,`6872-JK-4`,`6872-JK-9`,`6872-JK-14`,`6872-JK-19`,`6872-JK-24`) %>% 
  dplyr::rename("NOD_Sample4_t1"="6872-JK-4","NOD_Sample4_t2"="6872-JK-9","NOD_Sample4_t3"="6872-JK-14","NOD_Sample4_t4"="6872-JK-19","NOD_Sample4_t5"="6872-JK-24")

# NOD Sample 5
nod_sample_5 <-gene_expected %>% dplyr::select(entrezgene_id,external_gene_name,`6872-JK-77`,`6872-JK-10`,`6872-JK-15`,`6872-JK-20`,`6872-JK-25`) %>% 
  dplyr::rename("NOD_Sample5_t1"="6872-JK-77","NOD_Sample5_t2"="6872-JK-10","NOD_Sample5_t3"="6872-JK-15","NOD_Sample5_t4"="6872-JK-20","NOD_Sample5_t5"="6872-JK-25")

# Combining Data sets
NOD_dataset<- nod_sample_1 %>% 
  full_join(nod_sample_2) %>% full_join(nod_sample_3) %>%
  full_join(nod_sample_4) %>% full_join(nod_sample_5)

#removing the NOD sample datasets for now
rm(nod_sample_1,nod_sample_2,nod_sample_3,nod_sample_4,nod_sample_5)

###
### NS SAMPLES
###

# NS Sample 1
ns_sample_1 <-gene_expected %>% dplyr::select(entrezgene_id,external_gene_name,`6872-JK-51`,`6872-JK-55`,`6872-JK-59`,`6872-JK-63`,`6872-JK-67`) %>% 
  dplyr::rename("NS_Sample1_t1"="6872-JK-51","NS_Sample1_t2"="6872-JK-55","NS_Sample1_t3"="6872-JK-59","NS_Sample1_t4"="6872-JK-63","NS_Sample1_t5"="6872-JK-67")

# NS Sample 2 == Seq 52 missing or the real column is 52_REPREP
ns_sample_2 <-gene_expected %>% dplyr::select(entrezgene_id,external_gene_name,`6872-JK-52-REPREP`,`6872-JK-56`,`6872-JK-60`,`6872-JK-64`,`6872-JK-68`) %>% 
  dplyr::rename("NS_Sample2_t1"="6872-JK-52-REPREP","NS_Sample2_t2"="6872-JK-56","NS_Sample2_t3"="6872-JK-60","NS_Sample2_t4"="6872-JK-64","NS_Sample2_t5"="6872-JK-68")

# NS Sample 3
ns_sample_3 <-gene_expected %>% dplyr::select(entrezgene_id,external_gene_name,`6872-JK-53`,`6872-JK-57`,`6872-JK-61`,`6872-JK-65`,`6872-JK-69`) %>% 
  dplyr::rename("NS_Sample3_t1"="6872-JK-53","NS_Sample3_t2"="6872-JK-57","NS_Sample3_t3"="6872-JK-61","NS_Sample3_t4"="6872-JK-65","NS_Sample3_t5"="6872-JK-69")

# NS Sample 4
ns_sample_4 <-gene_expected %>% dplyr::select(entrezgene_id,external_gene_name,`6872-JK-54`,`6872-JK-58`,`6872-JK-62`,`6872-JK-66`,`6872-JK-70`)%>% 
  dplyr::rename("NS_Sample4_t1"="6872-JK-54","NS_Sample4_t2"="6872-JK-58","NS_Sample4_t3"="6872-JK-62","NS_Sample4_t4"="6872-JK-66","NS_Sample4_t5"="6872-JK-70")

# NS Sample 5 = none

# Combining Data sets
NS_dataset<-ns_sample_1 %>% full_join(ns_sample_2) %>%
  full_join(ns_sample_3) %>%full_join(ns_sample_4)

#
rm(ns_sample_1,ns_sample_2,ns_sample_3,ns_sample_4)

# Free up some more space and RAM
rm(gene_expected)


# Loading & Tidying Data --------------------------------------------------

# DATASET W/ ALL SAMPLES WITH HLA MUTATION
mutation_dataset<- NOD_dataset %>% full_join(NS_dataset)

#Mutate duplicate column for Entrezgene IDs
mutation_dataset<-mutation_dataset %>% mutate(
  "Duplicate1"=duplicated(entrezgene_id),
  "Duplicate2"=duplicated(external_gene_name))

#Filtering out the duplcates
mutation_dataset<-mutation_dataset %>% 
  filter(Duplicate1==FALSE,Duplicate2==FALSE,entrezgene_id!=".") %>% 
  dplyr::select(-Duplicate1,-Duplicate2)

# GSVA requires a matrix so make sure there are no duplicates
sum(duplicated(mutation_dataset$entrezgene_id)) # = 0

# Matrix for Log FC & P-Val ---------------------------------------------------------

# IF YOU RUN PREVIOUS CHUNK, change gene_expected to filtered_gene_expected
mutation_matrix <- mutation_dataset %>%
  # GSVA can't the Ensembl IDs so we should drop this column as well as the means
  dplyr::select(-entrezgene_id) %>%
  # We need to store our gene identifiers as row names
  tibble::column_to_rownames("external_gene_name") %>%
  # Now we can convert our object into a matrix
  as.matrix()


# Creating Mutation Calculation Table ------------------------------------------------------------------

n_samples <- 45 ## number of samples
nod_Grp1 <- 25 ## number of samples in group 1
ns_Grp2 <- n_samples - nod_Grp1 ## number of samples in group 2

## build design matrix
gene_design <- cbind(sampleGroup1=1, sampleGroup2vs1=c(rep(0, nod_Grp1), rep(1, ns_Grp2)))

## fit linear model
fit_mut <- lmFit(mutation_matrix, gene_design)

## estimate moderated t-statistics
fit_mut <- eBayes(fit_mut)

## Making the results into a dataset, Mutation Calc
mut_calc<-topTable(fit_mut, coef="sampleGroup2vs1",number = 22238) #10 vs 27 looks weird

# Selecting the columns I need
mut_calc<- mut_calc %>% dplyr::select(logFC,P.Value) %>%
  rownames_to_column() %>% dplyr::rename("gene_ID"=rowname) %>% 
  filter(logFC>=1.5 & P.Value<=0.1)


# Spearman Correlation Analysis (NOD) -------------------------------------------

## gene_gene_df is from gene_matrix.R

# the gene sets of comparison
set1 <- gene_gene_df$geneID  # the 12 genes connected to at least 4 others
set2 <- mut_calc$gene_ID  # the genes from calc above

# Filter the data for the genes in set1 and set2
genes_set1 <- mutation_dataset %>% filter(entrezgene_id %in% set1)
genes_set2 <- mutation_dataset %>% filter(external_gene_name %in% set2)

# Select columns that start with "NOD_Sample"
genes_set1 <- genes_set1 %>% dplyr::select(entrezgene_id,external_gene_name, starts_with("NOD_Sample"))
genes_set2 <- genes_set2 %>% dplyr::select(entrezgene_id,external_gene_name, starts_with("NOD_Sample"))

# Remove the gene_id column to keep only the expression data
expr_set1 <- genes_set1 %>% dplyr::select(-entrezgene_id,-external_gene_name)
expr_set2 <- genes_set2 %>% dplyr::select(-entrezgene_id,-external_gene_name)

# Transpose the data so that genes are columns and samples are rows
expr_set1 <- t(expr_set1)
expr_set2 <- t(expr_set2)

# Calculate the Spearman correlation
# This will give you a correlation matrix
nod_correlation_matrix <- cor(expr_set1, expr_set2, method = "spearman")

# Rename the columns and rows for clarity
colnames(nod_correlation_matrix) <- set2
rownames(nod_correlation_matrix) <- set1

# Display the correlation matrix
print(nod_correlation_matrix)


# Strong Gene Pairs (NOD) -------------------------------------------------

# Find the positions of correlations greater than 0.5
positions <- which(nod_correlation_matrix > 0.6 | nod_correlation_matrix < -0.6, arr.ind = TRUE)

# Extract the gene names corresponding to these positions
nod_gene_pairs <- data.frame(
  gene_set1 = rownames(nod_correlation_matrix)[positions[, 1]],
  gene_set2 = colnames(nod_correlation_matrix)[positions[, 2]],
  correlation = nod_correlation_matrix[positions]
)

# Display the gene pairs with correlations greater than 0.5
print(nod_gene_pairs)

# Spearman Correlation Analysis (NS) -------------------------------------------

# the gene sets of comparison
set1 <- gene_gene_df$geneID  # the 12 genes connected to at least 4 others
set2 <- mut_calc$gene_ID  # the genes from calc above

# Filter the data for the genes in set1 and set2
genes_set1 <- mutation_dataset %>% filter(entrezgene_id %in% set1)
genes_set2 <- mutation_dataset %>% filter(external_gene_name %in% set2)

# Select columns that start with "NS_Sample"
genes_set1 <- genes_set1 %>% dplyr::select(entrezgene_id,external_gene_name, starts_with("NS_Sample"))
genes_set2 <- genes_set2 %>% dplyr::select(entrezgene_id,external_gene_name, starts_with("NS_Sample"))

# Remove the gene_id column to keep only the expression data
expr_set1 <- genes_set1 %>% dplyr::select(-entrezgene_id,-external_gene_name)
expr_set2 <- genes_set2 %>% dplyr::select(-entrezgene_id,-external_gene_name)

# Transpose the data so that genes are columns and samples are rows
expr_set1 <- t(expr_set1)
expr_set2 <- t(expr_set2)

# This will give you a correlation matrix
ns_correlation_matrix <- cor(expr_set1, expr_set2, method = "spearman")

# Rename the columns and rows for clarity
colnames(ns_correlation_matrix) <- set2
rownames(ns_correlation_matrix) <- set1

# Display the correlation matrix
print(ns_correlation_matrix)


# Strong Gene Pairs (NS) -------------------------------------------------

# Find the positions of correlations greater than 0.5
positions <- which(ns_correlation_matrix > 0.6 | ns_correlation_matrix < -0.6, arr.ind = TRUE)

# Extract the gene names corresponding to these positions
ns_gene_pairs <- data.frame(
  gene_set1 = rownames(ns_correlation_matrix)[positions[, 1]],
  gene_set2 = colnames(ns_correlation_matrix)[positions[, 2]],
  correlation = ns_correlation_matrix[positions]
)

# Display the gene pairs with correlations greater than 0.5
print(ns_gene_pairs)



# Gene Freq Dataset & Overlap --------------------------------------------------------------

# Seeing the overlap of gene with a correlation >0.5
intersect(nod_gene_pairs$gene_set2,ns_gene_pairs$gene_set2)

# Count the occurrences of each gene in the second column of NOD dataset
nod_gene_count <- nod_gene_pairs %>%
  count(gene_set2, name = "count") %>% filter(count>=4)

# Count the occurrences of each gene in the second column of NS dataset
ns_gene_count <- ns_gene_pairs %>%
  count(gene_set2, name = "count") %>% filter(count>=4)

# Seeing the overlap of gene with 4 or more connections
overlap<-intersect(nod_gene_count$gene_set2,ns_gene_count$gene_set2)



# Saving Overlapping Genes Identifier -------------------------------------

#
geneset<-mutation_dataset %>% filter(external_gene_name%in%overlap)
geneset<-geneset$entrezgene_id

#
#genes_of_interest<-c(geneset,set1)
genes_of_interest<-c(geneset,set1) #TESTING








# Plugging Results into PathfindR (NOT NEEDED, JUST TESTING) ---------------------------------------------------------------

#pathfindR::run_pathfindR(mut_calc,p_val_threshold = 0.5) #threshold should be 0.05
        # ## Testing input
        # The input looks OK
        # ## Processing input. Converting gene symbols,
        # if necessary (and if human gene symbols provided)
        # Number of genes provided in input: 10
        # Number of genes in input after p-value filtering: 10 
##
        # Error in input_processing(input, p_val_threshold, pin_path, convert2alias) : 
        #   None of the genes were in the PIN
        # Please check your gene symbols
              #THIS CODE NEEDS GENE SYMBOLS










