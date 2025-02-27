### Bootstrapping a Possible Machine Learning Function
## Antonio H.
# 11/01/2024



# To-Do List --------------------------------------------------------------

# Add calc matrix code [DONE] 
# Add an gene expression column in the gene_expected dataset [DONE]
# what is an ideal LogFC to say gene expression is high?
# 167 genes that were not in our 'genes of interest'
# flip X train around samples = rows & genes = columns (mutation_dataset)
# add a diabetic or not column 

# Set Up & Req Dataset (Working) -------------------------------------------------------

# Load required libraries
library(randomForest)
library(caret)
library(dplyr)
library(limma)
library(ggplot2)
library(readr)
library(readxl)
library(tidyverse)

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

# Tidying Dataset (Working = UNKNOWN) ---------------------------------------------------------

# Need the samples=rows and genes=col
test1<- NOD_dataset %>% 
  #select(-entrezgene_id) %>%
  pivot_wider(names_from = external_gene_name,
              values_from = starts_with("NOD_sample"))

# Need the samples=rows and genes=col
test2<- NS_dataset %>% 
  #select(-entrezgene_id) %>%
  pivot_wider(names_from = external_gene_name,
              values_from = starts_with("NS_sample"))

#####################    NEW CODE ABOVE     ################################

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

# Removing files once done with them
rm(NOD_dataset,NS_dataset)


# Calculation Matrix (Working but Skipping for now) ------------------------------------------------------

# IF YOU RUN PREVIOUS CHUNK, change gene_expected to filtered_gene_expected
mutation_matrix <- mutation_dataset %>%
  # GSVA can't the Ensembl IDs so we should drop this column as well as the means
  dplyr::select(-entrezgene_id) %>%
  # We need to store our gene identifiers as row names
  tibble::column_to_rownames("external_gene_name") %>%
  # Now we can convert our object into a matrix
  as.matrix()

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
   filter(logFC>=0 & P.Value<=0.1)

#
 gene_expected_filtered <-mutation_dataset %>% left_join(mut_calc, by = join_by(external_gene_name == gene_ID)) %>% filter(!is.na(logFC))

 # freeing up space
 #rm(gene_design,fit_mut,mutation_matrix,mutation_dataset)
 
# Prediction System (Working) -----------------------------------------

# Separate features (X) and target (y)
# Assuming 'Gene' column represents gene names and inflammation labels are in another column named 'inflammation_score'
X <- gene_expected_filtered %>% select(-external_gene_name,-entrezgene_id)  # Drop gene names and target columns
y <- gene_expected_filtered$logFC

# Split data into training and testing sets
set.seed(42)
train_index <- createDataPartition(y, p = 0.8, list = FALSE)
X_train <- X[train_index, ]
X_test <- X[-train_index, ]
y_train <- y[train_index]
y_test <- y[-train_index]

# Apply Random Forest for feature selection
rf <- randomForest(X_train, y_train, ntree = 100, importance = TRUE)
important_features <- varImp(rf, scale = FALSE)
selected_genes <- rownames(important_features)[important_features$Overall > mean(important_features$Overall)]

# Filter training and test sets to include only important features
X_train_selected <- X_train %>% select(all_of(selected_genes))
X_test_selected <- X_test %>% select(all_of(selected_genes))

# Train a Random Forest on the selected features
rf_selected <- randomForest(X_train_selected, y_train, ntree = 100)

# Predict and calculate inflammation score on test data
y_pred <- predict(rf_selected, X_test_selected)
mse_score <- mean((y_test - y_pred)^2)

cat("Inflammation Score (Mean Squared Error):", mse_score, "\n")


stop()

# RETIRED CODE ---------------------------------------------------------------

# Load required libraries
library(randomForest)
library(caret)
library(dplyr)

# Load datasets
gene_expected <- read.csv("path/to/gene_expected.csv")  # Update path as necessary
gene_interest <- read.csv("path/to/gene_interest.csv")

gene_interest <- read.csv("Data/Private Data (from Lab)/interest_tidy_genes.csv")

# Filter gene_expected to include only genes in gene_interest
genes_of_interest <- gene_interest$external_gene_name  
gene_expected_filtered <- gene_expected %>% filter(external_gene_name %in% genes_of_interest)

# Assuming gene_interest is used to create the inflammation score
# Here, I assume we have an inflammation label or score variable
# in gene_expected. If not, I will need to create this based on gene_interest.

# Separate features (X) and target (y)
# Assuming 'Gene' column represents gene names and target column is 'inflammation_score'
X <- gene_expected %>% select(-Gene, -inflammation_score)  # Exclude gene names and target column
y <- gene_expected$inflammation_score  # Replace with your actual target variable

# Split data into training and testing sets
set.seed(42)
train_index <- createDataPartition(y, p = 0.8, list = FALSE)
X_train <- X[train_index, ]
X_test <- X[-train_index, ]
y_train <- y[train_index]
y_test <- y[-train_index]

# Train a Random Forest model on all genes without feature selection
rf_model <- randomForest(X_train, y_train, ntree = 100)

# Predict inflammation score on the test data
y_pred <- predict(rf_model, X_test)
mse_score <- mean((y_test - y_pred)^2)

cat("Inflammation Score (Mean Squared Error):", mse_score, "\n")

