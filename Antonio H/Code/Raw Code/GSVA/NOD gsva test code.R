### NOD sample 1 GSVA Test Code
## Antonio Holmes
# 06/18/2024

#Load Data & Req Packages  -------------------------------------------------

#Running Data Tidying Script
source("/Users/antonioholmes/project/Code/Raw Code/data_tidying.R")

#Installing GSVA for this script
install.packages("BiocManager")
BiocManager::install("GSVA", version="3.19", force=TRUE)

#Loading the required packages
library(GSVA)
library(msigdbr)

# Load & Process Data -----------------------------------------------------

## USING THE TIDY DATASET (GENE_EXPECTED) FROM DATA_TIDYING CODE 

# GSVA requires a matrix so make sure there are no duplicates
sum(duplicated(NOD_dataset$entrezgene_id)) # = 34235

#Mutate duplicate column for Entrezgene IDs
NOD_dataset<-NOD_dataset %>% mutate("Duplicate"=duplicated(entrezgene_id))

#Filtering out the duplcates
NOD_dataset<-NOD_dataset %>% filter(Duplicate==FALSE)

# Checking once again
sum(duplicated(NOD_dataset$entrezgene_id)) # should = 0

# Matrix for GSVA ---------------------------------------------------------

# 
gene_expected_matrix <- NOD_dataset %>%
  # GSVA can't the Ensembl IDs so we should drop this column as well as the means
  dplyr::select(-Duplicate,-entrezgene_id) %>%
  # We need to store our gene identifiers as row names
  tibble::column_to_rownames("external_gene_name") %>%
  # Now we can convert our object into a matrix
  as.matrix()

# Results ---------------------------------------------------------

# Making hallmarks list
nod_sets_hallmark <- msigdbr(species="Mus musculus", category="H") %>% filter(entrez_gene%in%genes_of_interest)
nod_pwl_hallmark <- split(nod_sets_hallmark$gene_symbol,nod_sets_hallmark$gs_name)

# Param used in gsva()
nod_gsva_param<-GSVA::gsvaParam(gene_expected_matrix,nod_pwl_hallmark)

#GSVA Function Results
nod_gsva_results<-GSVA::gsva(nod_gsva_param)

# Print 10 rows (checking the data out)
head(nod_gsva_results[, 1:10])

# [WORK ON THIS, PURE TRASH] Plots (Visuals) ---------------------------------------------------------

# Create a data frame with the column names classified as "NOD" or "NS"
my_gene_col <- data.frame(group = ifelse(grepl("NOD_",
                                               colnames(nod_gsva_results)), 
                                         yes = "NOD", no = "NS"))

# Set row names of the metadata to match the column names of the gene matrix
rownames(my_gene_col) <- colnames(nod_gsva_results)

# Ensure the 'group' column is treated as a factor
my_gene_col$group <- as.factor(my_gene_col$group)

# Normalizing the data
#gene_normalized <- t(scale(t(gene_expected_matrix)))

# Plot the heatmap with reordered columns
pheatmap::pheatmap(nod_gsva_results,annotation_col = my_gene_col, cluster_cols = FALSE,   # Disable column clustering
                   cluster_rows = FALSE    # Disable row clustering (optional)
)


# # Multiple line graph
# x <- #timepoints
# y1 <- #NOD GSVA Results
# y2 <- #NOR GSVA Results
# y3 <- #NS GSVA Results
# 
# # Plot multiple lines using matplot
# matplot(x, cbind(y1, y2, y3), type = "l", lty = 1,
#         col = c("red", "blue", "green"), xlab = "X",
#         ylab = "Y", main = "Multiple Lines Plot")
# legend("topright", legend = c("NOD", "NOR", "NS"),
#        col = c("red", "blue", "green"),
#        lty = 1)