### Creating an Pathway Visualization
## Antonio Holmes
# 11/12/2024



# Library Packages --------------------------------------------------------

#Loading the required packages
library(GSVA)
library(msigdbr)

# Load & Process Data -----------------------------------------------------

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


# Matrix for Pathways -----------------------------------------------------

# 
mutation_matrix <- mutation_dataset %>%
  # GSVA can't the Entrez IDs so we should drop this column as well as the means
  dplyr::select(-entrezgene_id) %>%
  # We need to store our gene identifiers as row names
  tibble::column_to_rownames("external_gene_name") %>%
  # Now we can convert our object into a matrix
  as.matrix()


# Results ---------------------------------------------------------

# Making genes of interest list
genes_of_interest <- mutation_dataset$entrezgene_id

# Making hallmarks list
nod_sets_hallmark <- msigdbr(species="Mus musculus", category="H") %>% filter(entrez_gene%in%genes_of_interest)
nod_pwl_hallmark <- split(nod_sets_hallmark$gene_symbol,nod_sets_hallmark$gs_name)

# Param used in gsva()
nod_gsva_param<-GSVA::gsvaParam(mutation_matrix,nod_pwl_hallmark)

#GSVA Function Results
nod_gsva_results<-GSVA::gsva(nod_gsva_param)

# Print 10 rows (checking the data out)
head(nod_gsva_results[, 1:10])

# Heatmap ---------------------------------------------------------

# Create a data frame with the column names classified as "NOD" or "NS"
my_gene_col <- data.frame(group = ifelse(grepl("NOD_",
                                               colnames(nod_gsva_results)), 
                                         yes = "NOD", no = "NS"))

# Set row names of the metadata to match the column names of the gene matrix
rownames(my_gene_col) <- colnames(nod_gsva_results)

# Ensure the 'group' column is treated as a factor
my_gene_col$group <- as.factor(my_gene_col$group)

# Normalizing the data
#gene_normalized <- t(scale(t(mutation_matrix)))

# Define the annotation data frame to highlight "Hallmark_Inflammatory_Response"
row_annotation <- data.frame(Pathway = ifelse(rownames(nod_gsva_results) == "HALLMARK_INFLAMMATORY_RESPONSE", 
                                              "Inflammatory Response", "Normal"))
rownames(row_annotation) <- rownames(nod_gsva_results)

# Set colors for row annotation to make it stand out
ann_colors <- list(Pathway = c(`Inflammatory Response` = "red", Normal = "grey"))

# Plot the heatmap with row annotation to highlight the pathway
pheatmap::pheatmap(nod_gsva_results, 
                   annotation_col = my_gene_col, 
                   annotation_row = row_annotation, 
                   annotation_colors = ann_colors,
                   cluster_cols = FALSE,    # Disable column clustering
                   cluster_rows = FALSE)    # Disable row clustering







