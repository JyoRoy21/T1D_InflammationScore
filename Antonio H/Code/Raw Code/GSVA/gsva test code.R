### GSVA Test Code
## Antonio Holmes
# 05/29/2024

# Install & Load Packages -------------------------------------------------

# # GSVA Package
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("GSVA")
# 
# # DPLYR Package
# install.packages("dplyr")
# 
# #TIDYVERSE Package
# install.packages("tidyverse")

# Load necessary libraries
library(GSVA)
library(Biobase)
library(GSEABase)
library(tidyverse,dplyr)
library(readr)
library(readxl)
library(msigdbr)



# Retrieving Data (TEMPORARY) ----------------------------------------------

# Load gene expression data 
gene_expected_count <- read_csv("Data/gene_expected_count.annot.csv") # all in one dataset

#Filtering only the genes of interest
gene_expected<-gene_expected_count %>% 
  filter(entrezgene_id%in%genes_of_interest) # add the rest

#Free up some space and RAM
rm(gene_expected_count)


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


# NOR Sample 1
nor_sample_1 <-gene_expected %>% dplyr::select(entrezgene_id,external_gene_name,`6872-JK-26`,`6872-JK-31`,`6872-JK-36`,`6872-JK-41`,`6872-JK-46`) %>% 
  dplyr::rename("NOR_Sample1_t1"="6872-JK-26","NOR_Sample1_t2"="6872-JK-31","NOR_Sample1_t3"="6872-JK-36","NOR_Sample1_t4"="6872-JK-41","NOR_Sample1_t5"="6872-JK-46")

# NOR Sample 2 == Seq 42 is missing from Gene_Expected
nor_sample_2 <-gene_expected %>% dplyr::select(entrezgene_id,external_gene_name,`6872-JK-27`,`6872-JK-32`,`6872-JK-37`,#`6872-JK-42`,
                                               `6872-JK-47`) %>% 
  dplyr::rename("NOR_Sample2_t1"="6872-JK-27","NOR_Sample2_t2"="6872-JK-32","NOR_Sample2_t3"="6872-JK-37","NOR_Sample2_t5"="6872-JK-47")


# NOR Sample 3 == Seq 28 missing from Gene_Expected
nor_sample_3 <-gene_expected %>% dplyr::select(entrezgene_id,external_gene_name,#`6872-JK-28`,
                                               `6872-JK-33`,`6872-JK-38`,`6872-JK-43`,`6872-JK-48`)%>% 
  dplyr::rename("NOR_Sample3_t2"="6872-JK-33","NOR_Sample3_t3"="6872-JK-38","NOR_Sample3_t4"="6872-JK-43","NOR_Sample3_t5"="6872-JK-48")


# NOR Sample 4
nor_sample_4 <-gene_expected %>% dplyr::select(entrezgene_id,external_gene_name,`6872-JK-29`,`6872-JK-34`,`6872-JK-39`,`6872-JK-44`,`6872-JK-49`) %>% 
  dplyr::rename("NOR_Sample4_t1"="6872-JK-29","NOR_Sample4_t2"="6872-JK-34","NOR_Sample4_t3"="6872-JK-39","NOR_Sample4_t4"="6872-JK-44","NOR_Sample4_t5"="6872-JK-49")

# NOR Sample 5
nor_sample_5 <-gene_expected %>% dplyr::select(entrezgene_id,external_gene_name,`6872-JK-30`,`6872-JK-35`,`6872-JK-40`,`6872-JK-45`,`6872-JK-50`) %>% 
  dplyr::rename("NOR_Sample5_t1"="6872-JK-30","NOR_Sample5_t2"="6872-JK-35","NOR_Sample5_t3"="6872-JK-40","NOR_Sample5_t4"="6872-JK-45","NOR_Sample5_t5"="6872-JK-50")

# Combining Data sets
NOR_dataset<- nor_sample_1 %>% 
  full_join(nor_sample_2) %>% full_join(nor_sample_3) %>%
  full_join(nor_sample_4) %>% full_join(nor_sample_5)

#removing the NOD sample datasets for now
rm(nor_sample_1,nor_sample_2,nor_sample_3,nor_sample_4,nor_sample_5)

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

# Full Dataset
gene_expected_count<-NOD_dataset %>% full_join(NOR_dataset) %>% full_join(NS_dataset)

# Load & Process Data -----------------------------------------------------

## USING THE TIDY DATASET (GENE_EXPECTED) FROM DATA_TIDYING CODE 

# GSVA requires a matrix so make sure there are no duplicates
sum(duplicated(gene_expected_count$entrezgene_id)) # = 34235

#Mutate duplicate column for Entrezgene IDs
gene_expected_count<-gene_expected_count %>% mutate("Duplicate1"=duplicated(entrezgene_id),
                                                    "Duplicate2"=duplicated(external_gene_name))

#Filtering out the duplcates
gene_expected_count<-gene_expected_count %>% filter(Duplicate1==FALSE,
                                                    Duplicate2==FALSE)

# Checking once again
sum(duplicated(gene_expected_count$entrezgene_id)) # should = 0

# MAKE THIS THE 18 + NEW 9 GENES
#gene_expected_count<-gene_expected_count %>% filter(external_gene_name%in%c("Ptgs2","Nos1","Nos2","Sod1","Sod2","Sod3","Cat","Gpx1","Gpx2","Gsr"))

# Matrix for GSVA ---------------------------------------------------------

# IF YOU RUN PREVIOUS CHUNK, change gene_expected to filtered_gene_expected
gene_expected_matrix <- gene_expected_count %>%
  # GSVA can't the Ensembl IDs so we should drop this column as well as the means
  dplyr::select(-Duplicate1,-Duplicate2,-entrezgene_id) %>%
  # We need to store our gene identifiers as row names
  tibble::column_to_rownames("external_gene_name") %>%
  # Now we can convert our object into a matrix
  as.matrix()

# Results ---------------------------------------------------------

# May need to run GS Algorithm BS script from the Gene Set & Pathways section

# Making hallmarks list
gsva_sets_hallmark <- msigdbr(species="Mus musculus", category="H") %>% filter(entrez_gene%in%genes_of_interest)

# How to split the data to get what I am looking for
gsva_pwl_hallmark <- split(gsva_sets_hallmark$gene_symbol,gsva_sets_hallmark$gs_cat) # A total score
#gsva_pwl_hallmark <- split(gsva_sets_hallmark$gene_symbol,gsva_sets_hallmark$gs_name) #shows by HALLMARK
#gsva_pwl_hallmark <- split(gsva_sets_hallmark$gene_symbol,gsva_sets_hallmark$gene_symbol) #shows by GENE

# Param used in gsva()
gsva_param<-GSVA::gsvaParam(gene_expected_matrix,gsva_pwl_hallmark)

#GSVA Function Results
gsva_results<-GSVA::gsva(gsva_param)


# Print 10 rows (checking the data out)
head(gsva_results[, 1:10])

# # GSVA over Time Plots (Visuals) ---------------------------------------------------------
# 
# # Open a pdf file
# #pdf("GSVA Score vs Timepoints.pdf") 
# 
# 
# ### Sample 1 Plot ###
# 
# 
# # Multiple line graph *(ONLY INCLUDES SAMPLE 1 OF NOD, NOR, AND NS)
# x <- c(1,2,3,4,5)
# y1 <- gsva_results[,1:5]
# y2 <- gsva_results[,26:30]
# y3 <- gsva_results[,49:53]
# 
# # Plot multiple lines using matplot
# matplot(x, cbind(y1, y2, y3), type = "l", lty = 1,
#         col = c("red", "blue", "green"), xlab = "Timepoints",
#         ylab = "GSVA Score", main = "Sample 1 Line Plot")
# 
# legend("topleft", legend = c("NOD", "NOR", "NS"),
#        col = c("red", "blue", "green"),
#        lty = 1, xpd = TRUE,cex = .7)
# 
# 
# ### Sample 2 Plot ###
# 
# 
# # Multiple line graph *(ONLY INCLUDES SAMPLE 2 OF NOD, NOR, AND NS)
# y1 <- gsva_results[,6:10]
# y2 <- gsva_results[,c(31,32,33,34,NA)]
# y3 <- gsva_results[,54:58]
# 
# # Plot multiple lines using matplot
# matplot(x, cbind(y1, y2, y3), type = "l", lty = 1,
#         col = c("red", "blue", "green"), xlab = "Timepoints",
#         ylab = "GSVA Score", main = "Sample 2 Line Plot")
# 
# legend("topright", legend = c("NOD", "NOR", "NS"),
#        col = c("red", "blue", "green"),
#        lty = 1,xpd = TRUE,cex = .6)
# 
# 
# ### Sample 3 Plot ###
# 
# 
# # Multiple line graph *(ONLY INCLUDES SAMPLE 3 OF NOD, NOR, AND NS)
# y1 <- gsva_results[,11:15]
# y2 <- gsva_results[,35:38]
# y3 <- gsva_results[,59:63]
# 
# # Plot multiple lines using matplot
# matplot(x, cbind(y1, y2, y3), type = "l", lty = 1,
#         col = c("red", "blue", "green"), xlab = "Timepoints",
#         ylab = "GSVA Score", main = "Sample 3 Line Plot")
# 
# legend("topleft", legend = c("NOD", "NOR", "NS"),
#        col = c("red", "blue", "green"),
#        lty = 1, xpd = TRUE,cex = .7)
# 
# 
# ### Sample 4 Plot ###
# 
# 
# # Multiple line graph *(ONLY INCLUDES SAMPLE 4 OF NOD, NOR, AND NS)
# y1 <- gsva_results[,16:20]
# y2 <- gsva_results[,39:43]
# y3 <- gsva_results[,64:68]
# 
# # Plot multiple lines using matplot
# matplot(x, cbind(y1, y2, y3), type = "l", lty = 1,
#         col = c("red", "blue", "green"), xlab = "Timepoints",
#         ylab = "GSVA Score", main = "Sample 4 Line Plot")
# 
# legend("top", legend = c("NOD", "NOR", "NS"),
#        col = c("red", "blue", "green"),
#        lty = 1, xpd = TRUE,cex = .7)
# 
# 
# ### Sample 5 Plot ###
# 
# 
# # Multiple line graph *(ONLY INCLUDES SAMPLE 4 OF NOD, NOR, AND NS)
# y1 <- gsva_results[,21:25]
# y2 <- gsva_results[,40:44]
# y3 <- NA
# 
# # Plot multiple lines using matplot
# matplot(x, cbind(y1, y2, y3), type = "l", lty = 1,
#         col = c("red", "blue", "green"), xlab = "Timepoints",
#         ylab = "GSVA Score", main = "Sample 5 Line Plot")
# 
# legend("topleft", legend = c("NOD", "NOR", "NS"),
#        col = c("red", "blue", "green"),
#        lty = 1, xpd = TRUE,cex = .7)
# 
#  ### All Samples Together ###
# x <- c(1,2,3,4,5)
# # NOD
# y1 <- gsva_results[,1:5]
# y2 <- gsva_results[,6:10]
# y3 <- gsva_results[,11:15]
# y4 <- gsva_results[,16:20]
# y5 <- gsva_results[,21:25]
# # NOR
# y6 <- gsva_results[,26:30]
# y7 <- gsva_results[,c(31,32,33,34,NA)]
# y8 <- gsva_results[,35:38]
# y9 <- gsva_results[,39:43]
# y10 <- gsva_results[,40:44]
# # NS
# y11 <- gsva_results[,49:53]
# y12 <- gsva_results[,54:58]
# y13 <- gsva_results[,59:63]
# y14 <- gsva_results[,64:68]
# y15 <- NA
# 
# 
# # Plot multiple lines using matplot
# matplot(x, cbind(y1, y2, y3, y4, y5,
#                  y6, y7, y8, y9, y10,
#                  y11, y12, y13, y14, y15), type = "l", lty = 1,
#         col = c("red","red","red","red","red",
#                 "blue","blue","blue","blue","blue",
#                 "green","green","green","green","green"), xlab = "Timepoints",
#         ylab = "Inflammation Score", main = "All Samples Line Plot")
# 
# legend("topleft", legend = c("NOD", "NOR", "NS"),
#        col = c("red", "blue", "green"),
#        lty = 1, xpd = TRUE,cex = .6)
# 
# rm(y1, y2, y3, y4, y5,y6, y7, y8, y9, y10, y11, y12, y13, y14, y15)
# 
# # Close the pdf file
# #dev.off() 
# 
# 

# Average Mean w/ SD (All Samples) ----------------------------------------

### All Samples Together ###
x <- c(1,2,3,4,5)

# Calculate means and standard deviations

mean_NOD <- rowMeans(cbind(gsva_results[,1:5], gsva_results[,6:10], gsva_results[,11:15], gsva_results[,16:20], gsva_results[,21:25]), na.rm=TRUE)

mean_NOR <- rowMeans(cbind(gsva_results[,26:30], gsva_results[,c(31,32,33,34,NA)], gsva_results[,35:38], gsva_results[,39:43], gsva_results[,40:44]), na.rm=TRUE)

mean_NS <- rowMeans(cbind(gsva_results[,49:53], gsva_results[,54:58], gsva_results[,59:63], gsva_results[,64:68], NA), na.rm=TRUE)

#
sd_NOD <- apply(cbind(gsva_results[,1:5], gsva_results[,6:10], gsva_results[,11:15], gsva_results[,16:20], gsva_results[,21:25]), 1, sd, na.rm=TRUE)
sd_NOR <- apply(cbind(gsva_results[,26:30], gsva_results[,c(31,32,33,34,NA)], gsva_results[,35:38], gsva_results[,39:43], gsva_results[,40:44]), 1, sd, na.rm=TRUE)
sd_NS <- apply(cbind(gsva_results[,49:53], gsva_results[,54:58], gsva_results[,59:63], gsva_results[,64:68], NA), 1, sd, na.rm=TRUE)

# Create a data frame for ggplot
df <- data.frame(
  Timepoint = rep(x, 3),
  Mean = c(mean_NOD, mean_NOR,
           mean_NS),
  SD = c(sd_NOD, sd_NOR, sd_NS),
  Group = rep(c("NOD", "NOR", "NS"), each = length(x))
)

# Plot with ggplot2
avg_inflammation_scale <- ggplot(df, aes(x = Timepoint, y = Mean, color = Group)) +
  geom_line(linewidth = 1) +
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width = 0.2, linewidth = 0.8) +
  scale_color_manual(values = c("NOD" = "red", "NOR" = "white",
                                "NS" = "blue")) +
  labs(
    title = "Inflammation Scores Across Timepoints",
    subtitle = "All Samples",
    x = "Timepoints",
    y = "Inflammation Score",
    color = "Group"
  ) +
  theme_minimal(base_size = 15) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "top"
  )

# Print the plot
print(avg_inflammation_scale)


