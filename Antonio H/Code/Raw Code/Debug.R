


# Why it does not matter if there are more than one row of genes ----------

#https://www.bioconductor.org/packages/release/bioc/vignettes/GSVA/inst/doc/GSVA.html


# GSVA Code that works (EX CODE) ----------------------------------------------------

BiocManager::install("limma")

library(limma)

p <- 10 ## number of genes
n <- 30 ## number of samples
nGrp1 <- 15 ## number of samples in group 1
nGrp2 <- n - nGrp1 ## number of samples in group 2

## consider three disjoint gene sets
geneSets <- list(set1=paste("g", 1:3, sep=""),
                 set2=paste("g", 4:6, sep=""),
                 set3=paste("g", 7:10, sep=""))

## sample data from a normal distribution with mean 0 and st.dev. 1
y <- matrix(rnorm(n*p), nrow=p, ncol=n,
            dimnames=list(paste("g", 1:p, sep="") , paste("s", 1:n, sep="")))

## genes in set1 are expressed at higher levels in the last 'nGrp1+1' to 'n' samples
y[geneSets$set1, (nGrp1+1):n] <- y[geneSets$set1, (nGrp1+1):n] + 2

## build GSVA parameter object
gsvapar <- gsvaParam(y, geneSets, maxDiff=TRUE)

## estimate GSVA enrichment scores for the three sets
gsva_es <- gsva(gsvapar)

## fit the same linear model now to the GSVA enrichment scores
fit <- lmFit(gsva_es, design)

## estimate moderated t-statistics
fit <- eBayes(fit)

## set1 is differentially expressed
topTable(fit, coef="sampleGroup2vs1")

# Making Hallmarks List like GeneSets (FUNCTIONING)-------------------------------------

sample_test<-sets_hallmark %>% filter(entrez_gene%in%genes)

test_hallmark <- split(sample_test$gene_symbol,sample_test$gs_name)

gsvapa_test <- gsvaParam(gene_expected_matrix, test_hallmark, maxDiff=FALSE)

gsva_test <- gsva(gsvapa_test)




# Log FC & P-Val (EX CODE) -----------------------------------------

## sample data from a normal distribution with mean 0 and st.dev. 1
y <- matrix(rnorm(n*p), nrow=p, ncol=n,
            dimnames=list(paste("g", 1:p, sep="") , paste("s", 1:n, sep="")))

## genes in set1 are expressed at higher levels in the last 'nGrp1+1' to 'n' samples
y[geneSets$set1, (nGrp1+1):n] <- y[geneSets$set1, (nGrp1+1):n] + 2

## build design matrix
design <- cbind(sampleGroup1=1, sampleGroup2vs1=c(rep(0, nGrp1), rep(1, nGrp2)))

## fit linear model
fit <- lmFit(y, design)

## estimate moderated t-statistics
fit <- eBayes(fit)

## genes in set1 are differentially expressed
topTable(fit, coef="sampleGroup2vs1")


# Making Log FC & P-Val (FUNCTIONING W/IN CORRECT SCRIPT) ---------------------------------

## I TRANSFERED THIS OVER TO DD DENESET ALGORITHM, I NEED TO MAKE A NEW MATRIX

# DATASET W/ ALL SAMPLES WITH HLA MUTATION
mutation_dataset<- NOD_dataset %>% full_join(NS_dataset)

n_samples <- 45 ## number of samples
nod_Grp1 <- 25 ## number of samples in group 1
ns_Grp2 <- n_samples - nod_Grp1 ## number of samples in group 2


## the matrix needed for lmFit()
gene_expected_matrix

## not sure if I need this yet
#gene_expected_matrix[hallmarks_list$set1, (nGrp1+1):n] <- gene_expected_matrix[hallmarks_list$set1, (nGrp1+1):n] + 2

## build design matrix
gene_design <- cbind(sampleGroup1=1, sampleGroup2vs1=c(rep(0, nod_Grp1), rep(1, ns_Grp2)))

## fit linear model
#fit <- 
  lmFit(gene_expected_matrix, gene_design)














# Interaction Graph w/ pathways (FUNCTIONING) -------------------------------------

install.packages(c("igraph", "ggplot2"))
BiocManager::install("clusterProfiler")
BiocManager::install("ReactomePA",force = TRUE)
library(igraph)
library(clusterProfiler)
library(ReactomePA)
library(ggplot2)
library(msigdbr)


debug_sets_hallmark <- msigdbr(species="Mus musculus", category="H") %>%
  mutate(hallmark_true= str_detect(sets_hallmark$gs_name,
                                   pattern = "HALLMARK")) %>% 
    filter(hallmark_true==TRUE) %>% 
  filter(entrez_gene%in%NOD_dataset$entrezgene_id)

# Example gene list
genes <- debug_sets_hallmark$entrez_gene

# Pathway enrichment analysis
enriched_pathways <- enrichPathway(gene = genes, organism = "mouse")

# Extract pathway and gene information
pathways <- as.data.frame(enriched_pathways)

# CAN I REUSE THE CODE BELOW PIPED INTO THE PATHFINDR FUNX

`# Create an edge list for the graph
edge_list <- data.frame(
  from = rep(pathways$Description, each = length(genes)),
  to = rep(genes, times = nrow(pathways)))

# Create the graph
g <- graph_from_data_frame(edge_list, directed = FALSE)

plot(g, vertex.label.cex = 0.8, vertex.size = 5)



  


# Interaction Graph w/out pathways (FUNCTIONING) ----------------------

#
gene_data <- data.frame(
  GeneIDs = pathways$geneID)

#
result <- data.frame(from = character(), to = character(), stringsAsFactors = FALSE)

# Loop through each row in the dataset
for (i in 1:nrow(pathways)) {
  # Split the GeneIDs string into individual gene IDs
  gene_ids <- unlist(strsplit(gene_data$GeneIDs[i], "/"))
  
  # Create the from and to columns
  from_to_pairs <- data.frame(
    from = gene_ids[-length(gene_ids)], # All but the last element
    to = gene_ids[-1]                   # All but the first element
  )
  
  # Add the from-to pairs to the result data frame
  result <- rbind(result, from_to_pairs)

}

# Create a graph object
gene_network <- graph_from_data_frame(d = result, directed = FALSE)


# Plot the gene network
plot(gene_network, vertex.label = V(gene_network)$name, vertex.size = 3, main = "Gene Interaction Network")



# Exploring the results dataset with gene connection ----------------------------------------

genes_of_interest<-c("14776","14775","14782","12359","18125","18126","19225","20655","20656","20657")

result<-result %>% dplyr::filter(to%in%genes_of_interest|
                           from%in%genes_of_interest)

# Create a graph object
gene_network <- graph_from_data_frame(d = result, directed = FALSE)


# Plot the gene network
plot(gene_network, vertex.label = V(gene_network)$name, vertex.size = 3, main = "Gene Interaction Network")





# Retrieving all 27 genes -------------------------------------------------

# Retrieving the mouse data
#genesets_hallmark<-msigdbr(species="Mus musculus", category="H") # only contains 19 out of 27
genesets_all<-msigdbr(species = "Mus musculus") # only contains 26 out of 27

#
genes_of_interest<-gene_expected_count$entrezgene_id

# Only the 27 genes we are interest in
pathway2 <- genesets_all %>% filter(entrez_gene%in%genes_of_interest)

# Transforming the data frame
pathway2 <- pathway2 %>% dplyr::select(gs_name, entrez_gene) %>% 
  mutate(presence = 1) %>%
  dplyr::rename(pathway = gs_name) %>% 
  pivot_wider(names_from = pathway,values_from = presence,values_fill = list(presence = 0))

# Calculate the degree centrality
network_matrix <- pathway2 %>%
  rowwise() %>%
  mutate(degree_centrality = sum(c_across(-entrez_gene))) %>%
  ungroup() 

#everything below network_matrix in gene-gene matrix should work

# A better Avg graph ------------------------------------------------------

# Load necessary libraries
library(ggplot2)
library(dplyr)

# Create a data frame for ggplot
df <- data.frame(
  Timepoint = rep(x, 3),
  Mean = c(mean_NOD, mean_NOR, mean_NS),
  SD = c(sd_NOD, sd_NOR, sd_NS),
  Group = rep(c("NOD", "NOR", "NS"), each = length(x))
)

# Plot with ggplot2
avg_inflammation_scale <- ggplot(df, aes(x = Timepoint, y = Mean, color = Group)) +
  geom_line(size = 1) +
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width = 0.2, linewidth = 0.8) +
  scale_color_manual(values = c("NOD" = "red", "NOR" = "blue", "NS" = "green")) +
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

# To-Do List --------------------------------------------------------------

#
# change correlation cutoff to >0.6 and <-0.6
#
# count gene pair connections (how freq is a gene is geneset 2 popping up)
#
# Add new genes to IS






