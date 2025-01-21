### Interactive Graph, not a matrix but i wanted to save it in the matrix folder for quick access
## Antonio Holmes
# 06/07/2024

# Install & Load Packages ------------------------------------------------------------------

install.packages(c("igraph", "ggplot2"))
BiocManager::install("clusterProfiler")

BiocManager::install("org.Mm.eg.db",force = TRUE) #mouse
#BiocManager::install("org.Hs.eg.db",force = TRUE) #human

BiocManager::install("ReactomePA",force = TRUE)

library(igraph)
library(clusterProfiler)
library(ReactomePA)
library(org.Mm.eg.db)
#library(org.Hs.eg.db)
library(ggplot2)
library(msigdbr)
library(dplyr)
library(stringr)


# Gene Set & Pathways ------------------------------------------------------------------

# Retrieving the mouse data
genesets_hallmark<-msigdbr(species="Mus musculus", category="H")

# Pre-defined (27) gene list
genesets_hallmark<- genesets_hallmark %>% filter(entrez_gene%in%gene_expected_count$entrezgene_id) %>% distinct(entrez_gene)

# Making the Entre Gene ID independent values
genes <- genesets_hallmark$entrez_gene

# Pathway enrichment analysis
enriched_pathways <- enrichPathway(gene = genes, organism = "mouse")

# Extract pathway and gene information
pathways <- as.data.frame(enriched_pathways)

#### CAN I REUSE THE CODE BELOW PIPED INTO THE PATHFINDR FUNX

`# Create an edge list for the graph
edge_list <- data.frame(
  from = rep(pathways$Description, each = length(genes)),
  to = rep(genes, times = nrow(pathways)))

# Excluding NON-hallmark pathways [IS THIS STILL NEEDED]
exclude_pathway <- c("Biological oxidations","Arachidonic acid metabolism",
                    "Platelet homeostasis","Paracetamol ADME",
                    "Fatty acid metabolism")

# Seeing which pathways we do not want is into Edge List
 edge_list %>% distinct(from) %>% filter(from %in% exclude_pathway)

# # [IS THIS STILL NEEDED]
 edge_list<-edge_list %>% dplyr::filter(!(from %in% exclude_pathway)) %>% distinct()

# Gene-Pathway Interaction Graph ------------------------------------------------------------------

# Create the graph
graph <- graph_from_data_frame(edge_list, directed = FALSE)

# Setting the Layout for the graph
layout<-layout_with_graphopt(graph) #drl w/ vertex size 3 looks the best
                                #graphopt w/ vertex size 1 is next best

# The graph
plot(graph, layout = layout,
     vertex.label.cex = 0.8, vertex.size = 3,
     vertex.label.dist = 0.5)



# Gene-Gene Interaction Graph ---------------------------------------------

# Dataset with only the gene IDS
gene_data <- data.frame(
  GeneIDs = pathways$geneID)

# Creating a dataframe with 2 columns, from and to
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



