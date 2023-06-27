library(devtools)
install_github("wjawaid/enrichR")
library(ggplot2)
library(dplyr)
library(enrichR)

analyze_enrichment <- function(data_file, dbs) {
  data <- read.csv(data_file, header = TRUE)
  
  unique_clusters <- unique(data$cluster)
  
  for (cluster in unique_clusters) {
    cluster_data <- data[data$cluster == cluster, ]
    
    up.idx <- which(cluster_data$p_val < 0.05 & cluster_data$avg_log2FC > 0.2)
    dn.idx <- which(cluster_data$p_val < 0.05 & cluster_data$avg_log2FC < 0.2)
    
    up.genes <- cluster_data[up.idx, ]$gene
    dn.genes <- cluster_data[dn.idx, ]$gene
    
    for (db in dbs) {
      enriched_pw_up <- enrichr(genes = up.genes, databases = db)
      write.csv(enriched_pw_up[[1]], file = paste0(db, "_", cluster, "_up.csv"))
      
      enriched_pw_dn <- enrichr(genes = dn.genes, databases = db)
      write.csv(enriched_pw_dn[[1]], file = paste0(db, "_", cluster, "_dn.csv"))
    }
  }
}

# List of CSV files to process
file_list <- c("CD4+ NKT-like Cells.csv",
               "CD8+ NKT-like Cells.csv",
               "Naive CD4+ T Cells.csv",
               "Memory CD4+ T Cells.csv",
               "Classical Monocytes.csv",
               "Myeloid Dendritic Cells.csv",
               "Naive B Cells.csv")

# Run analysis on each file
for (file in file_list) {
  analyze_enrichment(file, c("Reactome_2022", "GO_Biological_Process_2021", 
                             "GO_Cellular_Component_2021", "GO_Molecular_Function_2021"))
}
