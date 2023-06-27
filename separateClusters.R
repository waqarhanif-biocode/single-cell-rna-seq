# Read the data from the file
data <- read.csv("marker_genes_combined.csv", header = TRUE)

# Group the data by the 'cluster' column
grouped <- split(data, data$cluster)

# Save each group as a separate CSV file
for (i in seq_along(grouped)) {
  cluster <- names(grouped[i])
  filename <- paste0(cluster, ".csv")
  write.csv(grouped[[i]], file = filename, row.names = FALSE)
}
