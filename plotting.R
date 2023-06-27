# Load packages
library(ggplot2)

# Function to generate plot and save results
generate_combined_plot <- function(up_file, down_file, output_file) {
  up_data <- read.csv(up_file)
  up_data <- up_data[order(up_data$P.value),]
  up_data <- up_data[1:10,]
  
  down_data <- read.csv(down_file)
  down_data <- down_data[order(down_data$P.value),]
  down_data <- down_data[1:10,]
  
  # Combine up and down data
  combined_data <- rbind(up_data, down_data)
  combined_data$is_upregulated <- ifelse(combined_data$Term %in% up_data$Term, TRUE, FALSE)
  combined_data <- combined_data[combined_data$P.value < 0.05, ]
  # Shorten labels
  combined_data$Term <- substr(combined_data$Term, 1, 50)
  # Get the cluster name from the file name
  cluster_name <- gsub("_up.csv", "", basename(up_file))
  
  # Horizontal version
  p <- ggplot(combined_data, aes(x = Term, y = P.value)) +
    geom_segment(aes(x = Term, xend = Term, y = 0, yend = P.value, color = is_upregulated), size = 1.2) +
    geom_point(aes(color = is_upregulated), size = 4, alpha = 1) +
    theme_light() +
    coord_flip() +
    xlab("Pathway") +
    ylab("P-value") +
    ggtitle(cluster_name) +
    scale_x_discrete(position = "top") +
    scale_color_manual(values = c("#439A97", "skyblue")) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.border = element_blank(),
      axis.ticks.y = element_blank(),
      legend.position = "top"
    )
  
  ggsave(output_file, plot = p, dpi = 600, width = 6.5, height = 4.5)
}

# Directory where the files are located
directory <- getwd()

# List all the up and down files
up_files <- list.files(directory, pattern = "_up.csv", full.names = TRUE)
down_files <- list.files(directory, pattern = "_dn.csv", full.names = TRUE)

# Iterate over the files and generate the combined plot
for (i in seq_along(up_files)) {
  up_file <- up_files[i]
  down_file <- down_files[i]
  
  # Extract the cluster name from the file name
  cluster <- gsub("_up.csv", "", basename(up_file))
  
  # Specify the output file name
  output_file <- paste0("combined_biological-process_", cluster, ".png")
  
  # Generate the combined plot and save the results
  generate_combined_plot(up_file, down_file, output_file)
}
