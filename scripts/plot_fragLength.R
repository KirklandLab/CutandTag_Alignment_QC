#!/usr/bin/env Rscript

# Load required libraries
library(dplyr)
library(ggplot2)
library(viridis)
library(ggpubr)

# Capture command-line arguments
args <- commandArgs(trailingOnly = TRUE)
input_files <- args[-length(args)]  # All but last argument are input files
output_dir <- args[length(args)]    # Last argument is the output directory

# Ensure the output directory exists
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Extract sample names and unique histone names
sampleList <- gsub("_fragmentLen.txt", "", basename(input_files))
histList <- unique(sapply(strsplit(sampleList, "_"), `[`, 2))

# Initialize an empty data frame to store fragment length results
fragLen <- data.frame()

# Loop through each file to collect fragment length data
for (file_path in input_files) {
  # Read the fragment length data from each file
  fragData <- read.table(file_path, header = FALSE, col.names = c("fragLen", "fragCount"))
  
  # Extract histone and replicate information from file name
  sample_name <- gsub("_fragmentLen.txt", "", basename(file_path))
  histInfo <- strsplit(sample_name, "_")[[1]]
  histone <- histInfo[1]
  replicate <- histInfo[2]
  
  # Add calculated columns and combine with previous data
  fragData <- fragData %>%
    mutate(Weight = fragCount / sum(fragCount), Histone = histone, Replicate = replicate, sampleInfo = sample_name)
  fragLen <- rbind(fragLen, fragData)
}

# Convert columns to factors for ordered plotting
fragLen$sampleInfo <- factor(fragLen$sampleInfo, levels = sampleList)
fragLen$Histone <- factor(fragLen$Histone, levels = histList)

# Generate the fragment size density plot (violin plot)
fig5A <- ggplot(fragLen, aes(x = sampleInfo, y = fragLen, weight = Weight, fill = Histone)) +
  geom_violin(bw = 5) +
  scale_y_continuous(breaks = seq(0, 800, 50)) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
  theme_bw(base_size = 20) +
  ggpubr::rotate_x_text(angle = 20) +
  ylab("Fragment Length") +
  xlab("")

# Generate the fragment count plot
fig5B <- ggplot(fragLen, aes(x = fragLen, y = fragCount, color = Histone, group = sampleInfo, linetype = Replicate)) +
  geom_line(linewidth = 1) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma") +
  theme_bw(base_size = 20) +
  xlab("Fragment Length") +
  ylab("Count") +
  coord_cartesian(xlim = c(0, 500))

# Arrange the two plots side by side
final_plot <- ggarrange(fig5A, fig5B, ncol = 2)

# Save the arranged plot to a file in the specified output directory
output_file <- file.path(output_dir, "fragment_length_plot.png")
ggsave(output_file, final_plot, width = 14, height = 10)
