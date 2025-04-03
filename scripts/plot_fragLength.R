#!/usr/bin/env Rscript

# Load required libraries
library(dplyr)
library(ggplot2)
library(viridis)
library(ggpubr)

# -------------------------
# Capture command-line arguments
# -------------------------
args <- commandArgs(trailingOnly = TRUE)
input_files <- args[1:(length(args) - 2)]        # All but last 2 args are input files
metadata_file <- args[length(args) - 1]          # Second to last = samples.csv
output_dir <- args[length(args)]                 # Last arg = output directory

# -------------------------
# Ensure the output directory exists
# -------------------------
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# -------------------------
# Load sample metadata
# -------------------------
sample_metadata <- read.csv(metadata_file, stringsAsFactors = FALSE)

# Extract sample names from input file paths
sampleList <- gsub("_fragmentLen.txt", "", basename(input_files))

# Initialize data frame
fragLen <- data.frame()

# -------------------------
# Loop through fragment length files
# -------------------------
for (file_path in input_files) {
  sample_name <- gsub("_fragmentLen.txt", "", basename(file_path))
  
  # Lookup metadata
  histone <- sample_metadata$histone[sample_metadata$sample == sample_name]
  replicate <- sample_metadata$replicate[sample_metadata$sample == sample_name]
  
  # Read fragment length data
  fragData <- read.table(file_path, header = FALSE, col.names = c("fragLen", "fragCount"))
  
  # Annotate and add to combined dataframe
  fragData <- fragData %>%
    mutate(
      Weight = fragCount / sum(fragCount),
      Histone = histone,
      Replicate = replicate,
      sampleInfo = sample_name
    )
  
  fragLen <- rbind(fragLen, fragData)
}

# Convert to factors
fragLen$sampleInfo <- factor(fragLen$sampleInfo, levels = sampleList)
fragLen$Histone <- factor(fragLen$Histone, levels = unique(fragLen$Histone))
fragLen$Replicate <- as.factor(fragLen$Replicate)

# -------------------------
# Generate Plots
# -------------------------

# Violin plot (fragment size distribution)
fig5A <- ggplot(fragLen, aes(x = sampleInfo, y = fragLen, weight = Weight, fill = Histone)) +
  geom_violin(bw = 5) +
  scale_y_continuous(breaks = seq(0, 800, 50)) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
  theme_bw(base_size = 20) +
  ggpubr::rotate_x_text(angle = 20) +
  ylab("Fragment Length") +
  xlab("")

# Line plot (fragment counts)
fig5B <- ggplot(fragLen, aes(x = fragLen, y = fragCount, color = Histone, group = sampleInfo, linetype = Replicate)) +
  geom_line(linewidth = 1) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma") +
  theme_bw(base_size = 20) +
  xlab("Fragment Length") +
  ylab("Count") +
  coord_cartesian(xlim = c(0, 500))

# Combine and save
final_plot <- ggarrange(fig5A, fig5B, ncol = 2)
output_file <- file.path(output_dir, "fragment_length_plot.png")
ggsave(output_file, final_plot, width = 14, height = 10)
