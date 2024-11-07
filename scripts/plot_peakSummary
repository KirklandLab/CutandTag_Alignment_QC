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
sampleList <- gsub("_0.05_peaks.narrowPeak", "", basename(input_files))
histList <- unique(sapply(strsplit(sampleList, "_"), `[`, 1))

# Initialize an empty data frame to store peak information
peakData <- data.frame()

# Loop through each file to collect peak information
for (file_path in input_files) {
  # Read the peak data from each MACS2 output file
  peakInfo <- read.table(file_path, header = FALSE, col.names = c("chrom", "start", "end", "name", "score"))
  
  # Extract histone and replicate information from file name
  sample_name <- gsub("_0.05_peaks.narrowPeak", "", basename(file_path))
  histInfo <- strsplit(sample_name, "_")[[1]]
  histone <- histInfo[1]
  replicate <- histInfo[2]
  
  # Add calculated columns and combine with previous data
  peakInfo <- peakInfo %>%
    mutate(width = end - start, Histone = histone, Replicate = replicate, sampleInfo = sample_name)
  peakData <- rbind(peakData, peakInfo)
}

# Convert columns to factors for ordered plotting
peakData$sampleInfo <- factor(peakData$sampleInfo, levels = sampleList)
peakData$Histone <- factor(peakData$Histone, levels = histList)

# Plot: Number of Peaks (similar to fig7A)
fig7A <- peakData %>%
  group_by(Histone, Replicate) %>%
  summarise(peakN = n()) %>%
  ggplot(aes(x = Histone, y = peakN, fill = Histone)) +
  geom_boxplot() +
  geom_jitter(aes(color = Replicate), position = position_jitter(0.15)) +
  theme_bw(base_size = 18) +
  ylab("Number of Peaks") +
  xlab("")

# Plot: Peak Widths (similar to fig7B)
fig7B <- ggplot(peakData, aes(x = Histone, y = width, fill = Histone)) +
  geom_violin() +
  facet_grid(Replicate ~ .) +
  theme_bw(base_size = 18) +
  ylab("Width of Peaks") +
  xlab("")

# Arrange and save the final plot
final_plot <- ggarrange(fig7A, fig7B, ncol = 2)
output_file <- file.path(output_dir, "peak_summary_plot.png")
ggsave(output_file, final_plot, width = 14, height = 10)
