#!/usr/bin/env Rscript

# Load required libraries
library(ggplot2)
library(dplyr)
library(ggpubr)

# Define the path to the Bowtie2 summary files
projPath <- "results/alignment/sam/bowtie2_summary"

# List all summary files in the directory
summary_files <- list.files(projPath, pattern = "_bowtie2.txt", full.names = TRUE)

# Extract sample names from file names (e.g., "K27me3_rep1" from "K27me3_rep1_bowtie2.txt")
sampleList <- gsub("_bowtie2.txt", "", basename(summary_files))

# Extract unique histone names (e.g., "K27me3", "K4me3", "IgG") for ordering purposes
histList <- unique(sapply(strsplit(sampleList, "_"), `[`, 2))

# Initialize an empty data frame to store alignment results
alignResult <- data.frame()

# Loop through each file to collect alignment results
for (file_path in summary_files) {
  # Read Bowtie2 summary file
  alignRes <- read.table(file_path, header = FALSE, fill = TRUE, stringsAsFactors = FALSE)
  
  # Extract metrics from the summary file
  sequencingDepth <- as.numeric(gsub("[^0-9]", "", alignRes$V1[1]))  # Total reads
  mappedFragNum <- as.numeric(gsub("[^0-9]", "", alignRes$V1[4])) + 
    as.numeric(gsub("[^0-9]", "", alignRes$V1[5]))  # Mapped fragments
  alignmentRate <- as.numeric(gsub("[^0-9.]", "", alignRes$V1[6]))  # Alignment rate as a percentage
  unalignedRate <- 100 - alignmentRate  # Calculate unaligned rate
  
  # Extract histone and replicate information from file name
  sample_name <- gsub("_bowtie2.txt", "", basename(file_path))
  histInfo <- strsplit(sample_name, "_")[[1]]
  histone <- histInfo[2]
  replicate <- histInfo[1]
  
  # Append results to the data frame
  alignResult <- rbind(alignResult, data.frame(Histone = histone, Replicate = replicate, 
                                               SequencingDepth = sequencingDepth, 
                                               MappedFragNum = mappedFragNum, 
                                               AlignmentRate = alignmentRate, 
                                               UnalignedRate = unalignedRate))
}

# Convert Histone to a factor for ordered plotting based on histList
alignResult$Histone <- factor(alignResult$Histone, levels = histList)

# Plot 1: Sequencing Depth
fig1 <- ggplot(alignResult, aes(x = Histone, y = SequencingDepth / 1e6, fill = Histone)) +
  geom_boxplot() +
  geom_jitter(aes(color = Replicate), position = position_jitter(0.15)) +
  theme_bw(base_size = 14) +
  ylab("Sequencing Depth (Millions)") +
  xlab("") +
  ggtitle("Sequencing Depth per Histone")

# Plot 2: Mapped Fragments (Alignable Fragment Count)
fig2 <- ggplot(alignResult, aes(x = Histone, y = MappedFragNum / 1e6, fill = Histone)) +
  geom_boxplot() +
  geom_jitter(aes(color = Replicate), position = position_jitter(0.15)) +
  theme_bw(base_size = 14) +
  ylab("Mapped Fragments (Millions)") +
  xlab("") +
  ggtitle("Mapped Fragments per Histone")

# Plot 3: Alignment Rate
fig3 <- ggplot(alignResult, aes(x = Histone, y = AlignmentRate, fill = Histone)) +
  geom_boxplot() +
  geom_jitter(aes(color = Replicate), position = position_jitter(0.15)) +
  theme_bw(base_size = 14) +
  ylab("Alignment Rate (%)") +
  xlab("") +
  ggtitle("Alignment Rate per Histone")

# Plot 4: Unaligned Rate
fig4 <- ggplot(alignResult, aes(x = Histone, y = UnalignedRate, fill = Histone)) +
  geom_boxplot() +
  geom_jitter(aes(color = Replicate), position = position_jitter(0.15)) +
  theme_bw(base_size = 14) +
  ylab("Unaligned Rate (%)") +
  xlab("") +
  ggtitle("Unaligned Rate per Histone")

# Arrange the four plots into a single figure
final_plot <- ggarrange(fig1, fig2, fig3, fig4, ncol = 2, nrow = 2, common.legend = TRUE, legend = "bottom")

# Save the arranged plot to a file
ggsave("results/alignment/alignment_summary_plot.png", final_plot, width = 14, height = 10)
