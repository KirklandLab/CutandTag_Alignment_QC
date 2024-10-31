#!/usr/bin/env Rscript

# Load required libraries
library(ggplot2)
library(dplyr)
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
sampleList <- gsub("_bowtie2.txt", "", basename(input_files))
histList <- unique(sapply(strsplit(sampleList, "_"), `[`, 2))

# Initialize an empty data frame to store alignment results
alignResult <- data.frame()

# Loop through each file to collect alignment results
for (file_path in input_files) {
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
  histone <- histInfo[1]
  replicate <- histInfo[2]
  
  # Append results to the data frame
  alignResult <- rbind(alignResult, data.frame(Histone = histone, Replicate = replicate, 
                                               SequencingDepth = sequencingDepth, 
                                               MappedFragNum = mappedFragNum, 
                                               AlignmentRate = alignmentRate, 
                                               UnalignedRate = unalignedRate))
}

# Debugging: Print alignResult to check assignments
print(alignResult)

# Convert Histone to a factor for ordered plotting based on histList
alignResult$Histone <- factor(alignResult$Histone, levels = histList)

# Create individual plots for each metric
fig1 <- ggplot(alignResult, aes(x = Histone, y = SequencingDepth / 1e6, fill = Histone)) +
  geom_boxplot() +
  geom_jitter(aes(color = Replicate), position = position_jitter(0.15)) +
  theme_bw(base_size = 14) +
  ylab("Sequencing Depth (Millions)") +
  xlab("") +
  ggtitle("Sequencing Depth per Histone")

fig2 <- ggplot(alignResult, aes(x = Histone, y = MappedFragNum / 1e6, fill = Histone)) +
  geom_boxplot() +
  geom_jitter(aes(color = Replicate), position = position_jitter(0.15)) +
  theme_bw(base_size = 14) +
  ylab("Mapped Fragments (Millions)") +
  xlab("") +
  ggtitle("Mapped Fragments per Histone")

fig3 <- ggplot(alignResult, aes(x = Histone, y = AlignmentRate, fill = Histone)) +
  geom_boxplot() +
  geom_jitter(aes(color = Replicate), position = position_jitter(0.15)) +
  theme_bw(base_size = 14) +
  ylab("Alignment Rate (%)") +
  xlab("") +
  ggtitle("Alignment Rate per Histone")

fig4 <- ggplot(alignResult, aes(x = Histone, y = UnalignedRate, fill = Histone)) +
  geom_boxplot() +
  geom_jitter(aes(color = Replicate), position = position_jitter(0.15)) +
  theme_bw(base_size = 14) +
  ylab("Unaligned Rate (%)") +
  xlab("") +
  ggtitle("Unaligned Rate per Histone")

# Arrange the four plots into a single figure
final_plot <- ggarrange(fig1, fig2, fig3, fig4, ncol = 2, nrow = 2, common.legend = TRUE, legend = "bottom")

# Save the arranged plot to a file in the specified output directory
output_file <- file.path(output_dir, "alignment_summary_plot.png")
ggsave(output_file, final_plot, width = 14, height = 10)
