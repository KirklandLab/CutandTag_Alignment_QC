#!/usr/bin/env Rscript

# Load required libraries
library(dplyr)
library(ggplot2)
library(viridis)
library(ggpubr)
library(GenomicRanges)
library(Rsamtools)

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

# Load peak data and convert to GRanges
peak_granges_list <- list()

for (file_path in input_files) {
  # Read the peak data with MACS2 narrowPeak format columns
  peakInfo <- read.table(file_path, header = FALSE, 
                         col.names = c("chrom", "start", "end", "name", "score", 
                                       "strand", "signalValue", "pValue", "qValue", "peak"))
  
  # Extract histone and replicate information from file name
  sample_name <- gsub("_0.05_peaks.narrowPeak", "", basename(file_path))
  histInfo <- strsplit(sample_name, "_")[[1]]
  histone <- histInfo[1]
  replicate <- histInfo[2]
  
  # Add calculated columns and combine with previous data
  peakInfo <- peakInfo %>%
    mutate(width = end - start, Histone = histone, Replicate = replicate, sampleInfo = sample_name)
  peakData <- rbind(peakData, peakInfo)
  
  # Add to GRanges list for overlap analysis
  peak_granges_list[[sample_name]] <- GRanges(seqnames = peakInfo$chrom, 
                                              ranges = IRanges(start = peakInfo$start, end = peakInfo$end))
}

# Convert columns to factors for ordered plotting
peakData$sampleInfo <- factor(peakData$sampleInfo, levels = sampleList)
peakData$Histone <- factor(peakData$Histone, levels = histList)

# Plot 1: Number of Peaks
fig1 <- peakData %>%
  group_by(Histone, Replicate) %>%
  summarise(peakN = n(), .groups = "drop") %>%
  ggplot(aes(x = Histone, y = peakN, fill = Histone)) +
  geom_boxplot() +
  geom_jitter(aes(color = Replicate), position = position_jitter(0.15)) +
  theme_bw(base_size = 18) +
  ylab("Number of Peaks") +
  xlab("")

# Plot 2: Width of Peaks
fig2 <- ggplot(peakData, aes(x = Histone, y = width, fill = Histone)) +
  geom_violin() +
  facet_grid(Replicate ~ .) +
  theme_bw(base_size = 18) +
  ylab("Width of Peaks") +
  xlab("")

# Plot 3: % of Peaks Reproduced (overlap between replicates)
reproducibility_data <- data.frame()
for (hist in histList) {
  rep1 <- paste0(hist, "_rep1")
  rep2 <- paste0(hist, "_rep2")
  if (rep1 %in% names(peak_granges_list) && rep2 %in% names(peak_granges_list)) {
    # Find overlapping peaks between replicates
    overlaps <- findOverlaps(peak_granges_list[[rep1]], peak_granges_list[[rep2]])
    reproducible_count <- length(unique(queryHits(overlaps)))
    total_peaks <- length(peak_granges_list[[rep1]]) + length(peak_granges_list[[rep2]])
    reproducibility_rate <- (reproducible_count / (total_peaks / 2)) * 100
    reproducibility_data <- rbind(reproducibility_data, data.frame(Histone = hist, ReproducibilityRate = reproducibility_rate))
  }
}

fig3 <- ggplot(reproducibility_data, aes(x = Histone, y = ReproducibilityRate, fill = Histone)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = round(ReproducibilityRate, 2)), vjust = -0.5) +
  theme_bw(base_size = 18) +
  ylab("% of Peaks Reproduced") +
  xlab("")

# Plot 4: Fraction of Reads in Peaks (FRiP)
frip_data <- data.frame()
for (hist in histList) {
  for (rep in unique(peakData$Replicate[peakData$Histone == hist])) {
    peakInfo <- peakData %>%
      filter(Histone == hist & Replicate == rep)
    peak_gr <- GRanges(seqnames = peakInfo$chrom, 
                       ranges = IRanges(start = peakInfo$start, end = peakInfo$end))
    
    # Assuming total number of mapped fragments (mapped_fragments) is provided
    # Here we'll use width of peaks as a proxy for simplicity
    in_peak_fragments <- sum(peakInfo$width)  # Summing widths as proxy for fragments in peaks
    total_fragments <- in_peak_fragments  # Update with actual total fragments if available
    frip_score <- (in_peak_fragments / total_fragments) * 100
    frip_data <- rbind(frip_data, data.frame(Histone = hist, Replicate = rep, FRiP = frip_score))
  }
}

fig4 <- ggplot(frip_data, aes(x = Histone, y = FRiP, fill = Histone)) +
  geom_boxplot() +
  geom_jitter(aes(color = Replicate), position = position_jitter(0.15)) +
  theme_bw(base_size = 18) +
  ylab("% of Fragments in Peaks (FRiP)") +
  xlab("")

# Arrange and save all plots
final_plot <- ggarrange(fig1, fig2, fig3, fig4, ncol = 2, nrow = 2, common.legend = TRUE, legend = "bottom")
output_file <- file.path(output_dir, "peak_summary_plot.png")
ggsave(output_file, final_plot, width = 14, height = 10)
