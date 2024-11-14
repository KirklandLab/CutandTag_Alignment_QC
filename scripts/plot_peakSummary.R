#!/usr/bin/env Rscript

# Load required libraries
library(dplyr)
library(ggplot2)
library(viridis)
library(ggpubr)
library(GenomicRanges)
library(chromVAR)

# Capture command-line arguments
args <- commandArgs(trailingOnly = TRUE)
input_files <- args[-c(length(args) - 1, length(args))]  # All but last two arguments are input files
output_dir <- args[length(args) - 1]                    # Second-to-last argument is the output directory
bam_dir <- args[length(args)]                           # Last argument is the BAM directory

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
  # Get all replicates for this histone mark
  reps <- unique(peakData$Replicate[peakData$Histone == hist])
  
  # Only proceed if there are at least two replicates
  if (length(reps) >= 2) {
    # Generate all pairwise combinations of replicates
    rep_combinations <- combn(reps, 2)
    
    for (k in 1:ncol(rep_combinations)) {
      rep1_id <- rep_combinations[1, k]
      rep2_id <- rep_combinations[2, k]
      rep1 <- paste0(hist, "_", rep1_id)
      rep2 <- paste0(hist, "_", rep2_id)
      
      if (rep1 %in% names(peak_granges_list) && rep2 %in% names(peak_granges_list)) {
        # Find overlapping peaks between replicates
        overlaps <- findOverlaps(peak_granges_list[[rep1]], peak_granges_list[[rep2]])
        reproducible_count <- length(unique(queryHits(overlaps)))
        total_peaks <- (length(peak_granges_list[[rep1]]) + length(peak_granges_list[[rep2]])) / 2
        reproducibility_rate <- (reproducible_count / total_peaks) * 100
        
        # Store the reproducibility rate
        reproducibility_data <- rbind(reproducibility_data, 
                                      data.frame(Histone = hist, 
                                                 ReplicatePair = paste(rep1_id, rep2_id, sep = "-"),
                                                 ReproducibilityRate = reproducibility_rate))
      }
    }
  }
}

# Adjust your plotting code to accommodate the new data structure
fig3 <- ggplot(reproducibility_data, aes(x = Histone, y = ReproducibilityRate, fill = Histone)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_text(aes(label = round(ReproducibilityRate, 2), group = ReplicatePair), 
            vjust = -0.5, position = position_dodge(0.9)) +
  theme_bw(base_size = 18) +
  ylab("% of Peaks Reproduced") +
  xlab("") +
  facet_wrap(~ ReplicatePair)


# Plot 4: Fraction of Reads in Peaks (FRiP) using chromVAR
inPeakData <- data.frame()

for (hist in histList) {
  for (rep in unique(peakData$Replicate[peakData$Histone == hist])) {
    # Define the peak regions for the current histone and replicate
    peakInfo <- peakData %>%
      filter(Histone == hist & Replicate == rep)
    peak_gr <- GRanges(seqnames = peakInfo$chrom, 
                       ranges = IRanges(start = peakInfo$start, end = peakInfo$end),
                       strand = "*")
    
    # Path to the BAM file for the current histone and replicate
    bamFile <- file.path(bam_dir, paste0(hist, "_", rep, "_bowtie2.mapped.bam"))
    print(paste("Reading in file:", bamFile))
    
    # Use chromVAR to get fragment counts within peaks and total fragments
    fragment_counts <- getCounts(bamFile, peak_gr, paired = TRUE, by_rg = FALSE, format = "bam")
    
    # Check if fragment_counts has valid data
    if (length(fragment_counts$total) > 0) {
      # Total fragments (MappedFragNum) and fragments in peaks (inPeakN)
      inPeakN <- sum(counts(fragment_counts)[,1])  # Fragments in peaks
      MappedFragNum <- fragment_counts$total  # Total mapped fragments in BAM
      
      # Store results in a data frame
      inPeakData <- rbind(inPeakData, data.frame(Histone = hist, Replicate = rep, inPeakN = inPeakN, MappedFragNum = MappedFragNum))
    } else {
      print(paste("Warning: No data found in fragment counts for", hist, rep))
    }
  }
}

# Calculate FRiP scores
frip <- inPeakData %>% 
  mutate(FRiP = (inPeakN / MappedFragNum) * 100)

fig4 <- ggplot(frip, aes(x = Histone, y = FRiP, fill = Histone)) +
  geom_boxplot() +
  geom_jitter(aes(color = Replicate), position = position_jitter(0.15)) +
  theme_bw(base_size = 18) +
  ylab("% of Fragments in Peaks (FRiP)") +
  xlab("")

# Arrange and save all plots
final_plot <- ggarrange(fig1, fig2, fig3, fig4, ncol = 2, nrow = 2, common.legend = TRUE, legend = "bottom")
output_file <- file.path(output_dir, "peak_summary_plot.png")
ggsave(output_file, final_plot, width = 14, height = 10)
