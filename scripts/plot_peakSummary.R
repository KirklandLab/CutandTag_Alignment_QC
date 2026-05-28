#!/usr/bin/env Rscript

# Load required libraries
library(dplyr)
library(ggplot2)
library(viridis)
library(ggpubr)
library(GenomicRanges)

# -------------------------
# Capture command-line arguments
# -------------------------
args <- commandArgs(trailingOnly = TRUE)
peak_files <- strsplit(args[1], " ")[[1]]
frip_files <- strsplit(args[2], " ")[[1]]
metadata_file <- args[3]
output_dir <- args[4]

# -------------------------
# Ensure output directory exists
# -------------------------
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# -------------------------
# Load sample metadata
# -------------------------
sample_metadata <- read.csv(metadata_file, stringsAsFactors = FALSE)

# -------------------------
# Process peak files
# -------------------------
sampleList <- gsub("_0.05_peaks.narrowPeak", "", basename(peak_files))

peakData <- data.frame()
peak_granges_list <- list()

for (file_path in peak_files) {
  peakInfo <- read.table(file_path, header = FALSE,
                         col.names = c("chrom", "start", "end", "name", "score",
                                       "strand", "signalValue", "pValue", "qValue", "peak"))
  
  sample_name <- gsub("_0.05_peaks.narrowPeak", "", basename(file_path))
  histone <- sample_metadata$histone[sample_metadata$sample == sample_name]
  replicate <- sample_metadata$replicate[sample_metadata$sample == sample_name]
  
  peakInfo <- peakInfo %>%
    mutate(width = end - start,
           Histone = histone,
           Replicate = replicate,
           sampleInfo = sample_name)
  
  peakData <- rbind(peakData, peakInfo)
  
  peak_granges_list[[sample_name]] <- GRanges(seqnames = peakInfo$chrom,
                                              ranges = IRanges(start = peakInfo$start, end = peakInfo$end))
}

peakData$sampleInfo <- factor(peakData$sampleInfo, levels = sampleList)
peakData$Histone <- factor(peakData$Histone, levels = unique(peakData$Histone))
peakData$Replicate <- as.factor(peakData$Replicate)

# -------------------------
# Plot 1: Number of Peaks
# -------------------------
fig1 <- peakData %>%
  group_by(sampleInfo, Histone, Replicate) %>%
  summarise(peakN = n(), .groups = "drop") %>%
  ggplot(aes(x = Histone, y = peakN, fill = Histone)) +
  geom_boxplot() +
  geom_jitter(aes(color = Replicate), position = position_jitter(0.15)) +
  theme_bw(base_size = 18) +
  ylab("Number of Peaks") +
  xlab("") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

# -------------------------
# Plot 2: Width of Peaks
# -------------------------
fig2 <- ggplot(peakData, aes(x = Histone, y = width, fill = Histone)) +
  geom_violin() +
  facet_grid(Replicate ~ .) +
  theme_bw(base_size = 18) +
  ylab("Width of Peaks") +
  xlab("") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

# -------------------------
# Plot 3: Reproducibility
# -------------------------
# Compare all sample pairs within each histone group.

reproducibility_data <- data.frame()

for (hist in unique(as.character(peakData$Histone))) {
  
  # Get all samples for this histone group from the metadata.
  hist_samples <- sample_metadata$sample[sample_metadata$histone == hist]
  
  # Keep only samples that actually have peak GRanges loaded.
  hist_samples <- hist_samples[hist_samples %in% names(peak_granges_list)]
  
  # Need at least two samples to calculate pairwise reproducibility.
  if (length(hist_samples) >= 2) {
    
    sample_combinations <- combn(hist_samples, 2, simplify = FALSE)
    
    for (pair in sample_combinations) {
      sample1 <- pair[1]
      sample2 <- pair[2]
      
      peaks1 <- peak_granges_list[[sample1]]
      peaks2 <- peak_granges_list[[sample2]]
      
      # Skip if one sample has no peaks.
      if (length(peaks1) == 0 || length(peaks2) == 0) {
        next
      }
      
      overlaps <- findOverlaps(peaks1, peaks2)
      
      reproducible_count <- length(unique(queryHits(overlaps)))
      total_peaks <- (length(peaks1) + length(peaks2)) / 2
      reproducibility_rate <- (reproducible_count / total_peaks) * 100
      
      rep1_id <- sample_metadata$replicate[sample_metadata$sample == sample1]
      rep2_id <- sample_metadata$replicate[sample_metadata$sample == sample2]
      
      reproducibility_data <- rbind(
        reproducibility_data,
        data.frame(
          Histone = hist,
          SamplePair = paste(sample1, sample2, sep = " vs "),
          ReplicatePair = paste(rep1_id, rep2_id, sep = "-"),
          ReproducibilityRate = reproducibility_rate,
          stringsAsFactors = FALSE
        )
      )
    }
  }
}

if (nrow(reproducibility_data) > 0) {
  fig3 <- ggplot(
    reproducibility_data,
    aes(x = Histone, y = ReproducibilityRate, fill = Histone)
  ) +
    geom_col(position = position_dodge()) +
    geom_text(
      aes(label = round(ReproducibilityRate, 1)),
      vjust = -0.5,
      size = 3
    ) +
    facet_wrap(~ SamplePair) +
    theme_bw(base_size = 14) +
    ylab("% of Peaks Reproduced") +
    xlab("") +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      strip.text = element_text(size = 8)
    )
} else {
  fig3 <- ggplot() +
    geom_blank() +
    theme_void() +
    ggtitle("No Reproducibility Data Available")
}

# -------------------------
# Plot 4: FRiP Scores
# -------------------------
frip_data <- do.call(rbind, lapply(frip_files, function(file_path) {
  read.table(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
}))

if (!"Sample" %in% colnames(frip_data)) {
  stop("Error: 'Sample' column missing in FRiP data.")
}

# Attach histone and replicate info to FRiP data
frip_data$Histone <- sample_metadata$histone[match(frip_data$Sample, sample_metadata$sample)]
frip_data$Replicate <- as.factor(sample_metadata$replicate[match(frip_data$Sample, sample_metadata$sample)])

fig4 <- ggplot(frip_data, aes(x = Histone, y = FRiP, fill = Histone)) +
  geom_boxplot() +
  geom_jitter(aes(color = Replicate), position = position_jitter(0.15)) +
  theme_bw(base_size = 18) +
  ylab("% of Fragments in Peaks (FRiP)") +
  xlab("") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

# -------------------------
# Arrange and Save
# -------------------------
final_plot <- ggarrange(fig1, fig2, fig3, fig4,
                        ncol = 2, nrow = 2,
                        common.legend = TRUE, legend = "bottom")

output_file <- file.path(output_dir, "peak_summary_plot.png")
ggsave(output_file, final_plot, width = 14, height = 10)
