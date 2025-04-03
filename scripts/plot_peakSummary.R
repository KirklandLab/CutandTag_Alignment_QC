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
  group_by(Histone, Replicate) %>%
  summarise(peakN = n(), .groups = "drop") %>%
  ggplot(aes(x = Histone, y = peakN, fill = Histone)) +
  geom_boxplot() +
  geom_jitter(aes(color = Replicate), position = position_jitter(0.15)) +
  theme_bw(base_size = 18) +
  ylab("Number of Peaks") +
  xlab("")

# -------------------------
# Plot 2: Width of Peaks
# -------------------------
fig2 <- ggplot(peakData, aes(x = Histone, y = width, fill = Histone)) +
  geom_violin() +
  facet_grid(Replicate ~ .) +
  theme_bw(base_size = 18) +
  ylab("Width of Peaks") +
  xlab("")

# -------------------------
# Plot 3: Reproducibility
# -------------------------
reproducibility_data <- data.frame()

for (hist in unique(peakData$Histone)) {
  reps <- unique(peakData$Replicate[peakData$Histone == hist])
  if (length(reps) >= 2) {
    rep_combinations <- combn(reps, 2)
    for (k in 1:ncol(rep_combinations)) {
      rep1_id <- rep_combinations[1, k]
      rep2_id <- rep_combinations[2, k]
      rep1 <- sample_metadata$sample[sample_metadata$histone == hist & sample_metadata$replicate == rep1_id]
      rep2 <- sample_metadata$sample[sample_metadata$histone == hist & sample_metadata$replicate == rep2_id]
      
      if (rep1 %in% names(peak_granges_list) && rep2 %in% names(peak_granges_list)) {
        overlaps <- findOverlaps(peak_granges_list[[rep1]], peak_granges_list[[rep2]])
        reproducible_count <- length(unique(queryHits(overlaps)))
        total_peaks <- (length(peak_granges_list[[rep1]]) + length(peak_granges_list[[rep2]])) / 2
        reproducibility_rate <- (reproducible_count / total_peaks) * 100
        
        reproducibility_data <- rbind(reproducibility_data,
                                      data.frame(Histone = hist,
                                                 ReplicatePair = paste(rep1_id, rep2_id, sep = "-"),
                                                 ReproducibilityRate = reproducibility_rate))
      }
    }
  }
}

if (nrow(reproducibility_data) > 0) {
  fig3 <- ggplot(reproducibility_data, aes(x = Histone, y = ReproducibilityRate, fill = Histone)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    geom_text(aes(label = round(ReproducibilityRate, 2), group = ReplicatePair),
              vjust = -0.5, position = position_dodge(0.9)) +
    theme_bw(base_size = 18) +
    ylab("% of Peaks Reproduced") +
    xlab("") +
    facet_wrap(~ ReplicatePair)
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
