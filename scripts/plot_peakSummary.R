#!/usr/bin/env Rscript

# Load required libraries
library(dplyr)
library(ggplot2)
library(viridis)
library(ggpubr)
library(GenomicRanges)

# Capture command-line arguments
args <- commandArgs(trailingOnly = TRUE)
peak_files <- strsplit(args[1], ",")[[1]]  # First argument: comma-separated list of narrowPeak files
frip_files <- strsplit(args[2], ",")[[1]]  # Second argument: comma-separated list of FRiP files
output_dir <- args[3]                      # Third argument: output directory

# Ensure the output directory exists
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# --- Process Peak Files ---
sampleList <- gsub("_0.05_peaks.narrowPeak", "", basename(peak_files))
histList <- unique(sapply(strsplit(sampleList, "_"), `[`, 1))

peakData <- data.frame()
peak_granges_list <- list()

for (file_path in peak_files) {
  peakInfo <- read.table(file_path, header = FALSE, 
                         col.names = c("chrom", "start", "end", "name", "score", 
                                       "strand", "signalValue", "pValue", "qValue", "peak"))
  
  sample_name <- gsub("_0.05_peaks.narrowPeak", "", basename(file_path))
  histInfo <- strsplit(sample_name, "_")[[1]]
  histone <- histInfo[1]
  replicate <- histInfo[2]
  
  peakInfo <- peakInfo %>%
    mutate(width = end - start, Histone = histone, Replicate = replicate, sampleInfo = sample_name)
  peakData <- rbind(peakData, peakInfo)
  
  peak_granges_list[[sample_name]] <- GRanges(seqnames = peakInfo$chrom, 
                                              ranges = IRanges(start = peakInfo$start, end = peakInfo$end))
}

peakData$sampleInfo <- factor(peakData$sampleInfo, levels = sampleList)
peakData$Histone <- factor(peakData$Histone, levels = histList)

# --- Plot 1: Number of Peaks ---
fig1 <- peakData %>%
  group_by(Histone, Replicate) %>%
  summarise(peakN = n(), .groups = "drop") %>%
  ggplot(aes(x = Histone, y = peakN, fill = Histone)) +
  geom_boxplot() +
  geom_jitter(aes(color = Replicate), position = position_jitter(0.15)) +
  theme_bw(base_size = 18) +
  ylab("Number of Peaks") +
  xlab("")

# --- Plot 2: Width of Peaks ---
fig2 <- ggplot(peakData, aes(x = Histone, y = width, fill = Histone)) +
  geom_violin() +
  facet_grid(Replicate ~ .) +
  theme_bw(base_size = 18) +
  ylab("Width of Peaks") +
  xlab("")

# --- Plot 3: % of Peaks Reproduced ---
reproducibility_data <- data.frame()
for (hist in histList) {
  reps <- unique(peakData$Replicate[peakData$Histone == hist])
  if (length(reps) >= 2) {
    rep_combinations <- combn(reps, 2)
    for (k in 1:ncol(rep_combinations)) {
      rep1_id <- rep_combinations[1, k]
      rep2_id <- rep_combinations[2, k]
      rep1 <- paste0(hist, "_", rep1_id)
      rep2 <- paste0(hist, "_", rep2_id)
      
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

fig3 <- if (nrow(reproducibility_data) > 0) {
  ggplot(reproducibility_data, aes(x = Histone, y = ReproducibilityRate, fill = Histone)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    geom_text(aes(label = round(ReproducibilityRate, 2), group = ReplicatePair), 
              vjust = -0.5, position = position_dodge(0.9)) +
    theme_bw(base_size = 18) +
    ylab("% of Peaks Reproduced") +
    xlab("") +
    facet_wrap(~ ReplicatePair)
} else {
  ggplot() +
    geom_blank() +
    theme_void() +
    ggtitle("No Reproducibility Data Available")
}

# --- Process FRiP Files ---
frip_data <- data.frame()

for (file_path in frip_files) {
  temp <- tryCatch({
    read.table(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  }, error = function(e) {
    message(paste("Error reading FRiP file:", file_path))
    return(NULL)
  })
  
  if (!is.null(temp)) {
    frip_data <- rbind(frip_data, temp)
  }
}

if ("Sample" %in% colnames(frip_data) && nrow(frip_data) > 0) {
  fig4 <- ggplot(frip_data, aes(x = Sample, y = FRiP, fill = Sample)) +
    geom_boxplot() +
    geom_jitter(position = position_jitter(0.15)) +
    theme_bw(base_size = 18) +
    ylab("% of Fragments in Peaks (FRiP)") +
    xlab("")
} else {
  fig4 <- ggplot() +
    geom_blank() +
    theme_void() +
    ggtitle("No FRiP Data Available")
}

# --- Arrange and Save Plots ---
final_plot <- ggarrange(fig1, fig2, fig3, fig4, ncol = 2, nrow = 2, common.legend = TRUE, legend = "bottom")
output_file <- file.path(output_dir, "peak_summary_plot.png")
ggsave(output_file, final_plot, width = 14, height = 10)
