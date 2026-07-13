#!/usr/bin/env Rscript

# Load required libraries
library(ggplot2)
library(dplyr)
library(ggpubr)

# -------------------------
# Capture command-line arguments
# -------------------------
args <- commandArgs(trailingOnly = TRUE)
input_files <- args[1:(length(args) - 2)]         # All but last 2 args are the summary files
metadata_file <- args[length(args) - 1]           # Second to last = samples.csv
output_dir <- args[length(args)]                  # Last arg = output directory

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

# Initialize an empty data frame to store alignment results
alignResult <- data.frame()

# -------------------------
# Loop through each Bowtie2 summary file
# -------------------------
for (file_path in input_files) {
  sample_name <- gsub("_bowtie2.txt", "", basename(file_path))

  histone <- sample_metadata$histone[sample_metadata$sample == sample_name]
  replicate <- sample_metadata$replicate[sample_metadata$sample == sample_name]

  if (length(histone) != 1 || length(replicate) != 1) {
    stop("Sample not found uniquely in metadata: ", sample_name)
  }

  lines <- readLines(file_path, warn = FALSE)

  get_first_number <- function(pattern) {
    line <- grep(pattern, lines, value = TRUE)
    if (length(line) == 0) return(NA_real_)
    as.numeric(sub("^\\s*([0-9,]+).*", "\\1", gsub(",", "", line[1])))
  }

  get_percent <- function(pattern) {
  line <- grep(pattern, lines, value = TRUE)
  if (length(line) == 0) return(NA_real_)
  match <- regmatches(line[1], regexpr("[0-9.]+(?=%)", line[1], perl = TRUE))
  if (length(match) == 0) return(NA_real_)
  as.numeric(match)
  }

  sequencingDepth <- get_first_number("reads; of these:")
  pairedReads <- get_first_number("were paired; of these:")
  concordant_zero <- get_first_number("aligned concordantly 0 times")
  concordant_one <- get_first_number("aligned concordantly exactly 1 time")
  concordant_multi <- get_first_number("aligned concordantly >1 times")
  alignmentRate <- get_percent("overall alignment rate")

  mappedFragNum <- concordant_one + concordant_multi

  uniqueMappedRate <- ifelse(
    pairedReads > 0,
    100 * concordant_one / pairedReads,
    NA_real_
  )

  multiMappedRate <- ifelse(
    pairedReads > 0,
    100 * concordant_multi / pairedReads,
    NA_real_
  )

  unalignedRate <- ifelse(
    pairedReads > 0,
    100 * concordant_zero / pairedReads,
    NA_real_
  )

  alignResult <- rbind(
    alignResult,
    data.frame(
      Sample = sample_name,
      Histone = histone,
      Replicate = replicate,
      SequencingDepth = sequencingDepth,
      PairedReads = pairedReads,
      ConcordantZero = concordant_zero,
      ConcordantOne = concordant_one,
      ConcordantMulti = concordant_multi,
      MappedFragNum = mappedFragNum,
      AlignmentRate = alignmentRate,
      UniqueMappedRate = uniqueMappedRate,
      MultiMappedRate = multiMappedRate,
      UnalignedRate = unalignedRate
    )
  )
}

# -------------------------
# Ensure consistent order
# -------------------------
alignResult$Histone <- factor(alignResult$Histone, levels = unique(alignResult$Histone))
alignResult$Replicate <- as.factor(alignResult$Replicate)

# -------------------------
# Create plots
# -------------------------

fig1 <- ggplot(alignResult, aes(x = Histone, y = SequencingDepth / 1e6, fill = Histone)) +
  geom_boxplot() +
  geom_jitter(aes(color = Replicate), width = 0.2, height = 0) +
  theme_bw(base_size = 14) +
  ylab("Sequencing Depth (Millions)") +
  xlab("") +
  ggtitle("Sequencing Depth per Histone") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

fig2 <- ggplot(alignResult, aes(x = Histone, y = MappedFragNum / 1e6, fill = Histone)) +
  geom_boxplot() +
  geom_jitter(aes(color = Replicate), width = 0.2) +
  theme_bw(base_size = 14) +
  ylab("Mapped Fragments (Millions)") +
  xlab("") +
  ggtitle("Mapped Fragments per Histone") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

fig3 <- ggplot(alignResult, aes(x = Histone, y = AlignmentRate, fill = Histone)) +
  geom_boxplot() +
  geom_jitter(aes(color = Replicate), width = 0.2) +
  theme_bw(base_size = 14) +
  ylab("Alignment Rate (%)") +
  xlab("") +
  ggtitle("Alignment Rate per Histone") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

fig4 <- ggplot(alignResult, aes(x = Histone, y = UniqueMappedRate, fill = Histone)) +
  geom_boxplot() +
  geom_jitter(aes(color = Replicate), width = 0.2) +
  theme_bw(base_size = 14) +
  ylab("Unique Concordant Alignment Rate (%)") +
  xlab("") +
  ggtitle("Uniquely Mapped Concordant Pairs per Histone") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

# -------------------------
# Arrange and save
# -------------------------

final_plot <- ggarrange(fig1, fig2, fig3, fig4, 
                        ncol = 2, nrow = 2, 
                        common.legend = TRUE, legend = "bottom")

output_file <- file.path(output_dir, "alignment_summary_plot.png")
ggsave(output_file, final_plot, width = 14, height = 10)
