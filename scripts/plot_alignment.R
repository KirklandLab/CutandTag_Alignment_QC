#!/usr/bin/env Rscript

# Load required libraries
library(ggplot2)
library(dplyr)
library(ggpubr)

# -------------------------
# Capture command-line arguments
# -------------------------
args <- commandArgs(trailingOnly = TRUE)
input_files <- args[1:(length(args) - 2)]
metadata_file <- args[length(args) - 1]
output_dir <- args[length(args)]

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
  
  # Lookup metadata
  histone <- sample_metadata$histone[sample_metadata$sample == sample_name]
  replicate <- sample_metadata$replicate[sample_metadata$sample == sample_name]
  
  # Read summary file
  alignRes <- read.table(file_path, header = FALSE, fill = TRUE, stringsAsFactors = FALSE)
  
  # Use keyword matching to extract data safely
  total_reads_line <- alignRes$V1[grepl("reads; of these:", alignRes$V1)]
  sequencingDepth <- as.numeric(gsub(" .*", "", total_reads_line))
  
  mapped1_line <- alignRes$V1[grepl("aligned concordantly exactly 1 time", alignRes$V1)]
  mapped1 <- as.numeric(gsub("[^0-9]", "", mapped1_line))
  
  mapped2_line <- alignRes$V1[grepl("aligned concordantly >1 times", alignRes$V1)]
  mapped2 <- as.numeric(gsub("[^0-9]", "", mapped2_line))
  
  alignmentRate_line <- alignRes$V1[grepl("overall alignment rate", alignRes$V1)]
  alignmentRate <- as.numeric(gsub("[^0-9.]", "", alignmentRate_line))
  
  # Handle missing values
  if (any(is.na(c(sequencingDepth, mapped1, mapped2, alignmentRate)))) {
    warning(paste("Skipping sample", sample_name, "due to missing values."))
    next
  }
  
  mappedFragNum <- mapped1 + mapped2
  unalignedRate <- 100 - alignmentRate
  
  # Append to results
  alignResult <- rbind(
    alignResult,
    data.frame(
      Sample = sample_name,
      Histone = histone,
      Replicate = replicate,
      SequencingDepth = sequencingDepth,
      MappedFragNum = mappedFragNum,
      AlignmentRate = alignmentRate,
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
  geom_jitter(aes(color = Replicate), width = 0.2) +
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

fig4 <- ggplot(alignResult, aes(x = Histone, y = UnalignedRate, fill = Histone)) +
  geom_boxplot() +
  geom_jitter(aes(color = Replicate), width = 0.2) +
  theme_bw(base_size = 14) +
  ylab("Unaligned Rate (%)") +
  xlab("") +
  ggtitle("Unaligned Rate per Histone") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

# -------------------------
# Arrange and save
# -------------------------

final_plot <- ggarrange(fig1, fig2, fig3, fig4,
                        ncol = 2, nrow = 2,
                        common.legend = TRUE, legend = "bottom")

output_file <- file.path(output_dir, "alignment_summary_plot.png")
ggsave(output_file, final_plot, width = 14, height = 10)
