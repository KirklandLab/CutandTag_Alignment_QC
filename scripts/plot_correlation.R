#!/usr/bin/env Rscript

# Load required libraries
library(dplyr)
library(corrplot)
library(Cairo)  # Add Cairo library

# Capture command-line arguments
args <- commandArgs(trailingOnly = TRUE)
input_files <- args[-length(args)]  # All but last argument are input files
output_dir <- args[length(args)]    # Last argument is the output directory

# Ensure the output directory exists
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Initialize variables
fragCount <- NULL

# Loop through each input file to collect fragment count data
for (file_path in input_files) {
  
  # Extract sample name from file name
  sample_name <- gsub("_bowtie2.fragmentsCount.bin500.bed", "", basename(file_path))
  
  # Read the fragment count data
  fragCountTmp <- read.table(file_path, header = FALSE)
  colnames(fragCountTmp) <- c("chrom", "bin", sample_name)
  
  # Combine with existing data
  if (is.null(fragCount)) {
    fragCount <- fragCountTmp
  } else {
    fragCount <- full_join(fragCount, fragCountTmp, by = c("chrom", "bin"))
  }
}

# Calculate correlation matrix
M <- cor(fragCount %>% select(-chrom, -bin) %>% log2(), use = "complete.obs")

# Generate correlation plot
output_file <- file.path(output_dir, "fragCount_correlation_plot.png")
Cairo::CairoPNG(filename = output_file, width = 800, height = 800)
corrplot(M, method = "color", outline = TRUE, addgrid.col = "darkgray", order = "hclust",
         addrect = 3, rect.col = "black", rect.lwd = 3, cl.pos = "b", tl.col = "indianred4",
         tl.cex = 1, cl.cex = 1, addCoef.col = "black", number.digits = 2, number.cex = 1,
         col = colorRampPalette(c("midnightblue", "white", "darkred"))(100))
dev.off()
