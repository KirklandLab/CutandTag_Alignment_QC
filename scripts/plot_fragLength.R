#!/usr/bin/env Rscript

# ------------------------------------------------------------
# Fragment length plotting script
# ------------------------------------------------------------
# Input:
#   One or more fragment length files:
#     {sample}_fragmentLen.raw.txt
#     {sample}_fragmentLen.final.txt
#
#   samples.csv
#   output directory
#
# Output:
#   fragment_length_plot.png
# ------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(viridis)
  library(ggpubr)
})

# -------------------------
# Capture command-line arguments
# -------------------------
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop(
    "Usage: Rscript plot_fragLength.R <fragLen files...> <samples.csv> <output_dir>\n",
    "Example: Rscript plot_fragLength.R sample1_fragmentLen.final.txt sample2_fragmentLen.final.txt config/samples.csv results/plots"
  )
}

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

required_cols <- c("sample", "histone", "replicate")
missing_cols <- setdiff(required_cols, colnames(sample_metadata))

if (length(missing_cols) > 0) {
  stop(
    "samples.csv is missing required columns: ",
    paste(missing_cols, collapse = ", ")
  )
}

# -------------------------
# Helper function to extract sample name
# -------------------------
extract_sample_name <- function(file_path) {
  sample_name <- basename(file_path)

  sample_name <- gsub(
    "_fragmentLen\\.(raw|final)\\.txt$",
    "",
    sample_name
  )

  return(sample_name)
}

sampleList <- vapply(input_files, extract_sample_name, character(1))

# -------------------------
# Read and combine fragment length files
# -------------------------
fragLen <- data.frame()

for (file_path in input_files) {
  sample_name <- extract_sample_name(file_path)

  if (!file.exists(file_path)) {
    stop("Fragment length file does not exist: ", file_path)
  }

  if (!(sample_name %in% sample_metadata$sample)) {
    stop(
      "Sample name extracted from fragment file does not match samples.csv: ",
      sample_name,
      "\nFile: ",
      file_path
    )
  }

  histone <- sample_metadata$histone[sample_metadata$sample == sample_name]
  replicate <- sample_metadata$replicate[sample_metadata$sample == sample_name]

  fragData <- read.table(
    file_path,
    header = FALSE,
    col.names = c("fragLen", "fragCount")
  )

  if (nrow(fragData) == 0) {
    warning("Fragment length file is empty: ", file_path)
    next
  }

  fragData <- fragData %>%
    mutate(
      fragLen = as.numeric(fragLen),
      fragCount = as.numeric(fragCount),
      Weight = fragCount / sum(fragCount),
      Histone = histone,
      Replicate = replicate,
      sampleInfo = sample_name
    )

  fragLen <- bind_rows(fragLen, fragData)
}

if (nrow(fragLen) == 0) {
  stop("No fragment length data was loaded.")
}

# -------------------------
# Factor ordering
# -------------------------
fragLen$sampleInfo <- factor(fragLen$sampleInfo, levels = sampleList)
fragLen$Histone <- factor(fragLen$Histone, levels = unique(fragLen$Histone))
fragLen$Replicate <- as.factor(fragLen$Replicate)

# -------------------------
# Violin plot: weighted fragment size distribution
# -------------------------
fig5A <- ggplot(
  fragLen,
  aes(
    x = sampleInfo,
    y = fragLen,
    weight = Weight,
    fill = Histone
  )
) +
  geom_violin(bw = 5, scale = "width", trim = TRUE) +
  scale_y_continuous(breaks = seq(0, 800, 50)) +
  scale_fill_viridis(
    discrete = TRUE,
    begin = 0.1,
    end = 0.9,
    option = "magma",
    alpha = 0.8
  ) +
  coord_cartesian(ylim = c(0, 800)) +
  theme_bw(base_size = 20) +
  ggpubr::rotate_x_text(angle = 20) +
  ylab("Fragment Length") +
  xlab("") +
  ggtitle("Fragment length distribution")

# -------------------------
# Line plot: fragment counts by length
# -------------------------
fig5B <- ggplot(
  fragLen,
  aes(
    x = fragLen,
    y = fragCount,
    color = Histone,
    group = sampleInfo,
    linetype = Replicate
  )
) +
  geom_line(linewidth = 1) +
  scale_color_viridis(
    discrete = TRUE,
    begin = 0.1,
    end = 0.9,
    option = "magma"
  ) +
  theme_bw(base_size = 20) +
  xlab("Fragment Length") +
  ylab("Count") +
  coord_cartesian(xlim = c(0, 500)) +
  ggtitle("Fragment length counts")

# -------------------------
# Combine and save
# -------------------------
final_plot <- ggarrange(
  fig5A,
  fig5B,
  ncol = 2,
  labels = c("A", "B")
)

output_file <- file.path(output_dir, "fragment_length_plot.png")

ggsave(
  output_file,
  final_plot,
  width = 14,
  height = 10,
  dpi = 300
)
