#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(viridis)
  library(ggpubr)
  library(GenomicRanges)
})

# ============================================================
# Adjustable plotting settings
# ============================================================

# Number of detailed sample-pair comparisons shown per PDF page.
pairs_per_page <- 6

# Number of columns used for detailed pairwise comparison pages.
detail_page_columns <- 2

# Maximum number of pairs shown in the compact four-panel overview.
# All pairs are still included in the detailed PDF pages.
overview_pair_limit <- 10

# ============================================================
# Capture and validate command-line arguments
# ============================================================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 4) {
  stop(
    paste0(
      "Usage: Rscript plot_peakSummary.R ",
      "\"<peak files>\" \"<FRiP files>\" ",
      "<samples.csv> <output_dir>"
    )
  )
}

peak_files <- strsplit(args[1], "\\s+")[[1]]
frip_files <- strsplit(args[2], "\\s+")[[1]]
metadata_file <- args[3]
output_dir <- args[4]

peak_files <- peak_files[peak_files != ""]
frip_files <- frip_files[frip_files != ""]

if (length(peak_files) == 0) {
  stop("No peak files were supplied.")
}

if (length(frip_files) == 0) {
  stop("No FRiP files were supplied.")
}

if (!file.exists(metadata_file)) {
  stop("Metadata file does not exist: ", metadata_file)
}

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# ============================================================
# Load and validate sample metadata
# ============================================================

sample_metadata <- read.csv(metadata_file, stringsAsFactors = FALSE)

required_metadata_columns <- c("sample","histone","replicate")

missing_metadata_columns <- setdiff(
  required_metadata_columns,
  colnames(sample_metadata)
)

if (length(missing_metadata_columns) > 0) {
  stop(
    "Metadata file is missing required columns: ",
    paste(missing_metadata_columns, collapse = ", ")
  )
}

if (anyDuplicated(sample_metadata$sample)) {
  duplicated_samples <- unique(
    sample_metadata$sample[
      duplicated(sample_metadata$sample)
    ]
  )

  stop(
    "Sample names must be unique in samples.csv. Duplicated samples: ",
    paste(duplicated_samples, collapse = ", ")
  )
}

histone_levels <- unique(sample_metadata$histone)

# ============================================================
# Helper functions
# ============================================================

extract_peak_sample_name <- function(file_path) {
  file_name <- basename(file_path)
  sample_name <- sub(
    "_[^_]+_peaks\\.narrowPeak$",
    "",
    file_name
  )
  if (identical(sample_name, file_name)) {
    stop(
      "Peak filename does not match the expected MACS2 format: ",
      file_name,
      "\nExpected format: <sample>_<qvalue>_peaks.narrowPeak"
    )
  }
  return(sample_name)
}

empty_peak_dataframe <- function() {
  data.frame(
    chrom = character(),
    start = numeric(),
    end = numeric(),
    name = character(),
    score = numeric(),
    strand = character(),
    signalValue = numeric(),
    pValue = numeric(),
    qValue = numeric(),
    peak = numeric(),
    stringsAsFactors = FALSE
  )
}

make_empty_plot <- function(title_text) {
  ggplot() +
    geom_blank() +
    theme_void() +
    ggtitle(title_text)
}

# ============================================================
# Read peak files
# ============================================================

sample_list <- vapply(
  peak_files,
  extract_peak_sample_name,
  character(1)
)

peak_data_list <- list()
peak_count_list <- list()
peak_granges_list <- list()

for (file_path in peak_files) {

  if (!file.exists(file_path)) {
    stop("Peak file does not exist: ", file_path)
  }

  sample_name <- extract_peak_sample_name(file_path)

  metadata_match <- which(sample_metadata$sample == sample_name)

  if (length(metadata_match) != 1) {
    stop(
      "Peak-file sample was not found uniquely in samples.csv: ",
      sample_name
    )
  }

  histone <- sample_metadata$histone[metadata_match]
  replicate <- sample_metadata$replicate[metadata_match]

  if (file.info(file_path)$size == 0) {

    peak_info <- empty_peak_dataframe()

  } else {

    peak_info <- read.table(
      file_path,
      header = FALSE,
      sep = "\t",
      quote = "",
      comment.char = "",
      fill = TRUE,
      stringsAsFactors = FALSE,
      col.names = c(
        "chrom",
        "start",
        "end",
        "name",
        "score",
        "strand",
        "signalValue",
        "pValue",
        "qValue",
        "peak"
      )
    )

    peak_info <- peak_info %>%
      mutate(
        start = as.numeric(start),
        end = as.numeric(end),
        score = as.numeric(score),
        signalValue = as.numeric(signalValue),
        pValue = as.numeric(pValue),
        qValue = as.numeric(qValue),
        peak = as.numeric(peak)
      )
  }

  peak_count_list[[sample_name]] <- data.frame(
    sampleInfo = sample_name,
    Histone = histone,
    Replicate = replicate,
    peakN = nrow(peak_info),
    stringsAsFactors = FALSE
  )

  if (nrow(peak_info) > 0) {

    peak_info <- peak_info %>%
      mutate(
        width = end - start,
        Histone = histone,
        Replicate = replicate,
        sampleInfo = sample_name
      )

    peak_data_list[[sample_name]] <- peak_info

    # narrowPeak coordinates are zero-based, half-open.
    # GRanges uses one-based, closed coordinates.
    peak_granges_list[[sample_name]] <- GRanges(
      seqnames = peak_info$chrom,
      ranges = IRanges(
        start = peak_info$start + 1,
        end = peak_info$end
      )
    )

  } else {

    peak_granges_list[[sample_name]] <- GRanges()
  }
}

peakData <- bind_rows(peak_data_list)
peak_count_data <- bind_rows(peak_count_list)

peak_count_data$sampleInfo <- factor(
  peak_count_data$sampleInfo,
  levels = sample_list
)

peak_count_data$Histone <- factor(
  peak_count_data$Histone,
  levels = histone_levels
)

peak_count_data$Replicate <- as.factor(
  peak_count_data$Replicate
)

if (nrow(peakData) > 0) {

  peakData$sampleInfo <- factor(
    peakData$sampleInfo,
    levels = sample_list
  )

  peakData$Histone <- factor(
    peakData$Histone,
    levels = histone_levels
  )

  peakData$Replicate <- as.factor(
    peakData$Replicate
  )
}

# ============================================================
# Plot 1: Number of peaks
# ============================================================

fig1 <- ggplot(peak_count_data, aes(x = Histone, y = peakN, fill = Histone)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color = Replicate), width = 0.15, height = 0) +
  theme_bw(base_size = 18) +
  labs(x = "", y = "Number of Peaks", title = "Peak counts") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

# ============================================================
# Plot 2: Width of peaks
# ============================================================

if (nrow(peakData) > 0) {

  fig2 <- ggplot(peakData, aes( x = Histone, y = width, fill = Histone)) +
    geom_violin(scale = "width", trim = TRUE) +
    facet_grid(Replicate ~ .) +
    theme_bw(base_size = 18) +
    labs(x = "", y = "Width of Peaks", title = "Peak-width distributions") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

} else {

  fig2 <- make_empty_plot(
    "No peaks available for width analysis"
  )
}

# ============================================================
# Calculate pairwise peak reproducibility
# ============================================================

# Peak overlap is calculated only between samples belonging to
# the same histone group.
#
# For each sample pair:
#
#   Direction 1:
#     % of sample 1 peaks overlapping at least one sample 2 peak
#
#   Direction 2:
#     % of sample 2 peaks overlapping at least one sample 1 peak
#
#   Mean bidirectional overlap:
#     mean of the two directional percentages
#
# Any overlap of at least one base pair is counted.

reproducibility_summary_rows <- list()
reproducibility_directional_rows <- list()

for (hist in histone_levels) {

  hist_samples <- sample_metadata$sample[
    sample_metadata$histone == hist
  ]

  hist_samples <- hist_samples[
    hist_samples %in% names(peak_granges_list)
  ]

  if (length(hist_samples) < 2) {
    next
  }

  sample_combinations <- combn(hist_samples, 2, simplify = FALSE)

  for (pair in sample_combinations) {

    sample1 <- pair[1]
    sample2 <- pair[2]

    peaks1 <- peak_granges_list[[sample1]]
    peaks2 <- peak_granges_list[[sample2]]

    # A directional overlap percentage cannot be meaningfully
    # calculated if either sample has no peaks.
    if (length(peaks1) == 0 || length(peaks2) == 0) {
      warning(
        "Skipping peak-overlap comparison because at least one sample has no peaks: ",
        sample1,
        " vs ",
        sample2
      )
      next
    }

    sample1_overlapping_sample2 <- sum(
      countOverlaps(peaks1, peaks2) > 0
    )

    sample2_overlapping_sample1 <- sum(
      countOverlaps(peaks2, peaks1) > 0
    )

    pct_sample1_in_sample2 <- (
      100 *
        sample1_overlapping_sample2 /
        length(peaks1)
    )

    pct_sample2_in_sample1 <- (
      100 *
        sample2_overlapping_sample1 /
        length(peaks2)
    )

    mean_bidirectional_overlap <- mean(
      c(
        pct_sample1_in_sample2,
        pct_sample2_in_sample1
      )
    )

    replicate1 <- sample_metadata$replicate[
      sample_metadata$sample == sample1
    ]

    replicate2 <- sample_metadata$replicate[
      sample_metadata$sample == sample2
    ]

    pair_id <- paste(
      hist,
      sample1,
      sample2,
      sep = "__"
    )

    reproducibility_summary_rows[[length(reproducibility_summary_rows) + 1]] <- data.frame(
      PairID = pair_id,
      Histone = hist,
      Sample1 = sample1,
      Sample2 = sample2,
      Replicate1 = replicate1,
      Replicate2 = replicate2,
      Sample1Peaks = length(peaks1),
      Sample2Peaks = length(peaks2),
      Sample1OverlappingSample2 = sample1_overlapping_sample2,
      Sample2OverlappingSample1 = sample2_overlapping_sample1,
      PctSample1InSample2 = pct_sample1_in_sample2,
      PctSample2InSample1 = pct_sample2_in_sample1,
      MeanBidirectionalOverlap = mean_bidirectional_overlap,
      stringsAsFactors = FALSE
    )

    reproducibility_directional_rows[[length(reproducibility_directional_rows) + 1]] <- data.frame(
      PairID = pair_id,
      Histone = hist,
      Sample1 = sample1,
      Sample2 = sample2,
      Direction = paste0(sample1," -> ", sample2),
      DirectionOrder = "Sample 1 to Sample 2",
      SourceSample = sample1,
      TargetSample = sample2,
      SourcePeaks = length(peaks1),
      OverlappingPeaks = sample1_overlapping_sample2,
      OverlapRate = pct_sample1_in_sample2,
      stringsAsFactors = FALSE
    )

    reproducibility_directional_rows[[length(reproducibility_directional_rows) + 1]] <- data.frame(
      PairID = pair_id,
      Histone = hist,
      Sample1 = sample1,
      Sample2 = sample2,
      Direction = paste0(sample2," -> ", sample1),
      DirectionOrder = "Sample 2 to Sample 1",
      SourceSample = sample2,
      TargetSample = sample1,
      SourcePeaks = length(peaks2),
      OverlappingPeaks = sample2_overlapping_sample1,
      OverlapRate = pct_sample2_in_sample1,
      stringsAsFactors = FALSE
    )
  }
}

if (length(reproducibility_summary_rows) > 0) {

  reproducibility_summary <- bind_rows(
    reproducibility_summary_rows
  ) %>%
    arrange(
      desc(MeanBidirectionalOverlap),
      Histone,
      Sample1,
      Sample2
    ) %>%
    mutate(
      Rank = row_number(),
      PairLabel = sprintf(
        "%02d. %s | %s vs %s | mean %.1f%%",
        Rank,
        Histone,
        Sample1,
        Sample2,
        MeanBidirectionalOverlap
      ),
      OverviewLabel = paste0(Histone,": ", Sample1," vs ", Sample2)
    )

  reproducibility_directional <- bind_rows(
    reproducibility_directional_rows
  ) %>%
    left_join(
      reproducibility_summary %>%
        select(
          PairID,
          Rank,
          PairLabel,
          MeanBidirectionalOverlap
        ),
      by = "PairID"
    ) %>%
    arrange(
      Rank,
      DirectionOrder
    )

} else {

  reproducibility_summary <- data.frame()
  reproducibility_directional <- data.frame()
}

# ============================================================
# Plot 3: Compact reproducibility overview
# ============================================================

if (nrow(reproducibility_summary) > 0) {

  number_to_show <- min(
    overview_pair_limit,
    nrow(reproducibility_summary)
  )

  overview_data <- reproducibility_summary %>%
    slice_head(
      n = number_to_show
    ) %>%
    mutate(
      OverviewLabel = factor(
        OverviewLabel,
        levels = rev(OverviewLabel)
      ),
      Histone = factor(
        Histone,
        levels = histone_levels
      )
    )

  if (nrow(reproducibility_summary) > overview_pair_limit) {

    reproducibility_title <- paste0("Top ", overview_pair_limit," pairwise peak overlaps")

    reproducibility_subtitle <- paste0("All ", nrow(reproducibility_summary)," ranked pairs are shown in the PDF report")

  } else {

    reproducibility_title <- paste0("Pairwise peak reproducibility")

    reproducibility_subtitle <- paste0("Mean bidirectional overlap; at least 1 bp required")
  }

  fig3 <- ggplot(overview_data,aes(x = OverviewLabel, y = MeanBidirectionalOverlap, fill = Histone)) +
    geom_col(width = 0.7) +
    geom_text(aes(label = sprintf("%.1f%%", MeanBidirectionalOverlap)), hjust = -0.15, size = 3) +
    coord_flip(ylim = c(0, 108),clip = "off") +
    theme_bw(base_size = 12) +
    theme(legend.position = "none", axis.text.y = element_text(size = 7), panel.grid.major.y = element_blank(), 
          plot.margin = margin(t = 5.5, r = 25, b = 5.5, l = 5.5)) +
    labs(x = "", y = "Mean bidirectional overlap (%)", title = reproducibility_title, subtitle = reproducibility_subtitle)

} else {

  fig3 <- make_empty_plot("No pairwise peak reproducibility data available")
}

# ============================================================
# Plot 4: FRiP scores
# ============================================================

frip_data_list <- list()

for (file_path in frip_files) {

  if (!file.exists(file_path)) {
    stop("FRiP file does not exist: ", file_path)
  }

  if (file.info(file_path)$size == 0) {
    warning("Skipping empty FRiP file: ", file_path)
    next
  }

  frip_tmp <- read.table(file_path, header = TRUE, sep = "\t", quote = "", comment.char = "", stringsAsFactors = FALSE)

  frip_data_list[[length(frip_data_list) + 1]] <- frip_tmp
}

frip_data <- bind_rows(frip_data_list)

if (nrow(frip_data) == 0) {

  fig4 <- make_empty_plot(
    "No FRiP data available"
  )

} else {

  required_frip_columns <- c(
    "Sample",
    "FRiP"
  )

  missing_frip_columns <- setdiff(
    required_frip_columns,
    colnames(frip_data)
  )

  if (length(missing_frip_columns) > 0) {
    stop(
      "FRiP data is missing required columns: ",
      paste(missing_frip_columns, collapse = ", ")
    )
  }

  unknown_frip_samples <- setdiff(
    frip_data$Sample,
    sample_metadata$sample
  )

  if (length(unknown_frip_samples) > 0) {
    stop(
      "FRiP samples were not found in samples.csv: ",
      paste(unknown_frip_samples, collapse = ", ")
    )
  }

  frip_data$FRiP <- as.numeric(
    frip_data$FRiP
  )

  frip_data$Histone <- sample_metadata$histone[
    match(
      frip_data$Sample,
      sample_metadata$sample
    )
  ]

  frip_data$Replicate <- sample_metadata$replicate[
    match(
      frip_data$Sample,
      sample_metadata$sample
    )
  ]

  frip_data$Histone <- factor(
    frip_data$Histone,
    levels = histone_levels
  )

  frip_data$Replicate <- as.factor(
    frip_data$Replicate
  )

  fig4 <- ggplot(frip_data,aes(x = Histone, y = FRiP, fill = Histone)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(aes(color = Replicate), width = 0.15, height = 0) +
    theme_bw(base_size = 18) +
    labs(x = "", y = "% of Fragments in Peaks (FRiP)", title = "FRiP scores") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
}

# ============================================================
# Combine four-panel overview
# ============================================================

final_plot <- ggarrange(fig1, fig2, fig3, fig4, ncol = 2, nrow = 2, common.legend = TRUE, legend = "bottom", labels = c("A","B","C","D"))

# ============================================================
# Save overview PNG
# ============================================================

summary_png <- file.path(output_dir,"peak_summary_plot.png")

ggsave(filename = summary_png, plot = final_plot, width = 14, height = 10, dpi = 300)

# ============================================================
# Create multi-page PDF report
# ============================================================

report_pdf <- file.path(output_dir,"peak_summary_report.pdf")

pdf(file = report_pdf, width = 14, height = 10, onefile = TRUE)

# Page 1: four-panel summary
print(final_plot)

# Remaining pages: ranked directional pairwise comparisons
if (nrow(reproducibility_summary) > 0) {

  total_pairs <- nrow(
    reproducibility_summary
  )

  number_of_detail_pages <- ceiling(
    total_pairs / pairs_per_page
  )

  for (
    detail_page_number in seq_len(
      number_of_detail_pages
    )
  ) {

    first_rank <- (
      (detail_page_number - 1) *
        pairs_per_page
    ) + 1

    last_rank <- min(
      detail_page_number *
        pairs_per_page,
      total_pairs
    )

    page_summary <- reproducibility_summary %>%
      filter(
        Rank >= first_rank,
        Rank <= last_rank
      )

    page_directional <- reproducibility_directional %>%
      filter(
        Rank >= first_rank,
        Rank <= last_rank
      ) %>%
      mutate(
        PairLabel = factor(
          PairLabel,
          levels = page_summary$PairLabel
        ),
        DirectionOrder = factor(
          DirectionOrder,
          levels = c(
            "Sample 1 to Sample 2",
            "Sample 2 to Sample 1"
          )
        )
      ) %>%
      arrange(
        Rank,
        DirectionOrder
      ) %>%
      mutate(
        Direction = factor(
          Direction,
          levels = unique(Direction)
        )
      )

    page_mean_data <- page_summary %>%
      mutate(
        PairLabel = factor(
          PairLabel,
          levels = page_summary$PairLabel
        )
      )

    detail_plot <- ggplot(page_directional,aes(x = Direction, y = OverlapRate, fill = DirectionOrder)) +
      geom_col(width = 0.65) +
      geom_text(aes(label = sprintf("%.1f%%", OverlapRate)), vjust = -0.4, size = 3.2) +
      geom_hline(data = page_mean_data, aes(yintercept = MeanBidirectionalOverlap), inherit.aes = FALSE, linetype = "dashed", linewidth = 0.6) +
      facet_wrap(~ PairLabel, ncol = detail_page_columns, scales = "free_x") +
      scale_fill_viridis_d(option = "magma", begin = 0.2, end = 0.8) +
      coord_cartesian(ylim = c(0, 108), clip = "off") +
      theme_bw(base_size = 12) +
      theme(axis.text.x = element_text(angle = 25, hjust = 1, size = 8), legend.position = "bottom", strip.text = element_text(size = 9, face = "bold"),
            panel.grid.major.x = element_blank(), plot.margin = margin(t = 10, r = 15, b = 10, l = 10)) +
      labs(x = "", y = paste0("Source-sample peaks overlapping ","the paired sample (%)"), fill = "Direction",
           title = paste0("Ranked pairwise peak reproducibility: ","comparisons ", first_rank,"-", last_rank," of ", total_pairs),
        subtitle = paste0("Pairs are ordered by mean bidirectional overlap ","from highest to lowest. ","Dashed lines indicate the pair mean. ","At least 1 bp of overlap is required."))

    print(detail_plot)
  }
}

dev.off()

message("Saved peak summary PNG: ", summary_png)
message("Saved multi-page peak summary PDF: ", report_pdf)
