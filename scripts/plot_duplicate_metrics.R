#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 6) {
  stop(
    "Usage: Rscript plot_duplicate_metrics.R ",
    "<dup_metrics|NA> <down_metrics|NA> <down_target|NA> ",
    "<samples.csv> <combined_summary_out> <plot_out>"
  )
}

dup_files_arg <- args[1]
down_files_arg <- args[2]
down_target_file <- args[3]
metadata_file <- args[4]
combined_summary_out <- args[5]
plot_out <- args[6]

metadata <- read.csv(metadata_file, stringsAsFactors = FALSE)

required_metadata_cols <- c("sample", "histone", "replicate")
missing_metadata_cols <- setdiff(required_metadata_cols, colnames(metadata))

if (length(missing_metadata_cols) > 0) {
  stop(
    "Metadata file is missing required columns: ",
    paste(missing_metadata_cols, collapse = ", ")
  )
}

# -------------------------
# Helper: format large numbers without scientific notation
# -------------------------
format_full_number <- function(x) {
  x_num <- suppressWarnings(as.numeric(x))

  out <- format(
    x_num,
    big.mark = ",",
    scientific = FALSE,
    trim = TRUE
  )

  out[is.na(x_num)] <- "NA"

  return(out)
}

# -------------------------
# Load downsampling metrics
# -------------------------
if (down_files_arg == "NA") {
  down_df <- NULL
} else {
  down_files <- strsplit(down_files_arg, " ")[[1]]
  down_files <- down_files[down_files != ""]

  down_df <- down_files |>
    purrr::map_dfr(read.delim, stringsAsFactors = FALSE)

  if (!"DownsamplingEnabled" %in% colnames(down_df)) {
    down_df$DownsamplingEnabled <- NA
  }
}

# -------------------------
# Load downsampling target
# -------------------------
if (down_target_file == "NA") {
  down_target <- data.frame(
    target_fragments = NA_real_,
    mode = "not_available",
    lowest_fragments = NA_real_,
    lowest_above_floor = NA_real_,
    floor = NA_real_,
    manual_target = NA_real_,
    downsampling_enabled = NA,
    note = "target_file_not_available"
  )
} else {
  down_target <- read.delim(down_target_file, stringsAsFactors = FALSE)

  # Backward compatibility if older target files are used.
  if (!"lowest_above_floor" %in% colnames(down_target)) {
    down_target$lowest_above_floor <- NA_real_
  }
  if (!"downsampling_enabled" %in% colnames(down_target)) {
    down_target$downsampling_enabled <- NA
  }
}

# -------------------------
# Load duplicate metrics or create no-cap placeholder
# -------------------------
if (dup_files_arg == "NA") {

  if (is.null(down_df)) {
    stop("Both duplicate metrics and downsampling metrics are NA. Nothing to summarize.")
  }

  dup_df <- down_df |>
    transmute(
      sample = Sample,
      max_dups = "off",
      total_fragments_before = FragmentsBeforeDownsample,
      unique_fragment_positions = NA_real_,
      duplicate_fragments_before = NA_real_,
      pct_duplicate_before = NA_real_,
      kept_fragments_after = FragmentsBeforeDownsample,
      removed_fragments = 0,
      duplicate_fragments_after = NA_real_,
      pct_duplicate_after = NA_real_
    )

} else {
  dup_files <- strsplit(dup_files_arg, " ")[[1]]
  dup_files <- dup_files[dup_files != ""]

  dup_df <- dup_files |>
    purrr::map_dfr(read.delim, stringsAsFactors = FALSE)
}

# -------------------------
# If downsampling metrics are absent, create no-downsampling placeholder
# -------------------------
if (is.null(down_df)) {
  down_df <- dup_df |>
    transmute(
      Sample = sample,
      Action = "not_downsampled",
      TargetMode = "disabled",
      DownsamplingEnabled = FALSE,
      TargetNote = "downsampling_disabled",
      TargetFragments = NA_real_,
      Seed = "NA",
      SubsampleFraction = 1,
      FragmentsBeforeDownsample = kept_fragments_after,
      FragmentsAfterDownsample = kept_fragments_after
    )
}

# -------------------------
# Combine duplicate, downsampling, and metadata
# -------------------------
combined_df <- dup_df |>
  left_join(down_df, by = c("sample" = "Sample")) |>
  left_join(metadata, by = "sample") |>
  mutate(
    final_analysis_fragments = FragmentsAfterDownsample,
    pct_fragments_retained_after_dup_cap =
      ifelse(
        total_fragments_before > 0,
        100 * kept_fragments_after / total_fragments_before,
        NA_real_
      ),
    pct_fragments_retained_final =
      ifelse(
        total_fragments_before > 0,
        100 * final_analysis_fragments / total_fragments_before,
        NA_real_
      )
  ) |>
  select(
    sample,
    histone,
    replicate,
    max_dups,
    total_fragments_before,
    unique_fragment_positions,
    duplicate_fragments_before,
    pct_duplicate_before,
    kept_fragments_after,
    removed_fragments,
    duplicate_fragments_after,
    pct_duplicate_after,
    pct_fragments_retained_after_dup_cap,
    Action,
    TargetMode,
    DownsamplingEnabled,
    TargetNote,
    TargetFragments,
    Seed,
    SubsampleFraction,
    FragmentsBeforeDownsample,
    FragmentsAfterDownsample,
    final_analysis_fragments,
    pct_fragments_retained_final,
    everything()
  )

write.table(
  combined_df,
  combined_summary_out,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# -------------------------
# Duplicate plot dataframe
# -------------------------
dup_plot_df <- combined_df |>
  select(
    sample,
    histone,
    replicate,
    pct_duplicate_before,
    pct_duplicate_after
  ) |>
  pivot_longer(
    cols = c(pct_duplicate_before, pct_duplicate_after),
    names_to = "metric",
    values_to = "percent_duplicate"
  ) |>
  mutate(
    metric = recode(
      metric,
      pct_duplicate_before = "Before duplicate cap",
      pct_duplicate_after = "After duplicate cap"
    ),
    metric = factor(
      metric,
      levels = c("Before duplicate cap", "After duplicate cap")
    )
  )

has_dup_data <- any(!is.na(dup_plot_df$percent_duplicate))

# -------------------------
# Fragment/read retention plot dataframe
# -------------------------
read_plot_df <- combined_df |>
  select(
    sample,
    histone,
    replicate,
    total_fragments_before,
    kept_fragments_after,
    final_analysis_fragments
  ) |>
  pivot_longer(
    cols = c(total_fragments_before, kept_fragments_after, final_analysis_fragments),
    names_to = "metric",
    values_to = "fragments"
  ) |>
  mutate(
    metric = recode(
      metric,
      total_fragments_before = "Input fragments",
      kept_fragments_after = "After duplicate cap / pre-downsample",
      final_analysis_fragments = "Final analysis fragments"
    ),
    metric = factor(
      metric,
      levels = c(
        "Input fragments",
        "After duplicate cap / pre-downsample",
        "Final analysis fragments"
      )
    )
  )

# -------------------------
# Removed-fragment labels
# -------------------------
# Restores labels such as -10,000 or -250,000 on the duplicate burden plot.
# Labels are only shown when duplicate capping actually removed fragments.
removed_label_df <- combined_df |>
  filter(!is.na(removed_fragments), removed_fragments > 0) |>
  mutate(
    label = paste0("-", format_full_number(removed_fragments)),
    label_y = pmax(pct_duplicate_before, pct_duplicate_after, na.rm = TRUE)
  )

# -------------------------
# Plot titles
# -------------------------
dup_cap_label <- unique(combined_df$max_dups)
dup_cap_label <- dup_cap_label[!is.na(dup_cap_label)]

if (length(dup_cap_label) == 1 && dup_cap_label == "off") {
  dup_cap_title <- "Duplicate capping disabled"
} else if (length(dup_cap_label) == 1) {
  dup_cap_title <- paste0(
    "Duplicate burden before and after duplicate capping (cap = ",
    dup_cap_label,
    ")"
  )
} else {
  dup_cap_title <- "Duplicates before and after duplicate capping"
}

downsample_target <- suppressWarnings(as.numeric(down_target$target_fragments[1]))
downsample_mode <- as.character(down_target$mode[1])
downsample_target_label <- format_full_number(downsample_target)

if (!is.na(downsample_mode) && downsample_mode %in% c("disabled", "no_downsampling")) {
  if (!is.na(downsample_target) && downsample_target > 0) {
    downsample_title <- paste0(
      "Fragments retained after processing; BigWig scaling target = ",
      downsample_target_label,
      " (mode = ",
      downsample_mode,
      ")"
    )
  } else {
    downsample_title <- "Fragments retained after processing; downsampling disabled"
  }

} else if (!is.na(downsample_target) && downsample_target > 0) {
  downsample_title <- paste0(
    "Fragments retained after processing and downsampling",
    " (target = ",
    downsample_target_label,
    "; mode = ",
    downsample_mode,
    ")"
  )

} else {
  downsample_title <- "Fragments retained after processing and downsampling"
}

# -------------------------
# Duplicate plot
# -------------------------
if (has_dup_data) {

  p1 <- ggplot(
    dup_plot_df,
    aes(x = sample, y = percent_duplicate, fill = metric)
  ) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    facet_wrap(~ histone, scales = "free_x") +
    theme_bw(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.major.x = element_blank()
    ) +
    labs(
      x = "Sample",
      y = "% duplicate fragments",
      fill = "",
      title = dup_cap_title
    )

  if (nrow(removed_label_df) > 0) {
    p1 <- p1 +
      geom_text(
        data = removed_label_df,
        aes(
          x = sample,
          y = label_y,
          label = label
        ),
        inherit.aes = FALSE,
        vjust = -0.4,
        size = 3
      )
  }

  max_dup_y <- max(dup_plot_df$percent_duplicate, na.rm = TRUE)

  if (is.finite(max_dup_y) && max_dup_y > 0) {
    p1 <- p1 + expand_limits(y = max_dup_y * 1.15)
  }

} else {

  p1 <- ggplot(
    combined_df,
    aes(x = sample, y = 0)
  ) +
    geom_col(width = 0.7) +
    facet_wrap(~ histone, scales = "free_x") +
    theme_bw(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.major.x = element_blank()
    ) +
    labs(
      x = "Sample",
      y = "% duplicate fragments",
      title = "Duplicate capping disabled"
    )
}

# -------------------------
# Fragment/read retention plot
# -------------------------
p2 <- ggplot(
  read_plot_df,
  aes(x = sample, y = fragments, fill = metric)
) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  facet_wrap(~ histone, scales = "free_x") +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank()
  ) +
  labs(
    x = "Sample",
    y = "Mapped fragments / read pairs",
    fill = "",
    title = downsample_title
  )

p <- p1 / p2

ggsave(plot_out, p, width = 14, height = 10, dpi = 300)
