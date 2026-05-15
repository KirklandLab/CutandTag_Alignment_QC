#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)

dup_files_arg <- args[1]
down_files_arg <- args[2]
down_target_file <- args[3]
metadata_file <- args[4]
combined_summary_out <- args[5]
plot_out <- args[6]

dup_files <- strsplit(dup_files_arg, " ")[[1]]
down_files <- strsplit(down_files_arg, " ")[[1]]

metadata <- read.csv(metadata_file, stringsAsFactors = FALSE)

dup_df <- dup_files |>
  purrr::map_dfr(read.delim, stringsAsFactors = FALSE)

down_df <- down_files |>
  purrr::map_dfr(read.delim, stringsAsFactors = FALSE)

down_target <- read.delim(down_target_file, stringsAsFactors = FALSE)

combined_df <- dup_df |>
  left_join(down_df, by = c("sample" = "Sample")) |>
  left_join(metadata, by = "sample") |>
  mutate(
    final_analysis_fragments = FragmentsAfterDownsample,
    pct_fragments_retained_after_dup_cap =
      ifelse(total_fragments_before > 0,
             100 * kept_fragments_after / total_fragments_before,
             NA_real_),
    pct_fragments_retained_final =
      ifelse(total_fragments_before > 0,
             100 * final_analysis_fragments / total_fragments_before,
             NA_real_)
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
    )
  )

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
      total_fragments_before = "Raw fragments",
      kept_fragments_after = "After duplicate cap",
      final_analysis_fragments = "Final analysis fragments"
    ),
    metric = factor(
      metric,
      levels = c("Raw fragments", "After duplicate cap", "Final analysis fragments")
    )
  )

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
    title = "Duplicate burden before and after duplicate capping"
  )

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
    title = "Fragments retained after duplicate capping and downsampling"
  )

p <- p1 / p2

ggsave(plot_out, p, width = 14, height = 10, dpi = 300)
