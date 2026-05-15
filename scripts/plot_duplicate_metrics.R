#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
})

args <- commandArgs(trailingOnly = TRUE)

dup_files_arg <- args[1]
down_files_arg <- args[2]
metadata_file <- args[3]
dup_summary_out <- args[4]
down_summary_out <- args[5]
plot_out <- args[6]

dup_files <- strsplit(dup_files_arg, " ")[[1]]
down_files <- strsplit(down_files_arg, " ")[[1]]

metadata <- read.csv(metadata_file, stringsAsFactors = FALSE)

dup_df <- dup_files |>
  purrr::map_dfr(read.delim, stringsAsFactors = FALSE)

down_df <- down_files |>
  purrr::map_dfr(read.delim, stringsAsFactors = FALSE)

dup_df <- dup_df |>
  left_join(metadata, by = c("sample" = "sample"))

down_df <- down_df |>
  left_join(metadata, by = c("Sample" = "sample"))

write.table(
  dup_df,
  dup_summary_out,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  down_df,
  down_summary_out,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

plot_df <- dup_df |>
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

p <- ggplot(
  plot_df,
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
    title = "PCR/coordinate duplicate burden before and after duplicate capping"
  )

ggsave(plot_out, p, width = 12, height = 7, dpi = 300)
