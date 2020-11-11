cat("make data tables")

library(tidyverse)

data_dir <- "data"
data_table_dir <- "data-table"

# Functions ===================================================================

source(file.path(data_dir, "read_data.R"))

format_percent <- function(p, l, h) {
  f <- function(x) glue::glue("{signif(x * 100, 2)}%")
  glue::glue("{f(p)} ({f(l)}, {f(h)})")
}

summarise_binary <- function(bin_vec) {
  n <- length(bin_vec)
  s <- sum(bin_vec)
  f <- n - s
  p <- s / n
  ci <- PropCIs::exactci(s, n, 0.95)$conf.int
  tibble(
    total = n,
    success = s,
    failure = f,
    prop = p,
    low = ci[[1]],
    high = ci[[2]],
    summary = glue::glue("{format_percent(prop, low, high)} [{s} / {n}]")
  )
}

merge_groups <- function(onset, group) {
  onset_chr <- as.character(onset)
  group_chr <- as.character(group)
  if_else(group_chr == "covid", onset_chr, group_chr) %>%
    factor(levels = c(levels(onset), levels(group)))
}

save_data <- function(data, name) {
  write_csv(
    data,
    file.path(data_table_dir, paste0(name, ".csv")),
  )
  data
}

# Script ======================================================================

data <- read_data("data")

# Counts ----------------------------------------------------------------------

data_no_svnt_extra <- data %>%
  filter(!assay %in% c("svnt-20", "svnt-25"))

# Assay samples vs unique individuals
s <- function(data) {
  data %>%
    filter(!is.na(symptom_onset_cat)) %>%
    summarise(
      n_indiv = length(unique(id)),
      n_samples = n(),
      summary = glue::glue("{n_samples} ({n_indiv})"),
      .groups = "drop"
    ) %>%
    select(-n_indiv, -n_samples)
}
assay_counts_onsets <- data_no_svnt_extra %>%
  group_by(assay, symptom_onset_cat) %>%
  s() %>%
  rename(subset = symptom_onset_cat)
assay_counts_group <- data_no_svnt_extra %>%
  group_by(assay, group) %>%
  s() %>%
  rename(subset = group)
assay_counts_overall <- data_no_svnt_extra %>%
  group_by(assay) %>%
  s() %>%
  mutate(subset = "combined")

assay_counts <- bind_rows(
  list(assay_counts_onsets, assay_counts_group, assay_counts_overall)
) %>%
  pivot_wider(names_from = "assay", values_from = "summary") %>%
  mutate(subset = factor(
    subset,
    c(
      "<7", "7-14", ">14", "covid",
      "population", "cross-reactive", "no infection",
      "combined"
    ),
    c(
      "<7", "7-14", ">14", "covid total",
      "population", "cross-reactive", "no covid total",
      "overall total"
    )
  )) %>%
  arrange(subset)

save_data(assay_counts, "assay-counts")

# MN validation ---------------------------------------------------------------

mn <- read_data("mn")

mn_summ_indiv_groups <- mn %>%
  group_by(assay, symptom_onset_cat, group) %>%
  summarise(summarise_binary(result == mn), .groups = "drop")

mn_summ_onset_averaged <- mn_summ_indiv_groups %>%
  filter(group == "covid") %>%
  group_by(group, assay) %>%
  summarise(
    symptom_onset_cat = factor("onset-averaged"),
    across(c(prop, low, high), mean),
    summary = format_percent(prop, low, high),
    .groups = "drop"
  )

mn_summ <- bind_rows(mn_summ_indiv_groups, mn_summ_onset_averaged)

mn_summ %>%
  mutate(group_lbl = merge_groups(symptom_onset_cat, group)) %>%
  select(assay, group_lbl, summary) %>%
  pivot_wider(names_from = "assay", values_from = "summary") %>%
  arrange(group_lbl) %>%
  save_data("mn-agreement")
