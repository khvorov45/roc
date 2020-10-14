cat("make data tables")

library(tidyverse)

data_dir <- "data"
data_table_dir <- "data-table"

# Functions ===================================================================

source(file.path(data_dir, "read_data.R"))

source(file.path(data_dir, "calc_result_one_threshold.R"))

save_data <- function(data, name) {
  write_csv(
    data,
    file.path(data_table_dir, paste0(name, ".csv")),
  )
}

# Script ======================================================================

data <- read_data("data") %>%
  filter(!assay %in% c("svnt-20", "svnt-25"))

onset_assay_counts <- data %>%
  count(assay, symptom_onset_cat) %>%
  pivot_wider(names_from = "assay", values_from = "n")

save_data(onset_assay_counts, "onset-assay-counts")

f <- function(x) paste(x, collapse = " ")
assay_discrepancies <- data %>%
  calc_result_one_threshold() %>%
  group_by(id, group) %>%
  summarise(tibble(
    pos = f(assay[result == "pos"]),
    neg = f(assay[result == "neg"])
  ), .groups = "drop") %>%
  count(pos, neg, group, name = "count") %>%
  arrange(group, desc(count))

save_data(assay_discrepancies, "assay-discrepancies")

# Assay samples vs unique individuals
s <- function(data) {
  data %>%
    filter(!is.na(symptom_onset_cat)) %>%
    summarise(
      n_indiv = length(unique(id)),
      n_samples = n(),
      summary = glue::glue("{n_samples}"),
      .groups = "drop"
    ) %>%
    select(-n_indiv, -n_samples)
}
assay_counts_onsets <- data %>%
  group_by(assay, symptom_onset_cat) %>%
  s() %>%
  rename(subset = symptom_onset_cat)
assay_counts_group <- data %>%
  group_by(assay, group) %>%
  s() %>%
  rename(subset = group)
assay_counts_overall <- data %>%
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
      "healthy", "non-covid", "no infection",
      "combined"
    ),
    c(
      "<7", "7-14", ">14", "covid total",
      "healthy", "non-covid infection", "no covid total",
      "overall total"
    )
  )) %>%
  arrange(subset)

save_data(assay_counts, "assay-counts")
