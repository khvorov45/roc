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

data <- read_data("data")

onset_assay_counts <- data %>%
  count(assay, symptom_onset_cat) %>%
  pivot_wider(names_from = "assay", values_from = "n")

save_data(onset_assay_counts, "onset-assay-counts")

f <- function(x) paste(x, collapse = " ")
assay_discrepancies <- data %>%
  calc_result_one_threshold() %>%
  group_by(sample_id, true_covid) %>%
  summarise(tibble(
    pos = f(assay[result == "pos"]),
    neg = f(assay[result == "neg"])
  ), .groups = "drop") %>%
  count(pos, neg, true_covid, name = "count") %>%
  arrange(true_covid, desc(count))

save_data(assay_discrepancies, "assay-discrepancies")
