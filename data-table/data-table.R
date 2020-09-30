cat("make data tables")

library(tidyverse)

data_dir <- "data"
data_table_dir <- "data-table"

# Functions ===================================================================

source(file.path(data_dir, "read_data.R"))

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
