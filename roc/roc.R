cat("investigate test characteristics")

library(tidyverse)

data_dir <- "data"
roc_dir <- "roc"

# Functions ===================================================================

source(file.path(data_dir, "read_data.R"))

determine_result <- function(negative, positive) {
  function(measurement) {
    case_when(
      measurement < negative ~ "neg",
      measurement <= positive ~ "eqiv",
      measurement > positive ~ "pos"
    )
  }
}

calc_test_char <- function(data) {
  data %>%
    group_by(true_covid) %>%
    summarise(
      positive = sum(result == "pos"),
      negative = sum(result == "neg"),
      total = n(),
      .groups = "drop"
    )
}

# Script ======================================================================

data <- read_data("data")

data %>% count(assay)

data %>%
  mutate(
    result = case_when(
      startsWith(assay, "euro") ~ determine_result(0.8, 1.1)(measurement),
      startsWith(assay, "wantai") ~ determine_result(0.9, 1.1)(measurement),
      assay == "svnt" ~ determine_result(0.2, 0.2)(measurement),
    )
  ) %>%
  group_by(assay) %>%
  group_modify(~ calc_test_char(.x))
