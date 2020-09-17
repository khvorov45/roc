cat("extract data for analysis")

library(tidyverse)

data_raw_dir <- "data-raw"
data_dir <- "data"

# Functions ===================================================================

read_raw <- function(name, ...) {
  readxl::read_excel(file.path(data_raw_dir, paste0(name, ".xlsx")), ...)
}

lengthen_measurement <- function(data, name) {
  data %>%
    pivot_longer(
      contains("measurement"),
      names_to = "assay", values_to = "measurement"
    ) %>%
    mutate(assay = str_replace(assay, "measurement", name))
}

save_data <- function(data, name) {
  write_csv(data, file.path(data_dir, paste0(name, ".csv")))
}

# Script ======================================================================

euro_ncp <- read_raw("euro-ncp", range = "A2:L319") %>%
  select(
    id,
    cohort = notes,
    measurement = Index,
    symptom_onset_days
  ) %>%
  mutate(assay = "euro_ncp")

euro_s1 <- read_raw("euro-s1", range = "A2:T308") %>%
  select(
    id,
    cohort = notes,
    measurement_igg = Index_igg,
    measurement_iga = Index_iga,
    measurement_iga_new = Index_iga_new,
    symptom_onset_days
  ) %>%
  lengthen_measurement("euro")

svnt <- read_raw("svnt", range = "A3:L444") %>%
  select(
    id,
    cohort = notes,
    measurement = `% Inhib`,
    symptom_onset_days
  ) %>%
  mutate(assay = "svnt")

wantai <- read_raw("wantai", range = "B3:O349") %>%
  select(
    id,
    cohort,
    measurement_tot = sco_tot,
    measurement_igm = sco_igm,
    symptom_onset_days
  ) %>%
  lengthen_measurement("wantai")

# Unite them
all_data <- bind_rows(list(euro_ncp, euro_s1, svnt, wantai)) %>%
  filter(!is.na(measurement))

# Look for multiple measurement from the same individual for the same assay
length(unique(all_data$id)) # Total unique individuals

# Number of individuals that provided more than 1 measurement per assay
all_data %>%
  count(id, assay) %>%
  filter(n > 1)

# Since the number of individuals that provided more than 1 measument per assay
# (more than 1 being only 2 for all of these cases)
# is much smaller than the total number of unique individuals, I'll ignore this
# and pretend that all observations are independent.

save_data(all_data, "data")
