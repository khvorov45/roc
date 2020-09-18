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

# Modify variables
unique(all_data$cohort)
unique(all_data$symptom_onset_days)
unique(all_data$assay)

all_data %>%
  mutate(
    true_covid = str_detect(tolower(cohort), "pcr pos"),
    # Categorise onset to <=8, 9-14, >=15
    symptom_onset_cat = suppressWarnings(as.integer(
      # How hard is it code a goddamn integer duration consistently jesus christ
      recode(
        symptom_onset_days,
        # Look at this beautiful inconsistent spacing - almost all possible
        # variations!
        "A= 0 to 3" = "1",
        "B= 4 to 8" = "6",
        "C=9 to 14" = "12",
        "D = 15 to 20" = "18",
        "E= 21 to 30" = "25"
      ) %>%
        str_replace("[>|+]", "") %>%
        str_replace("More than ", "")
    )) %>%
      cut(c(-Inf, 8, 14, Inf)) %>%
      recode("(-Inf,8]" = "<=8", "(8,14]" = "9-14", "(14, Inf]" = ">=15")
  ) %>%
  select(-cohort, -symptom_onset_days)

save_data(all_data, "data")
