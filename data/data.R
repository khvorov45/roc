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

svnt_more <- read_raw("svnt-more") %>%
  # All id-measurement pairs are present in svnt, it's just measurement 2 that
  # I need
  filter(!is.na(measurement_2)) %>%
  select(id, measurement = measurement_2, cohort) %>%
  mutate(
    symptom_onset_days = NA_character_,
    id = as.character(id),
    assay = "svnt"
  )

svnt_final <- bind_rows(svnt, svnt_more)

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
all_data <- bind_rows(list(euro_ncp, euro_s1, svnt_final, wantai)) %>%
  filter(!is.na(measurement))

# Modify variables
unique(all_data$cohort)
unique(all_data$symptom_onset_days)
unique(all_data$assay)

all_data_mod <- all_data %>%
  mutate(
    group = case_when(
      str_detect(tolower(cohort), "pcr pos") ~ "covid",
      str_detect(tolower(cohort), "healthy control") ~ "healthy",
      TRUE ~ "non-covid"
    ),
    true_covid = group == "covid"
  ) %>%
  # Remove old iga for the non-covids since we are not interested in it
  filter(!(assay == "euro_iga" & !true_covid)) %>%
  mutate(assay = recode(assay, "euro_iga_new" = "euro_iga")) %>%
  # Remove the unneeded variables
  select(-cohort)

# Replace onset days for some covids with the newly provided
onset <- read_raw("onset") %>%
  select(id, symptom_onset_days = `Serum Days since Sx onset`) %>%
  mutate(
    symptom_onset_days = if_else(
      symptom_onset_days == "pending",
      35L,
      suppressWarnings(as.integer(symptom_onset_days))
    )
  ) %>%
  filter(!is.na(symptom_onset_days))

new_onset <- inner_join(onset, select(all_data_mod, -symptom_onset_days), "id")
old_onset <- all_data_mod %>% filter(!id %in% new_onset$id)

# Gotta fix the old onset
unique(old_onset$symptom_onset_days)

old_onset_fixed <- old_onset %>%
  mutate(
    symptom_onset_days = suppressWarnings(symptom_onset_days %>%
      str_replace("[+|>]", "") %>%
      as.integer(.))
  )

all_data_new_onset <- bind_rows(new_onset, old_onset_fixed)

# Create onset categories
all_data_onset_cats <- all_data_new_onset %>%
  mutate(
    symptom_onset_cat = case_when(
      symptom_onset_days < 7 ~ "<7",
      symptom_onset_days <= 14 ~ "7-14",
      symptom_onset_days > 14 ~ ">14",
      !true_covid ~ "no infection",
      is.na(symptom_onset_days) ~ NA_character_
    )
  ) %>%
  select(-symptom_onset_days)

# The "Lab ID" in raw data is SOMETIMES an individual id and SOMETIMES a sample
# id. There is NO systematic way to distinguish them, so I'll just trust
# Suellen when she tells me which ids belong to the same individual
all_data_fixids <- all_data_onset_cats %>%
  mutate(
    og_id = id,
    id = case_when(
      # Just random ids that are the same actually
      id %in% c("20510199", "20510596", "20111141") ~ "20510199",
      id %in% c(
        "20509672", "20509673", "20509674", "20509686", "20509687", "20509688"
      ) ~ "20509672",
      TRUE ~ id
    ) %>%
      # Handle CVD010 CVD010A situation
      str_replace("(\\d)[[:alpha:]]$", "\\1")
  ) %>%
  group_by(id, assay) %>%
  mutate(sample_id = paste(id, row_number(), sep = "-")) %>%
  ungroup()

# Look for multiple measurement from the same individual for the same assay
length(unique(all_data_fixids$id)) # Total unique individuals

# Number of individuals that provided more than 1 measurement per assay
all_data_fixids %>%
  count(id, assay, name = "n_samples") %>%
  filter(n_samples > 1) %>%
  group_by(id, n_samples) %>%
  summarise(n_assays = paste(assay, collapse = " "), .groups = "drop") %>%
  print(n = 50)

# Keep the FIRST observation per individual per assay.
# There is no date associated with repeat observations (of course there isn't)
# so just assume that the ones that appear first in the data are chronologicaly
# first (that's how they were recorded, so it should work)

all_data_firsts <- all_data_fixids %>%
  group_by(id, assay) %>%
  filter(row_number() == 1) %>%
  ungroup()

# Whoever doesn't have an onset category won't be used in the analysis
all_data_no_missing_onset <- all_data_firsts %>%
  filter(!is.na(symptom_onset_cat))

save_data(all_data_no_missing_onset, "data")
