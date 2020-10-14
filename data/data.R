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
      c(contains("measurement"), contains("result")),
      names_to = c(".value", "assay"),
      names_pattern = "([^_]*)_(.*)"
    ) %>%
    mutate(assay = paste(name, assay, sep = "_"))
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
    symptom_onset_days,
    result = Interpretation
  ) %>%
  mutate(assay = "euro_ncp")

euro_s1 <- read_raw("euro-s1", range = "A2:T308") %>%
  select(
    id,
    cohort = notes,
    measurement_igg = Index_igg,
    measurement_iga = Index_iga,
    measurement_iga_new = Index_iga_new,
    result_igg = Interpretation_igg,
    result_iga = Interpretation_iga,
    result_iga_new = Interpretation_iga_new,
    symptom_onset_days
  ) %>%
  lengthen_measurement("euro")

svnt <- read_raw("svnt", range = "A3:L444") %>%
  select(
    id,
    cohort = notes,
    measurement = `% Inhib`,
    result = `sVNT Int.`,
    symptom_onset_days
  ) %>%
  mutate(assay = "svnt")

wantai <- read_raw("wantai", range = "B3:O349") %>%
  select(
    id,
    cohort,
    measurement_tot = sco_tot,
    measurement_igm = sco_igm,
    result_tot = interpret_tot,
    result_igm = interpret_igm,
    symptom_onset_days
  ) %>%
  lengthen_measurement("wantai")

# Unite them
all_data <- bind_rows(list(euro_ncp, euro_s1, svnt, wantai)) %>%
  filter(!is.na(measurement))

# There shouldn't be any missing results
all_data %>%
  filter(is.na(result))

# Modify variables
unique(all_data$cohort)
unique(all_data$symptom_onset_days)
unique(all_data$assay)
unique(all_data$result)

all_data_mod <- all_data %>%
  mutate(
    group = case_when(
      str_detect(tolower(cohort), "pcr pos") ~ "covid",
      str_detect(tolower(cohort), "healthy control") ~ "healthy",
      str_detect(tolower(cohort), "corona") ~ "corona",
      str_detect(tolower(cohort), "negative control") ~ "neg-control",
      str_detect(tolower(cohort), "neg control") ~ "neg-control",
      TRUE ~ "non-covid"
    ),
    result = tolower(result) %>%
      str_replace_all("^equ$", "equiv") %>%
      str_replace_all("/equ$", "/equiv") %>%
      str_replace_all("/equ/", "/equiv/") %>%
      str_replace_all("^eqiv$", "equiv") %>%
      str_replace_all("equiv", "pos") %>%
      str_replace_all("p0s", "pos") %>%
      str_replace_all("pospos", "pos/pos") %>%
      str_replace_all("([[:alpha:]]{3})x2", "\\1/\\1") %>%
      map(., ~ str_split(.x, "/")) %>%
      map(flatten) %>%
      map_chr(
        .,
        function(res_vec) {
          if (sum(res_vec == "pos") > sum(res_vec == "neg")) "pos" else "neg"
        }
      ),
  ) %>%
  # Remove old iga for the non-covids since we are not interested in it
  filter(!(assay == "euro_iga" & group != "covid")) %>%
  mutate(assay = recode(assay, "euro_iga_new" = "euro_iga"))

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

all_data_new_onset <- bind_rows(new_onset, old_onset_fixed) %>%
  # Can't use covids with missing onset
  filter(!(group == "covid" & is.na(symptom_onset_days)))

# Create onset categories
all_data_onset_cats <- all_data_new_onset %>%
  mutate(
    symptom_onset_cat = case_when(
      symptom_onset_days < 7 ~ "<7",
      symptom_onset_days <= 14 ~ "7-14",
      symptom_onset_days > 14 ~ ">14",
      group != "covid" ~ "no infection",
    )
  )

# Shouldn't be any missing symptom onset category
all_data_onset_cats %>%
  filter(is.na(symptom_onset_cat))

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
  ungroup() %>%
  select(id, sample_id, everything()) %>%
  select(-og_id)

# Keep the first observation for covids going by symptom onset to avoid
# potentially correlated observations
all_data_covid_firsts <- all_data_fixids %>%
  group_by(id, assay, symptom_onset_cat) %>%
  filter(group != "covid" | symptom_onset_days == min(symptom_onset_days)) %>%
  ungroup()

# For id's that appear in seasonal corona and non-covid, keep only the
# seasonal corona
all_data_only_seas_corona <- all_data_covid_firsts %>%
  group_by(id, assay, symptom_onset_cat) %>%
  filter(n() == 1 | (!"corona" %in% group | group == "corona")) %>%
  filter(n() == 1 | (!"non-covid" %in% group | group == "non-covid")) %>%
  ungroup() %>%
  mutate(
    group = recode(group, "corona" = "non-covid", "neg-control" = "non-covid")
  )

# Look for multiple measurement from the same individual for the same assay
length(unique(all_data_only_seas_corona$id)) # Total unique individuals

# Number of individuals that provided more than 1 measurement per assay
all_data_only_seas_corona %>%
  count(id, assay, symptom_onset_cat, group, name = "n_samples") %>%
  filter(n_samples > 1) %>%
  group_by(id, n_samples, group) %>%
  summarise(n_assays = paste(assay, collapse = " "), .groups = "drop") %>%
  print(n = 50)

# At this point just keep the first observation for everyone

all_data_one_ind <- all_data_only_seas_corona %>%
  group_by(id, assay, symptom_onset_cat) %>%
  filter(row_number() == 1) %>%
  ungroup()

# Add the two extra assays - svnt 1st observation at 20 and svnt 1st observation
# at 25
calc_svnt_extra <- function(data, threshold) {
  data %>%
    filter(assay == "svnt") %>%
    mutate(
      result = if_else(measurement < threshold, "neg", "pos"),
      assay = paste0("svnt-", threshold)
    )
}
svnt_20 <- calc_svnt_extra(all_data_one_ind, 20)
svnt_25 <- calc_svnt_extra(all_data_one_ind, 25)

all_data_svnt_extra <- bind_rows(list(all_data_one_ind, svnt_20, svnt_25))

all_data_final <- all_data_svnt_extra %>%
  select(id, group, symptom_onset_cat, assay, measurement, result)

save_data(all_data_final, "data")
