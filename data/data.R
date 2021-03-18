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

fix_result <- function(result) {
  tolower(result) %>%
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
    )
}

fix_days_onset <- function(onset) {
  onset %>%
    str_replace("[+|>]", "") %>%
    as.integer()
}

categorise_onset <- function(onset, group) {
  case_when(
    onset < 7 ~ "<7",
    onset <= 14 ~ "7-14",
    onset > 14 ~ ">14",
    group != "covid" ~ "no infection",
  )
}

fix_group <- function(group) {
  group <- tolower(group)
  case_when(
    str_detect(group, "pcr pos") ~ "covid",
    str_detect(group, "pos pcr") ~ "covid",
    str_detect(group, "healthy control") ~ "population",
    str_detect(group, "negative control") ~ "population",
    str_detect(group, "neg control") ~ "population",
    str_detect(group, "neg pop") ~ "population",
    str_detect(group, "corona") ~ "cross-reactive",
    str_detect(group, "non-covid") ~ "cross-reactive",
    str_detect(group, "cross reac") ~ "cross-reactive",
    str_detect(group, "mers") ~ "cross-reactive",
    TRUE ~ NA_character_
  )
}

# Script ======================================================================

euro_ncp <- read_raw("euro-ncp", range = "A2:L319") %>%
  select(
    id,
    cohort = notes,
    measurement = Index,
    symptom_onset_days,
    result = Interpretation,
    diag = Diagnosis,
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
    symptom_onset_days,
    diag = Diagnosis,
  ) %>%
  lengthen_measurement("euro")

svnt <- read_raw("svnt", range = "A3:L444") %>%
  select(
    id,
    cohort = notes,
    measurement = `% Inhib`,
    result = `sVNT Int.`,
    symptom_onset_days,
    diag = Diagnosis,
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
    symptom_onset_days,
    diag,
  ) %>%
  lengthen_measurement("wantai")

more_results <- read_raw("more-results", range = "A4:AB77") %>%
  select(
    id = ID,
    symptom_onset_days = `Days post symptom onset`,
    measurement_euro_iga = euro_iga_Index,
    result_euro_iga = euro_iga_result,
    measurement_euro_igg = euro_igg_Index,
    result_euro_igg = euro_igg_result,
    measurement_euro_ncp = euro_ncp_Index,
    result_euro_ncp = euro_ncp_result,
    measurement_svnt = `svnt_% Inhib`,
    result_svnt_1 = svnt_first_result,
    result_svnt_2 = svnt_add_result,
  ) %>%
  mutate(
    result_svnt = if_else(is.na(result_svnt_2), result_svnt_1, result_svnt_2)
  ) %>%
  select(-result_svnt_1, -result_svnt_2) %>%
  lengthen_measurement("") %>%
  mutate(
    cohort = "pcr pos",
    diag = "SARS-CoV-2",
    assay = str_replace(assay, "^_", "")
  )

# No missing data is expected in extra results
any(!complete.cases(more_results))

# Unite them
all_data <- bind_rows(list(euro_ncp, euro_s1, svnt, wantai, more_results)) %>%
  filter(!is.na(measurement))

# There shouldn't be any missing results
all_data %>%
  filter(is.na(result))

# Modify variables
unique(all_data$cohort)
unique(all_data$symptom_onset_days)
unique(all_data$assay)
unique(all_data$result)
unique(all_data$diag)

all_data_mod <- all_data %>%
  mutate(
    # Group from cohort
    group = fix_group(cohort),
    # Result
    result_og = result,
    result = fix_result(result),
  ) %>%
  # Remove old iga for the non-covids since we are not interested in it
  filter(!(assay == "euro_iga" & group != "covid")) %>%
  mutate(assay = recode(assay, "euro_iga_new" = "euro_iga"))

count(all_data_mod, result_og, result) %>% print(n = 50)
count(all_data_mod, group, cohort) %>% print(n = 50)

# Shouldn't be any missing groups
all_data_mod %>%
  filter(is.na(group))

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
    symptom_onset_og = symptom_onset_days,
    symptom_onset_days = suppressWarnings(fix_days_onset(symptom_onset_days))
  )

count(old_onset_fixed, symptom_onset_og, symptom_onset_days) %>%
  print(n = 100)

all_data_new_onset <- bind_rows(
  new_onset, select(old_onset_fixed, -symptom_onset_og)
) %>%
  # Can't use covids with missing onset
  filter(!(group == "covid" & is.na(symptom_onset_days)))

# Create onset categories
all_data_onset_cats <- all_data_new_onset %>%
  mutate(symptom_onset_cat = categorise_onset(symptom_onset_days, group))

count(all_data_onset_cats, symptom_onset_days, symptom_onset_cat) %>%
  print(n = 50)

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

# Look for multiple measurement from the same individual for the same assay
length(unique(all_data_covid_firsts$id)) # Total unique individuals

# Number of individuals that provided more than 1 measurement per assay
all_data_covid_firsts %>%
  count(id, assay, symptom_onset_cat, group, name = "n_samples") %>%
  filter(n_samples > 1) %>%
  group_by(id, n_samples, group) %>%
  summarise(assays = paste(assay, collapse = " "), .groups = "drop") %>%
  # pull(id) -> bad_ids
  print(n = 50)

# walk(bad_ids, ~ all_data_covid_firsts %>% filter(id == .x) %>% print())

# At this point just keep the first observation for everyone

all_data_one_ind <- all_data_covid_firsts %>%
  group_by(id, assay, symptom_onset_cat, group) %>%
  filter(row_number() == 1) %>%
  ungroup()

# Shouldn't be any bad ids
all_data_one_ind %>%
  count(id, assay, symptom_onset_cat, group) %>%
  filter(n > 1)

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
  select(
    id, group, symptom_onset_cat, symptom_onset_days,
    assay, measurement, result
  )

save_data(all_data_final, "data")

# MN data =====================================================================

wantai_mn <- read_raw("wantai-mn", range = "A2:R117") %>%
  select(
    id,
    mn = result, wantai_tot = total_result, wantai_igm = igm_result,
    days_onset = `onset days`, group_og = cohort,
  ) %>%
  filter(!is.na(mn)) %>%
  pivot_longer(
    c(wantai_tot, wantai_igm),
    names_to = "assay", values_to = "result"
  )

euro_svnt_mn <- read_raw("euro-svnt-mn", range = "A2:U112") %>%
  select(
    id,
    mn = mn_result, euro_iga = iga_res, euro_igg = igg_result,
    euro_ncp = ncp_result, `svnt-20` = svnt_result,
    svnt = svnt_result_2, svnt_inh,
    days_onset = `onset days`, group_og = Group,
  ) %>%
  filter(!is.na(mn)) %>%
  mutate(
    svnt = if_else(is.na(svnt), `svnt-20`, svnt),
    `svnt-25` = if_else(svnt_inh < 25, "neg", "pos")
  ) %>%
  select(-svnt_inh) %>%
  pivot_longer(
    c(contains("euro"), contains("svnt")),
    names_to = "assay", values_to = "result"
  )

mn_all <- bind_rows(wantai_mn, euro_svnt_mn) %>%
  filter(!is.na(result)) %>%
  mutate(
    group_fixed = fix_group(group_og),
    result = fix_result(result),
    days_onset_fixed = fix_days_onset(days_onset),
    symptom_onset_cat = categorise_onset(days_onset_fixed, group_fixed),
  )

# Check fixes
unique(mn_all$mn)
mn_all %>% count(group_og, group_fixed)
mn_all %>%
  count(days_onset, days_onset_fixed, symptom_onset_cat, group_fixed) %>%
  print(n = 100)

mn_final <- mn_all %>%
  select(id, group = group_fixed, mn, symptom_onset_cat, assay, result)

save_data(mn_final, "mn")

# Assay name table ------------------------------------------------------------

assays <- tribble(
  ~assay, ~short, ~long, ~measure_name,
  "euro_iga", "E-S1-IgA", "Euroimmun S1 IgA", "OD/CO",
  "euro_igg", "E-S1-IgG", "Euroimmun S1 IgG", "OD/CO",
  "euro_ncp", "E-NCP-IgG", "Euroimmun NCP IgG", "OD/CO",
  "svnt", "sVNT", "GenScript sVNT (repeats for 18-22)", "% Inhibition",
  "svnt-20", "sVNT-20", "GenScript sVNT (first result c/o 20)", "% Inhibition",
  "svnt-25", "sVNT-25", "GenScript sVNT (first result c/o 25)", "% Inhibition",
  "wantai_igm", "W-IgM", "Wantai IgM", "OD/CO",
  "wantai_tot", "W-T", "Wantai total Ab", "OD/CO",
)

setdiff(assays$assay, mn_final$assay)
setdiff(mn_final$assay, assays$assay)

setdiff(assays$assay, all_data_final$assay)
setdiff(all_data_final$assay, assays$assay)

save_data(assays, "assay")
