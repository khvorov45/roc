cat("investigate test characteristics")

library(tidyverse)

data_dir <- "data"
roc_dir <- "roc"

# Functions ===================================================================

source(file.path(data_dir, "read_data.R"))

save_plot <- function(plot, name, ...) {
  ggdark::ggsave_dark(
    file.path(roc_dir, paste0(name, ".png")),
    plot,
    units = "cm",
    ...
  )
}

save_data <- function(data, name) {
  write_csv(data, file.path(roc_dir, paste0(name, ".csv")))
}

# Script ======================================================================

data <- read_data("data")

# Sensitivity and specificity -------------------------------------------------

count_results <- function(result) {
  tibble(
    positive = sum(result == "pos"),
    negative = sum(result == "neg"),
    total = length(result)
  )
}
calc_prop_ci <- function(success, total) {
  res <- PropCIs::exactci(success, total, 0.95)
  tibble(
    success = success,
    total = total,
    low = res$conf.int[[1]],
    point = success / total,
    high = res$conf.int[[2]]
  )
}

# Work out the sensitivity for each onset category
sens_onsets <- data %>%
  filter(group == "covid") %>%
  group_by(assay, symptom_onset_cat) %>%
  summarise(count_results(result), .groups = "keep") %>%
  summarise(calc_prop_ci(positive, total), .groups = "drop")
# Averaged sensitivity
sens_av <- sens_onsets %>%
  group_by(assay) %>%
  summarise(
    symptom_onset_cat = "averaged",
    across(c(low, point, high), mean),
    .groups = "drop"
  )
# All of the sensitivity
sens <- bind_rows(sens_onsets, sens_av) %>%
  mutate(
    symptom_onset_cat = factor(
      symptom_onset_cat,
      levels = c(levels(data$symptom_onset_cat), "averaged")
    )
  )

# Work out the specificity for healthy and non-covid
spec <- data %>%
  filter(group != "covid") %>%
  group_by(assay, group) %>%
  summarise(count_results(result), .groups = "keep") %>%
  summarise(calc_prop_ci(negative, total), .groups = "drop")

# Predictive values -----------------------------------------------------------

calc_pop_pred_vals <- function(sens, spec, prev) {

  # Estimated probabilities
  disease_and_positive <- sens * prev
  disease_and_negative <- (1 - sens) * prev
  healthy_and_positive <- (1 - spec) * (1 - prev)
  healthy_and_negative <- spec * (1 - prev)

  # Estimated predicted values
  ppv <- disease_and_positive / (disease_and_positive + healthy_and_positive)
  npv <- healthy_and_negative / (healthy_and_negative + disease_and_negative)

  bind_rows(
    as_tibble(as.list(ppv)) %>% mutate(char = "ppv"),
    as_tibble(as.list(npv)) %>% mutate(char = "npv")
  ) %>%
    mutate(prev = prev)
}

prevs <- c(0.001, 0.005, 0.01, 0.05, 0.1, 0.2)

pred_vals <- sens %>%
  inner_join(
    spec %>%
      filter(group == "healthy") %>%
      select(assay, spec_low = low, spec_point = point, spec_high = high),
    by = "assay"
  ) %>%
  group_by(assay, symptom_onset_cat) %>%
  summarise(
    map_dfr(prevs, ~ calc_pop_pred_vals(
      c("point" = point, "low" = low, "high" = high),
      c("point" = spec_point, "low" = spec_low, "high" = spec_high),
      .x
    )),
    .groups = "drop"
  )

# Tables out of results -------------------------------------------------------

summ_testchar <- function(point, low, high, success, total) {
  f <- function(x) paste0(round(x * 100, 1), "%")
  summary <- glue::glue(
    "{f(point)} [{f(low)} - {f(high)}]"
  )
  if (missing(success)) {
    return(summary)
  }
  if_else(
    is.na(success),
    summary,
    glue::glue("{summary} ({success} / {total})")
  )
}

sens_table <- sens %>%
  mutate(summary = summ_testchar(point, low, high, success, total)) %>%
  select(assay, symptom_onset_cat, summary) %>%
  pivot_wider(names_from = "assay", values_from = "summary")

save_data(sens_table, "assay-comp-sens")

spec_table <- spec %>%
  mutate(summary = summ_testchar(point, low, high, success, total)) %>%
  select(assay, group, summary) %>%
  pivot_wider(names_from = "assay", values_from = "summary")

save_data(spec_table, "assay-comp-spec")

pred_vals_table <- pred_vals %>%
  mutate(summary = summ_testchar(point, low, high)) %>%
  select(assay, symptom_onset_cat, char, prev, summary) %>%
  pivot_wider(names_from = "assay", values_from = "summary")

save_data(pred_vals_table, "assay-comp-predvals")

# Plots of results ------------------------------------------------------------

plot_testchar <- function(data, ylab) {
  data %>%
    ggplot(aes(assay, point)) +
    ggdark::dark_theme_bw(verbose = FALSE) +
    theme(
      strip.background = element_blank(),
      panel.spacing.x = unit(0, "null"),
      axis.text.x = element_text(angle = 30, hjust = 1)
    ) +
    scale_y_continuous(ylab, labels = scales::percent_format(1)) +
    xlab("Assay") +
    geom_pointrange(aes(ymin = low, ymax = high))
}

sens_plot <- plot_testchar(sens, "Sensitivity") +
  facet_wrap(~symptom_onset_cat, nrow = 1)
save_plot(sens_plot, "assay-comp-sens", width = 20, height = 8)

spec_plot <- plot_testchar(spec, "Specificity") +
  facet_wrap(~group, nrow = 1)
save_plot(spec_plot, "assay-comp-spec", width = 12, height = 8)

predvals_plot <- pred_vals %>%
  filter(symptom_onset_cat == "averaged") %>%
  plot_testchar("Estimate") +
  facet_grid(char ~ prev, scales = "free_y")
save_plot(predvals_plot, "assay-comp-predvals", width = 25, height = 15)
