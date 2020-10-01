cat("investigate test characteristics")

library(tidyverse)
library(furrr)

plan(multiprocess)

data_dir <- "data"
roc_dir <- "roc"

# Functions ===================================================================

source(file.path(data_dir, "read_data.R"))

source(file.path(data_dir, "calc_result_one_threshold.R"))

count_results <- function(data) {
  data %>%
    group_by(true_covid) %>%
    summarise(
      positive = sum(result == "pos"),
      negative = sum(result == "neg"),
      total = n(),
      .groups = "drop"
    )
}

count_true_outcome <- function(data) {
  data %>%
    group_by(result) %>%
    summarise(
      positive = sum(true_covid),
      negative = sum(!true_covid),
      total = n(),
      .groups = "drop"
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

calc_test_char <- function(data) {
  data <- count_results(data)
  sens <- data %>%
    summarise(calc_prop_ci(positive[true_covid], total[true_covid])) %>%
    mutate(char = "Sensitivity")
  spec <- data %>%
    summarise(calc_prop_ci(negative[!true_covid], total[!true_covid])) %>%
    mutate(char = "Specificity")
  bind_rows(sens, spec)
}

one_threshold <- function(min_euro,
                          min_wantai,
                          min_svnt,
                          data) {
  data %>%
    calc_result_one_threshold(min_euro, min_wantai, min_svnt) %>%
    group_by(assay) %>%
    group_modify(~ calc_test_char(.x)) %>%
    ungroup() %>%
    mutate(
      threshold = case_when(
        startsWith(assay, "euro") ~ min_euro,
        startsWith(assay, "wantai") ~ min_wantai,
        assay == "svnt" ~ min_svnt,
      )
    )
}

filter_one_onset <- function(onset, thresholds, data) {
  data %>%
    filter((true_covid & symptom_onset_cat == onset) | !true_covid) %>%
    pmap_dfr(thresholds, one_threshold, .) %>%
    mutate(onset = onset)
}

gen_beta_samples <- function(n_samples, success, total) {
  rbeta(n_samples, success + 1, total - success + 1)
}

quantile_to_df <- function(samples) {
  qs <- quantile(samples, c(0.025, 0.975))
  tibble(low = qs[[1]], high = qs[[2]])
}

calc_pop_pred_vals <- function(results, n_samples = 1e5, prev = 0.1) {
  # Sensitivity
  sens <- with(results, point[char == "Sensitivity"])
  sens_low <- with(results, low[char == "Sensitivity"])
  sens_high <- with(results, high[char == "Sensitivity"])
  # sens_samples <- with(
  #   results,
  #   gen_beta_samples(
  #     n_samples, success[char == "Sensitivity"], total[char == "Sensitivity"]
  #   )
  # )
  # Specificity
  spec <- with(results, point[char == "Specificity"])
  spec_low <- with(results, low[char == "Specificity"])
  spec_high <- with(results, high[char == "Specificity"])
  # spec_samples <- with(
  #   results,
  #   gen_beta_samples(
  #     n_samples, success[char == "Specificity"], total[char == "Specificity"]
  #   )
  # )
  # Estimated probabilities
  disease_and_positive <- sens * prev
  disease_and_positive_low <- sens_low * prev
  disease_and_positive_high <- sens_high * prev
  # disease_and_positive_samples <- sens_samples * prev
  disease_and_negative <- (1 - sens) * prev
  disease_and_negative_low <- (1 - sens_low) * prev
  disease_and_negative_high <- (1 - sens_high) * prev
  # disease_and_negative_samples <- (1 - sens_samples) * prev
  healthy_and_positive <- (1 - spec) * (1 - prev)
  healthy_and_positive_low <- (1 - spec_low) * (1 - prev)
  healthy_and_positive_high <- (1 - spec_high) * (1 - prev)
  # healthy_and_positive_samples <- (1 - spec_samples) * (1 - prev)
  healthy_and_negative <- spec * (1 - prev)
  healthy_and_negative_low <- spec_low * (1 - prev)
  healthy_and_negative_high <- spec_high * (1 - prev)
  # healthy_and_negative_samples <- spec_samples * (1 - prev)
  # Estimated predicted values
  ppv <- disease_and_positive / (disease_and_positive + healthy_and_positive)
  ppv_low <- disease_and_positive_low /
    (disease_and_positive_low + healthy_and_positive_low)
  ppv_high <- disease_and_positive_high /
    (disease_and_positive_high + healthy_and_positive_high)
  # ppv_samples <- disease_and_positive_samples /
  #  (disease_and_positive_samples + healthy_and_positive_samples)
  npv <- healthy_and_negative / (healthy_and_negative + disease_and_negative)
  npv_low <- healthy_and_negative_low /
    (healthy_and_negative_low + disease_and_negative_low)
  npv_high <- healthy_and_negative_high /
    (healthy_and_negative_high + disease_and_negative_high)
  # npv_samples <- healthy_and_negative_samples /
  #  (healthy_and_negative_samples + disease_and_negative_samples)
  # Combine into a df
  # ppv_est <- quantile_to_df(ppv_samples) %>% mutate(point = ppv, char = "PPV")
  # npv_est <- quantile_to_df(npv_samples) %>% mutate(point = npv, char = "NPV")
  ppv_est <- tibble(low = ppv_low, point = ppv, high = ppv_high, char = "PPV")
  npv_est <- tibble(low = npv_low, point = npv, high = npv_high, char = "NPV")
  bind_rows(ppv_est, npv_est)
}

one_prevalence <- function(prev, results, n_samples = 1e5) {
  results %>%
    group_by(assay, threshold, onset) %>%
    group_modify(~ calc_pop_pred_vals(.x, n_samples, prev)) %>%
    mutate(prev = prev) %>%
    ungroup()
}

filter_std_thresholds <- function(all_results,
                                  euro = 0.8, wantai = 0.9, svnt = 20) {
  all_results %>%
    filter(
      (startsWith(assay, "euro") & threshold == euro)
      | (startsWith(assay, "wantai") & threshold == wantai)
      | (assay == "svnt" & threshold == svnt)
    )
}

plot_calcs <- function(data, assay = "") {
  data_mod <- data %>%
    mutate(color = if_else(
      (startsWith(assay, "euro") & threshold == 0.8)
      | (startsWith(assay, "wantai") & threshold == 0.9)
      | (assay == "svnt" & threshold == 20), "red", "black",
    ))
  plot <- data_mod %>%
    ggplot(aes(threshold, point)) +
    ggdark::dark_theme_bw(verbose = FALSE) +
    theme(
      strip.background = element_blank(),
      panel.spacing = unit(0, "null"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.minor.x = element_blank()
    ) +
    facet_grid(char ~ onset, scales = "free_y") +
    scale_x_continuous(
      "Threshold",
      labels = scales::number_format(0.01)
    ) +
    scale_y_continuous("Estimate", labels = scales::percent_format(1)) +
    scale_color_identity() +
    geom_ribbon(aes(ymin = low, ymax = high), alpha = 0.4) +
    geom_line() +
    geom_point(aes(col = color))
  attr(plot, "assay") <- assay
  plot
}

plot_predvals <- function(data, assay = "") {
  plot <- data %>%
    ggplot(aes(prev, point)) +
    ggdark::dark_theme_bw(verbose = FALSE) +
    theme(
      strip.background = element_blank(),
      panel.spacing = unit(0, "null"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.minor.x = element_blank()
    ) +
    facet_grid(char ~ onset, scales = "free_y") +
    scale_x_log10(
      "Prevalence",
      breaks = unique(data$prev),
      labels = scales::percent_format(0.01)
    ) +
    scale_y_continuous("Estimate", labels = scales::percent_format(0.1)) +
    geom_pointrange(aes(ymin = low, ymax = high))
  attr(plot, "assay") <- assay
  plot
}

plot_assay_comp <- function(data) {
  data %>%
    ggplot(aes(assay, point)) +
    ggdark::dark_theme_bw(verbose = FALSE) +
    theme(
      strip.background = element_blank(),
      axis.text.x = element_text(angle = 30, hjust = 1),
      panel.spacing = unit(0, "null")
    ) +
    scale_y_continuous("Estimate", labels = scales::percent_format(1)) +
    geom_pointrange(aes(ymin = low, ymax = high))
}

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

n_to_test <- 15

thresholds <- tibble(
  min_euro = c(
    seq(0, 0.5, length.out = n_to_test),
    seq(1, 10, length.out = n_to_test),
    0.8
  ),
  min_wantai = c(
    seq(-0.05, 0.1, length.out = n_to_test),
    seq(0.2, 10, length.out = n_to_test),
    0.9
  ),
  min_svnt = c(
    seq(0, 15, length.out = n_to_test),
    seq(20, 100, length.out = n_to_test),
    20
  )
)

onsets <- na.omit(unique(data$symptom_onset_cat))
onsets <- onsets[onsets != "no infection"]

indiv_onsets <- future_map_dfr(onsets, filter_one_onset, thresholds, data)

# Not filtering by onset at all doesn't make sense because we want to compare
# assays and sample compositions are different - wantai is mostly late infection
# but svnt is mostly early infection, for example
averaged_onsets <- indiv_onsets %>%
  group_by(assay, char, threshold) %>%
  summarise(
    onset = "Averaged",
    low = mean(low),
    point = mean(point),
    high = mean(high),
    .groups = "drop"
  )

all_results <- bind_rows(indiv_onsets, averaged_onsets) %>%
  mutate(onset = factor(onset, levels = c(levels(onsets), "Averaged"))) %>%
  distinct(assay, onset, threshold, char, .keep_all = TRUE)

save_data(all_results, "roc")

# Reshape to only see the test characteristics
test_chars_only <- all_results %>%
  pivot_longer(c(low, point, high), names_to = "bound", values_to = "est") %>%
  select(-success, -total) %>%
  pivot_wider(names_from = "char", values_from = "est")

# Calculate ROC AUC

calc_roc_auc <- function(sens, spec) {
  DescTools::AUC(c(0, 1 - spec, 1), c(0, sens, 1))
}

aucs <- all_results %>%
  group_by(assay, onset) %>%
  summarise(
    point = calc_roc_auc(
      point[char == "Sensitivity"], point[char == "Specificity"]
    ),
    low = calc_roc_auc(
      low[char == "Sensitivity"], low[char == "Specificity"]
    ),
    high = calc_roc_auc(
      high[char == "Sensitivity"], high[char == "Specificity"]
    ),
    .groups = "drop"
  )

save_data(aucs, "aucs")

# Calculate predictive values at a range of prevalences
prevs <- c(0.001, 0.005, 0.01, 0.05, 0.1, 0.2)
pop_pred_vals <- future_map_dfr(prevs, one_prevalence, all_results)

save_data(pop_pred_vals, "pred-vals-pop")

# Plot all results
plots <- all_results %>%
  group_by(assay) %>%
  group_map(~ plot_calcs(.x, paste(.y$assay)))

walk(plots, ~ save_plot(.x, attr(.x, "assay"), width = 25, height = 15))

# Plot all results as a roc curve
test_chars_only_plot_mod <- test_chars_only %>%
  # Find thresholds to label
  mutate(
    label_threshold = (startsWith(assay, "euro") & threshold == 0.8)
    | (startsWith(assay, "wantai") & threshold == 0.9)
    | (assay == "svnt" & threshold == 20)
  ) %>%
  # Need to add the extremes
  group_by(assay, onset, bound) %>%
  group_modify(
    ~ bind_rows(
      .x,
      tribble(
        ~threshold, ~Sensitivity, ~Specificity,
        Inf, 0, 1,
        -Inf, 1, 0
      )
    )
  ) %>%
  arrange(threshold)

roc_plot <- test_chars_only_plot_mod %>%
  ggplot(aes(1 - Specificity, Sensitivity, group = bound)) +
  ggdark::dark_theme_bw(verbose = FALSE) +
  theme(
    strip.background = element_blank(),
    axis.text.x = element_text(angle = 30, hjust = 1)
  ) +
  facet_grid(assay ~ onset) +
  scale_x_continuous(labels = scales::percent_format(1)) +
  geom_abline(slope = 1, intercept = 0, lty = "11") +
  geom_path() +
  geom_point(data = filter(test_chars_only_plot_mod, bound == "point")) +
  geom_point(
    data = filter(test_chars_only_plot_mod, label_threshold, bound == "point"),
    col = "red"
  )

save_plot(
  roc_plot, "roc",
  width = 15, height = 15
)

# Filter down to the standard thresholds
std_threshold_results <- all_results %>% filter_std_thresholds()
std_threshold_predvals <- pop_pred_vals %>% filter_std_thresholds()

# Plot assay comparison
assay_comp_plot <- std_threshold_results %>%
  filter(char %in% c("Sensitivity", "Specificity")) %>%
  plot_assay_comp() +
  xlab("Assay at standard threshold") +
  facet_grid(char ~ onset, scales = "free_y")

save_plot(assay_comp_plot, "assay-comp", width = 20, height = 15)

# Plots assay comparision for ROC AUC
aucs_plot <- aucs %>%
  ggplot(aes(assay, point)) +
  ggdark::dark_theme_bw(verbose = FALSE) +
  theme(
    strip.background = element_blank(),
    axis.text.x = element_text(angle = 30, hjust = 1)
  ) +
  scale_x_discrete("Assay") +
  scale_y_continuous("ROC AUC") +
  facet_wrap(~onset) +
  geom_hline(yintercept = 0.5, lty = "11") +
  geom_hline(yintercept = 1, lty = "13") +
  geom_pointrange(aes(ymin = low, ymax = high))

save_plot(aucs_plot, "aucs", width = 15, height = 15)

# Plot assay comparison for predictive values
assay_comp_predvals <- std_threshold_predvals %>%
  filter(onset == "Averaged") %>%
  mutate(
    prev_lab = factor(
      prev,
      levels = prevs,
      labels = glue::glue("Prevalence {prevs * 100}%")
    )
  ) %>%
  plot_assay_comp() +
  xlab("Assay at standard threshold for unknown (averaged) onset") +
  facet_grid(char ~ prev_lab, scales = "free_y")

save_plot(
  assay_comp_predvals, "assay-comp-predvals",
  width = length(prevs) * 5, height = 15
)

# Plot one assay's predvals at different prevalences
assay_predvals <- std_threshold_predvals %>%
  group_by(assay) %>%
  group_map(~ plot_predvals(.x, .y$assay))

walk(
  assay_predvals,
  ~ save_plot(.x, paste0(attr(.x, "assay"), "-predvals"),
    width = 25, height = 15
  )
)

# Table of results at standard thresholds

f <- function(x) paste0(round(x * 100, 1), "%")
std_threshold_table <- std_threshold_results %>%
  mutate(
    summary = glue::glue(
      "{f(point)} [{f(low)} - {f(high)}]"
    ),
    summary = if_else(
      is.na(success),
      summary,
      glue::glue("{summary} ({success} / {total})")
    )
  ) %>%
  select(-success, -total, -low, -point, -high, -threshold) %>%
  pivot_wider(names_from = "assay", values_from = "summary") %>%
  select(
    onset, char, euro_ncp, euro_igg, euro_iga, svnt, wantai_tot, wantai_igm
  )

save_data(std_threshold_table, "assay-comp")

std_threshold_table_predvals <- std_threshold_predvals %>%
  mutate(
    summary = glue::glue(
      "{f(point)} [{f(low)} - {f(high)}]"
    )
  ) %>%
  select(-threshold, -low, -point, -high) %>%
  pivot_wider(names_from = "assay", values_from = "summary") %>%
  select(
    prev, onset, char,
    euro_ncp, euro_igg, euro_iga, svnt, wantai_tot, wantai_igm
  )

save_data(std_threshold_table_predvals, "assay-comp-predvals")

# AUCs table
f <- function(x) round(x, 3)
aucs_table <- aucs %>%
  mutate(
    summary = glue::glue("{f(point)} [{f(low)}, {f(high)}]")
  ) %>%
  select(-point, -low, -high) %>%
  pivot_wider(names_from = "assay", values_from = "summary") %>%
  select(
    onset,
    euro_ncp, euro_igg, euro_iga, svnt, wantai_tot, wantai_igm
  )

save_data(aucs_table, "assay-comp-auc")
