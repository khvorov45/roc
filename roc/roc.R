cat("investigate test characteristics")

library(tidyverse)
library(furrr)

plan(multiprocess)

data_dir <- "data"
roc_dir <- "roc"

# Functions ===================================================================

source(file.path(data_dir, "read_data.R"))

determine_result <- function(measurement, threshold, equivocal_range = 0) {
  if (equivocal_range > 0) {
    case_when(
      measurement < threshold ~ "neg",
      measurement <= threshold + equivocal_range ~ "eqiv",
      measurement > threshold + equivocal_range ~ "pos"
    )
  } else {
    if_else(
      measurement < threshold, "neg", "pos"
    )
  }
}

calc_result_one_threshold <- function(data,
                                      euro = 0.8,
                                      wantai = 0.9,
                                      svnt = 0.2) {
  data %>%
    mutate(
      result = case_when(
        startsWith(assay, "euro") ~
        determine_result(measurement, euro),
        startsWith(assay, "wantai") ~
        determine_result(measurement, wantai),
        assay == "svnt" ~ determine_result(measurement, svnt),
      )
    )
}

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

calc_pred_vals <- function(data) {
  data <- count_true_outcome(data)
  ppv <- tibble()
  if ("pos" %in% data$result) {
    ppv <- data %>%
      summarise(
        calc_prop_ci(positive[result == "pos"], total[result == "pos"])
      ) %>%
      mutate(char = "PPV")
  }
  npv <- tibble()
  if ("neg" %in% data$result) {
    npv <- data %>%
      summarise(
        calc_prop_ci(negative[result == "neg"], total[result == "neg"])
      ) %>%
      mutate(char = "NPV")
  }
  bind_rows(ppv, npv)
}

calc_both <- function(data) {
  bind_rows(calc_test_char(data), calc_pred_vals(data))
}

calc_both_one_threshold <- function(min_euro,
                                    min_wantai,
                                    min_svnt,
                                    data) {
  data %>%
    calc_result_one_threshold(min_euro, min_wantai, min_svnt) %>%
    group_by(assay) %>%
    group_modify(~ calc_both(.x)) %>%
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
    pmap_dfr(thresholds, calc_both_one_threshold, .) %>%
    mutate(onset = onset)
}

gen_beta_samples <- function(n_samples, success, total) {
  rbeta(n_samples, success + 1, total - success + 1)
}

quantile_to_df <- function(samples) {
  qs <- quantile(samples, c(0.025, 0.5, 0.975))
  tibble(low = qs[[1]], point = qs[[2]], high = qs[[3]])
}

calc_pop_pred_vals <- function(results, n_samples = 1e5, prev = 0.1) {
  sens_samples <- with(
    results,
    gen_beta_samples(
      n_samples, success[char == "Sensitivity"], total[char == "Sensitivity"]
    )
  )
  spec_samples <- with(
    results,
    gen_beta_samples(
      n_samples, success[char == "Specificity"], total[char == "Specificity"]
    )
  )
  disease_and_positive <- sens_samples * prev
  disease_and_negative <- (1 - sens_samples) * prev
  healthy_and_positive <- (1 - spec_samples) * (1 - prev)
  healthy_and_negative <- spec_samples * (1 - prev)
  ppv_samples <-
    disease_and_positive / (disease_and_positive + healthy_and_positive)
  npv_samples <-
    healthy_and_negative / (healthy_and_negative + disease_and_negative)
  ppv_est <- quantile_to_df(ppv_samples) %>% mutate(char = "PPV")
  npv_est <- quantile_to_df(npv_samples) %>% mutate(char = "NPV")
  bind_rows(ppv_est, npv_est)
}

one_prevalence <- function(prev, results, n_samples = 1e5) {
  results %>%
    group_by(assay, threshold, onset) %>%
    group_modify(~ calc_pop_pred_vals(.x, n_samples, prev)) %>%
    mutate(prev = prev)
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
    ggplot(aes(threshold, point, col = color)) +
    ggdark::dark_theme_bw(verbose = FALSE) +
    theme(
      strip.background = element_blank(),
      panel.spacing = unit(0, "null"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.minor.x = element_blank()
    ) +
    facet_grid(onset ~ char) +
    scale_x_continuous(
      "Threshold",
      breaks = unique(data_mod$threshold[data_mod$color != "red"]),
      labels = scales::number_format(0.01)
    ) +
    scale_y_continuous("Estimate", labels = scales::percent_format(1)) +
    scale_color_identity() +
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

n_to_test <- 10

thresholds <- tibble(
  min_euro = c(seq(0.5, 8, length.out = n_to_test), 0.8),
  min_wantai = c(seq(0.5, 7, length.out = n_to_test), 0.9),
  min_svnt = c(seq(18, 50, length.out = n_to_test), 20)
)

onsets <- na.omit(unique(data$symptom_onset_cat))
onsets <- onsets[onsets != "no infection"]

indiv_onsets <- future_map_dfr(onsets, filter_one_onset, thresholds, data)
any_onset <- future_pmap_dfr(thresholds, calc_both_one_threshold, data)

all_results <- bind_rows(indiv_onsets, mutate(any_onset, onset = "Any")) %>%
  mutate(onset = factor(onset, levels = c(levels(onsets), "Any")))

save_data(all_results, "roc")

# Calculate predictive values at a range of prevalences
prevs <- c(0.01, 0.05, 0.1, 0.2)
pop_pred_vals <- future_map_dfr(prevs, one_prevalence, all_results)

save_data(pop_pred_vals, "pred-vals-pop")

# Plot all results
plots <- all_results %>%
  group_by(assay) %>%
  group_map(~ plot_calcs(.x, paste(.y$assay)))

walk(plots, ~ save_plot(.x, attr(.x, "assay"), width = 20, height = 20))

# Filter down to the standard thresholds
std_threshold_results <- all_results %>% filter_std_thresholds()
std_threshold_predvals <- pop_pred_vals %>% filter_std_thresholds()

# Plot assay comparison
assay_comp_plot <- std_threshold_results %>%
  plot_assay_comp() +
  xlab("Assay at any threshold") +
  facet_grid(char ~ onset, scales = "free_y")

save_plot(assay_comp_plot, "assay-comp", width = 20, height = 20)

# Plot assay comparison for predictive values
assay_comp_predvals <- std_threshold_predvals %>%
  filter(onset == "Any") %>%
  plot_assay_comp() +
  xlab("Assay at any threshold for any onset") +
  facet_grid(char ~ prev, scales = "free_y")

save_plot(assay_comp_predvals, "assay-comp-predvals", width = 20, height = 15)

# Table of results at standard thresholds

f <- scales::percent_format(0.1)
std_threshold_table <- std_threshold_results %>%
  mutate(
    summary = glue::glue(
      "{f(point)} [{f(low)} - {f(high)}] ({success} / {total})"
    )
  ) %>%
  select(-success, -total, -low, -point, -high, -threshold) %>%
  pivot_wider(names_from = "assay", values_from = "summary") %>%
  select(
    onset, char, euro_ncp, euro_igg, euro_iga, svnt, wantai_tot, wantai_igm
  )

save_data(std_threshold_table, "assay-comp")
