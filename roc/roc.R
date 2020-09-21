cat("investigate test characteristics")

library(tidyverse)

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
  ppv <- data %>%
    summarise(
      calc_prop_ci(positive[result == "pos"], total[result == "pos"])
    ) %>%
    mutate(char = "PPV")
  npv <- data %>%
    summarise(
      calc_prop_ci(negative[result == "neg"], total[result == "neg"])
    ) %>%
    mutate(char = "NPV")
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

plot_calcs <- function(data, assay = "") {
  plot <- data %>%
    mutate(color = if_else(
      (startsWith(assay, "euro") & threshold == 0.8)
      | (startsWith(assay, "wantai") & threshold == 0.9)
      | (assay == "svnt" & threshold == 20), "red", "black",
    )) %>%
    ggplot(aes(threshold, point, col = color)) +
    ggdark::dark_theme_bw(verbose = FALSE) +
    theme(
      strip.background = element_blank(),
      panel.spacing = unit(0, "null"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.minor.x = element_blank()
    ) +
    facet_grid(onset ~ char) +
    scale_x_continuous("Threshold", breaks = unique(data$threshold)) +
    scale_y_continuous("Estimate", labels = scales::percent_format(1)) +
    scale_color_identity() +
    geom_pointrange(aes(ymin = low, ymax = high))
  attr(plot, "assay") <- assay
  plot
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

n_to_test <- 11

thresholds <- tibble(
  min_euro = seq(0.5, 1.5, length.out = n_to_test),
  min_wantai = seq(0.5, 1.5, length.out = n_to_test),
  min_svnt = seq(18, 28, length.out = n_to_test)
)

onsets <- na.omit(unique(data$symptom_onset_cat))
onsets <- onsets[onsets != "no infection"]

indiv_onsets <- map_dfr(onsets, filter_one_onset, thresholds, data)
any_onset <- pmap_dfr(thresholds, calc_both_one_threshold, data)

all_results <- bind_rows(indiv_onsets, mutate(any_onset, onset = "Any")) %>%
  mutate(onset = factor(onset, levels = c("<=8", "9-14", ">=15", "Any")))

save_data(all_results, "roc")

plots <- all_results %>%
  group_by(assay) %>%
  group_map(~ plot_calcs(.x, paste(.y$assay)))

walk(plots, ~ save_plot(.x, attr(.x, "assay"), width = 20, height = 20))

assay_comp_plot <- all_results %>%
  filter(
    (startsWith(assay, "euro") & threshold == 0.8)
    | (startsWith(assay, "wantai") & threshold == 0.9)
    | (assay == "svnt" & threshold == 20)
  ) %>%
  ggplot(aes(assay, point)) +
  ggdark::dark_theme_bw(verbose = FALSE) +
  theme(
    strip.background = element_blank(),
    axis.text.x = element_text(angle = 30, hjust = 1),
    panel.spacing.y = unit(0, "null")
  ) +
  facet_grid(onset ~ char, scales = "free_y") +
  scale_x_discrete("Assay at standard threshold") +
  scale_y_continuous("Estimate", labels = scales::percent_format(1)) +
  geom_pointrange(aes(ymin = low, ymax = high))

save_plot(assay_comp_plot, "assay-comp", width = 20, height = 20)
