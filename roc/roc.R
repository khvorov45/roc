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

calc_result_one_threshold <- function(data,
                                      min_euro = 0.8,
                                      min_wantai = 0.9,
                                      min_svnt = 0.2) {
  data %>%
    mutate(
      result = case_when(
        startsWith(assay, "euro") ~
        determine_result(min_euro, min_euro + 0.3)(measurement),
        startsWith(assay, "wantai") ~
        determine_result(min_wantai, min_wantai + 0.2)(measurement),
        assay == "svnt" ~ determine_result(min_svnt, min_svnt)(measurement),
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

plot_calcs <- function(data, name = "") {
  plot <- data %>%
    ggplot(aes(threshold, point)) +
    ggdark::dark_theme_bw(verbose = FALSE) +
    theme(
      strip.background = element_blank()
    ) +
    facet_wrap(~char, scales = "free") +
    scale_x_continuous("Threshold") +
    scale_y_continuous("Estimate", labels = scales::percent_format(1)) +
    geom_pointrange(aes(ymin = low, ymax = high))
  attr(plot, "name") <- name
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

roc <- pmap_dfr(thresholds, calc_both_one_threshold, data)

save_data(roc, "roc")

plots <- roc %>%
  group_by(assay) %>%
  group_map(~ plot_calcs(.x, .y$assay))

walk(plots, ~ save_plot(.x, attr(.x, "name"), width = 15, height = 15))

assay_comp_plot <- roc %>%
  filter(
    (startsWith(assay, "euro") & threshold == 0.8)
    | (startsWith(assay, "wantai") & threshold == 0.9)
    | (assay == "svnt" & threshold == 20)
  ) %>%
  ggplot(aes(assay, point)) +
  ggdark::dark_theme_bw(verbose = FALSE) +
    theme(
      strip.background = element_blank(),
      axis.text.x = element_text(angle = 30, hjust = 1)
    ) +
  facet_wrap(~char, scales = "free_y") +
  scale_x_discrete("Assay at standard threshold") +
  scale_y_continuous("Estimate", labels = scales::percent_format(1)) +
  geom_pointrange(aes(ymin = low, ymax = high))

save_plot(assay_comp_plot, "assay-comp", width = 15, height = 15)
