cat("investigate test characteristics")

library(tidyverse)

data_dir <- "data"
roc_dir <- "roc"

# Functions ===================================================================

source(file.path(data_dir, "read_data.R"))

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

calc_sens <- function(result) {
  count_results(result) %>%
    summarise(calc_prop_ci(positive, total))
}

calc_spec <- function(result) {
  count_results(result) %>%
    summarise(calc_prop_ci(negative, total))
}

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

plot_testchar <- function(data, ylab, xvar = "assay",
                          xlab = "Assay") {
  data %>%
    ggplot(aes(
      !!rlang::sym(xvar), point,
      ymin = low, ymax = high
    )) +
    ggdark::dark_theme_bw(verbose = FALSE) +
    theme(
      strip.background = element_blank(),
      panel.spacing.x = unit(0, "null"),
      axis.text.x = element_text(angle = 30, hjust = 1)
    ) +
    scale_color_discrete("Assay") +
    xlab(xlab) +
    scale_y_continuous(ylab, labels = scales::percent_format(1))
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

lighten_plot <- function(plot) {
  plot$theme <- ggdark::lighten_theme(plot$theme)
  plot
}

# Script ======================================================================

data <- read_data("data") %>% inner_join(read_data("assay"), by = "assay")

# Prevalences to evaluate predictive values at
prevs <- c(0.001, 0.0011, 0.0032, 0.005, 0.01, 0.05, 0.1, 0.2)

# Sensitivity and specificity -------------------------------------------------

# Work out the sensitivity for each onset category
sens_onsets <- data %>%
  filter(group == "covid") %>%
  group_by(assay, short, long, symptom_onset_cat) %>%
  summarise(calc_sens(result), .groups = "drop")
# Averaged sensitivity
sens_av <- sens_onsets %>%
  group_by(assay, short, long) %>%
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
  group_by(assay, short, long, group) %>%
  summarise(calc_spec(result), .groups = "drop")

# Predictive values -----------------------------------------------------------

pred_vals <- sens %>%
  inner_join(
    spec %>%
      filter(group == "population") %>%
      select(assay, spec_low = low, spec_point = point, spec_high = high),
    by = "assay"
  ) %>%
  group_by(assay, short, long, symptom_onset_cat) %>%
  summarise(
    map_dfr(prevs, ~ calc_pop_pred_vals(
      c("point" = point, "low" = low, "high" = high),
      c("point" = spec_point, "low" = spec_low, "high" = spec_high),
      .x
    )),
    .groups = "drop"
  )

# Tables out of results -------------------------------------------------------

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

sens_plot <- plot_testchar(sens, "Sensitivity", "short") +
  facet_wrap(
    ~symptom_onset_cat,
    nrow = 1, labeller = as_labeller(tools::toTitleCase)
  ) +
  geom_pointrange(aes(col = long))
save_plot(sens_plot, "assay-comp-sens", width = 28, height = 8)

spec_plot <- plot_testchar(spec, "Specificity", "short") +
  facet_wrap(~group, nrow = 1, labeller = as_labeller(tools::toTitleCase)) +
  geom_pointrange(aes(col = long))
save_plot(spec_plot, "assay-comp-spec", width = 18, height = 8)

# Sens/spec together
together_mod <- list(
  coord_cartesian(ylim = c(0, 1)),
  theme(legend.position = "none")
)
common_legend <- ggpubr::get_legend(lighten_plot(sens_plot))
sens_spec_together <- gridExtra::grid.arrange(
  lighten_plot(sens_plot) + together_mod,
  common_legend,
  lighten_plot(spec_plot) + together_mod,
  layout_matrix = rbind(
    rep(1, 13),
    c(rep(2, 6), rep(3, 7)) # Width-specific alignment
  )
)
save_plot(sens_spec_together, "sens-spec-together", width = 25, height = 18)

predvals_plot <- pred_vals %>%
  plot_testchar("Estimate", "short") +
  facet_grid(
    symptom_onset_cat + char ~ prev,
    switch = "y",
    labeller = labeller(
      char = toupper,
      prev = function(prev) glue::glue("Prevalence {as.numeric(prev) * 100}%"),
      symptom_onset_cat = tools::toTitleCase
    )
  ) +
  theme(
    strip.placement = "outside", axis.title.y = element_blank(),
    legend.position = "bottom"
  ) +
  geom_pointrange(aes(col = long))
save_plot(
  predvals_plot, "assay-comp-predvals",
  width = length(prevs) * 5, height = 25
)

# Try various thresholds for svnt ---------------------------------------------

data_svnt <- data %>%
  filter(assay == "svnt")

thresholds <- 15:30

calc_result <- function(threshold, measurement) {
  if_else(measurement < threshold, "neg", "pos")
}
sens_theshold <- function(threshold, data_svnt) {
  sens_onsets <- data_svnt %>%
    mutate(result = calc_result(threshold, measurement)) %>%
    filter(group == "covid") %>%
    group_by(symptom_onset_cat) %>%
    summarise(calc_sens(result), .groups = "drop")
  sens_av <- sens_onsets %>%
    summarise(across(c(low, point, high), mean)) %>%
    mutate(symptom_onset_cat = "averaged", )
  bind_rows(sens_onsets, sens_av) %>%
    mutate(
      threshold = threshold,
      symptom_onset_cat = factor(
        symptom_onset_cat,
        levels = c(levels(data$symptom_onset_cat), "averaged")
      )
    )
}
spec_threshold <- function(threshold, data_svnt) {
  data_svnt %>%
    mutate(result = calc_result(threshold, measurement)) %>%
    filter(group != "covid") %>%
    group_by(group) %>%
    summarise(calc_spec(result), .groups = "drop") %>%
    mutate(threshold = threshold)
}

sens_res <- map_dfr(thresholds, sens_theshold, data_svnt)

spec_res <- map_dfr(thresholds, spec_threshold, data_svnt)

pred_vals_res <- sens_res %>%
  inner_join(
    spec_res %>%
      filter(group == "population") %>%
      select(threshold, spec_low = low, spec_point = point, spec_high = high),
    by = "threshold"
  ) %>%
  group_by(threshold, symptom_onset_cat) %>%
  summarise(
    map_dfr(prevs, ~ calc_pop_pred_vals(
      c("point" = point, "low" = low, "high" = high),
      c("point" = spec_point, "low" = spec_low, "high" = spec_high),
      .x
    )),
    .groups = "drop"
  )

# Tables of thresholds --------------------------------------------------------

sens_res_table <- sens_res %>%
  mutate(summary = summ_testchar(point, low, high, success, total)) %>%
  select(symptom_onset_cat, threshold, summary) %>%
  pivot_wider(names_from = "symptom_onset_cat", values_from = "summary")

save_data(sens_res_table, "svnt-sens")

spec_res_table <- spec_res %>%
  group_by(group, threshold) %>%
  mutate(summary = summ_testchar(point, low, high, success, total)) %>%
  select(group, threshold, summary) %>%
  pivot_wider(names_from = "group", values_from = "summary")

save_data(spec_res_table, "svnt-spec")

predvals_res_table <- pred_vals_res %>%
  mutate(summary = summ_testchar(point, low, high)) %>%
  select(char, prev, threshold, symptom_onset_cat, summary) %>%
  pivot_wider(names_from = "symptom_onset_cat", values_from = "summary")

save_data(predvals_res_table, "svnt-predvals")

# Plot the various thresholds -------------------------------------------------

special_pointranges <- function(data) {
  list(
    geom_pointrange(
      data = data %>% filter(threshold == 20),
      col = "darkorange"
    ),
    geom_pointrange(
      data = data %>% filter(threshold == 25),
      col = "blue"
    )
  )
}

sens_res_plot <- plot_testchar(
  sens_res, "Sensitivity", "threshold", "Threshold"
) +
  facet_wrap(~symptom_onset_cat, nrow = 1) +
  scale_x_continuous(breaks = thresholds) +
  theme(panel.grid.minor.x = element_blank()) +
  geom_pointrange() +
  special_pointranges(sens_res)

save_plot(sens_res_plot, "svnt-sens", width = 26, height = 10)

spec_res_plot <- plot_testchar(
  spec_res, "Specificity", "threshold", "Treshold"
) +
  facet_wrap(~group, nrow = 1) +
  scale_x_continuous(breaks = thresholds) +
  geom_pointrange() +
  special_pointranges(spec_res)

save_plot(spec_res_plot, "svnt-spec", width = 15, height = 10)

pred_vals_res_plot <- pred_vals_res %>%
  plot_testchar("Estimate", "threshold", "Threshold") +
  facet_grid(
    symptom_onset_cat + char ~ prev,
    scales = "free_y", switch = "y",
    labeller = labeller(
      char = toupper,
      prev = function(prev) glue::glue("Prevalence {as.numeric(prev) * 100}%"),
      symptom_onset_cat = tools::toTitleCase
    )
  ) +
  theme(strip.placement = "outside", axis.title.y = element_blank()) +
  geom_pointrange() +
  special_pointranges(pred_vals_res)

save_plot(
  pred_vals_res_plot, "svnt-predvals",
  width = length(prevs) * 5, height = 25
)
