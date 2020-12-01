cat("summarise data")

library(tidyverse)

# Functions ===================================================================

source("data/read_data.R")

create_group_lbl <- function(group, symptom_onset_cat) {
  if_else(
    group == "covid", as.character(symptom_onset_cat), as.character(group)
  ) %>%
    factor(c(levels(symptom_onset_cat), levels(group)))
}

save_plot <- function(plot, name, ...) {
  ggdark::ggsave_dark(
    file.path("data-summary", paste0(name, ".png")),
    units = "cm",
    ...
  )
}

format_percent <- function(p, l, h) {
  f <- function(x) glue::glue("{signif(x * 100, 2)}%")
  glue::glue("{f(p)} ({f(l)}, {f(h)})")
}

summarise_binary <- function(bin_vec) {
  n <- length(bin_vec)
  s <- sum(bin_vec)
  f <- n - s
  p <- s / n
  ci <- PropCIs::exactci(s, n, 0.95)$conf.int
  tibble(
    total = n,
    success = s,
    failure = f,
    prop = p,
    low = ci[[1]],
    high = ci[[2]],
    summary = glue::glue("{format_percent(prop, low, high)} [{s} / {n}]")
  )
}

merge_groups <- function(onset, group) {
  onset_chr <- as.character(onset)
  group_chr <- as.character(group)
  if_else(group_chr == "covid", onset_chr, group_chr) %>%
    factor(levels = c(levels(onset), levels(group)))
}

save_data <- function(data, name) {
  write_csv(
    data,
    file.path("data-summary", paste0(name, ".csv")),
  )
  data
}

# Script ======================================================================

data <- read_data("data") %>%
  inner_join(read_data("assay"), by = "assay") %>%
  filter(!assay %in% c("svnt-20", "svnt-25")) %>%
  mutate(
    group_lbl = create_group_lbl(group, symptom_onset_cat),
    long = recode(
      long,
      "GenScript sVNT (repeats for 18-22)" = "GenScript sVNT (first measure)"
    )
  )

# Boxplots of measurements ----------------------------------------------------

data_boxplot <- data %>%
  filter(!is.na(symptom_onset_cat)) %>%
  mutate(
    min_threshold = case_when(
      startsWith(as.character(assay), "euro") ~ 0.8,
      startsWith(as.character(assay), "wantai") ~ 0.9,
      assay == "svnt" ~ 20,
    ),
    max_threshold = case_when(
      startsWith(as.character(assay), "euro") ~ 1.1,
      startsWith(as.character(assay), "wantai") ~ 1.1,
      assay == "svnt" ~ 21, # Need to make it visible
    )
  )

one_boxplot <- function(data, key) {
  data %>%
    ggplot(aes(group_lbl, measurement)) +
    theme_bw() +
    facet_wrap(~long, scales = "free_y") +
    theme(
      strip.background = element_blank(),
      axis.text.x = element_text(angle = 30, hjust = 1),
      axis.title.x = element_blank()
    ) +
    scale_x_discrete(
      expand = expansion(0), labels = as_labeller(tools::toTitleCase)
    ) +
    scale_y_continuous(key$measure_name) +
    geom_jitter(
      alpha = 0.4,
      shape = 18,
      position = position_jitter(0.1, 0, 1)
    ) +
    geom_boxplot(
      col = "blue",
      fill = NA,
      outlier.color = NA,
      outlier.fill = NA
    ) +
    geom_rect(
      aes(xmin = 0.5, xmax = 5.5, ymin = min_threshold, ymax = max_threshold),
      alpha = 0.01
    )
}
boxplots <- data_boxplot %>%
  group_by(assay, measure_name) %>%
  group_map(one_boxplot) %>%
  ggpubr::ggarrange(plotlist = ., nrow = 2, ncol = 3)

save_plot(boxplots, "boxplots", width = 20, height = 15)

# Heatmap of results ----------------------------------------------------------

data_heat <- data %>%
  group_by(id, group, symptom_onset_cat) %>%
  mutate(
    id_deid = paste(cur_group_id()),
    discordant = length(unique(result)) > 1, total = length(result)
  ) %>%
  ungroup() %>%
  filter(discordant) %>%
  mutate(
    id = reorder(id, total),
    id_deid = reorder(id_deid, total),
    result = recode(result, "pos" = "Positive", "neg" = "Negative")
  )

counts <- data_heat %>%
  select(id, group_lbl, measurement, assay) %>%
  pivot_wider(names_from = "assay", values_from = "measurement") %>%
  count(group_lbl)

gen_heatmap <- function(data, id_var) {
  data_heat %>%
    ggplot(aes(assay, !!rlang::sym(id_var))) +
    theme_bw() +
    theme(
      legend.position = "bottom",
    ) +
    scale_fill_manual("Test result", values = c("#88b1ff", "#ffb380")) +
    scale_y_discrete("Sample id", expand = expansion()) +
    scale_x_discrete("Assay", expand = expansion()) +
    geom_tile(aes(fill = result)) +
    geom_text(aes(label = signif(measurement, 2)), col = "black") +
    facet_wrap(
      ~group_lbl,
      ncol = 1, scales = "free_y", strip.position = "right"
    )
}

heatmap <- gen_heatmap(data_heat, "id")
heatmap_deid <- gen_heatmap(data_heat, "id_deid")

# Adjust the facet size like a psychopath

save_heatmap <- function(plot, name) {
  png(
    file.path("data-summary", paste0(name, ".png")),
    width = 20, height = 30, units = "cm", res = 150
  )
  gt <- ggplot_gtable(ggplot_build(plot))
  for (i in seq_along(counts$n)) {
    gt$heights[7 + (i - 1) * 4] <- unit(counts$n[[i]], "null")
  }
  grid::grid.draw(gt)
  dev.off()
}

save_heatmap(heatmap, "heatmap-discordant")
save_heatmap(heatmap_deid, "heatmap-discordant-deid")

# Counts ----------------------------------------------------------------------

data_no_svnt_extra <- data %>%
  filter(!assay %in% c("svnt-20", "svnt-25"))

# Assay samples vs unique individuals
s <- function(data) {
  data %>%
    filter(!is.na(symptom_onset_cat)) %>%
    summarise(
      n_indiv = length(unique(id)),
      n_samples = n(),
      summary = glue::glue("{n_samples} ({n_indiv})"),
      .groups = "drop"
    ) %>%
    select(-n_indiv, -n_samples)
}
assay_counts_onsets <- data_no_svnt_extra %>%
  group_by(assay, symptom_onset_cat) %>%
  s() %>%
  rename(subset = symptom_onset_cat)
assay_counts_group <- data_no_svnt_extra %>%
  group_by(assay, group) %>%
  s() %>%
  rename(subset = group)
assay_counts_overall <- data_no_svnt_extra %>%
  group_by(assay) %>%
  s() %>%
  mutate(subset = "combined")

assay_counts <- bind_rows(
  list(assay_counts_onsets, assay_counts_group, assay_counts_overall)
) %>%
  pivot_wider(names_from = "assay", values_from = "summary") %>%
  mutate(subset = factor(
    subset,
    c(
      "<7", "7-14", ">14", "covid",
      "population", "cross-reactive", "no infection",
      "combined"
    ),
    c(
      "<7", "7-14", ">14", "covid total",
      "population", "cross-reactive", "no covid total",
      "overall total"
    )
  )) %>%
  arrange(subset)

save_data(assay_counts, "assay-counts")

# MN validation ---------------------------------------------------------------

mn <- read_data("mn") %>% inner_join(read_data("assay"), by = "assay")

mn_summ_indiv_groups <- mn %>%
  group_by(assay, short, long, symptom_onset_cat, group) %>%
  summarise(summarise_binary(result == mn), .groups = "drop")

mn_summ_onset_averaged <- mn_summ_indiv_groups %>%
  filter(group == "covid") %>%
  group_by(group, assay, short, long) %>%
  summarise(
    symptom_onset_cat = factor("onset-averaged"),
    across(c(prop, low, high), mean),
    summary = format_percent(prop, low, high),
    .groups = "drop"
  )

mn_summ <- bind_rows(mn_summ_indiv_groups, mn_summ_onset_averaged) %>%
  mutate(group_lbl = merge_groups(symptom_onset_cat, group))

mn_summ %>%
  select(assay, group_lbl, summary) %>%
  pivot_wider(names_from = "assay", values_from = "summary") %>%
  arrange(group_lbl) %>%
  save_data("mn-agreement")

mn_summ_plot <- mn_summ %>%
  ggplot(aes(short, group_lbl)) +
  ggdark::dark_theme_bw(verbose = FALSE) +
  theme(
    panel.border = element_blank(),
    axis.ticks = element_blank()
  ) +
  scale_fill_viridis_c("Agreement", labels = scales::percent_format()) +
  scale_x_discrete("Assay") +
  scale_y_discrete("Group", labels = as_labeller(tools::toTitleCase)) +
  geom_tile(aes(fill = prop)) +
  geom_text(aes(label = paste0(round(prop * 100), "%")))

save_plot(mn_summ_plot, "mn-agreement", width = 18, height = 10)
