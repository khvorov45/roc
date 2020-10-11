cat("plot data")

library(tidyverse)

data_dir <- "data"
data_plot_dir <- "data-plot"

# Functions ===================================================================

source(file.path(data_dir, "read_data.R"))

source(file.path(data_dir, "calc_result_one_threshold.R"))

save_plot <- function(plot, name, ...) {
  ggdark::ggsave_dark(
    file.path(data_plot_dir, paste0(name, ".png")),
    units = "cm",
    ...
  )
}

# Script ======================================================================

data <- read_data("data")

boxplots <- data %>%
  filter(!is.na(symptom_onset_cat)) %>%
  mutate(
    min_threshold = case_when(
      startsWith(assay, "euro") ~ 0.8,
      startsWith(assay, "wantai") ~ 0.9,
      assay == "svnt" ~ 20,
    ),
    max_threshold = case_when(
      startsWith(assay, "euro") ~ 1.1,
      startsWith(assay, "wantai") ~ 1.1,
      assay == "svnt" ~ 21, # Need to make it visible
    )
  ) %>%
  ggplot(aes(symptom_onset_cat, measurement)) +
  ggdark::dark_theme_bw(verbose = FALSE) +
  facet_wrap(~assay, scales = "free_y") +
  theme(
    strip.background = element_blank()
  ) +
  scale_x_discrete("Symptom onset", expand = expansion(0)) +
  scale_y_continuous("Measure") +
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
    aes(xmin = 0.5, xmax = 4.5, ymin = min_threshold, ymax = max_threshold),
    alpha = 0.01
  )

save_plot(boxplots, "boxplots", width = 20, height = 15)

# Heatmap of results
data_heat <- data %>%
  calc_result_one_threshold() %>%
  group_by(sample_id) %>%
  mutate(discordant = length(unique(result)) > 1, total = length(result)) %>%
  ungroup() %>%
  filter(discordant) %>%
  mutate(
    sample_id = reorder(sample_id, total),
    group_lbl = if_else(true_covid, as.character(symptom_onset_cat), group) %>%
      factor(c("<7", "7-14", ">14", "non-covid", "healthy")),
    result = recode(result, "pos" = "Positive", "neg" = "Negative")
  )

counts <- data_heat %>%
  select(sample_id, group_lbl, measurement, assay) %>%
  pivot_wider(names_from = "assay", values_from = "measurement") %>%
  count(group_lbl)

heatmap <- data_heat %>%
  ggplot(aes(assay, sample_id)) +
  theme_bw() +
  theme(
    legend.position = "bottom",
  ) +
  scale_fill_manual("Test result", values = c("#d9e2f3", "#fbe4d5")) +
  scale_y_discrete("Sample id", expand = expansion()) +
  scale_x_discrete("Assay", expand = expansion()) +
  geom_tile(aes(fill = result)) +
  geom_text(aes(label = signif(measurement, 2)), col = "black") +
  facet_wrap(
    ~group_lbl,
    ncol = 1, scales = "free_y", strip.position = "right"
  )

# Adjust the facet size like a psychopath
gt <- ggplot_gtable(ggplot_build(heatmap))
for (i in seq_along(counts$n)) {
  gt$heights[7 + (i - 1) * 4] <- unit(counts$n[[i]], "null")
}

png(
  file.path(data_plot_dir, "heatmap-discordant.png"),
  width = 20, height = 30, units = "cm", res = 100
)
grid::grid.draw(gt)
dev.off()
