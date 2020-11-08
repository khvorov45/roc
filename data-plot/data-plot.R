cat("plot data")

library(tidyverse)

data_dir <- "data"
data_plot_dir <- "data-plot"

# Functions ===================================================================

source(file.path(data_dir, "read_data.R"))

create_group_lbl <- function(group, symptom_onset_cat) {
  if_else(
    group == "covid", as.character(symptom_onset_cat), as.character(group)
  ) %>%
    factor(c(levels(symptom_onset_cat), levels(group)))
}

save_plot <- function(plot, name, ...) {
  ggdark::ggsave_dark(
    file.path(data_plot_dir, paste0(name, ".png")),
    units = "cm",
    ...
  )
}

# Script ======================================================================

data <- read_data("data") %>%
  filter(!assay %in% c("svnt-20", "svnt-25")) %>%
  mutate(group_lbl = create_group_lbl(group, symptom_onset_cat))

boxplots <- data %>%
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
  ) %>%
  ggplot(aes(group_lbl, measurement)) +
  ggdark::dark_theme_bw(verbose = FALSE) +
  facet_wrap(~assay, scales = "free_y") +
  theme(
    strip.background = element_blank(),
    axis.text.x = element_text(angle = 30, hjust = 1)
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
    aes(xmin = 0.5, xmax = 5.5, ymin = min_threshold, ymax = max_threshold),
    alpha = 0.01
  )

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
    file.path(data_plot_dir, paste0(name, ".png")),
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
