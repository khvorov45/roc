cat("plot data")

library(tidyverse)

data_dir <- "data"
data_plot_dir <- "data-plot"

# Functions ===================================================================

source(file.path(data_dir, "read_data.R"))

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
  filter(!is.na(symptom_onset_cat), abs(measurement) < 100) %>%
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
