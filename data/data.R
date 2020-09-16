cat("extract data for analysis")

library(tidyverse)

data_raw_dir <- "data-raw"
data_dir <- "data"

# Functions ===================================================================

save_data <- function(data, name) {
  write_csv(data, file.path(data_dir, paste0(name, ".csv")))
}

# Script ======================================================================
