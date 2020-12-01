read_data <- function(name) {
  col_types_general <- list(
    assay = col_factor(c(
      "euro_iga", "euro_igg", "euro_ncp",
      "svnt", "svnt-20", "svnt-25",
      "wantai_igm", "wantai_tot"
    ))
  )
  col_types_specific <- list()
  if (name %in% c("data", "mn")) {
    col_types_specific <- list(
      id = col_character(),
      symptom_onset_cat = col_factor(c("<7", "7-14", ">14", "no infection")),
      group = col_factor(c("covid", "population", "cross-reactive"))
    )
  }
  col_types <- do.call(cols, c(col_types_general, col_types_specific))
  data <- read_csv(
    file.path("data", paste0(name, ".csv")),
    col_types = col_types
  )
  if (name == "assay") {
    data <- data %>%
      mutate(across(c(short, long), fct_reorder, as.integer(assay)))
  }
  data
}
