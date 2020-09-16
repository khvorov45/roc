read_data <- function(name) {
  read_csv(
    file.path("data", paste0(name, ".csv")),
    col_types = cols(
      subject_id = col_character(),
      days_post_symptom_onset = col_integer()
    )
  )
}
