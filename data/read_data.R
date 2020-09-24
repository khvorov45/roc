read_data <- function(name) {
  read_csv(
    file.path("data", paste0(name, ".csv")),
    col_types = cols(
      id = col_character(),
      symptom_onset_cat = col_factor(c("<7", "7-14", ">14", "no infection"))
    )
  )
}
