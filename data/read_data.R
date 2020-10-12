read_data <- function(name) {
  read_csv(
    file.path("data", paste0(name, ".csv")),
    col_types = cols(
      id = col_character(),
      assay = col_factor(c(
        "euro_iga", "euro_igg", "euro_ncp", "svnt", "wantai_igm", "wantai_tot"
      )),
      symptom_onset_cat = col_factor(c("<7", "7-14", ">14", "no infection"))
    )
  )
}
