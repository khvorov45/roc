determine_result <- function(measurement, threshold, equivocal_range = 0) {
  if (equivocal_range > 0) {
    case_when(
      measurement < threshold ~ "neg",
      measurement <= threshold + equivocal_range ~ "eqiv",
      measurement > threshold + equivocal_range ~ "pos"
    )
  } else {
    if_else(
      measurement < threshold, "neg", "pos"
    )
  }
}

calc_result_one_threshold <- function(data,
                                      euro = 0.8,
                                      wantai = 0.9,
                                      svnt = 20) {
  data %>%
    mutate(
      result = case_when(
        startsWith(as.character(assay), "euro")
        ~ determine_result(measurement, euro),
        startsWith(as.character(assay), "wantai")
        ~ determine_result(measurement, wantai),
        assay == "svnt" ~ determine_result(measurement, svnt),
      )
    )
}
