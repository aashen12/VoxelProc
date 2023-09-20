#' @title Incorporating feature engineering with clinical data
#'
#' @export

mergeData <- function(feat_eng_df, id_mapping = NULL, test_clinical_data) {

  # note: still need to generalize this to more colnames.
  if (!is.null(id_mapping)) {
  colnames(id_mapping) <- c("PatID", "pid")
  df_merged <- left_join(feat_eng_df, id_mapping, by = "pid") %>%
    dplyr::relocate(PatID, .after = "pid") %>%
    dplyr::select(-pid)
  result <- left_join(df_merged, test_clinical_data, by = "PatID") %>%
    na.omit()
  }

  else {
  result <- left_join(feat_eng_df, test_clinical_data, by = "pid") %>%
    na.omit()
  }

  return(result)

}
