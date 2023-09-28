#' @title Incorporating feature engineering with clinical data
#'
#' @description Combines the feature engineering results from \code{computeRegionEntropy()}
#' with clinical data to form a new data set.
#'
#' @details The function takes in three parameters: \code{feat_eng_df} is an output
#' of \code{computeRegionEntropy()}. \code{clinical_data} is the clinical data
#' that one wants to merge the entropy data with. \code{id_mapping} is an optional
#' argument that allows the practitioner to map different forms of patient
#' identification with each other. If this is \code{null}, then it will be assumed
#' that the 'pid' column in \code{feat_eng_df} is the same identification as what
#' is used in \code{clinical_data}. If \code{id_mapping} is not \code{NULL},
#' then 'pid' will be used as the common identification between \code{feat_eng_df}
#' and \code{clinical_data}.
#'
#' @param feat_eng_df A \code{tibble} that takes the same form as an output of
#' \code{computeRegionEntropy()}.
#' @param test_clinical_data A \code{dataframe} that includes the clinical data
#' of the patients.
#' @param match A \code{character} indicating which column name corresponds to
#' the identificaton seen in \code{feat_eng_df} and \code{id_mapping}
#' @param id_mapping A \code{dataframe} that has at least the columns 'PatID' and
#' match, the latter of which is an argument in the function. If 'match' is not
#' specified, then the function will assume that 'ScrID' is the second column name.
#'
#' @returns A \code{tibble} that outputs the original clinical dataframe with an
#' added entropy column.
#'
#' @import tidyr
#'
#' @export

mergeData <- function(feat_eng_df, clinical_data, match = "ScrID", id_mapping = NULL) {

  # note: still need to generalize this to more colnames.
  if (!is.null(id_mapping)) {
  id_mapping <- id_mapping %>% rename("pid" = match)
  df_merged <- left_join(feat_eng_df, id_mapping, by = "pid") %>%
    dplyr::relocate(PatID, .after = "pid") %>%
    dplyr::select(-pid, -SubjID, -DataDir)
  # why does this line exist?
  result <- left_join(clinical_data, df_merged, by = "PatID") %>%
    tidyr::drop_na(PatID)
  }

  else {
  colnames(clinical_data)[1] <- "pid"
  result <- left_join(clinical_data, feat_eng_df, by = "pid") %>%
    tidyr::drop_na(pid)
  }

  return(result)

}
