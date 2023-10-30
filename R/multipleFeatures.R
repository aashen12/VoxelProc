#' @import rlang
#'
#' @export
computeFeatures <- function(voxel_df,
                            methods = c("computeRegionEntropy", "computeRegionProp"),
                            base = 1,
                            epsilon = 0.05,
                            data_df = NULL) {
  i <- 1
  for (fun in methods) {
    if (i == 1) {
      if (fun == "computeRegionProp") {
        df <- eval(parse(text = paste(fun,
                                      "(voxel_df",
                                      ",",
                                      as.character(base),
                                      ",",
                                      as.character(epsilon),
                                      ",",
                                      "data_df)")))
        feat_df <- df
        i <- i + 1
      }
      else{
        df <- eval(parse(text = paste(fun, "(voxel_df, data_df)")))
        feat_df <- df
        i <- i + 1
      }
    }
    else {
      if (fun == "computeRegionProp") {
        df <- eval(parse(text = paste(fun,
                                      "(voxel_df",
                                      ",",
                                      as.character(base),
                                      ",",
                                      as.character(epsilon),
                                      ",",
                                      "data_df)")))
        desired_columns <- df[, !names(df) %in% c("pid", "cluster")]
        feat_df <- cbind(feat_df, desired_columns)
      }
      else{
        df <- eval(parse(text = paste(fun, "(voxel_df, data_df)")))
        desired_columns <- df[, !names(df) %in% c("pid", "cluster")]
        feat_df <- cbind(feat_df, desired_columns)
      }
    }
  }
  feat_df <- feat_df %>% tibble()
  return(feat_df)
}
