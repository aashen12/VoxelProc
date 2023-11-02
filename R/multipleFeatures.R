#' @title Compute Multiple Features
#'
#' @description Computes a specified list of features on a single long-format
#' \code{voxel_df}, stratifying results according to cluster if \code{data_df} is
#' specified.
#'
#' @details For each function that is included in the string-vector \code{methods},
#' \code{computeFeatures()} will run that function on a long-format \code{voxel_df}.
#' If \code{computeRegionProp()} is invoked, \code{base} and \code{epsilon} will
#' determine the epsilon-ball around which the proportion of MRI voxel values will be
#' calculated. \code{data_df} is meant to provide information about the cluster values,
#' and if this is not specified in the arguments of \code{computeFeatures()}, then
#' it will be assumed that all of the MRI voxel values come from the same cluster.
#'
#' @param voxel_df A \code{dataframe} in the following format: the first column
#' is assumed to be a column labeled "PID," which contains patient
#' identification. The second, third, and fourth columns are 'x', 'y', and 'z'
#' columns. The fifth column contains the value of the voxel that corresponds
#' with the given x-y-z coordinate.
#' @param methods A \code{character} vector that indicates which feature engineering
#' functions should be run on \code{voxel_df}. By defaul this is set to
#' \code{c("computeRegionEntropy", "computeRegionProp")}. When entering \code{methods},
#' type the function name as it appears when typing \code{VoxelProc::} into console.
#' @param base A \code{numeric} that indicates where the \code{epsilon}-ball in
#' the calculations should be centered, if \code{computeRegionProp()} is specified
#' in \code{methods}.
#' @param epsilon A \code{numeric} that indicates how wide the ball should be in
#' which we are calculating the proportion of MRI voxel values, if \code{computeRegionProp()}
#' is specified in \code{methods}.
#' @param data_df A \code{tibble} that is meant to be outputted from \code{dataDrivenClusters()}.
#' This tibble should be in the following format: the first three columns are
#' assumed to be 'x', 'y', and 'z' coordinates. The fourth and fifth columns are
#' UMAP coordinates, and the fifth column contains cluster labels. If this is
#' equal to \code{NULL}, \code{computeFeatures()} will assume that all data comes
#' from the same cluster.
#'
#' @import rlang
#'
#' @return A tibble with the following components:
#' \describe{
#' \item{pid}{a column indicating patient identification}
#' \item{feature_c}{columns indicating the features that have been extracted from
#' the data, with "c" indicating which cluster they came from.}
#' }
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
