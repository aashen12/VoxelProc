#' @title Computing Region Proportions
#'
#' @description Computes the proportion of MRI voxel values that lie within an epsilon
#' ball of a "base value," as is specified by the practitioner.
#'
#' @details \code{computeRegionProp()} takes in four arguments: a \code{tibble}
#' voxel_df and a \code{tibble} data_df, alongside two \code{numerics}, \code{epsilon} and \code{base}.
#' \code{voxel_df} is assumed to be in long format, with each individual pixel
#' being represented on a single row, so that the first columns is 'pid,' columns
#' 2 to 4 are 'x', 'y', 'z', and column 5 is 'value'. \code{data_df} is assumed
#' to be the output of \code{dataDrivenClusters()}. \code{computeRegionProp()}
#' computes the proportion of voxel values that lie between \code{base}-\code{epsilon} and
#' \code{base}+\code{epsilon}, and then displays those proportions as decimals.
#' These proportions are stratified on each individual by clusters.
#'
#' @param voxel_df A \code{dataframe} in the following format: the first column
#' is assumed to be a column labeled "PID," which contains patient
#' identification. The second, third, and fourth columns are 'x', 'y', and 'z'
#' columns. The fifth column contains the value of the voxel that corresponds
#' with the given x-y-z coordinate.
#' @param base A \code{numeric} that indicates where the \code{epsilon}-ball in
#' the calculations should be centered.
#' @param epsilon A \code{numeric} that indicates how wide the ball should be in
#' which we are calculating the proportion of MRI voxel values.
#' @param data_df A \code{tibble} that is meant to be outputted from \code{dataDrivenClusters()}.
#' This tibble should be in the following format: the first three columns are
#' assumed to be 'x', 'y', and 'z' coordinates. The fourth and fifth columns are
#' UMAP coordinates, and the fifth column contains cluster labels. If this is
#' equal to \code{NULL}, \code{computeRegionProp()} will assume that all data comes
#' from the same cluster.
#'
#' @return A tibble with the following components:
#' \describe{
#' \item{pid}{a column indicating patient identification}
#' \item{proportion_c}{columns indicating the calculated proportion of the patient for
#' the cluster that is specified.}
#' }
#'
#' @export

computeRegionProp <- function(voxel_df, base = 1, epsilon = 0.05, data_df = NULL) {

  if (!is.null(data_df)) {
    combine_df <- cbind(voxel_df, data_df[c("U1", "U2", "cluster")]) %>% tibble()
    result <- combine_df %>%
      mutate_at(vars("cluster"), factor) %>%
      group_by(pid, cluster) %>%
      summarize(proportion = sum(value >= base - epsilon & value <= base + epsilon)/length(value)) %>%
      pivot_wider(names_from = "cluster",
                  values_from = "proportion",
                  names_prefix = "proportion_c")
  }


  else {
    result <- voxel_df %>%
      group_by(pid) %>%
      summarize(proportion = sum(value >= base - epsilon & value <= base + epsilon)/length(value))
  }

  return(result)
}
