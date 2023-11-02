#' @title Computing Region Entropy
#'
#' @description Computes the entropy of the MRI voxel values for each patient
#' stratified on individual clusters.
#'
#' @details \code{computeRegionEntropy()} takes in two arguments: a \code{tibble}
#' voxel_df and a \code{tibble} data_df. \code{voxel_df} is assumed to be in long
#' format, with each individual pixel being represented on a single row, so that
#' the first columns is 'pid,' columns 2 to 4 are 'x', 'y', 'z', and column 5
#' is 'value'. \code{data_df} is assumed to be the output of \code{dataDrivenClusters()}.
#' \code{computeRegionEntropy()} combines \code{voxel_df} and \code{data_df} into
#' a single \code{tibble}, and then pivots this combined data with respect to
#' the columns 'pid' and 'cluster'. We then find the entropy of each pid-cluster
#' pair. If \code{data_df} is \code{null}, then the function assumes that all
#' the observations come from the same cluster.
#'
#' @param voxel_df A \code{dataframe} in the following format: the first column
#' is assumed to be a column labeled "PID," which contains patient
#' identification. The second, third, and fourth columns are 'x', 'y', and 'z'
#' columns. The fifth column contains the value of the voxel that corresponds
#' with the given x-y-z coordinate.
#' @param data_df A \code{tibble} that is meant to be outputted from \code{dataDrivenClusters()}.
#' This tibble should be in the following format: the first three columns are
#' assumed to be 'x', 'y', and 'z' coordinates. The fourth and fifth columns are
#' UMAP coordinates, and the fifth column contains cluster labels. If this is
#' equal to \code{NULL}, \code{computeRegionEntropy()} will assume that all data comes
#' from the same cluster.
#'
#' @import entropy
#'
#' @return A tibble with the following components:
#' \describe{
#' \item{pid}{a column indicating patient identification}
#' \item{entropy_c}{columns indicating the calculated entropy of the patient for
#' the cluster that is specified.}
#' }
#'
#' @export


computeRegionEntropy <- function(voxel_df, data_df = NULL) {

  if (!is.null(data_df)) {
  combine_df <- cbind(voxel_df, data_df[, 4:6]) %>% tibble()
  result <- combine_df %>%
    mutate_at(vars("cluster"), factor) %>%
    group_by(pid, cluster) %>%
    mutate(nvox = n()) %>%
    ungroup() %>%
    mutate(avg_vox = mean(nvox)) %>%
    group_by(pid, cluster) %>%
    summarise(entropy = entropy(discretize(value, numBins = floor(sqrt(avg_vox[1]))))
    ) %>%
    pivot_wider(names_from = "cluster", values_from = "entropy", names_prefix = "entropy_c")
  }


  else {
  result <- voxel_df %>%
    group_by(pid) %>%
    mutate(nvox = n()) %>%
    summarise(entropy = entropy(discretize(value, numBins = floor(sqrt(nvox[1])))))
  }

  return(result)
}
