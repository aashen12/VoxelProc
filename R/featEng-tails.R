#' @title computeTailMeans()
#'
#' @description Computes the means of the two-sided alpha-percentile tails of
#' each individual patient stratified on each individual cluster.
#'
#' @details \code{computeTailMeans()} finds the mean of the alpha-percentile tails
#' of each individual patient through the following steps: first, the function
#' extracts the 'PID,', 'x', 'y', 'z', 'value,' and 'cluster' columns from
#' \code{voxel_df} and \code{data_df}. We run \code{cbind} to combine these specific
#' columns. Following this, we isolate each individual cluster label and then
#' extract each individual PID. With each individual patient, we then find the
#' upper and lower alpha-percentile cutoffs, before finding the means of the
#' values that lie above and below those cutoffs, respectively. We then collect
#' these results and put them into a dataframe. If \code{data_df} is not specified,
#' then the function will assume that all of the data comes from a single cluster.
#'
#' @param voxel_df A \code{dataframe} in the following format: the first column
#' is assumed to be an index column counting from 1 to \code{nrow(voxel_df)},
#' and the second column is a column labeled "PID," which contains patient
#' identification. The third, fourth, and fifth columns are 'x', 'y', and 'z'
#' columns. The sixth column contains the value of the voxel that corresponds
#' with the given x-y-z coordinate.
#' @param data_df A \code{tibble} that is meant to be outputted from \code{dataDrivenClusters()}.
#' This tibble should be in the following format: the first three columns are
#' assumed to be 'x', 'y', and 'z' coordinates. The fourth and fifth columns are
#' UMAP coordinates, and the sixth column contains cluster labels. If this is
#' equal to \code{NULL}, \code{computeTailMeans()} will assume that all data comes
#' from the same cluster.
#' @param alpha A \code{double} that indicates what quantiles should be extracted
#' from the patient data. The lower the number, the tighter the tail will be.
#'
#' @return A tibble with the following components:
#' \describe{
#' \item{pid}{a column indicating patient identification}
#' \item{cluster}{column indicating which cluster the given feature was extracted
#' from}
#' \item{upper_quantile}{the upper alpha percentile of the voxel values in the
#' given cluster}
#' \item{lower_quantile}{the lower alpha percentile of the voxel values in the
#' given cluster}
#' \item{mean_upper}{the mean of the values above the upper quantile}
#' \item{mean_lower}{the mean of the values beneath the lower quantile}
#' }
#'
#' @export

computeTailMeans <- function(voxel_df, alpha = 0.05, data_df = NULL) {

  # assuming long format voxel_df in which column 1 is index, column 2 is PID,
  # column 3-5 is x,y,z, and column 6 is value. Assuming data_df to be output
  # from dataDrivenClusters(), in which columns 1-3 are x,y,z, column 4-5 are
  # UMAP coords, and column 6 is cluster labels.

  if (alpha != 0 & alpha != 1) {

    if (!is.null(data_df)) {

      # cbind voxel_df and data_df
      combine_df <- cbind(voxel_df, data_df[c("U1", "U2", "cluster")])

      # removing x, y, z columns since those are redundant after the xyz column.
      combine_df <- combine_df[c("pid", "value", "cluster")]

      # obtaining quantiles with cutoff 1-alpha and alpha
      quantile_summary <- combine_df %>%
        group_by(pid, cluster) %>%
        summarise(upper_quantile = quantile(value, 1-alpha), lower_quantile = quantile(value, alpha))

      # obtaining mean of values above and below the cutoffs
      mean_summary <- combine_df %>%
        group_by(pid, cluster) %>%
        summarise(mean_upper = mean(value[value > quantile(value, 1-alpha)]),
                  mean_lower = mean(value[value < quantile(value, alpha)]))

      # extracting the relevant columns from mean_summary
      mean_cols <- mean_summary[c("mean_upper", "mean_lower")]

      # appending together the correct columns in the right order
      result <- cbind(quantile_summary, mean_cols)[c("pid",
                                                      "cluster",
                                                      "upper_quantile",
                                                      "lower_quantile",
                                                      "mean_upper",
                                                      "mean_lower")] %>%
      pivot_wider(names_from = "cluster",
                  values_from = c("upper_quantile", "lower_quantile", "mean_upper", "mean_lower")) # pivoting the tibble into wide format
  }

    else {

      new_df <- voxel_df[c("pid", "value")] # without data_df, simply follow the same procedure as what is done above, but don't take into consideration a possible cluster label

      upper_summary <- cbind(new_df %>%
        group_by(pid) %>%
        top_frac(alpha) %>%
        summarise(mean_upper = mean(value)),
        new_df %>%
          group_by(pid) %>%
          summarise(upper_quantile = quantile(value, 1-alpha))
      )

      lower_summary <- cbind(new_df %>%
        group_by(pid) %>%
        top_frac(-alpha) %>%
        summarise(mean_lower = mean(value)),
        new_df %>%
          group_by(pid) %>%
          summarise(lower_quantile = quantile(value, alpha))
      )

      # putting all results into dataframe
      result <- as_tibble(cbind(upper_summary, lower_summary)[c("pid",
                                                      "upper_quantile",
                                                      "lower_quantile",
                                                      "mean_upper",
                                                      "mean_lower")])
      }
  }
  else {
    stop("alpha cannot be 0 or 1")
  }
  return(result)
}
