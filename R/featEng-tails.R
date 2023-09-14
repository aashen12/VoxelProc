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
#' UMAP coordinates, and the fifth column contains cluster labels. If this is
#' equal to \code{null}, \code{computeTailMeans()} will assume that all data comes
#' from the same cluster.
#' @param alpha A \code{double} that indicates what quantiles should be extracted
#' from the patient data. The lower the number, the tighter the tail will be.
#'
#' @export

computeTailMeans <- function(voxel_df, data_df = NULL, alpha = 0.05) {

  # assuming long format voxel_df in which column 1 is index, column 2 is PID,
  # column 3-5 is x,y,z, and column 6 is value. Assuming data_df to be output
  # from dataDrivenClusters(), in which columns 1-3 are x,y,z, column 4-5 are
  # UMAP coords, and column 6 is cluster labels.

  if (alpha != 0 & alpha != 1) {

    if (!is.null(data_df)) {

      # cbind voxel_df and data_df
      combine_df <- cbind(voxel_df, data_df[, 4:6])

      #removing x, y, z columns since those are redundant after the xyz column.
      combine_df <- combine_df[c("pid", "value", "cluster")]

      num_clusts <- length(unique(combine_df["cluster"]))

      iterated_df <- data.frame(pid = NA,
                                mean_upper = NA,
                                upper_quantile = NA,
                                mean_lower = NA,
                                lower_quantile = NA,
                                cluster = NA)

        upper_summary <- cbind(combine_df %>%
          group_by(pid, cluster) %>%
          summarise(upper_quantile = quantile(value, 1-alpha)),
          combine_df %>%
            group_by(pid, cluster) %>%
            top_frac(alpha) %>%
            summarise(mean_upper = mean(value))
        )

        lower_summary <- cbind(combine_df %>%
          group_by(pid, cluster) %>%
          summarise(lower_quantile = quantile(value, alpha)),
          combine_df %>%
            group_by(pid, cluster) %>%
            top_frac(-alpha) %>%
            summarise(mean_upper = mean(value))
        )

        result <- cbind(upper_summary, lower_summary)[c("pid",
                                                        "cluster",
                                                        "mean_upper",
                                                        "upper_quantile",
                                                        "mean_lower",
                                                        "lower_quantile")]
  }

    else {

      new_df <- voxel_df[c("pid", "value")]

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
                                                      "mean_upper",
                                                      "upper_quantile",
                                                      "mean_lower",
                                                      "lower_quantile")])
      }
  }
  else {
    stop("alpha cannot be 0 or 1")
  }
  return(result)
}
