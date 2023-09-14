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
      combine_df <- cbind(voxel_df, data_df)

      # extracting the columns in following order: PID, x, y, z, Value, Cluster
      combind_df <- combine_df[, c(2, 3, 4, 5, 6, 12)]

      # adding additional xyz identifier column
      combine_df <- combine_df %>% mutate(xyz = paste0(combine_df[, 1],", ",combine_df[, 2],", ",combine_df[, 3]))

      #removing x, y, z columns since those are redundant after the xyz column.
      combine_df <- combine_df[, c(1, 5, 6, 7)]

      num_clusts <- length(unique(combine_df[, 2]))

      # initializing empty vectors
      clust_vec <- rep(NA, num_clusts)
        PID_vec <-
        lower_mean_vec <-
        upper_mean_vec <-
        lower_quantile_vec <-
        upper_quantile_vec <- rep(NA, length(unique(combine_df[,1])))

      # cycling through each cluster label
      for (i in 1:num_clusts) {

        # isolating individual clusters
        filtered_df <- combine_df %>% filter(cluster == i)

          # for each cluster, isolate individual PID
          for(j in seq_along(unique(combine_df[,1]))) {
          ID <- unique(combine_df[,1])[j]

          # filter the dataframe so that we are only considering individual PIDs
          uniqueid_df <- filtered_df %>% filter(pid == ID)

          # find quantile cutoffs
          quantile_cutoff_lower <- quantile(uniqueid_df[, 3], alpha)
          quantile_cutoff_upper <- quantile(uniqueid_df[, 3], 1-alpha)

          # find values that fall beneath and above quantile cutoffs
          values_lower <- uniqueid_df[, 3][which(uniqueid_df[, 3] < quantile_cutoff_lower)]
          values_upper <- uniqueid_df[, 3][which(uniqueid_df[, 3] > quantile_cutoff_upper)]

          # finding averages
          mean_lower <- mean(values_lower)
          mean_upper <- mean(values_upper)

          # appending all values into a vector
          clust_vec[i] <- i
          PID_vec[j] <- ID
          lower_mean_vec[j] <- mean_lower
          upper_mean_vec[j] <- mean_upper
          lower_quantile_vec[j] <- quantile_cutoff_lower
          upper_quantile_vec[j] <- quantile_cutoff_upper

          }

        }

      # putting all results into a dataframe
      result <- data.frame(PID = PID_vec,
                       Cluster = clust_vec,
                       LowerMean = lower_mean_vec,
                       UpperMean = upper_mean_vec,
                       LowerQuantile = lower_quantile_vec,
                       UpperQuantile = upper_quantile_vec)

  }

    else {

      new_df <- voxel_df[, c(2:6)]

      PID_vec <-
        lower_mean_vec <-
        upper_mean_vec <-
        lower_quantile_vec <-
        upper_quantile_vec <- rep(NA, length(unique(new_df[,1])))

      for (i in seq_along(unique(new_df[,1]))) {

        ID <- unique(new_df[,1])[i]

        # filter the dataframe so that we are only considering individual PIDs
        uniqueid_df <- new_df %>% filter(pid == ID)

        # find quantile cutoffs
        quantile_cutoff_lower <- quantile(uniqueid_df[, 5], alpha)
        quantile_cutoff_upper <- quantile(uniqueid_df[, 5], 1-alpha)

        # find values that fall beneath and above quantile cutoffs
        values_lower <- uniqueid_df[, 5][which(uniqueid_df[, 5] < quantile_cutoff_lower)]
        values_upper <- uniqueid_df[, 5][which(uniqueid_df[, 5] > quantile_cutoff_upper)]

        # finding averages
        mean_lower <- mean(values_lower)
        mean_upper <- mean(values_upper)

        # appending all values into a vector
        PID_vec[i] <- ID
        lower_mean_vec[i] <- mean_lower
        upper_mean_vec[i] <- mean_upper
        lower_quantile_vec[i] <- quantile_cutoff_lower
        upper_quantile_vec[i] <- quantile_cutoff_upper

      }

      # putting all results into dataframe
      result <- data.frame(PID = PID_vec,
                          LowerMean = lower_mean_vec,
                          UpperMean = upper_mean_vec,
                          LowerQuantile = lower_quantile_vec,
                          UpperQuantile = upper_quantile_vec)
      }

  }

  else {

    stop("alpha cannot be 0 or 1")

  }

return(result)

}
