#' @title computeTailMeans()

computeTailMeans <- function(voxel_df, data_df, alpha) {

  # assuming long format voxel_df in which column 1 is index, column 2 is PID,
  # column 3-5 is x,y,z, and column 6 is value. Assuming data_df to be output
  # from dataDrivenClusters(), in which columns 1-3 are x,y,z, column 4-5 are
  # UMAP coords, and column 6 is cluster labels.

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
  clust_vec <- c()
  PID_vec <- c()
  lower_mean_vec <- c()
  upper_mean_vec <- c()

  # cycling through each cluster label
  for (i in 1:num_clusts) {

    # isolating individual clusters
    filtered_df <- combine_df %>% filter(Cluster == i)

    # for each cluster, isolate individual PID
    for(ID in unique(combine_df[, 1])) {

      # filter the dataframe so that we are only considering individual PIDs
      uniqueid_df <- filtered_df %>% filter(PID == ID)

      # find quantile cutoffs
      quantile_cutoff_lower <- quantile(unique_df[, 3], alpha)
      quantile_cutoff_upper <- quantile(unique_df[, 3], 1-alpha)

      # find values that fall beneath and above quantile cutoffs
      values_lower <- unique_df[, 3][which(unique_df[, 3] < quantile_cutoff_lower)]
      values_upper <- unique_df[, 3][which(unique_df[, 3] > quantile_cutoff_upper)]

      # finding averages
      mean_lower <- mean(values_lower)
      mean_upper <- mean(values_upper)

      clust_vec <- c(clust_vec, i)
      PID_vec <- c(PID_vec, ID)
      lower_mean_vec <- c(lower_mean_vec, mean_lower)
      upper_mean_vec <- c(upper_mean_vec, mean_upper)

    }
  }

  result <- data.frame(PID = PID_vec, Cluster = clust_vec, LowerMean = lower_mean_vec, UpperMean = upper_mean_vec)

  return(result)

}
