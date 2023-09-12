#' @title computeTailMeans()

computeTailMeans <- function(voxel_df, data_df, alpha) {

  #assuming long format data
  combine_df <- cbind(voxel_df, data_df)
  combind_df <- combine_df[, c(2, 3, 4, 5, 6, 12)]
  combine_df <- combine_df %>% mutate(xyz = paste0(combine_df[, 1],", ",combine_df[, 2],", ",combine_df[, 3]))
  combine_df <- combine_df[, c(1, 5, 6, 7)]

  num_clusts <- length(unique(combine_df[, 2]))

  clust_vec <- c()
  PID_vec <- c()
  lower_mean_vec <- c()
  upper_mean_vec <- c()

  for (i in 1:num_clusts) {
    filtered_df <- combine_df %>% filter(Cluster == i)

    for(ID in unique(combine_df[, 1])) {

      uniqueid_df <- filtered_df %>% filter(PID == ID)
      quantile_cutoff_lower <- quantile(unique_df[, 3], alpha)
      quantile_cutoff_upper <- quantile(unique_df[, 3], 1-alpha)
      values_lower <- unique_df[, 3][which(unique_df[, 3] < quantile_cutoff_lower)]
      values_upper <- unique_df[, 3][which(unique_df[, 3] > quantile_cutoff_upper)]
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
