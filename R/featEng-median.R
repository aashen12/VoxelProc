computeRegionMedian <- function(voxel_df, data_df = NULL) {

  if (!is.null(data_df)) {
    combine_df <- cbind(voxel_df, data_df[c("U1", "U2", "cluster")]) %>% tibble()
    result <- combine_df %>%
      mutate_at(vars("cluster"), factor) %>%
      group_by(pid, cluster) %>%
      summarize(median = median(value))
  }


  else {
    result <- voxel_df %>%
      group_by(pid) %>%
      summarize(median = median(value))
  }

  return(result)
}
