#' @export

computeRegionProp <- function(voxel_df, base = 1, epsilon = 0.05, data_df = NULL) {

  if (!is.null(data_df)) {
    combine_df <- cbind(voxel_df, data_df[c("U1", "U2", "cluster")]) %>% tibble()
    result <- combine_df %>%
      mutate_at(vars("cluster"), factor) %>%
      group_by(pid, cluster) %>%
      summarize(proportion = sum(value >= base - epsilon & value <= base + epsilon)/length(value))
  }


  else {
    result <- voxel_df %>%
      group_by(pid) %>%
      summarize(proportion = sum(value >= base - epsilon & value <= base + epsilon)/length(value))
  }

  return(result)
}
