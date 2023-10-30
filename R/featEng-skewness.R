#' @import moments
#'
#' @export

computeRegionSkewness <- function(voxel_df, data_df = NULL) {

  if (!is.null(data_df)) {
    combine_df <- cbind(voxel_df, data_df[c("U1", "U2", "cluster")]) %>% tibble()
    result <- combine_df %>%
      mutate_at(vars("cluster"), factor) %>%
      group_by(pid, cluster) %>%
      summarize(skewness = skewness(value)) %>%
      pivot_wider(names_from = "cluster",
                  values_from = "skewness",
                  names_prefix = "skewness_c")
  }


  else {
    result <- voxel_df %>%
      group_by(pid) %>%
      summarize(skewness = skewness(value))
  }

  return(result)
}
