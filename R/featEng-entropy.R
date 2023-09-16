#' @import entropy
#'
#' @export


computeRegionEntropy <- function(voxel_df, data_df, sep_clusters = TRUE) {

  combine_df <- cbind(voxel_df, data_df[, 4:6]) %>% tibble()

  if (sep_clusters){
  result <- combine_df %>%
    mutate_at(vars("cluster"), factor) %>%
    group_by(pid, cluster) %>%
    mutate(nvox = n()) %>%
    ungroup() %>%
    mutate(avg_vox = mean(nvox)) %>%
    group_by(pid, cluster) %>%
    summarise(entropy = entropy::entropy(entropy::discretize(value, numBins = floor(sqrt(avg_vox[1]))))
              )
  }

  else{
  result <- combine_df %>%
    mutate_at(vars("cluster"), factor) %>%
    group_by(pid) %>%
    mutate(nvox = n()) %>%
    summarise(entropy = entropy::entropy(entropy::discretize(value, numBins = floor(sqrt(nvox[1]))))
              )
  }
}
