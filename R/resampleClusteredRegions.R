#' @title resampleClusteredRegions
#'
#' @description Tests the stability of the output of \code{dataDrivenClusters()}
#' by resampling the voxel data set and using the adjusted Rand Index to compare
#' cluster values.
#'
#' @details \code{resampleClusteringRegions} randomly resamples \code{voxel_df}
#' according to the arguments \code{n_resamp} and \code{subsamp_prop}. Depending
#' on these parameters, a certain proportion of the patients in \code{voxel_df}
#' will be chosen, and \code{dataDrivenClusters()} will be ran on these subsets.
#' After running \code{dataDrivenClusters} on the subsets, the cluster
#' labels will be extracted for each resample. \code{resampleClusteringRegions()}
#' then uses the adjusted Rand Index to find how similar each pair of clusterings are to
#' each other. The function then averages the pairwise adjusted Rand indices as
#' a 'score' for how stable the clusterings are. In addition to this a matrix
#' plot is provided to visualize the pairwise rand indices.
#'
#'
#' @param voxel_df A \code{list} or \code{matrix} object that represents voxel data.
#' The first three columns should be coordinate data, while the rest of the
#' columns represent voxel data.
#' @param n_pca A \code{numeric} that indicates how many principal components
#' will be used in the pipeline. Default is set to 20.
#' @param n_umap A \code{numeric} that indicates how many umap components will
#' be used in the pipeline. Default is set to 2.
#' @param n_clust A \code{numeric} that indicates how many clusters the k-means
#' algorithm will search for in the pipeline. Default is set to 2.
#' @param region A \code{character} that indicates which region of the brain the
#' pipeline will be run on. Default is set to \code{NULL}.
#' @param n_resamp A \code{numeric} that indicates how many times
#' \code{resampleClusteredRegions()} should resample from \code{voxel_df}.
#' @param subsamp_prop A \code{numeric} indicating what proportion of \code{voxel_df}
#' will be resampled at each iteration of \code{resampleClusteredRegions}.
#'
#' @import mclust
#' @import ggplot2
#' @import reshape2
#' @import foreach
#'
#' @return A list with the following components:
#' \describe{
#' \item{Average}{the average of the Rand indices that have been recorded for
#' each pairwise resampling}
#' \item{Matrix}{a matrix plot that visualizes the pairwise rand indices of the
#' resamplings.}
#' }
#'
#' @export

resampleClusteredRegions <- function(voxel_df_long, n_pca = 20,
                                     n_umap = 2, n_clust = 2,
                                     n_resamp = 5,
                                     subsamp_prop = 0.8,
                                     region = NULL) {
  # failsafes
  if (n_umap < 2 || n_pca < 2){
    stop("n_umap and n_pca must be greater than 1")
  }

  else{

  if (subsamp_prop == 1 || subsamp_prop == 0){
    stop("subsamp_prop must be less than 1")
  }

  else{

  voxel_df <- voxel_df_long %>%
      tidyr::pivot_wider(names_from = "pid", values_from = "value")

  # find the number of columns needed to sample from subsamp_prop
  num <- floor((ncol(voxel_df) - 3)*subsamp_prop)

  # create empty data frame which clusters will be stored in later
  cluster_vec <- data.frame(matrix(NA,
                                   nrow = nrow(voxel_df),
                                   ncol = n_resamp))
  voxel_df_data <- select(voxel_df, -c("x", "y", "z"))

  # extract xyz coords
  xyz <- voxel_df[c("x", "y", "z")]

  # run parallelized stuff
  parloop <- foreach (i=1:n_resamp, .combine = "c") %dopar% {

    # sample voxel data from voxel_df
    subsamp <- voxel_df_data[, sample(ncol(voxel_df_data), size = num)]

    # bind xyz back to the subsample
    new_voxel_df <- cbind(xyz, subsamp)
    new_voxel_df <- new_voxel_df %>%
      tidyr::pivot_longer(cols = !(x | y | z), names_to = "pid", values_to = "value") %>%
      dplyr::relocate(pid, .before = "x")

    # run dataDrivenClusters() with specified values
    DDC <- suppressMessages(dataDrivenClusters(new_voxel_df, n_pca = n_pca, n_umap = n_umap,
                                               n_clust = n_clust, region = region))
    clusters <- DDC$data_df[, 6]

    # append to cluster_vec
    cluster_vec[, i] <- clusters
  }

  # run a fast for loop
  for (i in 1:n_resamp){
    cluster_vec[,i] <- parloop[i]
  }

  # build matrix of ARI values
  ARI <- matrix(NA, n_resamp, n_resamp)

  # reiterate to find pairwise ARIs
  for (i in 1:n_resamp) {
    for (j in c(1:n_resamp)[1:n_resamp != i]) {
      ARI_j <- adjustedRandIndex(cluster_vec[, i], cluster_vec[, j])
      ARI[i, j] <- ARI_j
    }
  }

  message("Matrix calculated")

  # turning the matrix upper-triangular
  ARI[upper.tri(ARI)] <- NA

  # setting diags to 1
  diag(ARI) <- NA

  # taking average
  avg <- round(mean(ARI[!is.na(ARI)]), digits = 2)

  # plotting matrix
  plot <- ggplot(melt(ARI), aes(x = Var1, y = Var2, fill = value)) +
    geom_tile(color = "white") +
    geom_text(aes(label = round(value, 2), fontface = "bold"), color = "black") +
    theme_bw() +
    scale_fill_gradient2(name = "ARI", limits = c(0, 1)) +
    labs(x = "Sample #", y = "Sample #",
          title = "ARI Matrix",
          subtitle = paste0("average = ", avg)) +
    theme(plot.title = element_text(face = "bold"))

  # appending to list
  result <- list(average = avg, matrix = plot, values = ARI)

  return(result)

    }
  }
}
