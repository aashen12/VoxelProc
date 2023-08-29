#' @title Data Driven Clusters
#'
#' @description \code{dataDrivenClusters()} applies a pipeline of dimensionality
#' reduction and cluster analysis to a given set of MRI voxel data. First, we
#' use principal component analysis to embed the data into some lower-dimensional
#' space. We then apply UMAP to that embedding and reduce it to yet another
#' smaller dimension. Lastly, k-means clustering is utilized to separate the
#' embedded voxel data into separate subgroups. Note that $\code{voxel_df}$
#' should be a n x p matrix, where n is the number of coordinates in that
#' particular region of the brain and p is the number of subjects.
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
#' pipeline will be run on. Default is set to \code{null}.
#'
#' @importFrom irlba prcomp_irlba
#' @importFrom umap umap
#' @importFrom tibble as_tibble
#' @importFrom stats kmeans
#' @import ggplot2
#' @import dplyr
#'
#' @return A tibble with the following components:
#' \describe{
#' \item{data_df}{a tibble that includes the original coordinate values,
#' the embedded UMAP coordinates, and the cluster value of the voxel data.}
#' \item{plot}{a scatterplot of the lower-dimensional UMAP embedding, colored
#' according to the k-means clustering algorithm.}
#' }
#' @export


dataDrivenClusters <- function(voxel_df, n_pca = 20, n_umap = 2, n_clust = 2,
                               region = NULL) {

  # failsafe in case n_pca or UMAP < 2
  if (n_pca < 2 || n_umap < 2){
    stop("Number of principal and UMAP components must be at least 2.")
  }

  else {
  # remove all columns that sum to 0
  voxel_df <- voxel_df[, colSums(voxel_df) != 0]
  xyz <- voxel_df[, 1:3]
  voxel_df <- voxel_df[, 4:ncol(voxel_df)]

  # scale the data
  voxel_df <- scale(voxel_df, center = TRUE, scale = TRUE)

  # perform PCA on data
  pca <- prcomp_irlba(voxel_df, n = n_pca, retx = TRUE)
  rot_data <- pca$x
  pcs <- pca$rotation

  # perform umap
  umap_sim <- umap(rot_data, n_components = n_umap, preserve.seed = TRUE)
  message("PCA and UMAP completed.")
  umap_coords <- umap_sim$layout
  umap_df <- as.data.frame(umap_coords)

  # apply kmeans to umap-reduced data
  km <- kmeans(umap_df, centers = n_clust, iter.max = 5000)
  message("Clustering completed.")

  # obtain the dataframe with xyz coords, UMAP coords, and cluster values
  data_df <- data.frame(cbind(xyz, umap_coords), cluster = km$cluster)
  colnames(data_df)[4:5] <- c("U1", "U2")
  data_df <- as_tibble(data_df)

  if (n_umap == 2) {
    plot <- data_df %>%
      ggplot(aes(x = U1, y = U2, color = factor(cluster))) +
      geom_point() +
      theme_bw() +
      labs(x = "UMAP 1", y = "UMAP 2", title = region, color = "Clusters")
    # For andy: import DNLI ggplot themes
    result <- list(data_df = data_df, plot = plot)
  }
  else {
    print("n_umap is greater than 2, so no visualization is possible")
    result <- list(data_df = data_df)
  }

  # the return object will consist just data_df
  return(result)
  }
}
