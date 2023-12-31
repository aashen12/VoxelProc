#' @title Data Driven Clusters
#'
#' @description \code{dataDrivenClusters()} applies a pipeline of dimensionality
#' reduction and cluster analysis in hopes of isolating signals of interest in
#' high-dimensional MRI voxel data.
#'
#' @param voxel_df_long A long \code{data.frame} or \code{matrix} object that
#' represents spatial patient-voxel data.
#' @param n_pca A \code{numeric} that indicates how many principal components
#' will be used in the pipeline. Default is set to 20.
#' @param n_umap A \code{numeric} that indicates how many umap components will
#' be used in the pipeline. Default is set to 2.
#' @param n_clust A \code{numeric} that indicates how many clusters the k-means
#' algorithm will search for in the pipeline. Default is set to 2.
#' @param region A \code{character} that indicates which region of the brain the
#' pipeline will be run on. Default is set to \code{NULL}.
#'
#' @details The pipeline uses two layers of dimensionality reduction in
#' order to isolate signals of interest in high-dimensional voxel data. The first
#' layer, PCA, will embed the data into a smaller but still non-visualizable
#' space. The second layer, UMAP, will embed this data into a space specified by
#' the user. The default for this is the visualizable 2D space. Lastly, k-means
#' is used to algorithmically determine where the various clusters lie. The
#' argument \code{voxel_df_long} should be an n by 5 data frame.
#' The first three columns should be coordinate data, while the rest of the
#' columns represent voxel data and patient ID values. The column names of the
#' coordinate values must be \code{x}, \code{y}, and \code{z}, verbatim,
#' representing the location of the given voxels in the brain.
#' The column name of the voxel values must be \code{value}, and the column
#' name of the patient IDs must be \code{pid}. The pipeline runs on the
#' voxel data, not on the coordinates.
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



dataDrivenClusters <- function(voxel_df_long, n_pca = 20, n_umap = 2, n_clust = 2,
                               region = NULL) {

  # failsafe in case n_pca or UMAP < 2
  if (n_pca < 2 || n_umap < 2){
    stop("Number of principal and UMAP components must be at least 2.")
  }

  if (ncol(voxel_df_long) != 5) {
    warning("Number of columns in input data does not equal the intended value (5).")
  }

  if (all(sort(names(voxel_df_long)) != c("pid", "value", "x", "y", "z"))) {
    warning("Column names in input data are not equal to the intended names.")
  }

  else {
  # remove all columns that sum to 0
  voxel_df <- voxel_df_long %>%
    tidyr::pivot_wider(names_from = "pid", values_from = "value")
  voxel_df <- voxel_df[, colSums(voxel_df) != 0]
  xyz <- voxel_df[c("x", "y", "z")]
  voxel_df <- select(voxel_df, -c("x", "y", "z"))

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

  # Adding region column
  if (!is.null(region)){
    region_col <- rep(region, nrow(voxel_df))
    data_df <- as_tibble(cbind(data_df, region_col))
  }
  else {
    data_df <- data_df
  }

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

  # the return object will consist of just data_df
  return(result)
  }
}
