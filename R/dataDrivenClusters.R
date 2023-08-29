#' @import irlba
#' @import umap


dataDrivenClusters <- function(voxel_df, n_pca = 20, n_umap = 2, n_clust = 2,
                               region = NULL) {

  # failsafe in case n_pca or UMAP < 2

  if (n_pca < 2 || n_umap < 2){

    stop("either n_pca < 2 or n_umap < 2, so no visualization is possible")

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
  umap_coords <- umap_sim$layout
  umap_df <- as.data.frame(umap_coords)

  # apply kmeans to umap-reduced data
  km <- kmeans(umap_df, centers = n_clust, iter.max = 5000)

  # obtain the dataframe with xyz coords, UMAP coords, and cluster values
  data_df <- cbind(cbind(xyz, umap_coords), km$cluster)

  # obtain a dataframe with umap coords and kmeans clusters to make plotting easier
  plot_df <- as.data.frame(cbind(umap_coords, km$cluster))
  colnames(plot_df) <- c("umap_coord_1", "umap_coord_2", "clusters")

  if (n_umap == 2) {
    plot <- plot_df %>%
      ggplot(aes(x = umap_coord_1, y = umap_coord_2, color = clusters)) +
      geom_point() +
      theme_bw() +
      labs(x = "UMAP 1", y = "UMAP 2", title = region)
    result <- list(data_df = data_df, plot = plot)
  }
  else {
    print("n_umap is greater than 2, so no visualization is possible")
    result <- list(data_df = data_df)
  }

  # the return object will consist of two dataframes, data_df and plot_df
  return(result)
  }
}
