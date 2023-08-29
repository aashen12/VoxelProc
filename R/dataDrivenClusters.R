#' @import irlba
#' @import umap


dataDrivenClusters <- function(voxel_df, n_pca = 20, n_umap = 2, n_clust = 2,
                               region = NULL) {

  #remove all columns that sum to 0

  voxel_df <- voxel_df[,colSums(voxel_df) != 0]
  xyz <- voxel_df[,1:3]
  voxel_df <- voxel_df[,4:ncol(voxel_df)]

  #scale the data

  voxel_df <- scale(voxel_df, center=TRUE, scale=TRUE)

  #perform PCA on data

  pca <- prcomp_irlba(voxel_df, n=n_pca, retx=TRUE)
  rot_data <- pca$x
  pcs <- pca$rotation

  #perform umap

  umap_sim <- umap(rot_data, n_components = n_umap, preserve.seed= TRUE)
  umap_coords <- umap_sim$layout
  umap_df <- as.data.frame(umap_coords)

  #apply kmeans to umap-reduced data

  km <- kmeans(umap_df, centers = n_clust, iter.max = 5000)

  #obtain the dataframe with xyz coords, UMAP coords, and cluster values

  data_df <- cbind(cbind(xyz, umap_coords), km$cluster)

  #obtain a dataframe with umap coords and kmeans clusters to make plotting easier

  plot_df <- as.data.frame(cbind(umap_coords, km$cluster))
  colnames(plot_df) <- c('umap coord 1', 'umap coord 2', 'clusters')
  plot <- ggplot(plot_df, aes(x = 'umap coord 1', y = 'umap coord 2', color = 'clusters')) + geom_point()

  return(data_df, plot)
}