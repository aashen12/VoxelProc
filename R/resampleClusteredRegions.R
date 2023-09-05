# This function resamples the data multiple times and returns an ARI matrix
# You will need to use your completed dataDrivenClusters code
#' @importFrom mclust adjustedRandIndex

resampleClusteredRegions <- function(voxel_df, n_pca = 20,
                                     n_umap = 2, n_clust = 2,
                                     n_resamp = 5,
                                     subsamp_prop = 0.8,
                                     region = NULL) {
  # n-resamp: number of subsamples to take
  # subsamp_prop: proportion of subjects to subsample from (floor this)
  # for ari, you will need mclust::adjustedRandIndex. please add it to namespace

  # find the number of columns needed to sample from subsamp_prop
  num <- floor(ncol(voxel_df[,4:ncol(voxel_df)])*subsamp_prop)

  # create empty data frame which clusters will be stored in later
  cluster_vec <- data.frame(matrix(NA,
                                   nrow = nrow(voxel_df),
                                   ncol = n_resamp))

  # run a for loop to resample and calculate clusters
  for (i in 1:n_resamp){

    # extract xyz coords
    xyz <- voxel_df[,1:3]

    # sample voxel data from voxel_df
    voxel_df <- voxel_df[,4:ncol(voxel_df)]
    subsamp <- voxel_df[,sample(ncol(voxel_df), size=num)]

    # bind xyz back to the subsample
    new_voxel_df <- cbind(xyz, subsamp)

    # run dataDrivenClusters() with specified values
    DDC <- dataDrivenClusters(new_voxel_df, n_pca = n_pca, n_umap = n_umap,
                              n_clust = n_clust, region = region)
    clusters <- DDC$data_df[,6]

    # append to cluster_vec
    cluster_vec[,i] <- clusters
  }

  ARI <- matrix(NA, n_resamp, n_resamp)

  for (i in 1:n_resamp){
    for (j in c(1:n_resamp)[1:n_resamp != i]) {
      ARI_j <- adjustedRandIndex(cluster_vec[,i], cluster_vec[,j])
      ARI[i,j] <- ARI_j
    }

  }


  ARI[lower.tri(ARI)] <- 0
  values <- ARI[col(ARI)!=row(ARI)]
  avg <- mean(values)
  result <- list(Average = avg, Matrix = ARI)

  return(result)

}
