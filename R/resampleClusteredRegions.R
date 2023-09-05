# This function resamples the data multiple times and returns an ARI matrix
# You will need to use your completed dataDrivenClusters code

resampleClusteredRegions <- function(voxel_df, n_pca = 20,
                                     n_umap = 2, n_clust = 2,
                                     n_resamp = 5,
                                     subsamp_prop = 0.8,
                                     region = NULL) {
  # n-resamp: number of subsamples to take
  # subsamp_prop: proportion of subjects to subsample from (floor this)
  # for ari, you will need mclust::adjustedRandIndex. please add it to namespace
}
