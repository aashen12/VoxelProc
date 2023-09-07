#' @import mclust
#' @import xpectr
#' @import ggplot2
#' @importFrom reshape2 melt
#' @export

resampleClusteredRegions <- function(voxel_df, n_pca = 20,
                                     n_umap = 2, n_clust = 2,
                                     n_resamp = 5,
                                     subsamp_prop = 0.8,
                                     region = NULL) {

  if (n_umap < 2 || n_pca < 2){

    stop("n_umap and n_pca must be greater than 1")

  }

  else{

  if (subsamp_prop == 1){

    stop("subsamp_prop must be less than 1")

  }

  else{

  # find the number of columns needed to sample from subsamp_prop
  num <- floor(ncol(voxel_df[,4:ncol(voxel_df)])*subsamp_prop)

  # create empty data frame which clusters will be stored in later
  cluster_vec <- data.frame(matrix(NA,
                                   nrow = nrow(voxel_df),
                                   ncol = n_resamp))

  voxel_df_data <- voxel_df[,4:ncol(voxel_df)]

  # extract xyz coords
  xyz <- voxel_df[,1:3]

  # run a for loop to resample and calculate clusters
  for (i in 1:n_resamp){

    # sample voxel data from voxel_df
    subsamp <- voxel_df_data[,sample(ncol(voxel_df_data), size=num)]

    # bind xyz back to the subsample
    new_voxel_df <- cbind(xyz, subsamp)

    # run dataDrivenClusters() with specified values
    DDC <- suppressMessages(dataDrivenClusters(new_voxel_df, n_pca = n_pca, n_umap = n_umap,
                              n_clust = n_clust, region = region))
    clusters <- DDC$data_df[,6]

    # append to cluster_vec
    cluster_vec[,i] <- clusters
  }

message("Cluster vectors calculated")

    # build matrix of ARI values
    ARI <- matrix(NA, n_resamp, n_resamp)

    # reiterate to find pairwise ARIs
  for (i in 1:n_resamp){
    for (j in c(1:n_resamp)[1:n_resamp != i]) {
      ARI_j <- adjustedRandIndex(cluster_vec[,i], cluster_vec[,j])
      ARI[i,j] <- ARI_j
    }

  }

message("Matrix calculated")

    # turning the matrix upper-triangular
    ARI[lower.tri(ARI)] <- 0

    # extracting individual values that are off-diagonal
    values <- ARI[col(ARI)!=row(ARI)]

    # setting diags to 1
    diag(ARI) <- 1

    # taking average
    avg <- mean(values)

    # plotting matrix
    plot <- ggplot(melt(ARI), aes(x = Var1, y = Var2, fill = value)) +
      geom_tile(color = "white") +
      geom_text(aes(label = round(value, 2)), color = "white") +
      theme_bw() +
      scale_fill_gradient2(name = "ARI", limits = c(0, 1)) +
      labs(x = "Sample #", y = "Sample #")

    # appending to list
    result <- list(Average = avg, Matrix = plot)

  return(result)

    }
  }
}
