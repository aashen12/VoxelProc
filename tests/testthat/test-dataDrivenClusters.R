library(MASS)
library(tibble)
set.seed(5284)

# define the (nonsensical but valid) data-generating process
cov <- diag(1000)
data <- mvrnorm(1000, rep(0, 1000), cov)
result <- dataDrivenClusters(data)

test_that(paste("dataDirvenClusters returns an error when n_umap < 2"),
          {
            expect_error(dataDrivenClusters(data, n_umap = 1))
          })

test_that(paste("dataDrivenClusters returns an error when n_pca < 2"),
          {
            expect_error(dataDrivenClusters(data, n_pca = 1))
          })

test_that(paste("dataDrivenClusters returns data_df in correct dimensions
                when regions is specified"),
          {
            result_with_region <- dataDrivenClusters(data, region = "Fill")
            num_cols <- dim(result_with_region$data_df)[2]
            expect_equal(num_cols, 7)
          })

test_that(paste("dataDrivenClusters returns data_df in correct dimension
                when regions are not specified"),
          {
            num_cols <- dim(result$data_df)[2]
            expect_equal(num_cols, 6)
          })

test_that(paste("dataDrivenClusters returns data_df as a tibble"),
          {
            expect_equal(is_tibble(result$data_df), TRUE)
          })

# this last test will be updated more, right now it is quite crude
test_that(paste("dataDrivenClusters returns the correct result given the
                specified pipeline"),
          {
            voxel_df <- data[,4:ncol(data)]
            voxel_df <- scale(voxel_df, center=TRUE, scale=TRUE)
            pca <- prcomp_irlba(voxel_df, n = 20, retx = TRUE)
            rot_data <- pca$x
            pcs <- pca$rotation
            umap_sim <- umap(rot_data, n_components = 2, preserve.seed = TRUE)
            umap_coords <- umap_sim$layout
            umap_df <- as.data.frame(umap_coords)
            colnames(umap_df) <- c("U1", "U2")

            expect_equal(sum(umap_df[1]), sum(result$data_df[,4]))
            expect_equal(sum(umap_df[2]), sum(result$data_df[,5]))
          })
