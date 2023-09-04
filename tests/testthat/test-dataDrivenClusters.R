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
test_that(paste("dataDrivenClusters returns the correct UMAP coordinates"),
          {
            set.seed(5284)
            result <- dataDrivenClusters(data)
            xyz <- data[,1:3]
            voxel_df <- data[,4:ncol(data)]
            voxel_df <- scale(voxel_df, center=TRUE, scale=TRUE)
            pca <- prcomp_irlba(voxel_df, n = 20, retx = TRUE)

            rot_data <- pca$x
            pcs <- pca$rotation
            umap_sim <- umap(rot_data, n_components = 2, preserve.seed = TRUE)
            umap_coords <- umap_sim$layout
            umap_df <- as.data.frame(umap_coords)
            km <- kmeans(umap_df, centers = 2, iter.max = 5000)
            sim_data_df <- data.frame(cbind(xyz, umap_coords), cluster = km$cluster)
            colnames(sim_data_df)[4:5] <- c("U1", "U2")
            sim_data_df <- as_tibble(sim_data_df)

            plot <- sim_data_df %>%
              ggplot(aes(x = U1, y = U2, color = factor(cluster))) +
              geom_point() +
              theme_bw() +
              labs(x = "UMAP 1", y = "UMAP 2", title = NULL, color = "Clusters")

            sim_result <- list(data_df = sim_data_df, plot = plot)
            test <- sim_result$data_df

            expect_equal(test[,4], result$data_df[,4])
            expect_equal(test[,5], result$data_df[,5])
          })
