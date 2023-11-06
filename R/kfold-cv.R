#' @import ISLR

crossValidation <- function(voxel_df, k = 10) {
  voxel_df_wide <- voxel_df %>%
    tidyr::pivot_wider(names_from = "pid", values_from = "value")

  cv.error <- c()
  for (i in 1:k) {
    glm.fit <- glm(dependentvariable ~ poly(independentvariable, i), data = data)
    cv.error <- c(cv.error, cv.glm(data, glm.fit, K = k)$delta[1])
  }

  return(cv.error)
}
