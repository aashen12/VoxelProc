#' @import ISLR
#' @import caret
#' @import survival
#' @import survcomp
#' @import glmnet

crossValidation <- function(clinical_data,
                            k = 5,
                            method = "coxph") {

  df <- clinical_data[-1]
  covariates <- clinical_data[, -which(names(df) %in% c("time", "status"))]
  proportion <- 1-1/k

  if (method == "coxph") {
    covar_names <- colnames(covariates)
    score <- c()
    for (j in 1:ncol(covariates)) {
      include <- covar_names[1:j]
      mini_df <- df[, which(names(df) %in% c(include, "time", "status"))]
      cvector <- c()
      for (i in 1:k) {
        inTrain <- createDataPartition(y = mini_df$time,
                                       p = proportion,
                                       list = FALSE)
        train <- mini_df[inTrain, ]
        test <- mini_df[-inTrain, ]
        coxph_model <- coxph(Surv(time, status) ~., data = train)
        test_pred <- coxph_model %>% predict(test)
        cindex <- concordance.index(test_pred, test$time, test$status)
        cvector <- c(cvector, cindex)
      }
      i_score <- mean(cvector)
      score <- c(score, i_score)
    }
    best_index <- which(score == min(score))
    best_covariates <- covar_names[1:best_index]
    best_df <- df[, which(names(df) %in% c(best_covariates, "time", "status"))]
    coxph_mod <- coxph(Surv(time, status) ~., data = best_df)
    result <- coxph_mod

  }
  else if (method == "LASSO") {
    df <- df %>% select(where(is.numeric))
    scale_data <- as.data.frame(scale(df, center = TRUE))
    covariates <- scale_data[, -which(names(scale_data) %in% c("status", "time"))] %>% as.matrix()
    status <- df["status"][[1]]
    time <- df["time"][[1]]
    cvfit <- cv.glmnet(covariates,
                       Surv(time, status),
                       nfolds = k,
                       alpha = 1,
                       family = "cox",
                       type.measure = "C")
    result <- cvfit
  }
  else if (method == "Ridge") {
    df <- df %>% select(where(is.numeric))
    scale_data <- scale(df, center = TRUE)
    covariates <- scale_data[, -which(names(scale_data) %in% c("status", "time"))] %>% as.matrix()
    status <- df["status"][[1]]
    time <- df["time"][[1]]
    cvfit <- cv.glmnet(covariates,
                       status,
                       nfolds = k,
                       alpha = 0,
                       family = "cox",
                       type.measure = "C")
    result <- cvfit
  }
  return(result)
}
