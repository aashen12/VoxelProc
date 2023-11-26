#' @import ISLR
#' @import caret
#' @import survival
#' @import survcomp
#' @import glmnet

crossValidation <- function(clinical_data,
                            k = 5,
                            n = 5,
                            method = "coxph",
                            id = "pid") {
  df <- clinical_data[, !names(clinical_data) %in% id]
  proportion <- 1-1/k
  if (method == "coxph") {
    scores <- rep(NA, n)
    for (j in 1:n) {
      cvector <- rep(NA, k)
      for (i in 1:k) {
        inTrain <- createDataPartition(y = df$time, #figure out how to create actual partitions
                                      p = proportion,
                                      list = FALSE)
        train <- df[inTrain, ]
        test <- df[-inTrain, ]
        coxph_model <- coxph(Surv(time, status) ~., data = train)
        test_pred <- coxph_model %>% predict(test)
        cindex <- concordance.index(test_pred, test$time, test$status)$c.index
        cvector[i] <- cindex
      }
      score <- mean(cvector)
      scores[j] <- score
    }
    result <- list(scores = scores,
                   mean = mean(scores),
                   sd = sd(scores))
  }
  else if (method == "lasso") {
    df <- df %>% select(where(is.numeric))
    status_time_col <- df[c("time", "status")]
    rest_of_df <- df[, -which(names(df) %in% c("status", "time"))]
    scale_data <- as.data.frame(scale(rest_of_df, center = TRUE))
    scale_data_full <- cbind(scale_data, status_time_col)
    for (j in 1:n) {
      for (i in 1:k) {
        inTrain <- createDataPartition(y = df$time, #figure out how to create actual partitions
                                       p = proportion,
                                       list = FALSE)
        train <- scale_data_full[inTrain, ]
        test <- scale_data_full[-inTrain, ]
        test_covariates <- test[, -which(colnames(test) %in% c("status", "time"))] %>% as.matrix()
        time <- train["time"][[1]]
        status <- train["status"][[1]]
        covariates <- train[, -which(colnames(train) %in% c("status", "time"))] %>% as.matrix()
        cvfit <- cv.glmnet(covariates,
                       Surv(time, status),
                       nfolds = k,
                       alpha = 1,
                       family = "cox",
                       type.measure = "C")
        lambda <- cvfit$lambda.min
        trained_model <- glmnet(covariates,
                                Surv(time, status),
                                family = "cox",
                                alpha = 1,
                                lambda = lambda)
        test_pred <- trained_model %>% predict(test_covariates)

      }
    }
  }
  else if (method == "Ridge") {
    df <- df %>% select(where(is.numeric))
    scale_data <- as.data.frame(scale(df, center = TRUE))
    covariates <- scale_data[, -which(names(scale_data) %in% c("status", "time"))] %>% as.matrix()
    status <- df["status"][[1]]
    time <- df["time"][[1]]
    cvfit <- cv.glmnet(covariates,
                       Surv(time, status),
                       nfolds = k,
                       alpha = 0,
                       family = "cox",
                       type.measure = "C")
    result <- list(score = mean(cvfit$cvm),
                   concordances = cvfit$cvm,
                   sd = sd(cvfit$cvm))
  }
  return(result)
}
