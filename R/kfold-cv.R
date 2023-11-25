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
    cvector <- rep(NA, k)
    for (i in 1:k) {
      inTrain <- createDataPartition(y = df$time,
                                    p = proportion,
                                    list = FALSE)
      train <- df[inTrain, ]
      test <- df[-inTrain, ]
      coxph_model <- coxph(Surv(time, status) ~., data = train)
      test_pred <- coxph_model %>% predict(test)
      cindex <- concordance.index(test_pred, test$time, test$status)$c.index
      cvector[i] <- cindex
    }
    result <- list(score = mean(cvector),
                   concordances = cvector,
                   sd = sd(cvector))
  }
  else if (method == "lasso") {
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
    result <- list(score = mean(cvfit$cvm),
                   concordances = cvfit$cvm,
                   sd = sd(cvfit$cvm))
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
