#' @import ISLR
#' @import caret
#' @import survival
#' @import survcomp

crossValidation <- function(clinical_data,
                            k = 5,
                            method = "coxph") {

  df <- clinical_data[-1]
  proportion <- 1-1/k

  if (method == "coxph") {
    cvector <- c()
    for (i in 1:k) {
      inTrain <- createDataPartition(y = df$time,
                                    p = proportion,
                                    list = FALSE)
      train <- df[inTrain, ]
      test <- df[-inTrain, ]
      coxph_model <- coxph(Surv(time, status) ~., data = train)
      test_pred <- coxph_model %>% predict(test)
      cindex <- concordance.index(test_pred, test$time, test$status)
      cvector <- c(cvector, cindex)
    }
    result <- cvector()
  }
  else if (method == "LASSO") {
    scale_data <- scale(df, center = TRUE)
    covariates.w.time <- scale_data[, -which(names(scale_data) %in% c("status"))]
    status <- scale_data["status"]
    cvfit <- cv.glmnet(covariates.w.time,
                       status,
                       nfolds = k,
                       alpha = 1,
                       family = "cox",
                       type.measure = "C")
    cverror <- cvfit$cvm
    result <- cverror
  }
  else if (method == "Ridge") {
    scale_data <- scale(df, center = TRUE)
    covariates.w.time <- scale_data[, -which(names(scale_data) %in% c("status"))]
    status <- scale_data["status"]
    cvfit <- cv.glmnet(covariates.w.time,
                       status,
                       nfolds = k,
                       alpha = 0,
                       family = "cox",
                       type.measure = "C")
    cverror <- cvfit$cvm
    result <- cverror
  }
  return(result)
}
