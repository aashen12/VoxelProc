#' @import ISLR
#' @import caret
#' @import survival
#' @import survcomp

crossValidation <- function(clinical_data,
                            k = 5,
                            method = "coxph") {

  df <- clinical_data[-1]
  proportion <- nrow(clinical_data)/k

  for (i in 1:k) {
    inTrain <- createDataPartition(y = clinical_data$time,
                                   p = proportion,
                                   list = FALSE)
    train <- df[inTrain, ]
    test <- df[-inTrain, ]
    coxph_model <- coxph(Surv(time, status) ~., data = train)
    test_pred <- coxph_model %>% predict(test)
    cindex <- concordance.index(test_pred, test$time, test$status)
  }

}
