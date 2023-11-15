#' @import ISLR

crossValidation <- function(clinical_data,
                            k = 5,
                            method = "coxph",
                            proportion = 0.8) {

  df <- clinical_data[-1]

  for (i in 1:k) {
    inTrain <- createDataPartition(y = clinical_data$time,
                                   p = proportion,
                                   list = FALSE)
    train <- df[inTrain, ]
    test <- df[-inTrain, ]
    coxph_model <- coxph(Surv(time, status) ~., data = train)
    test_pred <- coxph_model %>% predict(test)
    # do something for concordance index
  }

}
