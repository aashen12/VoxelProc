#' @import ISLR
#' @import caret
#' @import survival
#' @import survcomp
#' @import glmnet
#' @import fastDummies
#' @import purrr

crossValidation <- function(clinical_data,
                            k = 5,
                            n = 5,
                            method = "coxph",
                            id = "PatID") {
  if (method != "coxph" & method != "lasso" & method != "ridge") {
    stop("method must either be coxph, lasso, or ridge")
  }
  else{
    df <- clinical_data[, !names(clinical_data) %in% id] %>% na.omit()
    proportion <- 1-1/k
    if (method == "coxph") {
      has_character_column <- df %>%
        select_if(is.character) %>%
        map_lgl(~ any(!is.na(.)))
      if (any(has_character_column)) {
        df <- df %>% dummy_cols(remove_selected_columns = TRUE)
      }
      else {
        df <- df
      }
      cindex_list <- list()
      for (j in 1:n) {
        cvector <- rep(NA, k)
        lower <- rep(NA, k)
        upper <- rep(NA, k)
        for (i in 1:k) {
          inTrain <- createDataPartition(y = df$time, #figure out how to create actual partitions
                                        p = proportion,
                                        list = FALSE)
          train <- df[inTrain, ]
          test <- df[-inTrain, ]
          coxph_model <- coxph(Surv(time, status) ~., data = train)
          test_pred <- coxph_model %>% predict(test)
          concordance <- concordance.index(test_pred, test$time, test$status)
          cindex <- concordance$c.index
          lowerindex <- concordance$lower
          upperindex <- concordance$upper
          cvector[i] <- cindex
          lower[i] <- lowerindex
          upper[i] <- upperindex
        }
        clist <- list(c.indices = cvector,
                      lowers = lower,
                      uppers = upper,
                      sd = sd(cvector),
                      mean = mean(cvector))
        cindex_list[[j]] <- clist
      }
      result <- cindex_list
    }
    else if (method == "lasso") {
      has_character_column <- df %>%
        select_if(is.character) %>%
        map_lgl(~ any(!is.na(.)))
      if (any(has_character_column)) {
        df <- df %>% dummy_cols(remove_selected_columns = TRUE)
      }
      else {
        df <- df
      }
      status_time_col <- df[c("time", "status")]
      rest_of_df <- df[, -which(names(df) %in% c("status", "time"))]
      scale_data <- as.data.frame(scale(rest_of_df, center = TRUE))
      scale_data_full <- cbind(scale_data, status_time_col)
      cindex_list <- list()
      for (j in 1:n) {
        cvector <- rep(NA, k)
        lower <- rep(NA, k)
        upper <- rep(NA, k)
        for (i in 1:k) {
          inTrain <- createDataPartition(y = df$time, # figure out how to create actual partitions
                                        p = proportion,
                                        list = FALSE)
          train <- scale_data_full[inTrain, ]
          test <- scale_data_full[-inTrain, ]
          test_covariates <- test[, -which(colnames(test) %in% c("status", "time"))] %>% as.matrix()
          train_covariates <- train[, -which(colnames(train) %in% c("status", "time"))] %>% as.matrix()
          time <- train["time"][[1]]
          status <- train["status"][[1]]
          cvfit <- cv.glmnet(train_covariates,
                            Surv(time, status),
                            nfolds = k,
                            alpha = 1,
                            family = "cox",
                            type.measure = "C")
          lambda <- cvfit$lambda.min
          test_pred <- predict(cvfit, test_covariates, s = lambda)
          concordance <- concordance.index(test_pred, test$time, test$status)
          cindex <- concordance$c.index
          lowerindex <- concordance$lower
          upperindex <- concordance$upper
          cvector[i] <- cindex
          lower[i] <- lowerindex
          upper[i] <- upperindex
        }
        clist <- list(c.indices = cvector,
                      lowers = lower,
                      uppers = upper,
                      sd = sd(cvector),
                      mean = mean(cvector))
        cindex_list[[j]] <- clist
      }
      result <- cindex_list
    }
    else if (method == "ridge") {
      has_character_column <- df %>%
        select_if(is.character) %>%
        map_lgl(~ any(!is.na(.)))
      if (any(has_character_column)) {
        df <- df %>% dummy_cols(remove_selected_columns = TRUE)
      }
      else {
        df <- df
      }
      status_time_col <- df[c("time", "status")]
      rest_of_df <- df[, -which(names(df) %in% c("status", "time"))]
      scale_data <- as.data.frame(scale(rest_of_df, center = TRUE))
      scale_data_full <- cbind(scale_data, status_time_col)
      cindex_list <- list()
      for (j in 1:n) {
        cvector <- rep(NA, k)
        lower <- rep(NA, k)
        upper <- rep(NA, k)
        for (i in 1:k) {
          inTrain <- createDataPartition(y = df$time, #figure out how to create actual partitions
                                        p = proportion,
                                        list = FALSE)
          train <- scale_data_full[inTrain, ]
          test <- scale_data_full[-inTrain, ]
          test_covariates <- test[, -which(colnames(test) %in% c("status", "time"))] %>% as.matrix()
          train_covariates <- train[, -which(colnames(train) %in% c("status", "time"))] %>% as.matrix()
          time <- train["time"][[1]]
          status <- train["status"][[1]]
          cvfit <- cv.glmnet(train_covariates,
                            Surv(time, status),
                            nfolds = k,
                            alpha = 0,
                            family = "cox",
                            type.measure = "C")
          lambda <- cvfit$lambda.min
          test_pred <- predict(cvfit, test_covariates, s = lambda)
          concordance <- concordance.index(test_pred, test$time, test$status)
          cindex <- concordance$c.index
          lowerindex <- concordance$lower
          upperindex <- concordance$upper
          cvector[i] <- cindex
          lower[i] <- lowerindex
          upper[i] <- upperindex
        }
        clist <- list(c.indices = cvector,
                      lowers = lower,
                      uppers = upper,
                      sd = sd(cvector),
                      mean = mean(cvector))
        cindex_list[[j]] <- clist
      }
      result <- cindex_list
    }
    return(result)
  }
}
