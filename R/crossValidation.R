#' @title Cross Validation
#'
#' @description A cross validation procedure with a set amount of repeats towards
#' the goal of optimizing cox proportional hazards models, lasso models, and ridge
#' models in the context of features from clinical data.
#'
#' @param clinical_data A long-format \code{data.frame} or \code{tibble} that has
#' features represented by columns. The function additionally assumes the existence
#' of a patient ID column, marked by the \code{id} argument in \code{crossValidation()}.
#' @param k The number of folds in the k-fold cross validation procedure
#' @param n The number of repeats of the k-fold cross validation procedure
#' @param method A \code{character} indicating which method will be used in the
#' cross validation procedure. The default is set to 'coxph', indicating the
#' cox proportional hazards model. The other available options are ridge and lasso
#' regression.
#' @param id A \code{character} indicating the name of the patient ID column. The
#' default is assumed to be 'PatID'.
#' @param demo A \code{logical} that indicates whether or not \code{crossValidation()}
#' should print out a statement indicating whether each repeat uses disjoint intervals.
#' Included solely for debugging purposes.
#'
#' @details \code{crossValidation()} runs k-fold cross validation using the specified
#' model from the argument \code{method} on the clinical data that is provided. In
#' particular, for each repeat, the function splits the data into k folds. The model
#' is first fit according to the data in a single fold, and is then tested against the
#' remaining k-1 folds. The function uses concordance index as a measure of error between
#' the predicted and actual time-status pairs. We then repeat this process k-1 more times
#' for each of the other folds, and then repeat this overall process n-1 more times for
#' each repeat.
#'
#' @return A list with n elements, one for each repeat. Each element is another
#' list with the following components:
#' \describe{
#' \item{c.indices}{A collection of k concordance indices from that particular repeat}
#' \item{lowers}{A collection of lower bounds on the confidence intervals for the concordance indices}
#' \item{uppers}{A collection of upper bounds on the confidence intervals for the concordance indices}
#' \item{sd}{The standard deviation of c.indices}
#' \item{mean}{The mean of c.indices}
#' }
#'
#' @import ISLR
#' @import caret
#' @import survival
#' @import survcomp
#' @import glmnet
#' @import fastDummies
#' @import purrr
#' @import randomForestSRC
#'
#' @export

crossValidation <- function(clinical_data,
                            k = 5,
                            n = 5,
                            method = "coxph",
                            id = "PatID",
                            demo = FALSE) {
  if (method != "coxph" & method != "lasso" & method != "ridge" & method != "rf") {
    stop("method must either be coxph, lasso, ridge, or rf.")
  }
  else{
    num_rows_with_na <- as.character(sum(apply(is.na(clinical_data), 1, any)))
    num_rows <- as.character(nrow(clinical_data))
    message <- paste0("Removed ", num_rows_with_na, " out of ", num_rows, " rows.")
    message(message)
    df <- clinical_data[, !names(clinical_data) %in% id] %>% na.omit()
    proportion <- 1-1/k
    if (method == "coxph") {
      has_character_column <- df %>%
        select_if(is.character) %>%
        map_lgl(~ any(!is.na(.)))
      if (any(has_character_column)) {
        df <- df %>% dummy_cols(remove_selected_columns = TRUE)
        indicator <- c(1:nrow(df))
        newdf <- cbind(indicator, df)
      }
      else {
        df <- df
        indicator <- c(1:nrow(df))
        newdf <- cbind(indicator, df)
      }
      cindex_list <- list()
      for (j in 1:n) {
        cvector <- rep(NA, k)
        lower <- rep(NA, k)
        upper <- rep(NA, k)
        folds <- createFolds(newdf$indicator, k = k)
        if (demo == TRUE){
          if (all(sort(unlist(folds)) == c(1:nrow(newdf)))) {
            message("Folds are disjoint.")
          }
          else {
            stop("Folds are not disjoint.")
          }
        }
        for (i in 1:k) {
          newtrain <- newdf[-folds[[i]], ]
          newtest <- newdf[folds[[i]], ]
          train <- newtrain[, -which(names(newtrain) == "indicator")]
          test <- newtest[, -which(names(newtest) == "indicator")]
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
        message(paste0("Repeat ", j, " completed"))
      }
      result <- cindex_list
    }
    else if (method == "lasso") {
      has_character_column <- df %>%
        select_if(is.character) %>%
        map_lgl(~ any(!is.na(.)))
      if (any(has_character_column)) {
        df <- df %>% dummy_cols(remove_selected_columns = TRUE)
        indicator <- c(1:nrow(df))
        newdf <- cbind(indicator, df)
      }
      else {
        df <- df
        indicator <- c(1:nrow(df))
        newdf <- cbind(indicator, df)
      }
      status_time_col <- newdf[c("time", "status")]
      rest_of_df <- newdf[, -which(names(newdf) %in% c("status", "time"))]
      scale_data <- as.data.frame(scale(rest_of_df, center = TRUE))
      scale_data_full <- cbind(scale_data, status_time_col)
      cindex_list <- list()
      for (j in 1:n) {
        cvector <- rep(NA, k)
        lower <- rep(NA, k)
        upper <- rep(NA, k)
        folds <- createFolds(scale_data_full$indicator, k = k)
        if (demo == TRUE){
          if (all(sort(unlist(folds)) == c(1:nrow(newdf)))) {
            message("Folds are disjoint.")
          }
          else {
            stop("Folds are not disjoint.")
          }
        }
        for (i in 1:k) {
          newtrain <- scale_data_full[-folds[[i]], ]
          newtest <- scale_data_full[folds[[i]], ]
          train <- newtrain[, -which(names(newtrain) == "indicator")]
          test <- newtest[, -which(names(newtest) == "indicator")]
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
          lambda <- cvfit$lambda.1se
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
        message(paste0("Repeat ", j, " completed"))
      }
      result <- cindex_list
    }
    else if (method == "ridge") {
      has_character_column <- df %>%
        select_if(is.character) %>%
        map_lgl(~ any(!is.na(.)))
      if (any(has_character_column)) {
        df <- df %>% dummy_cols(remove_selected_columns = TRUE)
        indicator <- c(1:nrow(df))
        newdf <- cbind(indicator, df)
      }
      else {
        df <- df
        indicator <- c(1:nrow(df))
        newdf <- cbind(indicator, df)
      }
      status_time_col <- newdf[c("time", "status")]
      rest_of_df <- newdf[, -which(names(newdf) %in% c("status", "time"))]
      scale_data <- as.data.frame(scale(rest_of_df, center = TRUE))
      scale_data_full <- cbind(scale_data, status_time_col)
      cindex_list <- list()
      for (j in 1:n) {
        cvector <- rep(NA, k)
        lower <- rep(NA, k)
        upper <- rep(NA, k)
        folds <- createFolds(scale_data_full$indicator, k = k)
        if (demo == TRUE){
          if (all(sort(unlist(folds)) == c(1:nrow(newdf)))) {
            message("Folds are disjoint.")
          }
          else {
            stop("Folds are not disjoint.")
          }
        }
        for (i in 1:k) {
          newtrain <- scale_data_full[-folds[[i]], ]
          newtest <- scale_data_full[folds[[i]], ]
          train <- newtrain[, -which(names(newtrain) == "indicator")]
          test <- newtest[, -which(names(newtest) == "indicator")]
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
        message(paste0("Repeat ", j, " completed"))
      }
      result <- cindex_list
    }
    else if (method == "rf") {
      has_character_column <- df %>%
        select_if(is.character) %>%
        map_lgl(~ any(!is.na(.)))
      if (any(has_character_column)) {
        df <- df %>% dummy_cols(remove_selected_columns = TRUE)
        indicator <- c(1:nrow(df))
        newdf <- cbind(indicator, df)
      }
      else {
        df <- df
        indicator <- c(1:nrow(df))
        newdf <- cbind(indicator, df)
      }
      newdf$status <- as.numeric(as.factor(as.character(newdf$status))) - 1
      cindex_list <- list()
      for (j in 1:n) {
        cvector <- rep(NA, k)
        lower <- rep(NA, k)
        upper <- rep(NA, k)
        folds <- createFolds(newdf$indicator, k = k)
        if (demo == TRUE){
          if (all(sort(unlist(folds)) == c(1:nrow(newdf)))) {
            message("Folds are disjoint.")
          }
          else {
            stop("Folds are not disjoint.")
          }
        }
        for (i in 1:k) {
          newtrain <- newdf[-folds[[i]], ]
          newtest <- newdf[folds[[i]], ]
          train <- newtrain[, -which(names(newtrain) == "indicator")]
          test <- newtest[, -which(names(newtest) == "indicator")]
          test_covariates <- test[, -which(colnames(test) %in% c("status", "time"))] %>% as.data.frame()
          rfsrc_model <- rfsrc(Surv(time, status) ~., data = train)
          test_pred <- predict(rfsrc_model, test_covariates)
          predicted_values <- test_pred$predicted
          concordance <- concordance.index(predicted_values, test$time, test$status)
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
        message(paste0("Repeat ", j, " completed"))
      }
      result <- cindex_list
    }
    return(result)
  }
}
