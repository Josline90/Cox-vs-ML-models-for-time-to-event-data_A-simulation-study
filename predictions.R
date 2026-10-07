#----------------------------------------------------------------------#
#   A function to obtain predictions: risks scores and surv probs    #
#----------------------------------------------------------------------#

get_predictions <- function(fitted_model, newdata, times = NULL){
  model <- class(fitted_model)[1]
  if(model == "rfsrc"){
    rsf_predict <- predict(fitted_model, newdata = newdata)
    surv <- rsf_predict$survival 
    risk <- rsf_predict$predicted
    times <- fitted_model$time.interest
  }
  if(model == "deepsurv"){
    nn_predict_survprob_all <- predict(object = fitted_model, newdata = newdata, type = "survival")
    #deepsurv predicts survival using all unique times in train_data
    #i.e. both censoring and death times
    surv <- nn_predict_survprob_all
    risk <- predict(object = fitted_model, newdata = newdata, type = "risk")
    times <- as.numeric(colnames(surv))
  }
  if(model == "coxph"){
    cph_survfit <- summary(survfit(fitted_model, newdata = newdata), times = times)
    surv <- t(cph_survfit$surv)
    risk <- predict(fitted_model, newdata = newdata, type = "risk")
  }
  if(model == "xgb.Booster.surv"){
    x_test <- as.matrix(newdata %>% select(-time, -status))
    surv <- predict(object = fitted_model, newdata = x_test, type = "surv", times = times)
    colnames(surv) <- NULL
    risk <- predict(object = fitted_model, newdata = x_test, type = "risk") 
  }
  return(list(surv = surv, risk = risk, times = times))
}