#------------------------------------------------------------------------------------------#
#          A function to fit the models and calculate metrics : IBS/Cstats    #
#-------------------------------------------------------------------------------------------#

models_metric_IBS <- function(train_df, test_df){
  
  #fit models
  fmla <- Surv(time = time, event = status) ~ .
  
  # 1 RSF model: 'logrank'
  rsf_model <- rfsrc(formula = fmla, data = train_df, splitrule = "logrank", 
                     mtry = as.integer(sqrt(ncol(train_df[, !names(train_df) %in% c("time", "status")]))))
  rsf_predict <- get_predictions(fitted_model = rsf_model, newdata = test_df) # times = rsf_model$time.interest
  
  # 2. Cox PH model
  cph_model <- coxph(formula = fmla, data = train_df, x=T, y=T)
  cph_predict <- get_predictions(fitted_model = cph_model, newdata = test_df, times = rsf_model$time.interest) 
  
  # 3. Survival neural network model
  nn_model <- deepsurv(fmla, data = train_df, frac = 0)
  nn_predict <- get_predictions(fitted_model = nn_model, newdata = test_df)
  
  # 4. XGBoost for survival analysis
  # Prepare the training and validation data for XGBoost
  x_label <- ifelse(train_df$status == 1, train_df$time, -train_df$time)
  x_train <- as.matrix(train_df %>% select(-time, -status))
  Dtrain <- xgb.DMatrix(x_train, label = x_label)
  
  surv_xgboost_model <- xgb.train.surv(
    params = list(objective = "survival:cox", eval_metric = "cox-nloglik"),
    data = x_train, label = x_label, nrounds = 100)
  xgb_predict <- get_predictions(fitted_model = surv_xgboost_model, newdata = test_df, times = rsf_model$time.interest)
  
  times <- test_df$time 
  status <- test_df$status
  
  # 1. RSF model
  IBS_rsf <- my_sbrier(Surv(times, status), pred = t(rsf_predict$surv), btime = rsf_predict$times)[1,1] 
  
  # 2. Cox PH model
  IBS_cox <- my_sbrier(Surv(times, status), pred = t(cph_predict$surv), btime = cph_predict$times)[1,1] 
  
  # 3. Survival neural network model
  IBS_nn <- my_sbrier(Surv(times, status), pred = t(nn_predict$surv), btime = nn_predict$times)[1,1] 
  
  # 4. XGBoost for survival analysis
  IBS_xgb <- my_sbrier(Surv(times, status), pred = t(xgb_predict$surv), btime = xgb_predict$times)[1,1] 
  
  # }else{
  prob = .95
  z <- qnorm((1-prob)/2)
  
  # Harrell's Cstats
  
  # 1. Cox PH model
  cox_index.H <- concordance(object = Surv(time = times, event = status) ~ cph_predict$risk, data=test_df, reverse = TRUE)
  cox_c_CI.H <- cox_index.H$concordance + c(0, z, -z) * sqrt(cox_index.H$var)
  
  # 2. RSF model
  rsf_index.H <- concordance(object = Surv(time = times, event = status) ~ rsf_predict$risk, data=test_df, reverse = TRUE)
  rsf_c_CI.H <- rsf_index.H$concordance + c(0, z, -z) * sqrt(rsf_index.H$var)
  
  # 3. XGBoost for survival analysis
  xgb_index.H <- concordance(object = Surv(time = times, event = status) ~ xgb_predict$risk, data=test_df, reverse = TRUE)
  xgb_c_CI.H <- xgb_index.H$concordance + c(0, z, -z) * sqrt(xgb_index.H$var)
  
  # 4. Survival neural network model
  nn_index.H <- concordance(object = Surv(time = times, event = status) ~ nn_predict$risk, data = test_df, reverse = TRUE)
  nn_c_CI.H <- nn_index.H$concordance + c(0, z, -z) * sqrt(nn_index.H$var)
  
  
  #Uno's Cstats
  # 1. Cox PH model
  cox_index.U <- concordance(object = Surv(time = times, event = status) ~ cph_predict$risk, data=test_df, reverse = TRUE,timewt= "n/G2")
  cox_c_CI.U <- cox_index.U$concordance + c(0, z, -z) * sqrt(cox_index.U$var)
  
  # 2. RSF model
  rsf_index.U <- concordance(object = Surv(time = times, event = status) ~ rsf_predict$risk, data=test_df, reverse = TRUE,timewt= "n/G2")
  rsf_c_CI.U <- rsf_index.U$concordance + c(0, z, -z) * sqrt(rsf_index.U$var)
  
  # 3. XGBoost for survival analysis
  xgb_index.U <- concordance(object = Surv(time = times, event = status) ~ xgb_predict$risk, data=test_df, reverse = TRUE,timewt= "n/G2")
  xgb_c_CI.U <- xgb_index.U$concordance + c(0, z, -z) * sqrt(xgb_index.U$var)
  
  # 4. Survival neural network model
  nn_index.U <- concordance(object = Surv(time = times, event = status) ~ nn_predict$risk, data = test_df, reverse = TRUE,timewt= "n/G2")
  nn_c_CI.U <- nn_index.U$concordance + c(0, z, -z) * sqrt(nn_index.U$var)
  
  
  # Output
  models <- c("XGBoost", "RSF-logrank", "Cox PH", "NN")
  c_index.H <- c(xgb_c_CI.H[1], rsf_c_CI.H[1], cox_c_CI.H[1], nn_c_CI.H[1])
  lower_ci.H <- c(xgb_c_CI.H[2], rsf_c_CI.H[2], cox_c_CI.H[2], nn_c_CI.H[2])
  upper_ci.H <- c(xgb_c_CI.H[3], rsf_c_CI.H[3], cox_c_CI.H[3], nn_c_CI.H[3])
  
  c_index.U <- c(xgb_c_CI.U[1], rsf_c_CI.U[1], cox_c_CI.U[1], nn_c_CI.U[1])
  lower_ci.U <- c(xgb_c_CI.U[2], rsf_c_CI.U[2], cox_c_CI.U[2], nn_c_CI.U[2])
  upper_ci.U <- c(xgb_c_CI.U[3], rsf_c_CI.U[3], cox_c_CI.U[3], nn_c_CI.U[3])
  
  IBS_val <- c(IBS_xgb, IBS_rsf, IBS_cox, IBS_nn)
  
  d.frame <- data.frame(Model = models,
                        C_index_H = c_index.H, Lower_CI_H = lower_ci.H, Upper_CI_H = upper_ci.H,
                        C_index_U = c_index.U, Lower_CI_U = lower_ci.U, Upper_CI_U = upper_ci.U, 
                        IBScore = IBS_val)
  
  return(d.frame)
  
  
}