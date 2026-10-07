#------------------------------------------------------------------------------------------#
#          A function to fit correct Cox-assumption and calculate metrics : IBS/Cstats    #
#-------------------------------------------------------------------------------------------#
# used this to report Harrel's with counting format only
cox_split <- function(train_df, test_df){
  
  test.split <- survSplit(Surv(time, status)~., data = test_df, cut = 0.25, id="id", episode= "tgroup")
  train_dtimes <- sort(unique(with(train_df, time[status==1])))
  
  train.split <-survSplit(Surv(time, status)~., data = train_df, cut = 0.25, id = "id", episode= "tgroup")
  
  # ----------1. using tt
  test.split0 <- test.split
  test.split0$ttk <- test.split0$x3 * (test.split0$tstart >= 0.25)
  train.split0 <- train.split
  train.split0$ttk <- train.split0$x3 * (train.split0$tstart >= 0.25)
  cph_model.cstat0 <- coxph(Surv(tstart, time, status) ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10 + x11 + 
                              x12 + x13 + x14 + x15 + x16 + x17 + x18 + ttk, data = train.split0, x = TRUE, y = TRUE)
  cph_survfit0 <- summary(survfit(formula = cph_model.cstat0, newdata = test.split0), times = train_dtimes)
  cph_predict.surv0 <- t(cph_survfit0$surv)
  cph_predict.risk0 <- predict(cph_model.cstat0, newdata = test.split0, type = "risk")
  
  prob = .95
  z <- qnorm((1-prob)/2)
  
  # Harrel's with counting format + predictions on test.split + metric on test_df
  cox_index.Hctt <- concordance(object = Surv(time = test.split0$tstart, time2 = test.split0$time, event = test.split0$status) ~ cph_predict.risk0, data = test.split0, reverse = TRUE)
  cox_c_CI.Hctt <- cox_index.Hctt$concordance + c(0, z, -z) * sqrt(cox_index.Hctt$var)
  
  #using duration
  train.split$time_duration <- train.split$time- train.split$tstart
  test.split$time_duration <- test.split$time- test.split$tstart
  cph_model.cstat <- coxph(Surv(time_duration, status) ~ x1 + x2 + x3:strata(tgroup) + x4 + x5 + x6 + x7 + x8 + x9 + x10 + x11 + 
                             x12 + x13 + x14 + x15 + x16 + x17 + x18, data = train.split, x = TRUE, y = TRUE)
  
  cph_predict.risk <- predict(cph_model.cstat, newdata = test.split, type = "risk")
  
  # Harrel's with counting format + strata()
  cox_index.Hcstr <- concordance(object = Surv(time = tstart, time2 = time, event = status) ~ cph_predict.risk, data= test.split, reverse = TRUE)
  cox_c_CI.Hcstr <- cox_index.Hcstr$concordance + c(0, z, -z) * sqrt(cox_index.Hcstr$var)
  
  #-using duration (collapsed with weighted risks) + predictions on test.split + metric on test_df
  test.split$risk_score <- cph_predict.risk
  test.split$weighted_risk <- test.split$risk_score * test.split$time_duration
  collapsed_data <- aggregate(cbind(weighted_risk, time_duration) ~ id, data = test.split, FUN = sum)
  
  weighted_risk_score <- collapsed_data$weighted_risk / collapsed_data$time_duration
  test_df$time_duration <- collapsed_data$time_duration
  
  cox_index.H <- concordance(object = Surv(time = time_duration, event = status) ~ weighted_risk_score, data = test_df, reverse = TRUE)
  cox_c_CI.H <- cox_index.H$concordance + c(0, z, -z) * sqrt(cox_index.H$var)
  
  cox_index.U <- concordance(object = Surv(time = time_duration, event = status) ~ weighted_risk_score, data = test_df, reverse = TRUE, timewt= "n/G2")
  cox_c_CI.U <- cox_index.U$concordance + c(0, z, -z) * sqrt(cox_index.U$var)
  
  #combine surv probs
  surv_df0 <- as.data.frame(cph_predict.surv0)
  combine_test0 <- cbind(test.split0, surv_df0)
  
  columns_to_multiply0 <- grep("V", colnames(combine_test0), value = TRUE) #surv prob. columns
  
  # cum. surv. prob per ID account for time of death
  agg_probs0 <- combine_test0 %>%
    arrange(id, tstart) %>%
    group_by(id) %>%
    summarise(across(all_of(columns_to_multiply0),
                     ~ ifelse(any(status == 1), .[which(status == 1)[1]], #if event occurs, take prob up to first event
                              prod(.)), #otherwise take both (censored)
                     .names = "surv_prob_{col}")) %>% 
    ungroup() %>% 
    select(-id) %>% as.matrix()#as.data.frame()
  
  
  times <- test_df$time_duration
  status2 <- test_df$status
  IBS_cox <- my_sbrier(Surv(times, status2), pred = t(agg_probs0), btime = train_dtimes)[1,1] 
  #----
  # Output
  models <- c("Cox PH")
  c_index.Hctt <- c(cox_c_CI.Hctt[1])
  lower_ci.Hctt <- c(cox_c_CI.Hctt[2])
  upper_ci.Hctt <- c(cox_c_CI.Hctt[3])
  
  c_index.Hcstr <- c(cox_c_CI.Hcstr[1])
  lower_ci.Hcstr <- c(cox_c_CI.Hcstr[2])
  upper_ci.Hcstr <- c(cox_c_CI.Hcstr[3])
  
  c_index.H <- c(cox_c_CI.H[1])
  lower_ci.H <- c(cox_c_CI.H[2])
  upper_ci.H <- c(cox_c_CI.H[3])
  
  c_index.U <- c(cox_c_CI.U[1])
  lower_ci.U <- c(cox_c_CI.U[2])
  upper_ci.U <- c(cox_c_CI.U[3])
  
  IBS_val <- c(IBS_cox)
  
  d.frame <- data.frame(Model = models,
                        C_index_Hctt = c_index.Hctt, Lower_CI_Hctt = lower_ci.Hctt, Upper_CI_Hctt = upper_ci.Hctt,
                        C_index_Hcstr = c_index.Hcstr, Lower_CI_Hcstr = lower_ci.Hcstr, Upper_CI_Hcstr = upper_ci.Hcstr,
                        C_index_H = c_index.H, Lower_CI_H = lower_ci.H, Upper_CI_H = upper_ci.H,
                        C_index_U = c_index.U, Lower_CI_U = lower_ci.U, Upper_CI_U = upper_ci.U, 
                        IBScore = IBS_val
  )
  
  return(d.frame)
}