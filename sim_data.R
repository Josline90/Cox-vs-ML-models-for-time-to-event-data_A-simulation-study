## To ensure we have same main feature values for with/without noise variables
gen_covariates <- function(n){
  
  # generate main covariates (based on Riksstroke 2015-2020 dataset)
  age <- rtruncnorm(n = n, a = 18, b = 108, mean = 74.36, sd = 12.44)
  x1 <- (age - mean(age)) / sd(age)           #1. age
  x2 <- rbinom(n = n, size = 1, prob = 0.55)   #2. sex
  x3 <- rbinom(n = n, size = 1, prob = 0.27)   #3. AF
  x4 <- rbinom(n = n, size = 1, prob = 0.22)   #4. diabetes
  x5 <- rbinom(n = n, size = 1, prob = 0.27)   #5. previous_stroke_tia
  x6 <- rbinom(n = n, size = 1, prob = 0.14)   #6. smoking
  x7 <- rbinom(n = n, size = 1, prob = 0.62)   #7. hypertension
  x8 <- rbinom(n = n, size = 1, prob = 0.32)   #8. lipid
  x9 <- rbinom(n = n, size = 1, prob = 0.13)   #9. prior_anticoagulation
  x10 <- rbinom(n = n, size = 1, prob = 0.25)  #10. rankin_scale_prestroke
  
  # loss of consciousness
  loss <- apply(rmultinom(n = n, size = 1, prob = c(0.91, 0.07, 0.02)), 2, which.max)
  encoded_loss <- data.frame(model.matrix(~ factor(loss, levels = 1:3) - 1))
  x11 <- encoded_loss[, 2] #11. mild/drowsy
  x12 <- encoded_loss[, 3] #12. severe/unconscious #11-12. nihss_loss_of_consciousness
  x13 <- rbinom(n = n, size = 1, prob = 0.04) #13. thrombectomy
  x14 <- rbinom(n = n, size = 1, prob = 0.15)  #14. thrombolysis
  x15 <- rbinom(n = n, size = 1, prob = 0.23)  #15. wake_up_stroke
  x16 <- rbinom(n = n, size = 1, prob = 0.85)  #16. adm at stroke unit
  x17 <- rbinom(n = n, size = 1, prob = 0.39)  #17. stroke_alert
  x18 <- rbinom(n = n, size = 1, prob = 0.73)  #18. ambulance
  
  dat <- data.frame(x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14,
                    x15, x16, x17, x18)
  return(covariates = dat)
}

#-----------------------------------------------------------------------------#
##  A function to generate the dataset with predefined censoring proportion   #
#-----------------------------------------------------------------------------#
simulate_data <- function(noise, n, P, betas, lambda){
  
  # linear combination of covariates with coefficients (betas)
  LP <- as.matrix(P) %*% betas
  
  # Generate survival times
  #i). create hazard function for the mortality
  h_f <- lambda*exp(LP) #rate = lambda
  
  #ii). Create survival time based on Bender formula (2005) for times T ~ Exponential(baseline_hazard)
  u <- runif(n)
  time <- (-(log(u))) / h_f
  
  # Generate censoring times (unrelated to survival time)
  cens <-  rexp(n = n, rate = 2.5) # set for a specific cens level
  
  # Create event indicator 
  status <- ifelse(time <= cens, 1, 0) 
  time <- pmin(time, cens) # observed follow-up time
  
  # For 0% censoring level 
  #status <- rep(1, n)
  
  # return the data frame of the simulated dataset
  if(noise == TRUE){
    
    # e.g 500 noise variables
    n_df <- replicate(150, rnorm(n = n, mean = 0, sd = 1)) 
    b_df <- replicate(350, rbinom(n = n, size = 1, prob = 0.5)) 
    colnames(n_df) <- paste0("x", 19:168) 
    colnames(b_df) <- paste0("x", 169:518) 
    
    data_sim1 <- data.frame(time, status, P, n_df, b_df)
    return(Simulated_dataset = data_sim1)
  } else{
    data_sim2 <- data.frame(time, status, P)
    return(Simulated_dataset = data_sim2)
  }
}

#-----------------------------------------------------------------------------#
##  A function to generate the dataset;  for assessing nonlinearity misspecification  #
#-----------------------------------------------------------------------------#
simulate_data_nonlinear <- function(n, P, betas, beta_m, lambda){
  
  # mispecification
  x1 <- P$x1;  m <- x1*x1; m_b <- beta_m
  
  # linear combination of covariates with coefficients (betas)
  LP <- as.matrix(P) %*% betas + m*m_b
  # Generate survival times
  #i). create hazard function for the mortality
  h_f <- lambda*exp(LP) #rate = lambda
  
  #ii). Create survival time based on Bender formula (2005) for times T ~ Exponential(baseline_hazard)
  u <- runif(n)
  time <- (-(log(u))) / h_f
  
  # Generate censoring times (unrelated to survival time)
  cens <-  rexp(n = n, rate = 2.7) # 85% cens. can be changed
  
  # Create event indicator 
  status <- ifelse(time <= cens, 1, 0) 
  time <- pmin(time, cens) # observed follow-up time
  
  #status <- rep(1, n) # for 0% cens level
  
  x19 <- m
  
  data_sim <- data.frame(time, status, P, x19)
  return(Simulated_dataset = data_sim)
  
}

#-----------------------------------------------------------------------------#
##  A function to generate the dataset; for assessing mispecified PH assumption   #
#-----------------------------------------------------------------------------#

simulate_data_nonPH <- function(n, P, betas, k, beta_t, t0, lambda){ 
  
  # linear combination of covariates with coefficients (betas) bf varying covariate
  LP <- as.matrix(P) %*% betas
  # Generate survival times
  #i). create hazard function for the mortality
  h_f <- lambda*exp(LP) 
  #ii). Create survival time based on Bender formula (2005) for times T ~ Exponential(baseline_hazard)
  u <- runif(n)
  time1 <- (-(log(u))) / h_f
  
  # after varying covariate
  betas[k] <- beta_t #set the time varying covariate effect to positive
  LP_new <- as.matrix(P) %*% betas
  h_f2 <- lambda*exp(LP_new) 
  time2 <- (-(log(u))) / h_f2
  
  # Times
  time <- ifelse(time1 > t0, t0 + time2, time1)
  
  # Generate censoring times (unrelated to survival time)
  cens <-  rexp(n = n, rate = 3) 
  # Create event indicator 
  status <- ifelse(time <= cens, 1, 0)
  time <- pmin(time, cens) # observed follow-up time
  
  #status <- rep(1, n) # for 0% cens level
  
  data_sim <- data.frame(time, status, P)
  return(Simulated_dataset = data_sim)
}

#----------------------------------------------------------------------#
#   A function to set seed for split index: createDataPartition fxn    #
#----------------------------------------------------------------------#
splitData <- function(sim_data, split_size = 0.7, seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)  
  }
  x <- createDataPartition(sim_data$status, p = split_size, list = FALSE, times = 1)
  #x <- createDataPartition(sim_data$time, p = split_size, list = FALSE, times = 1) # for 0% cens level
  
  # train and test data
  test_data <- sim_data[-x, ]
  train_data <- sim_data[x, ]
  
  return(list(tr_data = train_data, te_data = test_data, x = x))
}
