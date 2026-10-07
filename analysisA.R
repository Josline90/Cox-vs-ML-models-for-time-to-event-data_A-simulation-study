rm(list = ls())
source('libraries.R')
source('sim_data.R')
source('predictions.R')
source('ibs.R')
source('metrics.R')

#---------------------------------------------------------------#
#                     Run the functions                         #
#---------------------------------------------------------------#
coef_riks <- c(0.810, 0.242, 0.311, 0.269, 0.086, 0.269, 0.027,
               -0.060, -0.037, 0.870, 0.968, 1.819, 0.188, -0.201, 
               -0.044, -0.088, 0.094, 0.440) # beta estimates from Cox model on Riksstroke 2015-2020 dataset

#beta_m <- 0.1 #0.1; 0.9 for non-linearity misspecification
inner_seed <- 13
set.seed(47)
n = 500 
valid_df <- 0
i <- 1
cens.p <- numeric()
data.samples <- list()
metric_results <- data.frame()
while (valid_df < 1000) { #number of iterations needed 
  cov <- gen_covariates(n = n)
  sample.data <- simulate_data(noise = FALSE, n = n, P = cov, betas = coef_riks, lambda = 0.125)
  #sample.data <- simulate_data_nonlinear(n = n, P = cov, betas = coef_riks, beta_m = beta_m, lambda = 0.125) #non-linear
  
  old_seed  <- .Random.seed
  data_full <- splitData(sim_data = sample.data, split_size = 0.7, seed = inner_seed)
  .Random.seed <- old_seed
  
  c.indexes <- tryCatch({
    models_metric_IBS(train_df = data_full[[1]], test_df = data_full[[2]])
  }, error = function(e) {
    cat("Error in models_metric on iteration", i, ":", e$message, "\n")
    return(NULL)
  })
  
  # If models_metric fails, skip this iteration
  if (is.null(c.indexes)) {
    i <- i + 1
    next
  }
  cens_p <- sum(sample.data$status == 0) / nrow(sample.data)
  
  cens.p[valid_df + 1] <- cens_p #only for valid sets
  
  data.samples[[valid_df + 1]] <- sample.data
  metric_results <- rbind(metric_results, data.frame(Iteration = i, Sample = n, c.indexes))
  valid_df <- valid_df + 1
  i <- i + 1
}


#summary
avg_metric <- metric_results %>%
  group_by(Model) %>%
  summarise(H_index = mean(C_index_H),
            H_lower_CI = mean(Lower_CI_H),
            H_upper_CI = mean(Upper_CI_H),
            U_index = mean(C_index_U),
            U_lower_CI = mean(Lower_CI_U),
            U_upper_CI = mean(Upper_CI_U),
            IBScore = mean(IBScore)) %>%
  as.data.frame()

print(avg_metric)

mean(cens.p)

