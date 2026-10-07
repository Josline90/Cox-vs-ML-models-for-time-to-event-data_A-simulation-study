rm(list = ls())
source('libraries.R')
source('sim_data.R')
source('predictions.R')
source('ibs.R')
source('metric_correctPH.R')

#---------------------------------------------------------------#
#                     Run the functions                         #
#---------------------------------------------------------------#
coef_X <- c(0.810, 0.242, 0, 0.269, 0.086, 0.269, 0.027,
            -0.060, -0.037, 0.870, 0.968, 1.819, 0.188, -0.201, 
            -0.044, -0.088, 0.094, 0.440) #beta_t => 0.311

k = 3 # covariate effect to vary
tvary <- 0.25 # median(time1) #time at which it varies
beta_t <- 3
inner_seed <- 13
set.seed(47)
n = 500 #sample size
valid_df <- 0
i <- 1
cens.p <- numeric()
data.samples <- list()
metric_results <- data.frame()

while (valid_df < 1000) { #number of iterations needed
  cov <- gen_covariates(n = n)
  sample.data <- simulate_data_nonPH(n = n, P = cov, betas = coef_X, k = k, beta_t = beta_t, t0 = tvary, lambda = 0.125)
  
  old_seed  <- .Random.seed
  data_full <- splitData(sim_data = sample.data, split_size = 0.7, seed = inner_seed)
  .Random.seed <- old_seed
  
  c.indexes <- tryCatch({
    cox_split(train_df = data_full[[1]], test_df = data_full[[2]])
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
  
  cens.p[valid_df + 1] <- cens_p # only for valid sets
  
  data.samples[[valid_df + 1]] <- sample.data
  metric_results <- rbind(metric_results, data.frame(Iteration = i, Sample = n, c.indexes))
  valid_df <- valid_df + 1
  i <- i + 1
  print(i)
  
}

# summary
avg_metric <- metric_results %>%
  group_by(Model) %>%
  summarise(Hctt_index = mean(C_index_Hctt),
            Hctt_lower_CI = mean(Lower_CI_Hctt),
            Hctt_upper_CI = mean(Upper_CI_Hctt), 
            Hcstr_index = mean(C_index_Hcstr),
            Hcstr_lower_CI = mean(Lower_CI_Hcstr),
            Hcstr_upper_CI = mean(Upper_CI_Hcstr), 
            H_index = mean(C_index_H),
            H_lower_CI = mean(Lower_CI_H),
            H_upper_CI = mean(Upper_CI_H),
            U_index = mean(C_index_U),
            U_lower_CI = mean(Lower_CI_U),
            U_upper_CI = mean(Upper_CI_U), IBScore = mean(IBScore)
  ) %>%
  as.data.frame()

print(avg_metric)

mean(cens.p)
