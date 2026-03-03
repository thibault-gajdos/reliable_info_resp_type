rm(list=ls(all=TRUE)) 

library(tidyverse)
library(cmdstanr)
library(posterior)
library(bayesplot)
library(brms)
library(kableExtra)


exp <- 2; # specify which experiment to run
targetvar <- 1; 
whichmodel <- 'log_trunc_simplified';
whichmodel <- 'log_basic';

#####################################################
##  FIT THE MODEL
####################################################
if (targetvar == 1){
  target <- 'choice';
  whichmodel <- paste0('choice_', whichmodel)
  load(paste0('../data/data_list_', target, '_exp', as.character(exp), '.rdata'))
  data_list = data_list_choice;
} 

data_list$grainsize = 5 ## specify grainsize for within chain parallelization

## Compile the model
model <- cmdstan_model(
    stan_file = paste0('../stan/', whichmodel, '.stan'),
    force_recompile = TRUE, ## necessary if you change the mode
    cpp_options = list(stan_opencl = FALSE, stan_threads = TRUE), ## within chain parallel
    stanc_options = list("O1"), ## fastest sampling
    compile_model_methods = TRUE ## necessary for loo moment matching
)


## Sampling
fit <- model$sample(
  data = data_list,
  seed = 4321,
  chains = 4,
  parallel_chains = 4,
  threads_per_chain = 4,
  iter_warmup = 2000,
  iter_sampling = 1000,
  max_treedepth = 12,
  adapt_delta = .9,
  save_warmup = FALSE
)


## Compute LOO
loo <- fit$loo(cores = 10, moment_match = TRUE)

## Save results
save(fit, file = paste0('../results/fits/exp', as.character(exp), '/fit_',whichmodel, '_exp', as.character(exp), '.rdata'))
save(loo, file = paste0('../results/loo/exp', as.character(exp), '/loo_',whichmodel, '_exp', as.character(exp), '.rdata'))



load(file = paste0('../results/fits/exp', as.character(exp), '/fit_',whichmodel, '_exp', as.character(exp), '.rdata'))

posterior_samples <- fit$draws()
posterior_df <- as_draws_df(posterior_samples)
selected_params <- posterior_df[, grepl("^mu", colnames(posterior_df)) & !grepl("^mu_pr", colnames(posterior_df))]
plot <- mcmc_pairs(selected_params,
                   off_diag_args = list(size = 0.1, alpha=0.1))
plot

cor(selected_params$mu_theta, selected_params$mu_psi)

selected_params <- posterior_df[, grepl("^mu", colnames(posterior_df)) & !grepl("^mu_pr", colnames(posterior_df))]
mean <- sapply(selected_params, function(var) mean(var))
#   mu_alpha    mu_beta  mu_lambda   mu_theta 
#   1.36944760 0.22714114 0.04713966 2.67189763 

selected_params <- posterior_df[, grepl("^sigma", colnames(posterior_df))]
plot <- mcmc_pairs(selected_params,
                   off_diag_args = list(size = 0.1, alpha=0.1))
sigma <- sapply(selected_params, function(var) mean(var))
# sigma_pr[1] sigma_pr[2] sigma_pr[3] sigma_pr[4] 
#  1.0394026   0.2222452   0.7783150   0.4049217 

selected_params <- posterior_df[, grepl("^params", colnames(posterior_df)) & grepl("1]", colnames(posterior_df))]
params1 <- sapply(selected_params, function(var) mean(var))
selected_params <- posterior_df[, grepl("^param_raw", colnames(posterior_df)) & grepl("1]", colnames(posterior_df))]
params_raw1 <- sapply(selected_params, function(var) mean(var))
c(mean[1], sd(sigma[1]*params_raw1))
c(mean(params1), sd(params1))
# 1.4507284 0.8636498

selected_params <- posterior_df[, grepl("^params", colnames(posterior_df)) & grepl("2]", colnames(posterior_df))]
params2 <- sapply(selected_params, function(var) mean(var))
selected_params <- posterior_df[, grepl("^param_raw", colnames(posterior_df)) & grepl("2]", colnames(posterior_df))]
params_raw2 <- sapply(selected_params, function(var) mean(var))
c(mean[2], sd(sigma[2]*params_raw2))
c(mean(params2), sd(params2))
# 0.2272730 0.1829485

selected_params <- posterior_df[, grepl("^params", colnames(posterior_df)) & grepl("3]", colnames(posterior_df))]
params3 <- sapply(selected_params, function(var) mean(var))
selected_params <- posterior_df[, grepl("^param_raw", colnames(posterior_df)) & grepl("3]", colnames(posterior_df))]
params_raw3 <- sapply(selected_params, function(var) mean(var))
c(mean[3], sd(sigma[3]*params_raw3))
c(mean(params3), sd(params3))
# [1] 0.07574997 0.07230431

selected_params <- posterior_df[, grepl("^params", colnames(posterior_df)) & grepl("4]", colnames(posterior_df))]
params4 <- sapply(selected_params, function(var) mean(var))
selected_params <- posterior_df[, grepl("^param_raw", colnames(posterior_df)) & grepl("4]", colnames(posterior_df))]
params_raw4 <- sapply(selected_params, function(var) mean(var))
c(mean[4], sd(sigma[4]*params_raw4))
c(mean(params4), sd(params4))
# 2.7565026 0.8262306

