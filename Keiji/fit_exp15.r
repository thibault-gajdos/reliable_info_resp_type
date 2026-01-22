rm(list=ls(all=TRUE)) 

library(tidyverse)
library(cmdstanr)
library(posterior)
library(bayesplot)
library(brms)
library(kableExtra)


exp <- 15; # specify which experiment to run
targetvar <- 1; 
whichmodel <- 'log_trunc_simplified';

#####################################################
##  FIT THE MODEL
####################################################
if (targetvar == 1){
  target <- 'normproba';
  whichmodel <- paste0('proba_', whichmodel)
  load(paste0('../data/data_list_', target, '_exp', as.character(exp), '.rdata'))
  data_list = data_list_normproba;
} else if (targetvar == 2){
  target <- 'proba';
  whichmodel <- paste0('proba_', whichmodel)
} else if (targetvar == 3){
  target <- 'simchoice';
  whichmodel <- paste0('choice_', whichmodel)
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





