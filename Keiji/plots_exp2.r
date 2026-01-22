rm(list=ls(all=TRUE))  


library(tidyverse)
library(cmdstanr)
library(posterior)
library(bayesplot)
library(brms)
library(kableExtra)

# ---------------------------------------------------
# DEFINE AND LOAD EXPERIMENT
# ---------------------------------------------------

exp <- '2'; # specify which experiment to run
targetvar <- 1; 
whichmodel <- 'log_trunc_simplified';

if (targetvar == 1){
  target <- 'choice'; xlabel = 'choice'; ylabel = "Prediction (Choice Proba)";
  whichmodel <- paste0('choice_', whichmodel)
  load(paste0('../data/data_list_', target, '_exp', as.character(exp), '.rdata'))
  data_list = data_list_choice;
  targetvar = data_list$choice;
}
N <- data_list$N

# ---------------------------------------------------
# PAIRS AND TRACES PLOTS
# ---------------------------------------------------

## LOAD MODELS FIT
load(paste0('../results/fits/exp',exp,'/fit_', whichmodel,'_exp',exp,'.rdata'))

## Pairs plots
# Extract the current fit
posterior_samples <- fit$draws()
posterior_df <- as_draws_df(posterior_samples)
selected_params <- posterior_df[, grepl("^mu", colnames(posterior_df)) & !grepl("^mu_pr", colnames(posterior_df))]
plot <- mcmc_pairs(selected_params,
                   off_diag_args = list(size = 0.1, alpha=0.1))
plot_file <- paste0('../results/plots/exp',exp, '/',target,'/pairs_plot_', whichmodel,"_exp", exp)
ggsave(plot, file = paste0(plot_file, ".jpeg"), width = 8, height = 8)



## Traces
selected_params <- posterior_df[, grepl("^mu", colnames(posterior_df))]
plot <- mcmc_trace(selected_params, size=0.1)
plot_file <- paste0('../results/plots/exp',exp, '/',target,'/trace_plot_', whichmodel,"_exp", exp)
ggsave(plot, file = paste0(plot_file, ".jpeg"), width = 10, height = 8)


# ---------------------------------------------------
# CHECK MODEL PREDICTION 
# ---------------------------------------------------
# indv_theta = selected_params <- posterior_df[, grepl("^params", colnames(posterior_df)) & grepl("4]", colnames(posterior_df))]
# n = 1; t = 1;
# theta = as.numeric(indv_theta[t,n])
# diff_evidence <- posterior_df[, grepl("^diff_evidence", colnames(posterior_df)) ]
# diff_evidence <- as.numeric(diff_evidence[1,1])
# p <- 1 / (1+exp(-theta*diff_evidence))
# proba <- posterior_df[, grepl("^pred_proba", colnames(posterior_df)) ]
# proba[1,1]

proba <- posterior_df[, grepl("^pred_proba", colnames(posterior_df)) ]
median <- sapply(proba, function(var) median(var))
dat = c(); sub = c(); tr = c(); j = 0; 
for (t in 1:data_list$T_max){
  for (s in 1:data_list$N){
    j = j + 1; 
    sub[j] = s;
    tr[j] = t;
    dat[j] = 2-targetvar[s,t]; ## 0 = red, 1 = blue
  }
}
d <- c(); 
d <- data.frame(median = median, data = dat, subj=sub, trial=tr) 
cor(d$median, d$data)

# Create a single plot with all subjects
plot <- ggplot(d, aes(x = factor(data), y = median, color = factor(subj))) +
  geom_violin() +       # Function curve
  facet_wrap(~ subj, ncol = 8) + # Create a plot for each subject
  theme_minimal() +
  labs(title = "Predictive check by Subject", x = xlabel, y = ylabel) +
  theme(legend.position = "none")

plot_file <- paste0('../results/plots/exp',exp, '/',target,'/pred_check_', whichmodel,"_exp", exp)
ggsave(plot, file = paste0(plot_file, ".jpeg"), width = 7, height = 8)



