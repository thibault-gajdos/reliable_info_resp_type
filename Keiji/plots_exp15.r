rm(list=ls(all=TRUE))  ## efface les données


library(tidyverse)
library(cmdstanr)
library(posterior)
library(bayesplot)
library(brms)
library(kableExtra)

# ---------------------------------------------------
# DEFINE AND LOAD EXPERIMENT
# ---------------------------------------------------

exp <- '15'; # specify which experiment to run
targetvar <- 1; 
whichModel <- 'log_trunc_simplified';

if (targetvar == 1){
  target <- 'normproba';  xlabel = 'Lo(normalised proba estimates)';
  whichmodel <- paste0('proba_', whichModel)
  load(paste0('../data/data_list_', target, '_exp', as.character(exp), '.rdata'))
  data_list = data_list_normproba; 
  targetvar = data_list$lo_proba;
} else if (targetvar == 2){
  target <- 'proba';
  whichmodel <- paste0('proba_', whichmodel)
} else if (targetvar == 3){
  target <- 'simchoice';
  whichmodel <- paste0('choice_', whichmodel)
}
N <- data_list$N

levels = c(0.5,0.55,0.60,0.65,0.7,0.75);
colorcode = list();
colorcode[[1]] = c(41, 201, 255);
colorcode[[2]] = c(198, 222, 155);
colorcode[[3]] = c(90, 200, 82);
colorcode[[4]] = c(8, 143, 68);
colorcode[[5]] = c(4, 98, 5);
colorcode[[6]] = c(0, 50, 0);



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

selected_params <- posterior_df[, grepl("^y_pred", colnames(posterior_df)) ]
median <- sapply(selected_params, function(var) median(var))
dat = c(); sub = c(); tr = c(); j = 0; 
for (t in 1:data_list$T_max){
  for (s in 1:data_list$N){
    j = j + 1; 
    sub[j] = s;
    tr[j] = t;
    dat[j] = targetvar[s,t]; 
  }
}
d <- c(); 
d <- data.frame(median = median, data = dat, subj=sub, trial=tr) 

# Create a single plot with all subjects
plot <- ggplot(d, aes(x = data, y = median, color = factor(subj))) +
  geom_point(size = 1) +       # Function curve
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +  # Diagonal y = x
  facet_wrap(~ subj, ncol = 5) + # Create a plot for each subject
  #scale_x_continuous(limits = c(0, 1)) +
  #scale_y_continuous(limits = c(0, 1)) +
  theme_minimal() +
  labs(title = "Predictive check by Subject", x = xlabel, y = "Prediction") +
  theme(legend.position = "none")

plot_file <- paste0('../results/plots/exp',exp, '/',target,'/pred_check_', whichmodel,"_exp", exp)
ggsave(plot, file = paste0(plot_file, ".jpeg"), width = 8, height = 8)



# ---------------------------------------------------
#  MODELS SUMMARIES 
# ---------------------------------------------------

## Parameters values
selected_params <- posterior_df[, grepl("^mu", colnames(posterior_df)) & !grepl("^mu_pr", colnames(posterior_df))]
summary_stats <- posterior_summary(selected_params)
summary_df <- as.data.frame(summary_stats)

## GROUP LEVEL
parameters <- summary_df
p = c(1:99)/100

if (whichModel == 'log_trunc_simplified'){
  f_trunc_simplified <- function(p, psi, l_obs, alpha, beta) {
    l = log(p/(1-p)) 
    ll = alpha*l*(psi/l_obs) + (1-alpha)*beta
    fp = exp(ll)/(1+exp(ll))
    return(fp)
  }
  alpha = parameters["mu_alpha", "Estimate"]
  beta = parameters["mu_beta", "Estimate"]
  psi = parameters["mu_psi", "Estimate"]
  l_obs = data_list$l_obs
  fp <- sapply(p, f_trunc_simplified, psi = psi, l_ob=l_obs, alpha =alpha, beta=beta)
}

df <- data.frame(p = p, fp = fp); 
df2 <- df %>% filter(p %in% levels);
df2$col <- vapply(colorcode, 
                  function(x) rgb(x[1], x[2], x[3], maxColorValue = 255),
                  character(1))

plot <- ggplot(df, aes(x = p, y = fp)) +
  geom_line(color = "blue", size = 1) +        # function curve
  geom_abline(intercept = 0, slope = 1,        # diagonal y = x
              color = "red", linetype = "dashed") +
  geom_point(data=df2, aes(x = p, y = fp, color = col),size=3) +
  scale_color_identity()+
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, 1)) +
  theme_minimal()
plot_file <- paste0('../results/plots/exp',exp, '/',target,'/func_', whichmodel,"_exp", exp)
ggsave(plot, file = paste0(plot_file, ".jpeg"), width = 5, height = 5)



