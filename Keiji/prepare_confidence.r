rm(list=ls(all=TRUE)) 

library(tidyverse)
library(cmdstanr)
library(posterior)
library(bayesplot)
library(brms)
library(kableExtra)


######################################
## prepare the data for fits etc.
#####################################

exp <- 12; # specify which experiment to run
data <- read.csv(paste0("../data/DATA_EXPChoiceReliab_Exp", as.character(exp),".csv"))

## Generate the variables proba_1,.., proba_6 and color_1,...,color_6
data <- data %>%
  separate(Sample_Reliability, into = paste0("proba_", 1:6), sep = "\\s+") %>%
  mutate(across(starts_with("proba_"), as.numeric)) %>%     
  mutate(color = str_extract_all(Sample_Color, "blue|red"))  %>%
  unnest_wider(color, names_sep = "_") 

write.csv(data, paste0("../data/data_reliability_exp", as.character(exp),".csv"))
save(data, file = paste0("../data/data_reliability_exp", as.character(exp),".rdata"))



##################################################
##        PREPARE THE DATA       
##################################################

load(file = paste0("../data/data_reliability_exp", as.character(exp),".rdata"))

##  if Response  ResponseButtonOrder= 1:  blue->1, red->0
##  if Response  ResponseButtonOrder= 0:  blue->0, red->1
##  we recode: blue ->1, red ->2
data <- data %>%
  mutate(choice = case_when(
    (ResponseButtonOrder == 1 & Response == 0) ~ 2, # red   
    (ResponseButtonOrder == 1 & Response == 1) ~ 1, # blue 
    (ResponseButtonOrder == 0 & Response == 0) ~ 1, # blue
    (ResponseButtonOrder == 0 & Response == 1) ~ 2  # red
  )) %>%
  mutate_at(vars(starts_with("color")), ~ ifelse(. == "blue", 1, 2)) %>%
  rowwise() %>%
  mutate(sample_number = sum(!is.na(c_across(starts_with("proba_"))))) %>%  
  ungroup() %>%
  mutate(conf = SliderResp_50to100/100) %>%
  group_by(ParticipantPrivateID) %>%
  mutate(conf_blue = case_when(
    (choice == 1) ~ conf,  
    (choice == 2) ~ 1-conf 
  )) %>%
  mutate(conf_blue_safe = case_when(
    (conf_blue == 0) ~ 0.001,
    (conf_blue == 1) ~ 0.999,
    TRUE ~ conf_blue
  )) %>%
  mutate(lo_conf = log(conf_blue_safe / (1 - conf_blue_safe))
         # log odds of confidence of blue being correct will be fitted 
  )%>%  
  mutate( norm_conf = ( (conf - min(conf)) / (max(conf) - min(conf)) * (0.999 - 0.5) + 0.5 )
          # normalised so the min conf value per ppt becomes 0.5 (approx. equal likelihood in the model probability scale) and the max value becomes 1 (approx. maximum likelihood) 
  ) %>%
  mutate( norm_conf_blue = case_when(
    # if the choice is blue (red), the subj report represents their belief about the proba of BLUE (RED) being correct
    (choice == 1) ~ norm_conf, # belief that BLUE is correct (50%-100%)
    (choice == 2) ~ 1-norm_conf # if choice is red, we take 1-conf so this represents their belief that BLUE is correct (0%-50%)
  )) %>%
  mutate(lo_normconf = log(norm_conf_blue / (1 - norm_conf_blue))
         # log odds of normalised confidence of blue being correct will be fitted 
  ) %>% 
  ungroup() 


N = length(unique(data$ParticipantPrivateID))
T_max = max(data$TrialNumber)
I_max <- max(data$sample_number) ## max number of samples/trial
## compute trials by subject
d <- data %>%
  group_by(ParticipantPrivateID) %>%
  summarise(t_subjs = n())
t_subjs <- d$t_subjs
subjs <- unique(data$ParticipantPrivateID)
l_obs = log( max(data$proba_1/100) / (1 - max(data$proba_1/100))); # the maxima of the observed reliability


## Initialize data arrays
choice  <- array(-1, c(N, T_max))
lo_conf  <- array(-1, c(N, T_max))
lo_normconf  <- array(-1, c(N, T_max))
posteriorBlue  <- array(-1, c(N, T_max))
color <- array( -1, c(N, T_max, I_max))
proba <- array(-1, c(N, T_max, I_max))
sample <- array(-1, c(N, T_max))
## fill the  arrays
for (n in 1:N) { ## loop through subjects
  t <- t_subjs[n] ## number of trials for subj i
  data_subj <- data %>% filter(ParticipantPrivateID == subjs[n])
  choice[n, 1:t] <- data_subj$choice 
  lo_conf[n, 1:t] <- data_subj$lo_conf 
  lo_normconf[n, 1:t] <- data_subj$lo_normconf 
  posteriorBlue[n, 1:t] <- data_subj$PosteriorBlue 
  for (k in 1:t) { ## loop through trials
    data_subj_t <- data_subj[k,]
    sample[n,k] <- data_subj_t$sample_number
    for (i in 1:data_subj_t$sample_number) {
      color_var <- paste0("color_", i)
      proba_var <- paste0("proba_", i)
      color[n, k, i] <- data_subj[[color_var]][k]
      proba[n, k, i] <- data_subj[[proba_var]][k]/100
    }
  }
}

data_list_conf <- list(
    N = N,
    T_max = T_max,
    I_max = I_max,
    Tsubj = t_subjs,
    color = color,
    proba = proba,
    lo_conf = lo_conf, ## 
    posteriorBlue = posteriorBlue,
    sample = sample,
    l_obs = l_obs
)

data_list_normconf <- list(
  N = N,
  T_max = T_max,
  I_max = I_max,
  Tsubj = t_subjs,
  color = color,
  proba = proba,
  lo_conf = lo_normconf, ## 
  posteriorBlue = posteriorBlue,
  sample = sample,
  l_obs = l_obs
)

data_list_choice <- list(
  N = N,
  T_max = T_max,
  I_max = I_max,
  Tsubj = t_subjs,
  color = color,
  proba = proba,
  choice = choice, ## 
  posteriorBlue = posteriorBlue,
  sample = sample,
  l_obs = l_obs
)


save(data_list_conf, file = paste0('../data/data_list_conf_exp', as.character(exp), '.rdata'))
save(data_list_normconf, file = paste0('../data/data_list_normconf_exp', as.character(exp), '.rdata'))
save(data_list_choice, file = paste0('../data/data_list_choice_exp', as.character(exp), '.rdata'))

