functions {
  vector clamp_vector(vector x, real lo, real hi) {
    vector[num_elements(x)] out;
    for (i in 1:num_elements(x)) {
      out[i] = fmin(fmax(x[i], lo), hi);
    }
    return out;
  }
  
  // Worker function for reduce_sum
  real partial_sum(array[] int slice_indices,
  int start, int end,
  vector mu_pr, vector sigma_pr,
  array[] int Tsubj,
  array[,] int sample,
  array[,,] int color,
  array[,,] real proba,
  array[,] int choice,
  matrix param_raw){
    real lp = 0;
    for (i in 1:size(slice_indices)) {
      int n = slice_indices[i];
      
      // Transform parameters
      vector[4] params;
      params[1] = mu_pr[1] + sigma_pr[1] * param_raw[n, 1]; // alpha
      params[2] = mu_pr[2] + sigma_pr[2] * param_raw[n, 2]; // beta       
      params[3] = Phi_approx(mu_pr[3] + sigma_pr[3] * param_raw[n, 3]); // lambda
      params[4] = Phi_approx(mu_pr[4] + sigma_pr[4] * param_raw[n, 4])*10; // theta

      for (t in 1:Tsubj[n]) {
        vector[2] evidence = rep_vector(0.0, 2);
        int sample_size = sample[n, t];
        
        for (s in 1:sample_size) {
          real p = proba[n, t, s];
          real l = logit(p);
          int color_val = color[n, t, s];
          real log_odds = params[1] * l +  params[2];
          evidence[color_val] += exp(params[3] * (s - sample_size)) * log_odds;
        }
        vector[2] evidence_safe = clamp_vector(evidence, -100, 100);
        lp += categorical_lpmf(choice[n, t] | softmax(params[4]*evidence_safe));

      }
    }
    
    return lp;
  }
  
  vector compute_evidence(int sample_size, array[] int color_data, array[] real proba_data, 
  real alpha, real beta, real lambda, real theta) {
    vector[2] evidence = rep_vector(0.0, 2);
    
    for (s in 1:sample_size) {
      real p = proba_data[s];
      real l = logit(p);
      real log_odds = alpha * l +  beta; 
      evidence[color_data[s]] += exp(lambda * (s - sample_size)) * log_odds;
    }
    return evidence;
  }
  
  real compute_log_lik(int sample_size, array[] int color_data, array[] real proba_data, 
  int choice, real alpha, real beta, real lambda, real theta) {
    vector[2] evidence = compute_evidence(sample_size, color_data, proba_data, alpha, beta, lambda, theta);
    vector[2] evidence_safe = clamp_vector(evidence, -100, 100);
    return categorical_lpmf(choice | softmax(theta*evidence_safe));
  }
}

data {
  int<lower=1> N;
  int<lower=1> T_max;
  int<lower=1> I_max;
  array[N] int<lower=1> Tsubj;
  array[N, T_max] int<lower=-1> sample;
  array[N, T_max, I_max] int<lower=-1, upper=2> color;
  array[N, T_max, I_max] real<lower=-1, upper=1> proba;
  array[N, T_max] int<lower=-1, upper=2> choice;
  
  //  grainsize for reduce_sum
  int<lower=5> grainsize;
  
}

parameters {
  vector[4] mu_pr;
  vector<lower=0>[4] sigma_pr;
  matrix[N, 4] param_raw;
}

model {
  mu_pr ~ std_normal();
  sigma_pr ~ std_normal();
  to_vector(param_raw) ~ std_normal();
  
  // Create array of indices for parallelization
  array[N] int indices;
  for (n in 1:N) indices[n] = n;
  
  // Use reduce_sum for parallelization
  target += reduce_sum(partial_sum, indices, grainsize,
  mu_pr, sigma_pr,
  Tsubj, sample, color, proba, choice, param_raw);
}

generated quantities {
  matrix[N, 4] params;
  array[N, T_max] real y_pred = rep_array(-1.0, N, T_max); 
  array[N, T_max] real pred_proba = rep_array(-1.0, N, T_max); 
  array[N, T_max] real diff_evidence = rep_array(-1.0, N, T_max); 
  vector[sum(Tsubj)] log_lik;
  real mu_alpha = mu_pr[1];
  real mu_beta = mu_pr[2];
  real mu_lambda = Phi_approx(mu_pr[3]);
  real mu_theta = Phi_approx(mu_pr[4])*10;

  int k = 0;
  for (n in 1:N) {
    params[n, 1] = mu_pr[1] + sigma_pr[1] * param_raw[n, 1];
    params[n, 2] = mu_pr[2] + sigma_pr[2] * param_raw[n, 2];
    params[n, 3] = Phi_approx(mu_pr[3] + sigma_pr[3] * param_raw[n, 3]);
    params[n, 4] = Phi_approx(mu_pr[4] + sigma_pr[4] * param_raw[n, 4])*10;

    for (t in 1:Tsubj[n]) {
      k += 1;
      
      array[I_max] int color_trial;
      array[I_max] real proba_trial;
      
      for (i in 1:I_max) {
        color_trial[i] = color[n, t, i];
        proba_trial[i] = proba[n, t, i];
      }
      
      vector[2] evidence = compute_evidence(sample[n, t], color_trial, proba_trial,
      params[n, 1], params[n, 2], params[n, 3], params[n, 4]);
      vector[2] evidence_safe = clamp_vector(evidence, -100, 100);
      vector[2] prob = softmax(params[n,4]*evidence_safe);
      y_pred[n, t] = categorical_rng(prob);
      log_lik[k] = compute_log_lik(sample[n, t], color_trial, proba_trial,
      choice[n, t],
      params[n, 1], params[n, 2], params[n, 3], params[n, 4]);
      
      diff_evidence[n, t] = evidence_safe[1]-evidence_safe[2];
      pred_proba[n, t] = prob[1];
    }
  }
}
