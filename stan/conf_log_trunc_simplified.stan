functions {
  vector clamp_vector(vector x, real lo, real hi) {
    vector[num_elements(x)] out;
    for (i in 1:num_elements(x)) {
      out[i] = fmin(fmax(x[i], lo), hi);
    }
    return out;
  }
  
  real clamp(real x, real lo) {
    real out;
      out = fmax(x, lo);
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
                 array[,] real lo_conf, // log odds of confidence scale (<=0: a high confidence in favour of red; >0: a high confidence in favour of blue)
                 matrix param_raw,
                 real l_obs){
    real lp = 0;
    for (i in 1:size(slice_indices)) {
      int n = slice_indices[i];
      
      // Transform parameters
      vector[6] params;
      params[1] = Phi_approx(mu_pr[1] + sigma_pr[1] * param_raw[n, 1]); // alpha
      params[2] = mu_pr[2] + sigma_pr[2] * param_raw[n, 2]; // beta 
      params[3] = Phi_approx(mu_pr[3] + sigma_pr[3] * param_raw[n, 3]); // lambda
      params[4] = Phi_approx(mu_pr[4] + sigma_pr[4] * param_raw[n, 4])*2; // sigma for logodds of estimate
      params[5] = Phi_approx(mu_pr[5] + sigma_pr[5] * param_raw[n, 5])*5; // psi

      for (t in 1:Tsubj[n]) {
        vector[2] evidence = rep_vector(0.0, 2);
        int sample_size = sample[n, t];

        for (s in 1:sample_size) {
          real p = proba[n, t, s];
          real l = logit(p);
          int color_val = color[n, t, s];
          // without actual truncation
          real log_odds = params[1] * l * (params[5] / l_obs) + (1-params[1])* params[2]; // divided by l_obs

          evidence[color_val] += exp(params[3] * (s - sample_size)) * log_odds;
        }
        
        vector[2] evidence_safe = clamp_vector(evidence, -100, 100);
        real lp_raw = normal_lpdf(lo_conf[n,t] | evidence_safe[1]-evidence_safe[2],  clamp(params[4], 0.0001));
        real lp_safe = clamp(lp_raw, -1000);

        lp += lp_safe;
      }
    }

    return lp;
  }

  vector compute_evidence(int sample_size, array[] int color_data, array[] real proba_data, 
                          real alpha, real beta, real lambda, real psi, real l_obs) {
    vector[2] evidence = rep_vector(0.0, 2);

    for (s in 1:sample_size) {
      real p = proba_data[s];
      real l = logit(p);
      real log_odds = alpha * l * (psi / l_obs) + (1-alpha) * beta; // divided by l_obs
      evidence[color_data[s]] += exp(lambda * (s - sample_size)) * log_odds;
    }
    return evidence;
  }

  real compute_log_lik(int sample_size, array[] int color_data, array[] real proba_data, 
                       real lo_conf, real alpha, real beta, real lambda, real sigma_est, real psi, real l_obs) {
    vector[2] evidence = compute_evidence(sample_size, color_data, proba_data, alpha, beta, lambda, psi, l_obs);
    vector[2] evidence_safe = clamp_vector(evidence, -100, 100);
    real lp_raw = normal_lpdf(lo_conf | evidence_safe[1]-evidence_safe[2],  clamp(sigma_est, 0.0001));
    real lp_safe = clamp(lp_raw, -1000);
    
    return lp_safe;
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
  array[N, T_max] real<lower=-7, upper=7> lo_conf;

  //  grainsize for reduce_sum
  int<lower=5> grainsize;
  
  // precomputed largest proba logit
  real<lower=0> l_obs;
}

parameters {
  vector[5] mu_pr;
  vector[5] sigma_pr;
  matrix[N, 5] param_raw;
}

model {
  mu_pr ~ std_normal();
  sigma_pr ~ std_normal();
  to_vector(param_raw) ~ std_normal();
  
  // Create array of indices for parallelization
  array[N] int indices;
  for (n in 1:N) indices[n] = n;
  
  // Use map_rect for parallelization
 target += reduce_sum(partial_sum, indices, grainsize,
                     mu_pr, sigma_pr, 
                     Tsubj, sample, color, proba, lo_conf, param_raw, l_obs);
}

generated quantities {
  matrix[N, 5] params;
  array[N, T_max] real y_pred = rep_array(-1.0, N, T_max);
  vector[sum(Tsubj)] log_lik;
  real mu_alpha = Phi_approx(mu_pr[1]);
  real mu_beta = mu_pr[2];
  real mu_lambda = Phi_approx(mu_pr[3]);
  real mu_sigma_est = Phi_approx(mu_pr[4])*2;
  real mu_psi = Phi_approx(mu_pr[5])*5;

  int k = 0;
  for (n in 1:N) {
    params[n, 1] = Phi_approx(mu_pr[1] + sigma_pr[1] * param_raw[n, 1]);
    params[n, 2] = mu_pr[2] + sigma_pr[2] * param_raw[n, 2];
    params[n, 3] = Phi_approx(mu_pr[3] + sigma_pr[3] * param_raw[n, 3]);
    params[n, 4] = Phi_approx(mu_pr[4] + sigma_pr[4] * param_raw[n, 4])*2;
    params[n, 5] = Phi_approx(mu_pr[5] + sigma_pr[5] * param_raw[n, 5])*5;

    for (t in 1:Tsubj[n]) {
      k += 1;

      array[I_max] int color_trial;
      array[I_max] real proba_trial;

      for (i in 1:I_max) {
        color_trial[i] = color[n, t, i];
        proba_trial[i] = proba[n, t, i];
      }

      vector[2] evidence = compute_evidence(sample[n, t], color_trial, proba_trial,
                                            params[n, 1], params[n, 2], params[n, 3], params[n, 5], l_obs);
      vector[2] evidence_safe = clamp_vector(evidence, -100, 100);
      y_pred[n, t] = normal_rng(evidence_safe[1]-evidence_safe[2], params[n, 4]);
      log_lik[k] = compute_log_lik(sample[n, t], color_trial, proba_trial,
                                   lo_conf[n, t],
                                   params[n, 1], params[n, 2], params[n, 3], params[n, 4], params[n, 5], l_obs);
                                   
    }
  }
}


