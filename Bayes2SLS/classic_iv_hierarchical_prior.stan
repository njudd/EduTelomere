data {
  int N;
  int PX; // dimension of exogenous covariates
  int PZ; // dimension of instruments
  int J; // number of previous studies
  vector[J] b_j; // previous estimates
  vector<lower = 0>[J] se_j; // standard error of previous estimates
  matrix[N, PX] X_exog; // exogenous covariates
  vector[N] X_endog; // engogenous covariates
  matrix[N, PZ] Z; // instruments
  vector[N] Y_outcome; // outcome variable
  int<lower=0,upper=1> run_estimation; // simulate (0) or estimate (1)
}
transformed data {
  matrix[N, 2] Y;
  Y[,1] = X_endog;
  Y[,2] = Y_outcome;
}
parameters {
  // intercepts
  vector[2] alpha;
  
  // hierarchical prior parameters
  vector[J] z;
  real beta_hat;
  real<lower = 0> sigma_beta;
  
  // first stage coefficients
  vector[PZ] Delta;
  vector[PX] Gamma;
  
  // second stage coefficients
  vector[PX] delta;
  real beta_ours;
  
  // covariance matrix parameters
  cholesky_factor_corr[2] L_Omega;
  vector<lower = 0>[2] tau;
}
transformed parameters {
  matrix[N, 2] mu; // the conditional means of the process
  vector[J] beta_j;
  
  // first stage
  mu[:,1] = rep_vector(alpha[1], N) + append_col(X_exog, Z)*append_row(Delta, Gamma);
  // second stage
  mu[:,2] = rep_vector(alpha[2], N) + append_col(X_endog,X_exog) * append_row(beta_ours, delta);
  
  // this the non-centered parameterization (if beta_j ~ normal(beta_hat, sigma_beta) then
  // beta_j = beta_hat + sigma_beta *z with z ~ normal(0, 1))
  
  for(j in 1:J) {
    beta_j[j] = beta_hat + sigma_beta * z[j];
  }
}
model {
  // intercept prior
  alpha ~ normal(0, 1);
  
  // first stage priors
  Delta ~ normal(0, 1);
  Gamma ~ normal(0, 1);
  
  // second stage priors
  delta ~ normal(0, 1);
  
  // hierarchical prior 
  
  z ~ normal(0, 1);
  
  for(j in 1:J) {
    b_j[j] ~ normal(beta_j[j], se_j[j]);
  }
  
  beta_ours ~ normal(beta_hat, sigma_beta);
  
  // Covarince matrix prior
  tau ~ cauchy(0, 2);
  L_Omega ~ lkj_corr_cholesky(2);

  // likelihood
  if(run_estimation ==1){
    for(n in 1:N) {
      Y[n] ~ multi_normal_cholesky(mu[n], diag_pre_multiply(tau, L_Omega));
    }
  }
}
generated quantities {
  vector[N] y_sim;
  vector[N] x_endog;
  
  for(n in 1:N) {
    vector[2] error;
    error = multi_normal_cholesky_rng(rep_vector(0.0, 2), diag_pre_multiply(tau, L_Omega));
    
    x_endog[n] = mu[n, 1] + error[1];
    y_sim[n] = alpha[2] + x_endog[n] * beta_ours + X_exog[n]* delta + error[2];
  }
}
