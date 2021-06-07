data {
  int<lower=0> N;
  vector[N] y;
  int K;
  int start;
  int Age_idx[N];
  vector[K] mu;
}

transformed data {
  int count_up[abs(start-K)];
  int count_down[start-1];
  for (k in 1:num_elements(count_down)) count_down[k] = start - k;
  for (k in 1:num_elements(count_up)) count_up[k] = start + k;
}

parameters {
  real<lower=0> delta_intercept;
  vector[K] delta_raw;
  real<lower=0> delta_sd;
  real<lower=0,upper=1> theta_intercept;
  vector[K] theta_raw;
  real<lower=0> theta_sd;
  real<lower=.25> sigma1;
  real<lower=.25> sigma2;
}

transformed parameters {
  vector[K] theta;
  vector[K] delta;
  vector[K] mu1;
  vector[K] mu2;
  
  theta[start] = logit(theta_intercept);
  delta[start] = delta_intercept;
  for (k in count_up) {
    theta[k] = theta[k-1] + theta_raw[k]*theta_sd;
    delta[k] = delta[k-1] + delta_raw[k]*delta_sd;
    
  }
  for (k in count_down) {
    theta[k] = theta[k+1] + theta_raw[k]*theta_sd;
    delta[k] = delta[k+1] + delta_raw[k]*delta_sd;
  }
  theta = inv_logit(theta);
  mu1 = mu - delta;
  mu2 = (mu - theta .* mu1) ./ (1-theta);
}

model {
  for (n in 1:N) {
    target += log_mix(theta[Age_idx[n]],
                        normal_lpdf(y[n] | mu1[Age_idx[n]], sigma1),
                        normal_lpdf(y[n] | mu2[Age_idx[n]], sigma2));
  }
  delta_intercept ~ normal(3.5,1); 
  theta_intercept ~ beta(10,20); 
  theta_raw ~ std_normal();
  delta_raw ~ std_normal();
  theta_sd ~ normal(0,.05);
  delta_sd ~ normal(0,.1);
  sigma1 ~ std_normal();
  sigma2 ~ std_normal();
}

generated quantities {
  vector[K] yhat;
  for (k in 1:K) yhat[k] = theta[k] > uniform_rng(0,1) ? normal_rng(mu1[k],sigma1) : normal_rng(mu2[k],sigma2);
}
