data {
  // Data for time course analysis
  int N_DAY;                       // number of data points
  vector[N_DAY] Y_DAY;             // outocome log10 viral load
  vector[N_DAY] X_DAY;             // day of measurement in time series
  int G;                           // number of subjects
  int gstart_DAY[G];               // first day in time series
  int gend_DAY[G];                 // last day in time series
  int N_onset;                     // number of subjects with symptom-onset data
  int idx_onset[N_onset];          // index of subjects with symptom-onset data
  vector[N_onset] onset;           // onset in days from first test
  int N_NegTests;                  // Total number negative test results
  int idx_NegTests[N_NegTests];    // position of negative test results in Y_DAY
  vector[N_DAY] PCR;               // Type of PCR system used for test
  int N_centres;                   // Number of test centre types
  matrix[N_DAY,N_centres] centre;  // Matrix with centre types for each test
  int N_centre1;                   // Number of test centre types for 1st test
  int centre1[G];                  // centre of each participants 1st test
  int N_ld_centre;                 // Number of test centre types with longest "stay"
  int ld_centre[G];                // each subjects centre types with longest "stay"
  real imputation_limit;           // upper limit for imputed viral loads for neg tests
  int num_basis_Age;               // numberf of bases for Age spline model
  matrix[num_basis_Age, G] B_Age;  // regression matrix for Age spline model
  vector[G] Age;                   // Cases ages
  int K_PGH;                       // Number of covariates
  matrix[G,K_PGH] X_PGH;           // covariates 
  vector[G] max_load;              // Maximum measured load per case
  int condition_on_data;           // 1 for parameter estimation, 0 for prior predictive plots
  // Data for culture positivity analysis
  int N_CP;                        // number of data points
  int Y_CP[N_CP];                  // Outcome (0,1)
  matrix[N_CP,1] X_CP;             // Predictor (log10 viral load)
}

parameters {
  // grand means for key model parameters
  real log_slope_up_mu;
  real log_slope_down_mu;
  real log_intercept_mu;
  
  // covariates for key model parameters
  row_vector[num_basis_Age] a_raw_intercept_Age;
  real<lower=0> tau_intercept_Age;
  real<lower=0> a0_intercept_Age; 
  vector[K_PGH] betaPGH_intercept;
  row_vector[num_basis_Age] a_raw_slope_up_Age;
  real<lower=0> tau_slope_up_Age;
  real a0_slope_up_Age; 
  vector[K_PGH] betaPGH_slope_up;
  row_vector[num_basis_Age] a_raw_slope_down_Age;
  real<lower=0> tau_slope_down_Age;
  real a0_slope_down_Age; 
  vector[K_PGH] betaPGH_slope_down;

  // individual random effects for key model parameters
  real<lower=0> intercept_sigma;
  vector[G] intercept_raw; 
  real<lower=0> slope_up_sigma;
  vector[G] slope_up_raw;
  real<lower=0> slope_down_sigma;
  vector[G] slope_down_raw;
  real<lower=0> sigma_sigma;
  vector[G] sigma_raw;
  
  // first centre random effects for key model parameters
  real<lower=0> int_centre1_sigma;
  vector[N_centre1] int_centre1_raw;
  real<lower=0> slope_down_ld_centre_sigma;
  vector[N_ld_centre] slope_down_ld_centre_raw;
  real<lower=0> slope_up_ld_centre_sigma;
  vector[N_ld_centre] slope_up_ld_centre_raw;
  real<lower=0> int_ld_centre_sigma;
  vector[N_ld_centre] int_ld_centre_raw;
  
  // smoothness parameter for logistic weight function
  real<lower=0> beta_sweight_mu;
  
  // impute negative tests
  vector<lower=-25,upper=imputation_limit>[N_NegTests] imp_neg;
  
  // shifts for individuals
  vector<lower=-10,upper=20>[G] b_shift;
  // random effects of cirst centre on shifts
  real<lower=0> shift_centre1_sigma;
  vector[N_centre1] shift_centre1_raw;
  
  // error variance
  real sigma_mu; // mean
  
  // centre random effects on viral load
  vector[N_centres] centre_raw;
  real<lower=0> centre_sigma;
  // PCR system effect on viral load
  real intercept_PCR;
  
  // culture positivity parameters
  real alpha_CP;
  vector[1] beta_CP;
  
  real<lower=0,upper=1> theta;
}

transformed parameters {
  // spline coefficients
  row_vector[num_basis_Age] a_slope_up_Age = a_raw_slope_up_Age*tau_slope_up_Age;
  row_vector[num_basis_Age] a_intercept_Age = a_raw_intercept_Age*tau_intercept_Age;
  row_vector[num_basis_Age] a_slope_down_Age = a_raw_slope_down_Age*tau_slope_down_Age;
  
  // centre random effects on measured loads
  vector[N_centres] int_centr = centre_raw * centre_sigma;
  // error variance
  vector[G] sigma = exp(sigma_mu + sigma_sigma*sigma_raw);
  // initialize vectors for individual level parameters
  vector<lower=0>[G] intercept;
  vector<lower=0>[G] slope_up;
  vector<upper=0>[G] slope_down;
  vector<lower=0>[G] time2peak;
  vector[G] shift;
  vector[N_onset] shifted_onset;
  { // calculate model parameters
     vector[N_centre1] shift_centre1 = shift_centre1_sigma * shift_centre1_raw;
     vector[N_centre1] int_centre1 = int_centre1_sigma * int_centre1_raw;
     vector[N_ld_centre] slope_up_ld_centre = slope_up_ld_centre_sigma * slope_up_ld_centre_raw;
     vector[N_ld_centre] int_ld_centre = int_ld_centre_sigma * int_ld_centre_raw;
     vector[N_ld_centre] slope_down_ld_centre = slope_down_ld_centre_sigma * slope_down_ld_centre_raw;
     vector[G] log_intercept = 
                 log_intercept_mu + // intercept (fixed effect [FE], group level mean)
                 intercept_sigma * intercept_raw + // individual level random effects
                 a0_intercept_Age * Age + to_vector(a_intercept_Age * B_Age) + // age with splines (FE)
                 X_PGH * betaPGH_intercept; // covariates (FE)
     vector[G] log_slope_up = 
                 log_slope_up_mu +
                 slope_up_sigma * slope_up_raw + 
                 a0_slope_up_Age * Age + to_vector(a_slope_up_Age * B_Age) + 
                 X_PGH * betaPGH_slope_up; 
     vector[G] log_slope_down = 
                 log_slope_down_mu + 
                 slope_down_sigma * slope_down_raw + 
                 a0_slope_down_Age * Age + to_vector(a_slope_down_Age * B_Age) + 
                 X_PGH * betaPGH_slope_down;
     // test-centre category based random effects
     for (g in 1:G) {
       shift[g] = b_shift[g] + shift_centre1[centre1[g]];
       log_intercept[g] += int_centre1[centre1[g]] + int_ld_centre[ld_centre[g]];
       log_slope_down[g] += slope_down_ld_centre[ld_centre[g]];
       log_slope_up[g] += slope_up_ld_centre[ld_centre[g]];
     }
     // apply link function
     intercept = exp(log_intercept);
     slope_up = exp(log_slope_up);
     slope_down = -exp(log_slope_down);
     shifted_onset = onset + shift[idx_onset]; // calculate temporal location of onset relative to peak load
  }
  time2peak = intercept ./ slope_up;
}

model {
  // DAY
  log_intercept_mu ~ normal(2.1,.1); 
  log_slope_up_mu ~ normal(.6,.25); 
  log_slope_down_mu ~ normal(-1.75,.5); 
  beta_sweight_mu ~ normal(10,1); 
  
  shift_centre1_sigma ~ std_normal();
  shift_centre1_raw ~ std_normal();
  
  intercept_sigma ~ gamma(4,20);
  intercept_raw ~ std_normal();
  
  slope_up_sigma ~ gamma(4,20);
  slope_up_raw ~ std_normal();
  
  slope_down_sigma ~ gamma(4,20);
  slope_down_raw ~ std_normal();
  
  sigma_mu ~ std_normal();
  sigma_sigma ~ std_normal();
  sigma_raw ~ std_normal();
  
  a0_intercept_Age ~ normal(0,.125);
  a_raw_intercept_Age ~ std_normal();
  tau_intercept_Age ~ normal(0, .5); 
  betaPGH_intercept ~ normal(0,.125);
  // one sd change in age changes slope by up to around 20%
  // pnorm(.2,0,.125) = .945, exp(.2) = 1.2
  a0_slope_up_Age ~ normal(0,.125);
  a_raw_slope_up_Age ~ std_normal(); 
  tau_slope_up_Age ~ normal(0, .75);  
  betaPGH_slope_up ~ normal(0,.125);
  a0_slope_down_Age ~ normal(0,.225);
  a_raw_slope_down_Age ~ std_normal(); 
  tau_slope_down_Age ~ normal(0, .75); 
  betaPGH_slope_down ~ normal(0,.225);
  
  // random centre effects
  centre_raw ~ std_normal();
  centre_sigma ~ std_normal();
  
  int_centre1_sigma ~ normal(0,.125);
  int_centre1_raw ~ std_normal();
  slope_down_ld_centre_sigma ~ normal(0,.125);
  slope_down_ld_centre_raw ~ std_normal();
  
  intercept_PCR ~ normal(0,1);
  
  if (condition_on_data == 1) {
    vector[N_DAY] Y_DAY_imputed = Y_DAY;
    for (j in 1:N_NegTests) {
      Y_DAY_imputed[idx_NegTests[j]] = imp_neg[j];
    }
    
    for (g in 1:G) {
      int N_g = gend_DAY[g] - gstart_DAY[g] + 1;
      vector[N_g] DAY_shifted1 = X_DAY[gstart_DAY[g]:gend_DAY[g]] + shift[g];
      vector[N_g] weight_down1 = inv_logit(DAY_shifted1 * beta_sweight_mu);
      vector[N_g] est_up1 = intercept[g] + slope_up[g]*DAY_shifted1;
      vector[N_g] est_down1 = intercept[g] + slope_down[g]*DAY_shifted1;
      vector[N_g] yhat1 = (1-weight_down1) .* est_up1 + weight_down1 .* est_down1;

      yhat1 += PCR[gstart_DAY[g]:gend_DAY[g]] * intercept_PCR;
      yhat1 += centre[gstart_DAY[g]:gend_DAY[g],] * int_centr;
      target += normal_lpdf(Y_DAY_imputed[gstart_DAY[g]:gend_DAY[g]] | yhat1, sigma[g]);
    }
    
    //shifted_onset ~ skew_normal(-3,10,12);
    for (k in 1:N_onset) {
     target += log_mix(theta,
     skew_normal_lpdf(shifted_onset[k] | -2.5,10,10),
     normal_lpdf(shifted_onset[k] | 0,50));
    }
    
    // culture positivity analysis
    alpha_CP ~ normal(0,20);
    beta_CP ~ normal(0,3);
    Y_CP ~ bernoulli_logit_glm(X_CP,alpha_CP, beta_CP);
  }
}

generated quantities {
  real slope_up_mu = exp(log_slope_up_mu);
  real slope_down_mu = exp(log_slope_down_mu);
  real intercept_mu = exp(log_intercept_mu);
  real time2peak_mu = intercept_mu/slope_up_mu;
}
