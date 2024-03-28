//// stan code for Baysian CAP regression 

data {
  int<lower=0> p;  // number of regions 
  int<lower=1> d;  // rank of the latent components 
  int<lower=1> q;  // number of covariates (excluding intercept)
  int<lower=0> N;  // number of data points  
  int<lower=0> n;  // number of subjects  
  array[N] int<lower=1, upper=n> id; // each data point has an associated grouping factor 
  array[N] vector[p] Y; // p dimensional signals   
  array[N] vector[q] X; // q dimensional covariates
  matrix[p,p] whitening_mat;  // S_bar^{-1/2} (a matrix used to whiten the signal)
  vector[n] t_vec;  // number of time points (a vector of length n, one number per subject)
  array[n] matrix[p,p] S;   // an array of subject-specific covariance estimates (for DfD computation)
  real<lower=0> eta;  //  random effect covariance LKJ hyperparameter, specifying the amount of expected correlations
  real<lower=0> beta0_sd;  // sd of the intercept of linear predictor
  real<lower=0> beta_sd;   // sd of the regression coefficients of linear predictor
  int<lower=0,upper=1> lp_re;  // indicator of whether random effect is included in the linear predictor
}


transformed data{
  array[N] vector[p] Y_s; 
  for (i in 1:N) Y_s[i] = whitening_mat*Y[i]; // scaled Y: scaled by the sqaured root of the marginal covariance (S_bar^{-1/2})
}


parameters {
  vector[p*d] U;     // to represent Gamma
  matrix[d,q] beta;  // main regression coefficients (excluding the intercept)
  vector[d] beta0;   // intercept
  matrix[d,n] z;     // matrix of subject-specific (d-dimensional) random effects 
  cholesky_factor_corr[d] L_Omega;  // random effect correlation matrix cholesky factor
  vector<lower=0>[d] Omega_sd;      // random effect covariance matrix sd vector
}

transformed parameters{ 
  matrix[p,d] Gamma;  // (p-by-d) orthonormal matrix Gamma 
  {
    matrix[p,d] U_mat;
    vector[d] eval;
    vector[d] eval_trans;
    matrix[d,d] evec;
    U_mat = to_matrix(U, p, d);
    eval = eigenvalues_sym(U_mat'*U_mat);
    for(k in 1:d){
      eval_trans[k] = 1/sqrt(eval[k]);
    }
    evec = eigenvectors_sym(U_mat'*U_mat);
    Gamma = U_mat*evec*diag_matrix(eval_trans)*evec'; 
  }
}


model{
  array[N] vector[d] y;  // the projected signals -- this is a local variable (within this block). 
  array[N] vector[d] lp; // (d-dimensional) linear predictor 
  
  // priors   
  to_vector(z) ~ normal(0, 1);  // latent random effect prior
  to_vector(beta)  ~ normal(0, beta_sd);    //normal(0, 2.5);    
  to_vector(beta0) ~ normal(0, beta0_sd);   //normal(0, 2.5);  
  Omega_sd ~  cauchy(0, 1);  
  L_Omega ~ lkj_corr_cholesky(eta); // eta =1 //lkj_corr_cholesky(1); //lkj_corr_cholesky(2); 
  U ~ normal(0,1);   //elementwise standard normal
  
  // projected signal specification 
  for (i in 1:N) {
    y[i] = to_vector(Gamma'*Y_s[i]);   
  }
  
  // linear predictor (latent log-variance) specification 
  for (i in 1:N){
    if(lp_re==0) lp[i] = beta0 + beta * X[i];  
    else lp[i] = beta0 + beta * X[i] + diag_pre_multiply(Omega_sd, L_Omega)* z[,id[i]];  
  }
  
  // heteroscedasticity model specification      
  for (i in 1:N) {
    y[i] ~ normal(0, exp(lp[i]/2)+0.01);  // small value (0.01) is added to avoid a degenerative distribution
  } 
} 
 
generated quantities {
  corr_matrix[d] Omega_cor = multiply_lower_tri_self_transpose(L_Omega);
  matrix[d,d] Omega = quad_form_diag(Omega_cor, Omega_sd); 
  vector[d] eval = eigenvalues_sym(tcrossprod(Gamma'*whitening_mat));  
  matrix[d,d] evec = eigenvectors_sym(tcrossprod(Gamma'*whitening_mat));  
  vector[d] beta0_org = beta0 - diagonal(evec*diag_matrix(log(eval))*evec'); 
  vector[N] log_lik; 
  vector[N] log_lik1; 
  vector[N] log_lik_null;  
  {
    array[N] vector[d] lp0; // (d-dimensional) linear predictor without random effects
    array[N] vector[d] lp;  // (d-dimensional) linear predictor inclduing random effects
    for(i in 1:N){
      lp0[i] = beta0 + beta * X[i];  
      if(lp_re==0) lp[i] = lp0[i];
      else lp[i] = lp0[i] + diag_pre_multiply(Omega_sd, L_Omega)* z[,id[i]]; 
      log_lik[i] = multi_normal_lpdf(to_vector(Gamma'*Y_s[i]) | rep_vector(0,d), diag_matrix(exp(lp[i])));
      log_lik1[i] = multi_normal_lpdf(to_vector(Gamma'*Y_s[i]) | rep_vector(0,d), diag_matrix(exp(lp0[i])));
      log_lik_null[i] = multi_normal_lpdf(to_vector(Gamma'*Y_s[i]) | rep_vector(0,d), diag_matrix(rep_vector(1,d)));  
    }
  }
  
  real DfD; // DfD computation to assess the diagonal assumption of the dimension-reduced covariances
  {
    array[n] matrix[d,d] Psi_hat;
    vector[n] DfD_vec;
    for (i in 1:n){
      Psi_hat[i] = Gamma'*S[i]*Gamma;
      DfD_vec[i] = t_vec[i]*(log_determinant(diag_matrix(diagonal(Psi_hat[i]))+0.01) - log_determinant(Psi_hat[i]+0.01));
    }
    DfD = mean(DfD_vec);   
  }
}

/// the end of the code 
