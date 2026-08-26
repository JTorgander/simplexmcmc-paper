functions {
    // Returning log density of centered Gumbel distribution
    real gumbel_centered_lpdf(vector y, int n_classes, vector x){
    vector[n_classes] u;
    u[1:n_classes-1] = y;
    u[n_classes] = 0;

    return lgamma(n_classes) + sum(x - u) - n_classes * log_sum_exp(x - u);
    }
  }

  data {
    int<lower=1> n_classes;                    // Number of categories
    int n_non_missing;                         // Sample size
    int n_missing;                             // Number of missing categories to be estimated
    int dim;                                   // Dimension
    array[n_classes] vector[dim] lambda;       // Mean prior mean
    vector[n_classes] alpha;
    real tau;                                  // GS temperature parameter
    matrix[n_non_missing, dim] y_non_missing;  // Data
    matrix[n_missing, dim] y_missing;          // Observations with missing classes
    real<lower=0> upsilon;
    real<lower=0> gamma;
    real<lower=0> eta;
    vector[n_missing + n_non_missing] classes;
    array[n_missing] int classes_missing; 
    array[n_non_missing] int classes_non_missing; 
  }

  parameters {
    array[n_classes] vector[dim] mu;
    simplex[n_classes] pi;                     // Probability vector
    array[n_classes] cholesky_factor_corr[dim] L_Omega;
    array[n_classes] vector[dim] log_sigma;    // Scale
    matrix[n_missing, n_classes-1] u;          // GS variables in unconstrained space
  }

  transformed parameters {
    // Transforming from GS to CG-distribution
    matrix[n_missing, n_classes] u_softmax;
    matrix[n_missing, n_classes] u_transformed;
    array[n_classes] matrix[dim, dim] L_Sigma;
    array[n_classes] vector[dim] sigma; 
    for (k in 1:n_classes){
      sigma[k] = exp(log_sigma[k]);
      L_Sigma[k] = diag_pre_multiply(sigma[k], L_Omega[k]);  // Lower-triangular
      
    }
    
    real u_max;
    real sum_y;
    vector[n_missing] sum_u;

    for (i in 1:n_missing){
     // Inverse GS -> CG transformation
     u_transformed[i, 1:n_classes-1] = u[i]/tau;
     u_transformed[i, n_classes] = 0;
     u_max = max(u_transformed[i]);
     sum_y = sum(exp(u_transformed[i] - u_max));
     sum_u[i] = log_sum_exp(u_transformed[i]);
     u_softmax[i, n_classes] = 0;
    for (k in 1:n_classes-1) {
      u_softmax[i, k] = exp(u_transformed[i, k] - u_max)/sum_y;
    }
    u_softmax[i, n_classes] = 1 - sum(u_softmax[i]);
    }
  }

 model {

    vector[n_classes] mixture;
    row_vector[dim] y_current;
    int c;
    vector[n_classes] log_pi = log(pi);
    
    // Computing likelihood for observations where classes are unknown
    for (i in 1:n_missing){
        y_current = y_missing[i];
        for (j in 1:n_classes-1){
          mixture[j] = u_transformed[i, j] - sum_u[i] + multi_normal_cholesky_lpdf(to_vector(y_current) | to_vector(mu[j]), L_Sigma[j]);
        }
        mixture[n_classes] =  (-1)*sum_u[i] + multi_normal_cholesky_lpdf(to_vector(y_current) | to_vector(mu[n_classes]), L_Sigma[n_classes]);
        target += log_sum_exp(mixture) + gumbel_centered_lpdf(to_vector(u[i]) | n_classes, log_pi);
    }
    // Computing likelihood for observations where classes are KNOWN
    for (i in 1:n_non_missing){
        y_current = y_non_missing[i];
        c = classes_non_missing[i];
        target +=  multi_normal_cholesky_lpdf(to_vector(y_current) | to_vector(mu[c]), L_Sigma[c]);
  }
    // Computing priors
    for (i in 1:n_classes){
        target += multi_normal_lpdf(to_vector(mu[i]) | to_vector(lambda[i]),  diag_matrix(rep_vector(eta, dim))) ;
        target += normal_lpdf(log_sigma[i] | 0, gamma);
        target += lkj_corr_cholesky_lpdf(L_Omega[i] | upsilon);
    
  }
   target += dirichlet_lpdf(pi| alpha);
 }

 generated quantities {

  vector[n_classes] mixture;
  vector[n_missing] c;
  row_vector[dim] y_current;
  matrix[n_missing, dim] y_pred;
  real z_max;
  int c_known;
  int c_max;
  real lp_temp;
  real z_temp;
  real log_dens = 0;
  real lh_approx = 0;
  real lh_target = 0;
  real lh_ratio = 0;
  real guide_lp = 0;
  
  // Performing simplex projection step
  for (i in 1:n_missing){
    z_max = 0;
    c_max = 0;
    y_current = y_missing[i];
    
    for (j in 1:(n_classes)){
      z_temp = u_softmax[i,j];
       if ( z_temp> z_max){
          z_max = z_temp;
          c_max = j;
      }
    }
  c[i] = c_max;  
  
  // Outputing log density terms for external MH correction
  lp_temp = gumbel_centered_lpdf(to_vector(u[i]) | n_classes, log(pi));
  log_dens += lp_temp;
  guide_lp += lp_temp;
  lp_temp = multi_normal_cholesky_lpdf(to_vector(y_current) | to_vector(mu[c_max]), L_Sigma[c_max] );
  lh_target = lp_temp;
  guide_lp += lp_temp;
  y_pred[i] = to_row_vector(multi_normal_cholesky_rng(to_vector(mu[c_max]), L_Sigma[c_max]));
  for (j in 1:(n_classes-1)){
       mixture[j] = u_transformed[i, j] - sum_u[i] + multi_normal_cholesky_lpdf(to_vector(y_current) | to_vector(mu[j]), L_Sigma[j]);
  }
  mixture[n_classes] =  (-1)*sum_u[i] + multi_normal_cholesky_lpdf(to_vector(y_current) | to_vector(mu[n_classes]),  L_Sigma[n_classes]);
  lh_approx = log_sum_exp(mixture);
  log_dens += lh_target;
  lh_ratio += lh_target - lh_approx;
  }
  // Separating into cases where classes are assumed known/unknown
  for (i in 1:n_non_missing){
    y_current = y_non_missing[i];
    c_known = classes_non_missing[i];
    guide_lp +=  multi_normal_cholesky_lpdf(to_vector(y_current) | to_vector(mu[c_known]), L_Sigma[c_known]);
  }
  for (j in 1:n_classes){
    guide_lp += multi_normal_lpdf(to_vector(mu[j]) | to_vector(lambda[j]),  diag_matrix(rep_vector(eta, dim))) ;
    guide_lp += lkj_corr_cholesky_lpdf(L_Omega[j] | upsilon);
    guide_lp += normal_lpdf(log_sigma[j] | 0, gamma);
  }
  
   guide_lp += dirichlet_lpdf(pi| alpha);
 }
