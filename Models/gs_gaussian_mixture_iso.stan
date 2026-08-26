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
    int<lower=1> n_classes;                   // Number of categories
    int n_non_missing;                        // Sample size
    int n_missing;                            // Number of missing categories to be estimated
    int dim;                                  // Dimension
    array[n_classes] vector[dim] lambda;      // Mean prior mean
    vector[n_classes] alpha;
    real tau;                                 // GS temperature parameter
    real<lower = 0> sigma;
    matrix[n_non_missing, dim] y_non_missing; // Data
    matrix[n_missing, dim] y_missing;         // Observations with missing classes
    real<lower=0> upsilon;
    real<lower=0> gamma;
    real<lower=0> eta;
    vector[n_missing + n_non_missing] classes;
    array[n_missing] int classes_missing; 
    simplex[n_classes] pi;
    array[n_non_missing] int classes_non_missing; 
  }

  transformed data{
    vector[n_classes] log_pi = log(pi);
  }
  parameters {
    array[n_classes] vector[dim] mu;
    matrix[n_missing, n_classes-1] u; // Cluster assignments to be estimated
  }

  transformed parameters {
    matrix[n_missing, n_classes] u_softmax;
    matrix[n_missing, n_classes] u_transformed;
    real u_max;
    real sum_y;
    vector[n_missing] sum_u;
    
    u_transformed[:, 1:n_classes-1] = u/tau;
    u_transformed[:, n_classes] = rep_vector(0, n_missing);
  
   for (i in 1:n_missing){
     // Inverse GS -> CG transformation
    u_softmax[i] = log_softmax(u_transformed[i]')';
    }
  }

 model {
    vector[n_classes] mixture;
    row_vector[dim] y_current;
    int c;
    for (i in 1:n_missing){
        y_current = y_missing[i];
      
        for (j in 1:n_classes-1){
          mixture[j] = u_softmax[i,j] + normal_lpdf(y_current | (mu[j])', sqrt(sigma)); 
        }
        mixture[n_classes] = u_softmax[i,n_classes] + normal_lpdf(y_current | (mu[n_classes])',  sqrt(sigma)); 
        target += log_sum_exp(mixture) + gumbel_centered_lpdf(to_vector(u[i]) | n_classes, log_pi);

    }

    // Computing likelihood for observations where classes are KNOWN
    for (i in 1:n_non_missing){
        y_current = y_non_missing[i];
        c = classes_non_missing[i];
        target +=  normal_lpdf(y_current | (mu[c])', sqrt(sigma));
  }
  // Computing priors
    for (i in 1:n_classes){
        target += normal_lpdf(mu[i] | (lambda[i])', eta) ;
  }
 }

 generated quantities {

  vector[n_classes] mixture;
  vector[n_missing] c;
  row_vector[dim] y_current;
  real z_max;
  int c_max;
  int c_known;
  real z_temp;
  real lh_approx = 0;
  real lh_target = 0;
  real lh_ratio = 0;
  vector[n_missing] log_lik;
  real lh_temp;
  real guide_lp = 0;
  
  // Performing simplex projection step
  for (i in 1:n_missing){
    z_max = 0;
    c_max = 0;
    y_current = y_missing[i];
    
    for (j in 1:(n_classes)){
      z_temp = exp(u_softmax[i,j]);
       if ( z_temp> z_max){
          z_max = z_temp;
          c_max = j;
      }
      mixture[j] = u_softmax[i, j] + normal_lpdf(y_current | (mu[j])', sqrt(sigma));
    }
    // Outputing log density terms for external MH correction
    lh_temp = normal_lpdf(y_current | (mu[c_max])', sqrt(sigma));
    lh_target += lh_temp;
    guide_lp += lh_temp;
    guide_lp += gumbel_centered_lpdf(to_vector(u[i]) | n_classes, log(pi));
    c[i] = c_max;  
    

   lh_temp = log_sum_exp(mixture);
   log_lik[i] = lh_temp;
   lh_approx += lh_temp;
   
  }
  lh_ratio = lh_target - lh_approx;
  
  for (i in 1:n_non_missing){
        y_current = y_non_missing[i];
        c_known = classes_non_missing[i];
        guide_lp +=  normal_lpdf(y_current | (mu[c_known])', sqrt(sigma));
  }
  for (j in 1:n_classes){
    guide_lp += normal_lpdf(mu[j] | (lambda[j])', eta) ;
  }
 }
