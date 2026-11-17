functions {

    real gumbel_centered_lpdf(real y, int n_classes, vector x){
    vector[n_classes] u;
    u[1] = y;
    u[n_classes] = 0;

    return lgamma(n_classes) + sum(x - u) - n_classes * log_sum_exp(x - u);
    }
  }

  data {
                    // Number of categories
    int N;                 // Number of missing categories to be estimated
  
  
    vector[2] alpha;
    real tau;                       // Temperature parameter
    vector[N] y;                 // Data
    int n_coefs;
    matrix[N, n_coefs] X;

     real lambda;              // slab variance
     real upsilon;
     real s;
     simplex[2] pi;
  }

  parameters {
  real beta0;
  vector[n_coefs] beta;
  vector[n_coefs] u;
  // simplex[2] pi;   // Probability vector
 
    //real<lower=0> gamma;
    real sigma_unconstr;

  }

  transformed parameters {
    // Transforming from GS to CG-distribution
    //vector[1] beta_filtered;
    matrix[n_coefs, 2] u_softmax;
    matrix[n_coefs, 2] u_transformed;
    vector[2] log_pi = log(pi);
    vector[n_coefs] beta_z;
    real sigma = exp(sigma_unconstr); // Cluster width
   

  
     // Inverse GS -> CG transformation
    u_transformed[:,1] = u/tau;
    u_transformed[:,2] = rep_vector(0, n_coefs);
    
   for (j in 1:n_coefs){
     
     u_softmax[j] = softmax(u_transformed[j]')';
     
     
   }
   
   beta_z = beta .* u_softmax[:, 1];
   
   
   // beta_filtered = beta * u_softmax[1];
  }


model {
  int n_classes = 2; 
  vector[2] mixture;
  real z;
  vector[N] mu;

 // for (i in 1:N){

   //   mu = beta0;
    //for (j in 1:n_coefs){
      
    //  z = u_softmax[j,1];
    //  mu += z * X[i, j] * beta[j];
      

  //  }
    
  //  target +=   normal_lpdf(y[i]| mu, sqrt(sigma));

  //}
  
  //mu = X * beta_z + rep_vector(beta0, N);
  
  // y ~ multi_normal(mu , sigma * identity_matrix(N));
  y ~ normal_id_glm(X, beta0, beta_z, sqrt(sigma));

  for (j in 1:n_coefs){
    target += gumbel_centered_lpdf(u[j] | 2, log_pi);
  }
  
 // target +=  multi_normal_lpdf(beta | rep_vector(0, n_coefs), lambda*identity_matrix(n_coefs));
  target +=  multi_normal_lpdf(beta | rep_vector(0, n_coefs), lambda*sigma*identity_matrix(n_coefs));
  target += normal_lpdf(beta0 | 0, sqrt(sigma*lambda));
 // target += dirichlet_lpdf(pi| alpha);
  target +=  scaled_inv_chi_square_lpdf(sigma| upsilon, s) + sigma_unconstr;
  
}


 generated quantities {

   vector[n_coefs] c;
   real z;
   vector[N] log_lik;
   int c_max;
   real lh_approx = 0;
   real lh_target = 0;
   real lh_ratio = 0;
   real guide_lp = 0;
   vector[N] mu_z;
   vector[N] mu_c;
   real lh_temp;
  
  vector[n_coefs] beta_c;
   
   for (j in 1:n_coefs){
      
    
       if ( u_softmax[j,1] > 0.5){
          
          c[j] = 1;
      } else {
        
        c[j] = 0;
      }
    
    }
    beta_c = beta.* c;
    mu_z  = X * beta_z + rep_vector(beta0, N);
    mu_c  = X * beta_c + rep_vector(beta0, N);
    
    //lh_target = multi_normal_lpdf(y | mu_c , sigma * identity_matrix(N));
    
    for (i in 1:N){
      
      lh_temp = normal_lpdf(y[i]| mu_c[i], sqrt(sigma));
      log_lik[i] = lh_temp;
      lh_target += lh_temp;
      lh_approx += normal_lpdf(y[i] | mu_z[i], sqrt(sigma));
      guide_lp += lh_temp;
    }
   
    
    
    lh_ratio = lh_target - lh_approx;
    
    for (j in 1:n_coefs){
    guide_lp += gumbel_centered_lpdf(u[j] | 2, log_pi);
  }
  
   guide_lp +=  multi_normal_lpdf(beta | rep_vector(0, n_coefs), lambda*sigma*identity_matrix(n_coefs));
  guide_lp += normal_lpdf(beta0 | 0, sqrt(sigma*lambda));
 // target += dirichlet_lpdf(pi| alpha);
  guide_lp +=  scaled_inv_chi_square_lpdf(sigma| upsilon, s) + sigma_unconstr;

  
 }