library(cmdstanr)
library(LaplacesDemon)
library(VGAM)

stan_file_dir <- "data {
  int K;
  vector[K] alpha;
}
parameters {
  simplex[K] z;
}
model {
  target += dirichlet_lpdf(z| alpha);
}

generated quantities {
  real log_dens = dirichlet_lpdf(z | alpha); // constant across draws for alpha=1
}
"
dir_model <- cmdstan_model(write_stan_file(stan_file_dir))


stan_file_ln <- "data {
  int<lower=2> K;
  vector[K-1] mu;
  real<lower=0> sigma;
}
parameters {
  vector[K-1] z;                         // ALR coords
}
transformed parameters {
  simplex[K] z_trans = softmax(append_row(z, 0)); // map to simplex
}
model {
  z ~ normal(mu, sigma);                 // logistic-normal prior via z
}
generated quantities {
  real lp_z  = normal_lpdf(z | mu, sigma); // log p(z) on R^(K-1)
  // For the density on the simplex: p(z_trans) = p(z) * |det d(ALR)/d(pi)|
  // log |J_ALR(z_trans)| = -sum(log(z_trans))
  real log_dens = lp_z- sum(log(z_trans));        // log p(pi) (logistic-normal on simplex)
}
"

ln_model <- cmdstan_model(write_stan_file(stan_file_ln))

stan_file_gs <- " functions {
  
  real gumbel_centered_lpdf(vector y, int K, vector x){
    vector[K] u;
    u[1:K-1] = y;
    u[K] = 0;
    
    return lgamma(K) + sum(x - u) - K * log_sum_exp(x - u);
  }
  
  real gumbel_softmax_lpdf(vector log_z, int K, vector log_pi, real tau){
    return lgamma(K) +(K-1)*log(tau) - K*log_sum_exp(log_pi - tau*log_z) + sum(log_pi - (tau + 1)*log_z);
  }
}

data {
  int<lower=1> K;  // Number of categories
  simplex[K] pi;   // Probability vector
  real tau;        // Temperature parameter
}

transformed data {
  
  vector[K] log_pi = log(pi);
}

parameters {
  vector[K-1] z;
}

transformed parameters {
  
  vector[K] z_trans;
  vector[K] log_z; 
  real u_max;
  real sum_y;
  vector[K] u_transformed;
  
  // Transforming CS -> GS
  u_transformed[1:K-1] = z/tau;
  u_transformed[K] = 0;
 // u_max = fmax(max(u_transformed),0);
  //sum_y = sum(exp(u_transformed - u_max));
  
  
  //z_trans[K] = 0;
  //for (k in 1:K-1) {
  //  z_trans[k] = exp(u_transformed[k] - u_max)/sum_y;
    
//  }
  
  //z_trans[K] = 1 - sum(z_trans);
  
  log_z = log_softmax(u_transformed);
  
  z_trans = exp(log_z);
  
}

model {
  
  target += gumbel_centered_lpdf(z | K, log_pi);
  
}

generated quantities {
  
  int c; 
  real z_max = 0;
  //vector[K] log_z = log(z_trans);
  real log_dens = gumbel_softmax_lpdf(log_z | K, log_pi, tau);
  for (i in 1:K){
    
    if (z_trans[i] > z_max){
      z_max = z_trans[i];
      c = i;
    }
  }
  
  
}

"
gs_model <- cmdstan_model(write_stan_file(stan_file_gs))
#gs_mod <- cmdstan_model("~/bridgestan/test_models/gs/gs.stan")