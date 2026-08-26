 // Stan model for SimplexMCMC sampling from a Gumbel-Softmax distribution with categorical distribution pi and temperature tau
  functions {
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
  log_z = log_softmax(u_transformed);
  z_trans = exp(log_z);
}

model {
  target += gumbel_centered_lpdf(z | K, log_pi);
}

generated quantities {
  
  int c; 
  real z_max = 0;
  real guide_lp = gumbel_centered_lpdf(z | K, log_pi);
  // Implementing projection step for SimplexMCMC
  for (i in 1:K){
    if (z_trans[i] > z_max){
      z_max = z_trans[i];
      c = i;
    }
  }
}
