 // Stan model for SimplexMCMC sampling from a Logistic Normal distribution with independent components 
data {
  int<lower=2> K;
  vector[K-1] mu;
  real<lower=0> sigma;
}
parameters {
  vector[K-1] z;                         // ALR coords
}
transformed parameters {
  simplex[K] z_trans = softmax(append_row(z, 0)); // Map to simplex
}
model {
  z ~ normal(mu, sigma);                 // Logistic-normal prior via z
}
generated quantities {
  real lp_z  = normal_lpdf(z | mu, sigma); // Log p(z) on R^(K-1)
  real guide_lp  = normal_lpdf(z | mu, sigma); 
}
