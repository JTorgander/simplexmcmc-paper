 // Stan model for sampling from a Logistic Normal distribution with a full covariance matrix
data {
  int<lower=2> K;
  vector[K-1] mu;
  vector[K-1] Sigma_diag;
  vector[(K-1)*(K-2) %/% 2]  Sigma_tri;
}

transformed data {
  matrix[K-1, K-1] Sigma;
  int count = 1;

  // 1. Fill the diagonal
  for (k in 1:(K-1)){
    Sigma[k, k] = Sigma_diag[k];
  }
  // 2. Fill the off-diagonals (Symmetric)
  for (i in 2:(K-1)) {
    for (j in 1:(i - 1)) {
      Sigma[i, j] = Sigma_tri[count];
      Sigma[j, i] = Sigma_tri[count];
      count += 1;
    }
  }
  matrix[K-1, K-1] L_Sigma = cholesky_decompose(Sigma);
}
parameters {
  vector[K-1] z;                         // ALR coords
}
transformed parameters {
  simplex[K] z_trans = softmax(append_row(z, 0)); // map to simplex
}
model {
  z ~ multi_normal_cholesky(mu, L_Sigma);                // logistic-normal prior via z
}
generated quantities {
  real lp_z  = multi_normal_lpdf(z | mu, Sigma); // log p(z) on R^(K-1)
  // For the density on the simplex: p(z_trans) = p(z) * |det d(ALR)/d(pi)|
  // log |J_ALR(z_trans)| = -sum(log(z_trans))
  real log_dens = lp_z- sum(log(z_trans));        // log p(pi) (logistic-normal on simplex)
}