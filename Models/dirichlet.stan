 // Stan model for sampling from a Dirichlet distribution using SimplexMCMC
data {
  int K;
  vector[K] alpha;
}
parameters {
 vector[K-1] z;
}

transformed parameters{
   simplex[K] z_trans = softmax(append_row(z, 0));
}
model {
  target += dirichlet_lpdf(z_trans| alpha)  + sum(log(z_trans)) ;
}

generated quantities {
  
  real guide_lp = dirichlet_lpdf(z_trans | alpha) + sum(log(z_trans)); 
}
