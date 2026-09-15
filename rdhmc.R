library(bridgestan)
library(cmdstanr)
library(parallel)
library(jsonlite)
#source("~/Project/hmc-sandbox/R/NUTS_helpers.R")

SEED <- 123
MAX_TREEDEPTH = 10

# Initializing Stan model and extracting initial state, adapted mass matrix and  
# step size, together with names of the constrained parameter components. 
initialize_model <- function(model_object, warm_up, n_chains, max_treedepth){
  # Initializing Cmdstan model
  fit <- model_object$model$cmdstan_model$sample(
    data = model_object$data,
    parallel_chains = n_chains,
    chains = n_chains,
    iter_warmup = warm_up,
    iter_sampling = 1,
    sig_figs = 18,
    max_treedepth = max_treedepth,
    save_warmup = TRUE
  )
  
  # Extracting mass matrix, initial parameter values and step size
  n_unconstr_params <- length(model_object$model$bs_model$param_unc_names())
  is_matrix = if_else(n_unconstr_params > 1, TRUE, FALSE)
  mass_matrix = fit$inv_metric(matrix = is_matrix)
  theta_0 <- fit$draws(format = "list")
  step_size <- fit$metadata()$step_size_adaptation
  
  # Extracting constrained parameter components
  all_params <-  fit$summary()$variable
  all_params <- str_replace_all(str_replace_all(all_params,"[\\[\\,]", "\\."), "\\]", "") #Ensuring that names match between cmdstan and bridgestan
  constr_params_idx <- which(all_params %in% model_object$param_names)
  theta_0 <- lapply(theta_0, function(x) x[constr_params_idx])

  return(list(mass_matrix= mass_matrix,
              theta_0 = theta_0,
              step_size = step_size,
              param_names = model_object$param_names,
              fit = fit
  ))
}

# Performs T leapfrog steps given intial state-momentum pair (q,p)
leapfrog_nsteps_diag <- function(model, q, p, eps, T, inv_mass_matrix) {

  q <- as.numeric(q)
  p <- as.numeric(p)
  inv_mass_matrix <- as.matrix(inv_mass_matrix)
  half_eps     <- 0.5 * eps
  eps_inv_mass <- eps * inv_mass_matrix 
 
  # Initial gradient (used for first half-step)
  grad <- model$log_density_gradient(q)
  g <-  as.numeric(grad$gradient)
  for (i in 1:T) {
    p <- p + half_eps * g
    q <- q + drop(eps_inv_mass %*% p)
    # Comuting gradients using Bridgestan
    grad <- model$log_density_gradient(q)
    g <- as.numeric(grad$gradient)
    p <-     p + half_eps * g
  }
  list(q = q, p = -p, log_dens = grad$val)
}

# Initializing RDHMC sampler
init_rdhmc <- function(model_file, data, warmup_init, n_chains, max_treedepth = MAX_TREEDEPTH, model_seed = SEED){

  # Defining Stan model objects 
  cmdstan_model <- cmdstan_model(model_file)
  bs_model <-  bridgestan::StanModel$new(model_file, toJSON(data, auto_unbox = TRUE, digits = 16), model_seed, warn = FALSE)
  model_object <- list(model = list(bs_model = bs_model, 
                                    cmdstan_model = cmdstan_model), 
                                    param_names = bs_model$param_names(), 
                                    data = data,
                                    model_seed = model_seed)
  # Initializing Stan model
  init <- initialize_model(model_object, warmup_init, n_chains, max_treedepth)
  model <- model_object$model$bs_model
  init$rng <-  model$new_rng(model_object$model_seed)
  
  # Extracting initial parameter vector and corresponding generated quantities from the intialized Stan model
  theta_0_unconstr <- lapply(init$theta_0, function(x) model$param_unconstrain(x))
  gq0 <-  lapply(theta_0_unconstr, function(x) model$param_constrain(x, include_tp = TRUE, include_gq = TRUE, rng=  init$rng))
  dim_constr <- length(init$theta_0[[1]])
  init$dim_unconstr <- length(theta_0_unconstr[[1]])
  
  # Ensuring that parameter names match between Cmdstan and Bridgestan models
  param_names <- model$param_names()
  param_names <- str_replace_all(str_replace(param_names, "\\.", "[" ), "\\.", "\\,")
  param_names[str_detect(param_names, "\\[")] <- str_c(param_names[str_detect(param_names, "\\[")], "]")
  init$param_names <- param_names
  
  # Extracting initial state 
  param_names_full <- model$param_names(include_tp = TRUE, include_gq = TRUE)
  init$log_dens_idx <- which(param_names_full == "guide_lp")
  init$z_idx <- which(str_detect(param_names_full, "z_trans\\."))

  init$c_idx <- which(str_detect(param_names_full, "^c(\\.\\d+)?$"))
  init$param_unc_names <- model$param_unc_names()
  
  # Ensuring name consitency for inital state
  param_names_full <- str_replace_all(str_replace(param_names_full, "\\.", "[" ), "\\.", "\\,")
  param_names_full[str_detect(param_names_full, "\\[")] <- str_c(param_names_full[str_detect(param_names_full, "\\[")], "]")
  init$param_names_full <- param_names_full
  
  # Extracting and inverting the estimated mass matrix
  inv_mass_matrix <- init$mass_matrix
  mass_matrix <- lapply(inv_mass_matrix, function(x) solve(x))
  lf_fit_window <- 10
  diag <- init$fit$sampler_diagnostics(format = "list", inc_warmup = TRUE)

  # Computing integration length based on number of leaprog steps of the initialized model
  n_leapfrogs <- lapply(diag, function(x) x$n_leapfrog__)
  L <- lapply(n_leapfrogs, function(x){
      floor(as.numeric(min(x[(length(x)-lf_fit_window):length(x)])))
  })
  
  # Preparing initialized chains for parallel sampling 
  chain_inits <- mapply(function(x,y,z,w,v,u,t) list(theta_0 = x, 
                                                       gq0 = y, 
                                                       mass_matrix = z,  
                                                       inv_mass_matrix = w, 
                                                       step_size = v,
                                                       n_leapfrogs = u,
                                                       L = t), 
                        theta_0_unconstr, 
                        gq0, 
                        mass_matrix,
                        inv_mass_matrix,
                        init$step_size, 
                        n_leapfrogs,
                        L, SIMPLIFY = FALSE)

  return(list(init = init, model_object = model_object, chain_inits = chain_inits))
}

# Helper function for producing a proposed sample using RDHMC based on a initialized stan model bs_model and data hmc_data, 
# with initial states chain_init starting from a previously accepted state z_prev.
rdhmc_proposal <- function(bs_model, hmc_data, init,  chain_init, z_prev, param_names){
  
  # Drawing momentum and computing kinetic energy
  d  <- init$dim_unconstr
  sd <- sqrt(diag(chain_init$mass_matrix))
  p  <- rnorm(d, mean = 0, sd = sd)
  K_prev <- 0.5*t(p)%*%chain_init$inv_mass_matrix%*%p
 
  # Performing L leapfrog steps
  L <- chain_init$L
  lf_step <- leapfrog_nsteps_diag(bs_model, z_prev, p, chain_init$step_size, L, chain_init$inv_mass_matrix)
  
  # Extracting proposed RDHMC state and computing Hamiltonian
  theta_prop <- lf_step$q
  p_prop <- lf_step$p
  theta_gq <- bs_model$param_constrain(theta_prop, include_tp = TRUE, include_gq = TRUE, rng=  init$rng)
  z_trans <- theta_gq[init$z_idx]
  c_prop <- theta_gq[init$c_idx]
  U_prop <- -theta_gq[init$log_dens_idx]
  K_prop <- drop(0.5*t(p_prop)%*%(chain_init$inv_mass_matrix)%*%p_prop)
  H_prop <- as.numeric(U_prop + K_prop)
  
  return(list(U_prop = U_prop, 
              K_prop = K_prop, 
              K_prev = K_prev,
              H_prop = H_prop,
              theta_prop = theta_prop,
              z_trans = z_trans,
              c_prop = c_prop))
}

# Main function conducting RDHMC sampler based on an initialized Stan model init and initialized chains chain_init
# Producing n_samples using RDHMC using a warmup period of warmup_sampler samples and fixed integration length L. 
rdhmc <- function(init, chain_init, n_samples, warmup_sampler, L = NA){
  
  # Extracting hyperparameters and model states from the initialized model object
  bs_model <- init$model_object$model$bs_model
  gq0 <- chain_init$gq0
  theta_0_unconstr <- chain_init$theta_0
  hmc_init <- init$init
  c_prev <- gq0[hmc_init$c_idx]
  U_prev <- -bs_model$log_density(theta_0_unconstr, propto = FALSE, jacobian = TRUE)
  theta_prev <- theta_0_unconstr
  lp_prev <- 0
  param_names <- hmc_init$param_names
  param_names_full <- hmc_init$param_names_full
  accept <- FALSE    
  accepts <- rep(0, n_samples)
  lp <- rep(0, n_samples)
  accept_probs <- rep(0, n_samples)
  n_classes <- length(hmc_init$c_idx)  
  
  # Initializing matrices for the categorical and continuous RDHMC samples respectively
  c_samples <- matrix(0, n_samples, n_classes)
  samples <- matrix(0, n_samples, length(param_names_full))
 
  # Using the integration length from the initialized model if no input is given
  if (!is.na(L)){
    chain_init$L <- L
  }

  # Performing RDHMC sampling
  for (i in 1:(n_samples + warmup_sampler)){
  
    # Proposal step
    hmc_step <- tryCatch(
      {rdhmc_proposal(bs_model, hmc_data, hmc_init, chain_init, theta_prev, param_names)}
      , error = function(err) {
        print(paste0("Energy error: Rejecting sample", err ))
        hmc_step <- list()
        hmc_step$U_prop <- U_prev
        hmc_step$theta_prop <- theta_prev
        hmc_step$K_prev <- 0
        hmc_step$H_prop <- Inf
        return(hmc_step)}  
    )
  
    # Computing acceptance probability
    H_prev <- as.numeric(U_prev + hmc_step$K_prev)
    U_prop <- hmc_step$U_prop
    theta_prop <- hmc_step$theta_prop
    lp_prop <- 1
    target_ratio <- -hmc_step$H_prop + H_prev 
    if (is.nan(target_ratio)){
      print("Energy error II: Rejecting sample")
      target_ratio <- -Inf
      
    } 
    c_prop <- hmc_step$c_prop
    accept_prob <- min(0, target_ratio)
    
    # Performing Metropolis-Hastings acceptance
    if (log(runif(1)) < accept_prob){
      theta_prev <- theta_prop
      c_prev <- c_prop
      lp_prev <- lp_prop
      accept <- TRUE
      U_prev <- U_prop
    }
    
    if(i%%10==0) print(i)
    
    # Storing samples after warmup period has concluded
    if (i > warmup_sampler){
      c_samples[i - warmup_sampler, ] <- c_prev
      accepts[i - warmup_sampler] <- accept
      lp[i - warmup_sampler] <- lp_prev
      accept_probs[i - warmup_sampler] <- accept_prob
      samples[i - warmup_sampler,] <- bs_model$param_constrain(theta_prev, include_tp = TRUE, include_gq = TRUE, rng=  hmc_init$rng)
    }
    accept <- FALSE
  }
  theta_samples <- as.data.frame(samples)
  names(theta_samples) <- param_names_full
  c_samples_df <- as.data.frame(c_samples)
  
  names(c_samples_df) <- str_c("c[", 1:ncol(c_samples_df), "]")
  
  return(list(theta_samples = theta_samples,
              class_samples = c_samples_df,
              samples = bind_cols(theta_samples, c_samples_df),
              lp = lp,
              accepts = accepts))
}

# Wrapper for the RDHMC sampling function. Takes in the Stan model file name model_file and data model_data.
# Producing n_samples samples using RDHMC, using n_chains chains and a warmup period for the main sampler and initialization sampler
# of warmup_sampler and warmup_init samples respectively. 
rdhmc_sample <- function(model_file, model_data, n_chains, n_samples, warmup_sampler, warmup_init, max_treedepth = MAX_TREEDEPTH,  L= NA, model_seed = SEED){

  # Initializing Stan model 
  hmc_init <- init_rdhmc(model_file, model_data, warmup_init, n_chains, max_treedepth, model_seed)
  hmc_init$n_samples <- n_samples
  hmc_init$warmup_sampler <- warmup_sampler
  hmc_init$warmup_init <- warmup_init
  hmc_init$n_chains <- n_chains
  
  # Performing RDHMC sampling, using parallel chains. 
  if (n_chains > 1){
    sampler_wrapper <- function(i){rdhmc(hmc_init, hmc_init$chain_inits[[i]], n_samples, warmup_sampler, L) }
    experiments <- mclapply(1:n_chains, sampler_wrapper, mc.cores = n_chains)
    experiments_list <- lapply(experiments, function(x)x$samples)
    samples_df <- bind_rows(experiments_list, .id = ".chain")
    
  } else {
    experiments <- rdhmc(hmc_init, hmc_init$chain_inits[[1]], n_samples, warmup_sampler, L)
    samples_df <- experiments$samples
    samples_df[".chain"] <- 1
  }
  
  return(list(experiments = experiments, init = hmc_init, samples_df = samples_df))
  
}
