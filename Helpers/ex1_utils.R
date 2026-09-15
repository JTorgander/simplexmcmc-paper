source("rdhmc.R")

library(cmdstanr)
library(tidyverse)
library(showtext)
library(posterior)





# Defining density and proposal functions for the candidate SimplexMCMC proposal distributions
unif_dens <- function(x, args){ ddirichlet(x,  alpha = rep(1, K))}
dir_guide <-function(x, args){ ddirichlet(x,  alpha = args$alpha)}
cat_dens <- function(x, args){
  return(1)
}

gs_prop <- function(args){ 
  g <- rgumbel(args$K, 0, 1)
  z<- exp((log(args$pi) + g)/args$tau)
  z <- z/sum(z)
  
  return(z)
  
}
unif_prop <- function(args){ rdirichlet(1,  alpha = rep(1, args$K))}
dir_prop <- function(args){rdirichlet(1,  alpha = args$alpha)}
unit_vectors <- function(n) {
  split(diag(n), col(diag(n)))
}
cat_prop <- function(args){sample(unit_vectors(args$K), 1, replace = FALSE, prob = args$pi)[[1]]}
ln_prop <- function(args){
  pi <- args$pi
  mu <- pi
  z <- sapply(mu, function(x) rnorm(1, x, sd = 1))
  z <- c(z, 0)
  z_max <- max(z)
  z_trans <- exp(z - z_max)
  return(z_trans/sum(z_trans))
}


# Defining functions for computing the normalizing lambda_c constants for each of 
# the candidate proposal distributions using numerical integration
lambda_c_dir <- function(alpha, c, rel.tol = 1e-10){
  idx <- setdiff(seq_along(alpha), c)

  f <- function(g){
    # Vectorized in g
    logprod <- rowSums(sapply(alpha[idx], function(a)
      pgamma(g, shape = a, rate = 1, log.p = TRUE)))  # log of product of CDFs
    logf <- dgamma(g, shape = alpha[c], rate = 1, log = TRUE) + logprod
    out <- exp(logf)
    out[!is.finite(out)] <- 0
    out
  }
  upper <- qgamma(1 - 1e-14, shape = alpha[c], rate = 1)  # finite, safe upper bound
  integrate(f, lower = 0, upper = upper, rel.tol = rel.tol)$value
}

get_lambda_dir <- function(alpha){
  p <- vapply(seq_along(alpha), \(c) lambda_c_dir(alpha, c), 1e-10)
  p / sum(p)  # Renormalize tiny drift (≤1e-10) from quadrature
}

lambda_c_ln <- function(mu, c, rel.tol = 1e-10){
  idx <- setdiff(seq_along(mu), c)
  f <- function(g){
    logprod <- rowSums(sapply(mu[idx], function(a)
      pnorm(g, mean = a, sd = 1, log.p = TRUE)))  # log of product of CDFs
    dnorm(g, mean = mu[c], sd = 1) * exp(logprod)
  }
  upper <- qnorm(1 - 1e-16, mean = mu[c], sd = 1)  # finite, safe upper bound
  integrate(f, lower = 0, upper = Inf, rel.tol = rel.tol)$value
}

get_lambda_ln<- function(mu){
  p <- vapply(seq_along(mu), \(c) lambda_c_ln(mu, c), 1e-16)
  logprod <- sum(sapply(mu, function(x) pnorm(0, mean = x, sd = 1, log.p = TRUE)  ))
  p <- c(p , exp(logprod))
  return(p)
}


# Helper functions for updating and computing the total variation distance between the empirical distribution
# of the produced samples and the true target categorical distribution
update_pi_pred <- function(samples, pi_pred, i){
  c <- samples[i]
  if (i ==1){
    pi_pred[c] <- 1
  } else{
    pi_pred <- (i-1)*pi_pred
    pi_pred[c] <-  pi_pred[c] + 1
    pi_pred <- pi_pred/i
  }
  return(pi_pred)
}

get_pi_bias <- function(pi, samples){
  K <- length(pi)
  pi_bias <- rep(0 ,length(samples))
  pi_pred <- rep(0, K)
  for (j in 1:length(samples)){
    pi_pred <- update_pi_pred(samples, pi_pred, j)
    pi_bias[j] <- 0.5*sum(abs(pi_pred - pi))
  }
  return(pi_bias)
}

# Main function running the SimplexMCMC sampler. The function takes in a target categorical distribution
# pi, a string method indicating the type of SimplexMCMC sampler used,  SimplexMCMC proposal function proposal,
# proposal distribution arguments prop_args, guide distribution density function lambda, together with the number
# of desired samples n_samples and warmup period of the intialization and sampling steps warmup_init respective warmup_sampler.
# Based on the value of the method argument, the resulting SimplexMCMC sampler either uses a RDHMC based sampler or 
# a proposal from the exact guide distribution. Returns a list containing samples, convergence rate measured by total variation distance, acceptance rate, number of leapfrog steps
# and (if applicable) the Stan model used for the HMC initalization.
simplex_sampler <- function(pi, method, proposal, prop_args, lambda, n_samples, warmup_init, warmup_sampler){

  # Initializing sampling arguments
  K <- length(pi)
  pi_bias <- rep(0 ,n_samples)
  pi_pred <- rep(0, K)
  lp_prev <- 0
  accepts <- rep(0, n_samples)
  accept_probs <- rep(0, n_samples)
  accept <- FALSE
  n_leapfrogs <- rep(0, n_samples)
  lp <- rep(0, n_samples)
  samples <- rep(0, n_samples)
  
  # Initializing sampler based on the target distribution to be fitted and type of SimplexMCMC sampler
  if (method == "RDHMC"){
    # Initializing RDHMC based sampler
    hmc_data <- prop_args$data
    hmc_init <- init_rdhmc(prop_args$model_file, hmc_data, warmup_init, 1)
    init <- hmc_init$init
    chain_init <- hmc_init$chain_inits[[1]]
    
    model_object <- hmc_init$model_object
    bs_model <- model_object$model$bs_model
    theta_0_unconstr <- chain_init$theta_0
    gq0 <- chain_init$gq0
    z_trans_prev <- gq0[init$z_idx]
    c_prev <- which.max(z_trans_prev)
    U_prev <- -gq0[init$log_dens_idx]
    z_prev <- theta_0_unconstr
    bs_model <- model_object$model$bs_model
    param_names <- init$param_names
    lf_fit_window <- 10
    n_leap_prop <- chain_init$L
    n_leap_prev <- n_leap_prop
  } else { 
    # Initializing exact proposal based sampler
    z_prev <- unif_prop(list(K = K))
    z_trans <- z_prev
    c_prev <- which.max(z_trans)
    lp_prev <- 1
    target_ratio <- 1
    hmc_init <- NA
    n_leap_prev <- NA
    U_prev <- 0
    U_prop <- 0
  }
  
  # Running sampler
  for (i in 1:(n_samples + warmup_sampler)){
    if (method == "RDHMC"){
      hmc_step <- tryCatch(
        {proposal(bs_model, hmc_data, init, chain_init, z_prev, param_names)}
        , error = function(err) {
          print(paste0("Energy error: Rejecting sample", err ))
          hmc_step <- list()
          hmc_step$U_prop <- U_prev
          hmc_step$theta_prop <- z_prev
          hmc_step$z_trans <- z_trans_prev
          hmc_step$K_prev <- 0
          hmc_step$H_prop <- Inf
          return(hmc_step)}  
      )
      H_prev <- as.numeric(U_prev + hmc_step$K_prev)
      U_prop <- hmc_step$U_prop
      z_trans <- hmc_step$z_trans
      z_prop <- hmc_step$theta_prop
      lp_prop <- 1
      target_ratio <- exp(-hmc_step$H_prop + H_prev) 
      if (is.nan(target_ratio)){
        print("Energy error II: Rejecting sample")
        target_ratio <- 0
      } 
    }
    else {
      z_prop <- proposal(prop_args)
      lp_prop <- 1
      n_leap_prop <- NA
      target_ratio <- 1
      z_trans <- z_prop # Enables transforms
    }
    # Acceptance step
    c_prop <- which.max(z_trans)
    r <- (pi[c_prop]/ pi[c_prev]) * (lambda[c_prev] / lambda[c_prop])
    accept_prob <- min(1, r*target_ratio)
   
    if (runif(1) < accept_prob){
      z_prev <- z_prop
      if (method == "RDHMC") z_trans_prev <- z_trans
      c_prev <- c_prop
      lp_prev <- lp_prop
      n_leap_prev <- n_leap_prop
      accept <- TRUE
      U_prev <- U_prop
    }
    # Storing samples and computing convergence metrics
    if (i > warmup_sampler){
      samples[i - warmup_sampler] <- c_prev
      n_leapfrogs[i - warmup_sampler] <- n_leap_prev
      pi_pred <- update_pi_pred(samples, pi_pred, i - warmup_sampler)
      pi_bias[i - warmup_sampler] <- 0.5*sum(abs(pi_pred - pi))
      accepts[i - warmup_sampler] <- accept
      lp[i - warmup_sampler] <- lp_prev
      accept_probs[i - warmup_sampler] <- accept_prob
    }
    accept <- FALSE
  }
  return(list(samples = samples, 
              pi_bias = pi_bias, 
              accepts = accepts, 
              n_leapfrogs = n_leapfrogs, 
              lp = lp, 
              accept_prob = accept_probs, 
              hmc_inits = hmc_init))
}

# Wrapper function for running the simplex sampling experiments. 
# Takes in a function model_getters  intializing each experiment, the number of experiments n_experiments, 
# the number of samples produced from each experiment experiment_sample_size, the target categorical distributions cat_dists, 
# warmup period of the intialization and sampling steps warmup_init respective warmup_sampler,
# and boolean save_runs, indicating if the runs should be saved to disk. Runs the corresponding SimplexMCMC experiments 
# and saves the resulting samples and metrics to disk.
get_simplex_experiments <- function(model_getters, n_experiments, experiment_sample_size, cat_dists, warmup_init, warmup_sampling, save_runs = TRUE){
  
  # Creating experiment folder
  if (save_runs){
    folder_name <-  paste0("Experiments/", Sys.Date())
    if (!dir.exists(folder_name)){
      dir.create(folder_name, recursive = TRUE)
    }
  }
  # Running experiments for each provided target distribution
  out <- data.frame()
  for (dist in cat_dists){
    pi <- dist$pi
    K <- length(pi)
    models <- model_getters(pi, K, dist$alpha0)
    
    for (model in models){
        print(paste0("Fitting model: ", model$model_name))
        if (model$method == "RDHMC"){
        model$proposal_args$model <- cmdstan_model(stan_file= model$proposal_args$model_file)
        }
        # Setting up sampler wrapper for the current experiment and model
        sampler_wrapper <- function(chain_ids){ simplex_sampler(pi = pi, 
                                                                method = model$method, 
                                                                proposal = model$proposal,
                                                                prop_args = model$proposal_args, 
                                                                lambda = model$lambda, 
                                                                n_samples = experiment_sample_size, 
                                                                warmup_init = warmup_init, 
                                                                warmup_sampler = warmup_sampling)}
        # Setting up parallel sampling 
        if (n_experiments > 1){
          experiments <- mclapply(1:n_experiments, sampler_wrapper, mc.cores = n_experiments)
          data <- lapply(experiments, function(x){ x[-length(x)]})
          model_init <- lapply(experiments, function(x){ x$hmc_inits})
          row <- bind_rows(data, .id = "chain")
        } else {
          experiments <- sampler_wrapper(1)
          model_init <- experiments$hmc_inits
          row <- as.data.frame(experiments[-length(experiments)])
          row["chain"] <- 1
        }
      row <- row %>% group_by(chain) %>%  mutate(name = model$model_name, K = K, pi = dist$name, iter = row_number()) %>% transmute(name, K, pi, chain = as.numeric(chain), iter, samples, pi_bias, accepts, accept_prob, n_leapfrogs)
      out <- rbind(out, row)
     
      # Saving results to disk
      if (save_runs){
        file_name <- paste0(folder_name, "/", model$model_name, "_", K, "_", dist$name, ".RData")
        print(str_c("Saving object: ", file_name))
        save(row, pi, model_init, file = file_name)
      }
    }
  }
  return(out)
}
# Helper function loading experiments saved in .rds files and 
# collecting sampling runs from each experiment into a data.frame
load_experiments <- function(folder_name){
    file_names <- list.files(folder_name, full.names = TRUE)
    lst <- list()
    i <- 1
    name_prev <- ""
    samples <- list()
    for (name in file_names){
      tmp_env <- new.env()
      model_name <- str_split(name, "/")[[1]][3]
      model_name <- str_extract(model_name, "^(.*?)(?=_\\d+\\.[^.]+$)")
      load(name, envir = tmp_env)
      
      # Convert all objects in the environment into a list
      data_list <- as.list(tmp_env)
      if(!is.null(data_list$row)){
        data <- data_list$row
      } else{
        data <- data_list$out %>% mutate(chain = as.numeric(chain))
      }
      lst <- append(lst,list(data))
    }
    experiments <- dplyr::bind_rows(lst)
  }



# Defining plotting function for producing convergence plots
font_add_google("Lato", "lato")
font_add_google("Lora", "lora")
showtext_auto()

get_conv_plots <- function(experiments){
  experiments %>% group_by(name, iter, pi, K)  %>% 
    mutate(bias_mean = mean(pi_bias))  %>% 
    ggplot(aes(x= iter, y = pi_bias, color = name, group = interaction(name,chain))) +
    geom_line(alpha = 0.05)  + geom_line(aes(x= iter, y = bias_mean, color = name)) +
    theme_minimal(base_family = "lora") +
    facet_grid(
      rows = vars(pi),
      cols = vars(K),
      labeller = labeller(K = function(x) paste0("K = ", x)),  scales = "free_y"
    ) +
    labs(x = "Iteration", y = "Total variation distance") +
    scale_color_manual(
      name = "Proposal",
      values = c(
        "Categorical (exact)" = "#000000",  # black
        "GS1"                 = "#E69F00",   # green
        "Dirichlet"           = "#CC79A7",  # blue
        "Logistic normal"     = "#19A842",# orange
        "Unif"                = "#0072B2"   # purple
      ),
      labels = c(
        "Cat (exact)" = "Categorical (exact)",
        "GS1" = expression("Gumbel-Softmax"),
        "Dirichlet (HMC-100)" = "Dirichlet",
        "Logitnormal" = "Logistic Normal",
        "Unif" = "Uniform"
      ),
    ) + 
    theme_minimal(base_size = 11) +
    theme(
      panel.grid.major = element_line(colour = "grey90", size = 0.3),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(colour = "grey80", fill = NA, size = 0.4),
      # axis.line = element_line(colour = "grey30", size = 0.4),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", colour = "grey10"),
      #legend.title = element_text(face = "bold"),
      legend.title = element_blank(),
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.text = element_text(size = 12),
      legend.margin = margin(t = 5),
      #legend.key = element_blank(),
      plot.margin = margin(5.5, 20, 5.5, 5.5),
      strip.text.x = element_text(face = "bold", size = 12),
      strip.text.y = element_text(face = "bold", size = 12),
      
      
      strip.placement = "outside"
      
    ) +
    scale_x_continuous(breaks = c(0, 500, 1000),
                       labels = c("0", "500", "1000")) + guides(
                         color = guide_legend(
                           override.aes = list(linewidth = 1.8, alpha = 1)
                         )
                       ) 
}
