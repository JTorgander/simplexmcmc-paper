
# Helper function for generating normal mixture data. The function generates n_samples samples 
# from a normal mixture distribution with mean vectors mu, covariance matrices Sigma and mixture
# probabilities pi. The function randomly assigns missing_data_prop percent of the data as missing 
# and standardizes data and mean vectors if standardize = TRUE.
get_mixture_data <- function(n_samples, dim, pi, mu, Sigma, missing_data_prop = 0, standardize = TRUE){
  
  n_classes <- length(pi)
  data_full <- matrix(0, nrow = n_samples, ncol = dim)
  
  # Generating class labels, indicating which distribution each sample is drawn from
  classes_full <- sample(1:n_classes, n_samples, replace = TRUE, prob = pi)
  
  # Converting to one-hot matrix format
  classes_one_hot <- matrix(0, nrow=n_samples, ncol=n_classes)
  for (i in 1:n_samples){
    c <- classes_full[i]
    classes_one_hot[i,c] <- 1
    data_full[i, ] <- MASS::mvrnorm(1, mu[[c]], Sigma[[c]])
  }
  
  # Standardizing data and mean vectors (optional) 
  mu_std <- mu
  if (standardize){
    y <- data_full
    y_mean <- apply(y, 2, mean)
    y_sd   <- apply(y, 2, sd)
  
    y_std <- sweep(y, 2, y_mean, `-`)
    data_full <- sweep(y_std, 2, y_sd, `/`)
    mu_std <- lapply(mu, function(x) (x - y_mean)/y_sd)
  }
  
  # Sampling missing value indicators
  missing_classes <- sample(c(TRUE,FALSE), n_samples, replace=TRUE, prob = c(missing_data_prop, 1- missing_data_prop))
  n_missing <- sum(missing_classes)
  n_non_missing <- n_samples - n_missing
  missing_idx <- (1:n_samples)[missing_classes]
  
  # Separating data into missing and non-missing
  data_missing <- data_full[missing_classes, ]
  data_non_missing <- data_full[!missing_classes, ]
  classes_non_missing <- classes_full[!missing_classes]
  classes_missing <- classes_full[missing_classes]
  
  return(list(y=data_full,
              n_classes = n_classes,
              n_samples = n_samples,
              dim = dim,
              classes = classes_full,
              classes_one_hot = classes_one_hot,
              
              mu = mu,
              mu_std = mu_std,
              pi = pi,
              sigma = sigma,
              missing_classes = missing_idx,
              n_missing = n_missing,
              n_non_missing = n_non_missing,
              y_missing = data_missing,
              y_non_missing = data_non_missing,
              classes_non_missing = classes_non_missing,
              classes_missing = classes_missing)
  )
}

# Helper function preparing the simulated data to a format recognized by
# the stan model used for the experiment
prepare_model_data <- function(data){
  
  data_prepared <- data
  data_prepared$y_missing <- data_prepared$y 
  data_prepared$y_non_missing <- list()
  data_prepared$n_missing <- data$n_missing + data$n_non_missing
  data_prepared$n_non_missing <- 0
  data_prepared$classes_non_missing <- integer(0)
  data_prepared$classes_missing <- data$classes

  return(data_prepared)
}

# Helper function simulating a correlation matrix from an RLKJ-distribution with parameter eta and p. 
simulate_corr_matrix <- function(eta, p){
  OmegaL <- rlkj_corr_cholesky(n = 1, eta = eta, p = p)
  t(OmegaL) %*% OmegaL
}

# Helper function, creating a evenly spaced grid with K points.
sample_grid <- function(K){

  points <- matrix(0, K, 2)
  start_idx <- 1
  end_idx <- 1
  k <- 1
  stop <- FALSE
  while(!stop){
    
    row1 <- as.matrix(expand.grid(-k:k, -k))
    row2 <- row1
    row2[, 2] <- k
    row3 <- as.matrix(expand.grid(-k, -(k-1):(k-1)))
    row4 <- row3
    row4[,1] <- k
    
    rows <- rbind(row1, row2, row3, row4)
    start_idx <- end_idx + 1
    end_idx <- start_idx + nrow(rows) - 1
    
    if (end_idx >= K){
      n_sample_points <- K - start_idx + 1
      sample_row_idx <- sample(1:nrow(rows), n_sample_points)
      rows <- rows[sample_row_idx, ]
      points[start_idx:K, ] <- rows
      stop <- TRUE
      
    }else {
      points[start_idx:end_idx,  ]  <- rows
    }
    k <- k+1
  }
  return(points)
}

# Generating N data points from a normal mixture distribution with K isotropic components
# with means determined by mu_grid. If mu_grid is not provided the means are created on an 
# evenly spaced grid. 
get_iso_grid_data <- function(N, K, mu_grid = NULL, standardize = TRUE ){
 
  # Generating grid of means if not provided
  if (any(is.null(mu_grid))){
    mu_grid <- sample_grid(K)
  }

  # Ordering grid to prevent identification proglems
  mu_grid_ordered <- mu_grid[order(mu_grid[,1], mu_grid[,2]), ]
  mu <- lapply(seq_len(nrow(mu_grid_ordered)), function(i) as.numeric(mu_grid_ordered[i, ]))

  # Creating diagonal covariance matrix with variance eta
  eta <- 0.1
  sigma <- lapply(1:K, function(x) diag(eta, dim, dim))
  Omega <- lapply(1:K, function(x) diag(1, dim, dim))
  Sigma <- mapply(function(x,y) x %*% y %*% x, sigma, Omega, SIMPLIFY = FALSE)
  
  # Simulating data and preparing for modeling
  pi <- rep(1/K, K)
  mixture_data  <- get_mixture_data(N, dim, pi, mu, Sigma, missing_data_prop = 0.1, standardize= standardize)
  model_data <- prepare_model_data(mixture_data)
  return(model_data)
}

# Extracing mean posterior samples from samples_df and assures that 
# the means are indexed/labeled the same over the chains. 
get_mu_samples <- function(samples_df){
  mu_samples <- samples_df %>% select(contains("mu["), `.chain`) %>%  
                pivot_longer(cols = contains("mu")) %>% 
                transmute(chain = `.chain`, class = str_extract(name, "\\d+"), dim = str_match(name, "\\d(?=\\]$)"), value) %>% 
                group_by(chain, class, dim) %>%              
                mutate(draw = row_number()) %>%              
                ungroup() %>% 
                pivot_wider(names_from = dim, values_from = value) %>% 
                unnest(cols = c(`1`, `2`)) %>% 
                select(-draw)

  # Correcting for potential label switching by ordering the labeling such that a_1 < b_1 implies that a < b
  mu_order <- mu_samples %>% group_by(class, chain) %>% summarise(x1 = mean(`1`)) %>% ungroup(class) %>% group_by(chain) %>% arrange(chain, x1)  %>% mutate(rank = dense_rank(x1)) %>% select(-x1)
  mu_samples_reordered <- mu_samples %>% inner_join(mu_order) %>% transmute(chain, class = factor(rank), x1 = `1`, x2 = `2`)
  
  return(list(mu_samples = mu_samples_reordered, mu_order = mu_order))
  
}

# Helper function extracting the categorical samples from samples_df, generating the posterior mode and ordering the categories
# according to the order used for the means in mu_order. Returns posterior modes, categorical samples and empirical proportions of 
# each category
get_class_preds <- function(samples_df, mu_order){
  class_samples <- samples_df %>% select(contains("c["), `.chain`) %>% mutate(iter = row_number()) %>% pivot_longer(cols = contains("c["), values_to = "class_pred") %>% transmute(iter, chain = `.chain`, obs_no = as.numeric(str_match(name, "\\d+")), class_pred = as.factor(class_pred))
  class_samples <- class_samples %>% inner_join(mu_order, by = c("class_pred" ="class", "chain" = "chain")) %>% transmute(iter, obs_no, class_pred = rank)
  class_preds <- class_samples %>% group_by(obs_no) %>% count(class_pred)  %>% filter(n == max(n))
  class_props <- class_samples %>% group_by(obs_no) %>% count(class_pred) %>% mutate(N = sum(n))
  
  return(list(class_preds = class_preds, class_samples = class_samples, class_props = class_props))
}

# Helper function transforming the mean vectors from wide to long table format 
mu_pivot_longer <- function(df){
  df %>% 
  pivot_longer(cols = contains("mu")) %>% 
    transmute(chain = `.chain`, class = str_extract(name, "\\d+"), dim = str_match(name, "\\d(?=\\]$)"), value) %>% 
    group_by(chain, class, dim) %>%              # duplicates are multiple draws
    mutate(draw = row_number()) %>%              # index each draw within (chain, class, dim)
    ungroup() %>% 
    pivot_wider(names_from = dim, values_from = value) %>% 
    unnest(cols = c(`1`, `2`)) %>% 
    select(-draw)
  
}

# Helper function ordering the mean samples
order_mu <- function(df){
  df %>% group_by(class, chain) %>% 
    summarise(x1 = round(mean(`1`),0), x2 = round(mean(`2`),0)) %>% 
    ungroup(class) %>% 
    group_by(chain) %>% 
    arrange(chain, x1, x2)  %>% 
    mutate(rank = row_number()) %>% 
    select(-c(x1,x2))
}

# Generating RDHMC evaluation metrics given a fitted model fit and data model_data.
# Returns a data frame row containing acceptance rate of the RDHMC sampler, average R-hat 
# and ESS (bulk and tail) of the mean vector samples, and log scores
get_rdhmc_metrics <- function(model_fit){
  
  mu_samples <- get_mu_samples(model_fit$samples_df)
  mu_samples_reordered <- mu_samples$mu_samples
  acceptance_rate <- mean(sapply(model_fit$experiments, function(x) mean(x$accepts)))
  
  
  param_comb <- expand_grid(c("x1","x2"), unique(mu_samples_reordered$class))
  param_comb <- split(param_comb, seq_len(nrow(param_comb)))
  
  get_rhat_ess <- function(data, var_comb){
    draws_matrix <- data[c("chain", "class", pull(var_comb[1]))] %>% filter(class == pull(var_comb[2])) %>% pivot_wider(names_from = chain, values_from = 3, values_fn = list) %>% unnest(everything()) %>% select(-class) %>% as.matrix()
    r_hat <- posterior::rhat(draws_matrix)
    ess_bulk <- posterior::ess_bulk(draws_matrix)
    ess_tail <- posterior::ess_tail(draws_matrix)
    return(c(r_hat, ess_bulk, ess_tail))
    
  }
  L <- sapply(model_fit$init$chain_inits, function(x) x$L)
  
  r_hat_ess <- sapply(param_comb, function(x) get_rhat_ess(mu_samples_reordered, x))
  mean_lf <- mean(L)
  n_grads <- mean_lf * 2
  r_hat_ess_mean <- rowMeans(r_hat_ess, na.rm = TRUE)
  r_hat_mean <- r_hat_ess_mean[1]
  ess_bulk_mean <- r_hat_ess_mean[2]
  ess_bulk_grad <- ess_bulk_mean/n_grads
  
  ess_tail_mean <- r_hat_ess_mean[3]
  
  ess_tail_grad <- ess_tail_mean/n_grads
  log_lik <- model_fit$samples_df %>% select(contains("log_lik")) %>% as.matrix()
  log_score <- sum(log(colMeans(exp(log_lik))))
  row <- list(acceptance_rate = acceptance_rate, 
              r_hat_mean = r_hat_mean,
              lf = mean_lf,
              ess_bulk_mean = ess_bulk_mean, 
              ess_tail_mean = ess_tail_mean, 
              ess_bulk_grad = ess_bulk_grad,
              ess_tail_grad = ess_tail_grad,
              log_score = log_score)
  
  return(row)
}

# Wrapper function for running experiments for the isotropic normal mixture models with N samples  K mixture components and Gumbel-Softmax 
# temperature tau. The sampler uses n_chains chains and a warmup and sampling period of iter_warmup and iter_sampling
# respectively. The fitted model is saved to disk in a .rds file 
run_iso_experiment <- function(model,data, tau, N, K,n_chains, iter_warmup, iter_sampling, method, max_treedepth){
  
  # Setting model hyperparameters
  n_mixture_components <- length(data$pi)
  data$tau <- tau
  data$lambda <- lapply(1:K, function(x) c(0,0))
  data$upsilon <- 4
  data$sigma <- 0.1
  data$gamma <- 1 #Prior variance 
  data$n_mixture_components <- n_mixture_components
  data$alpha <- rep(1, K)
  data$eta <- 1
  data$pi <- rep(1/K, K)

  # Running sampler 
  fit <- rdhmc_sample(model, data, n_chains = n_chains, n_samples = iter_sampling, warmup_sampler = floor(iter_warmup/2), warmup_init= floor(iter_warmup/2), max_treedepth = max_treedepth, L =NA)
  metrics <- get_rdhmc_metrics(fit)
    
  # Saving result to disk
  folder_name <-  paste0("./Experiments 2/", method ,"/", Sys.Date())
  if (!dir.exists(folder_name)){
    dir.create(folder_name, recursive = TRUE)
  }
  
  filename <- paste0(folder_name, "/", N, "_", K, "_",str_replace(as.character(data$tau), "\\.", ""), ".RData")
  model_filename <- paste0(folder_name, "/", N, "_", K, "_",str_replace(as.character(data$tau), "\\.", ""), ".rds")
  saveRDS(fit, model_filename)
}

# Reading isotropic normal mixture experiments from folder_name and returns a data.frame with
# RDHMC acceptance rate, average R-hat, leap frog steps, average ESS (bulk and tail),
# ESS/number of gradient evaluatios and log score
load_experiments <- function(folder_name){
  
  file_names <- list.files(folder_name, full.names = TRUE)
  lst <- list()
  i <- 1
  name_prev <- ""
  samples <- list()
  model_data <- data_25_5
  for (name in file_names){
    print(name)
    tmp_env <- new.env()
    model_name <- str_split(name, "/")[[1]][4]
    
    model_name <- str_extract(model_name, "^(.*?)(?=_\\d+\\.[^.]+$)")
    model_tau <- str_extract(name, "\\d+(?=\\.rds$)")
    N <- str_extract(model_name, "^\\d+")
    print(model_tau)
    # load(name, envir = tmp_env)
    model <- readRDS(name)
    
    row <- get_rdhmc_metrics(model)
    row$tau <- model_tau
    row$N <- N
    lst <- append(lst,list(row))
  }
  df <- dplyr::bind_rows(lst) %>% arrange(desc(N), desc(tau)) %>% transmute(N, tau, acceptance_rate = round(acceptance_rate, 4), r_hat_mean= round(r_hat_mean, 4), lf = round(lf, 0),  ess_bulk_mean = round(ess_bulk_mean,0), ess_bulk_grad = round(ess_bulk_grad,1), log_score)
  return(df)
}

