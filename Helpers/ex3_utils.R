library(posterior)
library(tidyverse)

# Helper function taking in an fitted RDHMC model model_fit and samples draws.
# Returns data.frames with ESS for the regression coefficient and variance samples.  
# ESS is computed using both bulk and tail variants, and normalized by sampler runtime
get_ess_rdhmc <- function(model_fit, draws){
  
  t <- mean(fit$init$init$fit$time()$chains$total)
  samples <- draws
  draws <- list(samples = samples)
  accept_rate <- mean(sapply(model_fit$experiments, function(x) mean(x$accepts)))
    
  ess_bulk_full <- draws$samples %>%  select(contains("sigma"), contains("beta_c") , beta0) %>% summarise_all(posterior::ess_bulk) 
  ess_tail_full <- draws$samples %>%  select(contains("sigma"), contains("beta_c"), beta0 ) %>% summarise_all(posterior::ess_tail) 
  
  ess_bulk_mean <- ess_bulk_full%>% as.numeric() %>% mean(na.rm = TRUE)
  ess_tail_mean <- ess_tail_full %>% as.numeric() %>% mean(na.rm = TRUE)
  
  ess_bulk_mean_s <- ess_bulk_mean/t
  ess_tail_mean_s <- ess_tail_mean/t
  
  beta_ess_bulk <- ess_bulk_full %>% select(contains("beta_c"))  %>% as.numeric() %>% mean(na.rm = TRUE)
  beta_ess_tail <-ess_tail_full %>% select(contains("beta_c"))  %>% as.numeric() %>% mean(na.rm = TRUE)
  
  beta_ess_bulk_s <- beta_ess_bulk/t
  beta_ess_tail_s <-  beta_ess_tail/t
  
  beta0_ess_bulk <- ess_bulk_full %>% select(contains("beta0"))  %>% as.numeric() %>% mean(na.rm = TRUE)
  beta0_ess_tail <-ess_tail_full %>% select(contains("beta0"))  %>% as.numeric() %>% mean(na.rm = TRUE)
  
  beta0_ess_bulk_s <- beta0_ess_bulk/t
  beta0_ess_tail_s <-  beta0_ess_tail/t
  
  
  sigma_ess_bulk <- ess_bulk_full %>% select(contains("sigma")) %>% as.numeric() %>% mean(na.rm = TRUE)
  sigma_ess_tail <- ess_tail_full  %>% select(contains("sigma")) %>% as.numeric() %>% mean(na.rm = TRUE)
  
  sigma_ess_bulk_s <- sigma_ess_bulk/t
  sigma_ess_tail_s <- sigma_ess_tail/t
  
  c_ess_bulk <-ess_bulk_full %>%   select(starts_with("c["))  %>% as.numeric() %>% mean(na.rm = TRUE)
  c_ess_tail <- ess_tail_full%>%  select(starts_with("c[")) %>% as.numeric() %>% mean(na.rm = TRUE)
  c_ess_bulk_s <- c_ess_bulk/t
  c_ess_tail_s <- c_ess_tail/t
  
  ess_mean <-  draws$samples %>%  select(starts_with("c["), contains("sigma") )
  ess <- data.frame(tau, t, accept_rate,ess_bulk_mean, ess_tail_mean, beta0_ess_bulk, beta0_ess_tail, beta_ess_bulk, beta_ess_tail, sigma_ess_bulk,  sigma_ess_tail, c_ess_bulk, c_ess_tail )
  ess_norm <- data.frame(tau, t, accept_rate, ess_bulk_mean_s, ess_tail_mean_s, beta0_ess_bulk_s, beta0_ess_tail_s,  beta_ess_bulk_s, beta_ess_tail_s,  sigma_ess_bulk_s, sigma_ess_tail_s, c_ess_bulk_s, c_ess_tail_s)
  
  names(ess) <- c("tau", "t", "accept_rate", "ess_bulk_mean", "ess_tail_mean","beta0_ess_bulk", "beta0_ess_tail", "beta_ess_bulk", "beta_ess_tail", "sigma_ess_bulk",  "sigma_ess_tail", "c_ess_bulk", "c_ess_tail"  )
  names(ess_norm) <- c("tau", "t", "accept_rate", "ess_bulk_mean_s", "ess_tail_mean_s", "beta0_ess_bulk_s", "beta0_ess_tail_s", "beta_ess_bulk_s", "beta_ess_tail_s",  "sigma_ess_bulk_s", "sigma_ess_tail_s", "c_ess_bulk_s", "c_ess_tail_s"  )
  
  return(list(ess=ess, ess_norm=ess_norm, ess_bulk_full=ess_bulk_full))
}

# Helper function taking in samples from an external Julia Gibbs sampler julia_samples 
# and corresponding recorded sampler run time (s) sapmler_runtime.
# Returns a list with ESS and ESS/sampler_runtime (bulk and tail) for the regression coefficients betea
# and variance parameter sigma.
get_ess_julia <- function(julia_samples, sampler_runtime){
  
  ess_gibbs_mean <- julia_samples  %>% summarise_all(posterior::ess_bulk) %>% as.numeric() %>% mean(na.rm =TRUE)
  ess_gibbs_mean_s <- ess_gibbs_mean/sampler_runtime
 
  ess_gibbs_beta_bulk <- beta_julia[,1] %>% summarise_all(posterior::ess_bulk) %>% as.numeric() %>% mean()
  ess_gibbs_beta_tail <- beta_julia[,1] %>% summarise_all(posterior::ess_bulk) %>% as.numeric() %>% mean()

  ess_gibbs_sigma_bulk <- sigma_julia %>% posterior::ess_bulk()
  ess_gibbs_sigma_tail <- sigma_julia %>% posterior::ess_tail()
  
  ess_gibbs_sigma_bulk_s <- ess_gibbs_sigma_bulk/sampler_runtime
  ess_gibbs_sigma_tail_s <- ess_gibbs_sigma_tail/sampler_runtime
  
  ess_gibbs_beta0_bulk <- beta_julia[,1] %>% posterior::ess_bulk()
  ess_gibbs_beta0_tail <- beta_julia[,1] %>% posterior::ess_tail()
  
  ess_gibbs_beta0_bulk_s <- ess_gibbs_beta0_bulk/sampler_runtime
  ess_gibbs_beta0_tail_s <- ess_gibbs_beta0_tail/sampler_runtime
  
  return(list(ess_gibbs_mean = ess_gibbs_mean,
              ess_gibbs_mean_s = ess_gibbs_mean_s,
              ess_gibbs_beta_bulk = ess_gibbs_beta_bulk,
              ess_gibbs_beta_tail = ess_gibbs_beta_tail,
              ess_gibbs_sigma_bulk = ess_gibbs_sigma_bulk,
              ess_gibbs_sigma_tail = ess_gibbs_sigma_tail,
              ess_gibbs_sigma_bulk_s = ess_gibbs_sigma_bulk_s,
              ess_gibbs_sigma_tail_s = ess_gibbs_sigma_tail_s,
              ess_gibbs_beta0_bulk = ess_gibbs_beta0_bulk,
              ess_gibbs_beta0_tail = ess_gibbs_beta0_tail,
              ess_gibbs_beta0_bulk_s = ess_gibbs_beta0_bulk_s,
              ess_gibbs_beta0_tail_s = ess_gibbs_beta0_tail_s)
         )
}

# Defining ggplot template for small text fonts
theme_small <- theme(
  axis.text.x  = element_text(size = 5),
  axis.text.y  = element_text(size = 5),
  axis.title.x = element_text(size = 6),
  axis.title.y = element_text(size = 6),
  
  strip.text   = element_text(size = 6),
  
  legend.title = element_text(size = 6),
  legend.text  = element_text(size = 5),
  
  legend.key.height = unit(0.35, "cm"),
  legend.key.width  = unit(0.35, "cm")
)
