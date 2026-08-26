library(tidyverse)
source("ex1_utils.R")
source("models.R")


# Loading predetermined target categorical distributions from disk
# Creating the environment variables dist1, dist1b, dist2, dist2b, dist3, dist3b, dist4, dist4b
load("Experiments/Experiment 1/dists.RData")
SEED <- 123
# Setting up experiments 
models_rdhmc <- function(pi, K, alpha0){
  
  list(
       gs1 = list(method = "RDHMC",
                  model_name = "GS1_RDHMC",
                  proposal = rdhmc_proposal, 
                  proposal_args = list(method = "RDHMC", 
                                       model_file = "./Models/gs.stan",
                                       model_seed = SEED,
                                       data = list(K = K, 
                                                   pi = as.numeric(pi), 
                                                   tau = 1)), 
                  lambda = pi),
     ln_rdhmc = list(method = "RDHMC",
                     model_name = "LN_RDHMC",
                     proposal = rdhmc_proposal, 
                     proposal_args = list(method = "RDHMC", 
                                          model_file = "./Models/ln.stan", 
                                          guide_dist = ln_prop, 
                                          guide_args = list(mu = pi), 
                                          data = list(mu = pi[1:(K-1)], 
                                                      K = K, 
                                                      sigma = 1)), 
                    lambda = get_lambda_ln(pi[1:(K-1)])),
    dir_unif = list(method = "RDHMC", 
                    model_name = "Unif_RDHMC",
                    proposal = rdhmc_proposal, 
                    proposal_args = list(method = "RDHMC", 
                                         model_file = "./Models/dirichlet.stan", 
                                         guide_dist = dir_guide, 
                                         model_seed = SEED,
                                         guide_args = list(alpha = rep(1, K)), 
                                         data = list(pi = as.numeric(pi), K = K, alpha =  as.numeric(rep(1,K)))), 
                    lambda = rep(1/K, K)),
    dir_adapt = list(method = "RDHMC", 
                     model_name = "Dirichlet_RDHMC",
                     proposal = rdhmc_proposal, 
                     proposal_args = list(method = "RDHMC", 
                                          model_file = "./Models/dirichlet.stan", 
                                          model_seed = SEED,
                                          guide_dist = dir_guide, 
                                          guide_args = list(alpha = pi), 
                                          data = list(pi = as.numeric(pi), K = K, alpha =  as.numeric(pi))), lambda = get_lambda_dir(pi)))
}

# Setting up experiment for exact sampling.
models_exact <- function(pi, K, alpha0){list(cat = list(method = "Exact", model_name = "Cat (exact)", proposal = cat_prop, proposal_args = list(pi = pi, K = K, guide_dist = cat_dens), lambda = pi),
                                              gs1 = list(method = "Exact", model_name = "GS1 (exact)", proposal = gs_prop, proposal_args = list( pi = pi, K = K, tau =1), lambda = pi) ,
                                              unif = list(method = "Exact", model_name = "Unif (exact)", proposal = unif_prop, proposal_args = list(K = K, guide_dist = unif_dens), lambda = rep(1/K, K)) ,
                                              dir = list(method = "Exact", model_name = "Dirchlet (exact)", proposal = dir_prop, proposal_args = list( alpha =  pi, guide_dist = dir_guide), lambda = get_lambda_dir(pi )),
                                              ln = list(method = "Exact", model_name = "Logit normal (exact)", proposal = ln_prop, proposal_args = list( pi =  pi[1:(K-1)]), lambda = get_lambda_ln(pi[1:(K-1)] )))
                                        }


# Running RDHMC experiments for a Gumbel-Softmax target distribution
experiments_rdhmc <- get_simplex_experiments(model_getters = models_rdhmc,
                                                n_experiments = 10, 
                                                experiment_sample_size = 1000,
                                                cat_dists = list(
                                                               pi_bal1 = dist1,
                                                               pi1_unbal = dist1b,
                                                                pi2_bal = dist2,
                                                                pi2_unbal = dist2b,
                                                                pi4_bal = dist4,
                                                                pi4_unbal = dist4b
                                                                ),
                                                 warmup_init = 2000, 
                                                 warmup_sampling = 2000, 
                                                 save_runs = TRUE)

# Producing reference solutions based on exact sampling
experiments_exact <- get_simplex_experiments(model_getters = models_exact,
                                             n_experiments = 1, 
                                             experiment_sample_size = 1000, 
                                             cat_dists =  list( pi_bal1 = dist1, 
                                                                pi1_unbal = dist1b,
                                                                pi2_bal = dist2 ,
                                                                pi2_unbal = dist2b,
                                                                pi4_bal = dist4,
                                                                pi4_unbal = dist4b),
                                             warmup_init = 2000, 
                                             warmup_sampling = 0, 
                                             save_runs = TRUE)


# Loading previous experiments
experiments <- load_experiments("Experiments/Experiment 1") 

# Separating metric table into RDHMC and exact sampling 
experiments_rdhmc <- experiments %>% filter((str_detect(name, "RDHMC"))) %>% 
                                     mutate(name =  str_replace(name, "Dirchlet", "Dirichlet"))
experiments_exact <- experiments %>% filter((str_detect(name, "exact") | name == "Cat (exact)"), !name %in% c("GS10", "GS01"), K %in% c(10, 100, 1000)) %>% 
                                     mutate(name =  str_replace(name, "Dirchlet", "Dirichlet"), name = str_replace(name, "\\s\\(exact\\)", ""), name=str_replace(name, "Cat", "Categorical (exact)"), name = str_replace(name, "Logit", "Logistic")) 


# Generating convergence plots for SimplexMCMC samplers with exact proposal
p1 <- get_conv_plots(experiments_exact)

#### Generating table with comparison metrics for different RDHMC proposals ####

# Computing ESS (bulk and tail) and R-hat for each RDHMC sampler
sampler_summary <- experiments_rdhmc %>% 
                    select(name, iter, K, pi, chain, samples) %>%  
                    group_by(name, K, pi) %>% 
                    nest(samples = c(iter, chain, samples)) %>%
                    mutate(sample_matrix = map(samples, ~ {.x %>%
                                                           select(iter, chain, samples) %>%
                                                           pivot_wider(
                                                             names_from  = chain,
                                                             values_from = samples,
                                                             values_fill = NA
                                                             )  %>%
                                                           arrange(iter) %>%
                                                           select(-iter) %>%
                                                           as.matrix()}),
                          ess_bulk = map(sample_matrix, ~ ess_bulk(.x)),
                          ess_tail  = map(sample_matrix, ~ ess_tail(.x)),
                          rhat  = map(sample_matrix, ~ rhat(.x))) %>% 
                   unnest(ess_bulk, ess_tail, rhat) %>% 
                   select(-c(samples, sample_matrix))

# Summarizing acceptance rate, number of leapfrog steps and mini
experiment_summary <- experiments_rdhmc %>%  group_by(name, K, pi, chain) %>% 
                        summarize(n_leapfrogs = sum(n_leapfrogs), 
                                  accept_rate = mean(accepts), 
                                  pi_bias = min(pi_bias)
                                  )

# Combining experiment and sampler summaries, computes ESS over number of gradient evaluations 
# and extracts the relevant experiments, corresponding to K = 100. 
metrics <- inner_join(sampler_summary, experiment_summary) %>% 
           mutate(ess_bulk_per_grad = ess_bulk/(n_leapfrogs*2)) %>% 
           ungroup(chain)  %>% 
           summarise_all(mean) %>% 
           select(-chain) %>% 
           filter(!str_detect(name, "exact"), !name %in% c("GS01", "GS10")) %>% 
           arrange(K, pi, pi_bias) %>% 
           mutate(name = str_replace(name, "\\s\\(HMC-100\\)", "") ) 

# Cleaning up metrics to produce the final metrics table
metrics_table <- metrics %>% 
                  filter(K == 100, !str_detect(name, "GS0\\d+")) %>%  
                  transmute(name, 
                            K, 
                            pi, 
                            tv = round(pi_bias,3), 
                            accept_rate = round(accept_rate, 3),
                            rhat = round(rhat,3), 
                            ess_bulk_per_grad = round(ess_bulk_per_grad,3)
                            ) %>% 
                  arrange(name, pi)
