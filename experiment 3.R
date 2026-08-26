library(spikeslab)
library(tidyverse)
library(scales)
library(gridExtra)
library(loo)
source("ex3_utils.R")
source("rdhmc.R")

# Loading diabetes data set
data("diabetesI", package = "spikeslab")

# Preparing data for modeling
data <- as.matrix(diabetesI)
N <- nrow(data)
K <- ncol(data)

means <- as.matrix(colMeans(diabetesI))
sds <- apply(data, 2, sd)
ones <- matrix(1, nrow = N, ncol = 1)
means_matrix <- ones %*% t(means)

data_std <- (data - means_matrix) %*% diag(1/sds)
data_std <- scale(data, center = TRUE, scale = TRUE)
X <- data_std[, 2:K]
y <- data_std[,1]

# Defining model hyperparameters
tau <- 0.0025
lambda <- 1
lambda_slab <- 1 # Slab variance
lambda_spike <- 1e-6 # Spike variance
upsilon <- 0.01
eta <- 0.1
beta0priorvar <-  lambda_slab
s <- 1
p <- 0.5
pi <- c(p, 1-p)
n_coefs <- ncol(X)

# Defining input data for the RDHMC model
model_data <- list(X =X,
                   y = y, 
                   N = N,
                   n_coefs = n_coefs,
                   lambda = lambda,
                   lambda_spike = lambda_spike,
                   lambda_slab = lambda_slab,
                   beta0priorvar = beta0priorvar,
                   eta = eta,
                   tau = tau,
                   alpha = c(2,2),
                   s = s,
                   pi = pi, 
                   upsilon = upsilon
)

# RDHMC sampling 
#t_start <- Sys.time()
#fit <- rdhmc_sample("Models/sns.stan", model_data, n_chains = 10, n_samples = 1000, warmup_sampler = 1000, warmup_init = 2000, max_treedepth = 13)
#t_end <- Sys.time()
fit <- readRDS("Experiments/Experiment 3/ex3_rdhmc_full_13.rds")

# Extracting RDHMC draws
draws <- fit$samples_df

# Reading samples externally produced by a Gibbs sampler implemented in Julia
beta_julia <- read_csv("Experiments/Experiment 3/beta_julia.csv", col_names = FALSE)
sigma_julia <- read_csv("Experiments/Experiment 3/sigma_julia.csv", col_names = FALSE)
z_julia <- read_csv("Experiments/Experiment 3/z_julia.csv", col_names = FALSE)

# Setting recorded run-time (s) for the Julia sampler
t2 <- 3522.3500730991364

# Computing RMSE for the first moment estimates between the RDHMC and Gibbs samples
beta_hmc <- draws %>% select(beta0, contains("beta_c")) %>% summarize_all(mean) %>% pivot_longer(cols = everything()) %>% mutate(i = as.numeric(str_extract(name, "\\d+"))) %>% group_by(i) %>% summarise(value = mean(value))
beta_gibbs <- beta_julia %>% summarize_all(mean) %>% pivot_longer(cols = everything()) %>% mutate(i = as.numeric(str_extract(name, "\\d+")) - 1) %>% group_by(i) %>% summarise(value = mean(value))
beta_diff <- inner_join(beta_hmc, beta_gibbs, by = c("i" = "i")) %>% transmute(component = i, diff = (value.x - value.y)^2) 

sigma_hmc <- draws%>% pull(sigma) %>% mean()
sigma_gibbs <- sigma_julia[[1]] %>% mean()
sigma_diff <- (sigma_hmc - sigma_gibbs)^2

rmse <- mean(c(beta_diff$diff, sigma_diff))

# Computing ESS of the continuous parameters from the RDHMC and Gibbs sampler
ess_obj_rdhmc <- get_ess_rdhmc(fit, draws)
ess_rdhmc <- ess_obj_rdhmc$ess
ess_rdhmc_norm <- ess_obj_rdhmc$ess_norm
julia_samples <- bind_cols(beta_julia, sigma_julia, z_julia)
julia_ess_metrics <- get_ess_julia(julia_samples, t2)

# Computing ELPD for both samplers using the "loo" package
set.seed(19890226)
log_lik_hmc <- draws %>% select(contains("log_lik")) %>% as_draws_array()
loo_hmc <- loo(log_lik_hmc)

n_samples <- nrow(beta_julia)
n <- length(y)
log_lik_gibbs <- matrix(0, nrow = n_samples, ncol = n)
for (i in 1:n_samples){
  mu_hat <- X %*% t(beta_julia[i, -1])
  log_lik_gibbs[i,] <- dnorm(y, mean = mu_hat + as.numeric( beta_julia[i, 1]), sd = as.numeric(sqrt(sigma_julia[[1]][i])), log = TRUE)
}
loo_gibbs <- loo(log_lik_gibbs)

# Computing posterior modes for the categorical samples from the RDHMC and Gibbs sampler
c_hmc <- draws %>% select(starts_with("c[")) %>% summarize_all(mean) %>% pivot_longer(cols = everything()) %>% mutate(i = as.numeric(str_extract(name, "\\d+"))) %>% group_by(i) %>% summarise(value = mean(value))
c_gibbs <-  z_julia %>% summarize_all(mean) %>% pivot_longer(cols = everything()) %>% mutate(i = as.numeric(str_extract(name, "\\d+"))) %>% group_by(i) %>% summarise(value = mean(value))
c_hmc["method"] <- "GS-HMC"
c_gibbs["method"] <- "SNS-Gibbs"

# Extracting and collecting the significant selection indicators
is_significant_gibbs <- c_gibbs %>% filter(value >=0.5) %>% select(i)
is_significant_hmc <- c_hmc %>% filter(value >=0.5) %>% select(i)
c_dist <- bind_rows(c_hmc, c_gibbs)

# Producing plot over posterior distribution of the 10 most significant selection indications
n <- 10
top_n_gibbs <- union(is_significant_gibbs$i, is_significant_hmc$i)
gibbs_order <- c_dist %>% filter(i %in% top_n_gibbs, method == "SNS-Gibbs") %>% arrange(desc(value)) %>% pull(i)

p3a <- c_dist %>% filter(i %in% top_n_gibbs) %>% group_by(method) %>% arrange(method, desc(value)) %>% ggplot(aes(x = factor(i, levels = gibbs_order), y = value, fill = factor(method, levels = c("SNS-Gibbs", "GS-HMC") ))) + 
  geom_bar(color = "black", linewidth = 0.2, stat = "identity", position = "dodge", linewidth = 2, width = 0.7, alpha = 0.8) + theme_bw() + 
  scale_fill_discrete( palette = scales::pal_brewer(type = "qual", palette =6) ) +  
  scale_y_continuous(labels = scales::percent_format()) + 
  labs(x = expression(I[i]), y = "Posterior probability", fill = "Method") + 
  theme(legend.position = "none")

# Adding small font theme
p3a2 <- p3a +theme_small +  
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1
    )
  )

# Extracting regression coefficient posterior samples corresponding to the 8 
# most significant selection indicators
m <- 8
best_beta <- gibbs_order[1:m]
best_beta_str <-  paste0("i = ", best_beta)
beta_i_hmc <- draws %>% select(contains("beta_c")) 
beta_0_hmc <-  draws %>% select(contains("beta0")) 
beta_i_hmc <- beta_i_hmc[, best_beta]

beta_0_gibbs <- beta_julia[, 1]
beta_i_gibbs <- beta_julia[, best_beta + 1]

names(beta_i_hmc) <- best_beta_str
names(beta_i_gibbs) <- best_beta_str

beta_i_hmc$method <- "RD-HMC"
beta_i_gibbs$method <- "SNS-Gibbs"
beta_i_samples <- bind_rows(beta_i_hmc, beta_i_gibbs) %>% pivot_longer(cols =1:m, names_to = "beta") %>% mutate(beta = factor(beta, levels = best_beta_str))


# Generating figure of posterior distributions 
p3b <- beta_i_samples %>% 
  ggplot(aes(x = value, fill = factor(method, levels = c("SNS-Gibbs", "RD-HMC") ), color = factor(method, levels = c("SNS-Gibbs", "RD-HMC") ))) +
  geom_density(alpha = 0.7, adjust = 1, color = NA) +
  theme_minimal() +
  scale_fill_discrete(palette = scales::pal_brewer(type = "qual", palette = 6))  +
  facet_wrap(vars(beta), nrow = 2, ncol = 4,  scales = "free_x") + labs(x = expression(beta[i]), y = "Posterior density", fill = "Method", color = "Method") + theme(legend.position = "bottom",  legend.justification = c(-0.2, 0), legend.margin = margin(0, 0, 0, 0)) +
  scale_x_continuous(
    breaks = pretty_breaks(n = 3),
    labels = label_number(accuracy = 0.1) 
  ) + theme_small


# Collecting to one figure
p3 <- grid.arrange(p3a2, p3b, nrow= 1, ncol=2, widths = c(1.3, 1.7))


# Creating metrics table 
metrics_hmc <- c("RD-HMC", paste0(is_significant_hmc$i, collapse = ","), round(loo_hmc$estimates[1,1], 3), round(loo_hmc$estimates[1,2], 3), round((ess_rdhmc$sigma_ess_bulk + ess_rdhmc$beta0_ess_bulk)/2), round((ess_rdhmc_norm$sigma_ess_bulk + ess_rdhmc_norm$beta0_ess_bulk)/2,3), round(rmse,6) )
metrics_gibbs <- c("SNS-Gibbs", paste0(is_significant_gibbs$i, collapse = ","), round(loo_gibbs$estimates[1,1], 3), round(loo_gibbs$estimates[1,2], 3), round((julia_ess_metrics$ess_gibbs_sigma_bulk + julia_ess_metrics$ess_gibbs_beta0_bulk)/2), round((julia_ess_metrics$ess_gibbs_sigma_bulk_s + julia_ess_metrics$ess_gibbs_beta0_bulk_s)/2,3), "-" )

df_metrics <- as.data.frame(rbind(metrics_gibbs,metrics_hmc))
names(df_metrics) <- c("# significant param Gibbs", "# significant param HMc",  "PSIS-LOO Abs diff", "PSIS-SE", "ESS", "ESS/s", "Beta RMSE")
