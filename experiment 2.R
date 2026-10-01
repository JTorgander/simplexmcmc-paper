source("rdhmc.R")
source("Helpers/ex2_utils.R")
library(tidyverse)
library(MASS, exclude = 'select')
library(LaplacesDemon)
library(nimble)

# Setting number of mixtures and mixture component dimensions 
set.seed(19890226)
K <- 4
dim <- 2

# Setting number of data points that should be generated.
n_obs <- 100

# Sampling mean vectors from a multivariate normal with diagonal covariance
s1 <- 200
mu <- lapply(1:K, function(x){ mvrnorm(1, mu = rep(0, dim), Sigma = s1* diag(1, dim,dim))})
mu <- mu[ order(sapply(mu, `[[`, 1)) ]

# Sampling mixture components covariances
eta <- 4
sigma <- lapply(1:K, function(x) diag(rhalfcauchy(dim, scale=eta)))
Omega  <- lapply(1:K, function(x) simulate_corr_matrix(1, dim) )
Sigma <- mapply(function(x,y) x %*% y %*% x, sigma, Omega, SIMPLIFY = FALSE)

# Sampling mixture probabilities
pi = rdirichlet(1, alpha = rep(1, K))

# Generating normal mixture data set and preparing for the Stan model
mixture_data  <- get_mixture_data(n_obs, dim, pi, mu, Sigma, missing_data_prop = 0.1, standardize= TRUE)
model_data <- prepare_model_data(mixture_data)

# Producing plot of original data set
df <- data.frame(cbind(model_data$y, model_data$classes))
names(df) <- c("x", "y", "class_actual")
p_raw <- df %>%  ggplot(aes(x = x, y = y)) + geom_point(colour = "grey30", size = 0.5) + theme_minimal(base_size = 8) + ylim(c(-2,2)) +  labs(x = NULL, y = NULL, title = "(a) Original data") + theme(axis.text.x = element_blank(), axis.text.y = element_text(size = 5))

# Setting model hyperparameters
model_data$tau <- 0.01
model_data$lambda <- lapply(1:K, function(x) c(0,0))
model_data$upsilon <- 4
model_data$gamma <- 1 #Prior variance 
n_mixture_components <- length(pi)
model_data$n_mixture_components <- n_mixture_components
model_data$alpha <- rep(1, K)
model_data$eta <- 1
model_data$sigma <- rep(0.01, n_mixture_components)

# Fitting RDHMC model on mixture data (uncomment second line to load the model used in the paper)
#fit <- rdhmc_sample("Models/gs_gaussian_mixture.stan", model_data, n_chains = 10, n_samples = 1000, warmup_sampler = 1000, warmup_init= 1000)
fit <- readRDS("Experiments/Experiment 2/ex2a_rdhmc.rds") # Reading previously fitted model

# Extracting experiment samples into one data frame
experiments_list <- lapply(fit$experiments, function(x)x$samples)
samples_df <- bind_rows(experiments_list, .id = ".chain")

# Extracting mean vector posterior distributions and creating plot
mu_samples_obj <- get_mu_samples(samples_df)
mu_samples <-  mu_samples_obj$mu_samples
mu_order <- mu_samples_obj$mu_order
p_mu <- ggplot(mu_samples, aes(x = x1, y = x2)) + geom_point(alpha = 0.1, aes(color = class), size = 0.1) + 
  scale_alpha(range = c(0.1, 0.6)) +
  theme_minimal() + 
  scale_color_discrete( palette = scales::pal_brewer(type = "qual", palette =2) ) + labs(x = NULL, y = NULL, color = "Class", title = "(d) Group mean posteriors") + lims(x = c(-2,1.5), y = c(-2,2)) + theme(axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 5), legend.position = "none")

# Extracting categorical samples and computing corresponding posterior modes and proportions
classes <- get_class_preds(samples_df, mu_order)
class_preds <- classes$class_preds
class_samples <- classes$class_samples
class_props <- classes$class_props

# Creating posterior plot of categorical variables corresponding to misclassified points
df_w_preds <- df %>% mutate(obs_no = row_number()) %>% inner_join(class_preds) %>% mutate(is_correct = class_actual == class_pred )
incorrect_obs <- df_w_preds %>% filter(!is_correct) %>% select(obs_no)
incorrect_class_dists <- class_props %>% filter(obs_no %in% incorrect_obs$obs_no) 
p_dist <- incorrect_class_dists %>% ggplot(aes(x = factor(obs_no), y = n/N, fill =factor(class_pred ))) + geom_col(position = "dodge", linewidth = 1) + scale_fill_discrete( palette = scales::pal_brewer(type = "qual", palette =2) ) + theme_minimal(base_size = 5) +  scale_y_continuous(labels = scales::percent_format())  + labs(title = "(c) Class distribution", y = NULL, x = "Obs #", fill = "Class") + theme(axis.title.x = element_text(size = 6))

# Creating plot of categorical variable posterior modes
p_preds <- df_w_preds %>% ggplot(aes(x = x, y = y)) + 
  geom_point(aes(color = factor(class_pred)), size = 0.5) + 
  geom_point(
    data = ~ dplyr::filter(.x, !is_correct),
    color = "black",     # choose any color you want
    size = 2,         # slightly larger if desired,
    shape = 3
  ) + geom_text(data = ~ dplyr::filter(.x, !is_correct), aes(label=obs_no), nudge_x = -0.3, nudge_y = 0, size = 2) + 
  theme_minimal(base_size = 11) + ylim(c(-2,2)) + 
  scale_color_discrete( palette = scales::pal_brewer(type = "qual", palette =2) ) + #+ theme(legend.position = "none") +
  labs(x = NULL, y = NULL, color = "Predicted class", title = "(b) Class predictions") + theme(legend.position = "none", axis.text.y = element_blank(),  axis.text.x = element_blank()) + lims(x = c(-2,1.5), y = c(-2,2)) 

# Extracting and creating plot of posterior predictive distribution
y_preds <- samples_df  %>%  select(contains("y_pred")) %>% mutate(iter = row_number()) %>% pivot_longer(cols = contains("y_pred")) %>% transmute(iter, obs_no = as.numeric(str_match(name, "\\d+")), dim = str_match(name, "\\d(?=\\]$)"), value)
set.seed(19890226) # Thinning sample to simplify computation
thin_idx <- sample(1:nrow(y_preds), 200000)
y_preds <- y_preds[thin_idx, ]
p_ypred <- y_preds %>% left_join(class_preds) %>%  pivot_wider(names_from = dim, values_from = value) %>% unnest() %>% ggplot(aes(x=`1`, y = `2`, color = factor(class_pred))) + geom_point(size = 0.5, alpha = 0.08) + 
  scale_color_discrete( palette = scales::pal_brewer(type = "qual", palette =2) ) + theme_minimal() + theme(legend.position = "none",  axis.text.x = element_text(size = 6), axis.text.y = element_blank()) + labs(x = NULL, y = NULL, title = "(e) Predictive distribution") + lims(x = c(-2,1.5), y = c(-2,2)) 

#### Experiment (b) Missing value imputation #####

# Rearranging data to include known class labels
model_data2 <- model_data
model_data2$y_missing <- mixture_data$y_missing
model_data2$y_non_missing <- mixture_data$y_non_missing
model_data2$n_missing <- mixture_data$n_missing
model_data2$n_non_missing <-  mixture_data$n_non_missing
model_data2$classes_missing <- mixture_data$classes_missing
model_data2$classes_non_missing <- mixture_data$classes_non_missing

# Fitting missing value imputation model
#fit2 <- rdhmc_sample("Models/gs_gaussian_mixture.stan", model_data2, n_chains = 10, n_samples = 1000, warmup_sampler = 1000, warmup_init= 1000, L =NA)
fit2 <- readRDS("Experiments/Experiment 2/ex2b_rdhmc.rds")

# Extracting samples and producing imputed missing values plot
draws2 <- fit2$samples_df
mu_imp <- get_mu_samples(draws2)
mu_order_imp <- mu_imp$mu_order
classes_imp <- get_class_preds(draws2, mu_order_imp)
class_preds_imp <- classes_imp$class_preds
class_sampels_imp <- classes_imp$class_samples
class_props_imp <- classes_imp$class_props
class_preds_imp$obs_no <- model_data2$missing_classes
df_w_preds2 <- df %>% mutate(obs_no = row_number()) %>% inner_join(class_preds_imp)
p_imp <- df %>%  ggplot(aes(x = x, y = y)) + geom_point(aes(colour = factor(class_actual)), size = 1, alpha = 0.2) + geom_point(data = df_w_preds2, aes(color = factor(class_actual)), size = 0.6)+  theme_minimal(base_size = 8) + ylim(c(-2,2)) + scale_color_discrete( palette = scales::pal_brewer(type = "qual", palette =2) ) + theme(legend.position = "none", axis.text.y  = element_blank(), axis.text.x  = element_text(size = 6)) + labs(x = NULL, y = NULL, title = "(f) Class imputation") + lims(x = c(-2,1.5), y = c(-2,2)) 


#### Experiment (c) temperature calibration for isotropic normal mixture model 

# Generating/loading data
#set.seed(19890226)
#data_25_5 <- get_iso_grid_data(25, 5, NULL, TRUE)
#data_100_9 <- get_iso_grid_data(100, 9, NULL, TRUE)
data_25_5 <- readRDS("Experiments/Experiment 2/ex2_25_5.rds")
data_100_9 <- readRDS("Experiments/Experiment 2/ex2_100_9.rds")

# Generating plots of original data
p1 <- ggplot(as.data.frame(data_100_9$y), aes(x = `V1`, y = `V2`)) + geom_point() + theme_minimal() + labs(x = "x", y = "y")
p2 <- ggplot(as.data.frame(data_25_5$y), aes(x = `V1`, y = `V2`)) + geom_point() + theme_minimal()  + labs(x = "x", y = "y")

# Setting up and running experiments
configs <- list(list(N = 25, K = 5, data = data_25_5), 
                list(N = 100, K = 9, data = data_100_9))

for (config in configs){
  for (tau in c(1, 0.5, 0.1, 0.05, 0.01, 0.005, 0.001)){
    run_iso_experiment("Models/gs_gaussian_mixture_iso.stan", 
                       config$data, 
                       tau = tau, 
                       N = config$N, 
                       K = config$K, 
                       chains = 1, 
                       iter_warmup = 2000, 
                       iter_sampling = 1000, 
                       method = "RDHMC", 
                       max_treedepth = 10)
  }
}

# Reading experiments from disk and generating table
set.seed(0226)
folder_name <- "Experiments 2/RDHMC/Paper"
df_experiments <- load_experiments(folder_name)