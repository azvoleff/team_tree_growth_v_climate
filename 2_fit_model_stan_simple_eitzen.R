library(rstan)
library(foreach)
library(doParallel)

model_file <- "simple_model_eitzen.stan" 
n_chains <- 8
n_iter <- 1000

load("model_data.RData")
load("init_data.RData")

init_data[[1]] <- init_data[[1]][!(names(init_data[[1]]) %in% c("b_ijk", "b_jk", "b_k"))]
model_data <- model_data[!(names(model_data) %in% c("n_plot", "n_site", 
                                                    "plot_ID", "site_ID", "WD", 
                                                    "spi"))]
model_data$max_obs_per_tree <- ncol(model_data$dbh)
model_data$n_miss <- nrow(model_data$miss_indices)
model_data$n_obs <- nrow(model_data$obs_indices)
# Stan can't hold matrices of ints, but can do lists of ints
model_data$obs_indices_tree <- model_data$obs_indices[, 1]
model_data$obs_indices_period <- model_data$obs_indices[, 2]
model_data$miss_indices_tree <- model_data$miss_indices[, 1]
model_data$miss_indices_period <- model_data$miss_indices[, 2]
model_data <- model_data[names(model_data) != "obs_indices"]
model_data <- model_data[names(model_data) != "miss_indices"]
# Use linear indexing to select observeds from dbh:
obs_linear_ind <- (model_data$obs_indices_period - 1) * nrow(model_data$dbh) + model_data$obs_indices_tree # From http://bit.ly/1rnKrC3
model_data$dbh_obs <- model_data$dbh[obs_linear_ind]
model_data <- model_data[names(model_data) != "dbh"]
# Now use linear indexing to select latent_dbhs as initialization values for 
# missing dbhs:
miss_linear_ind <- (model_data$miss_indices_period - 1) * nrow(init_data[[1]]$dbh_latent) + model_data$miss_indices_tree # From http://bit.ly/1rnKrC3
init_data[[1]]$dbh_miss <- init_data[[1]]$dbh_latent[miss_linear_ind]

# Stan doesn't allow NAs in input, so store 10 in the NAs of dbh_latent (these 
# cells will never be accessed anyways)
init_data[[1]]$dbh_latent[is.na(init_data[[1]]$dbh_latent)] <- 10

seed <- 1638
stan_fit_simple_eitzen_test <- stan(model_file, data=model_data, iter=10,
                                    chains=1, init=rep(init_data, 1),
                                    refresh=1, seed=seed)
print("finished running test stan model")
save(stan_fit_simple_eitzen_test, file="stan_fit_simple_eitzen_test.RData")

# Fit initial model, on a single CPU. Run only one iteration.
stan_fit_simple_eitzen_initial <- stan(model_file, data=model_data, iter=1,
                                       chains=1, init=init_data, chain_id=1)
print("finished running initial stan model")
save(stan_fit_simple_eitzen_initial, file="stan_fit_simple_eitzen_initial.RData")

# # Function to extract the last parameter estimates from a chain in a stanfit 
# # object to use as initialization values for a new chain.
# get_inits <- function(stan_fit_simple_eitzen, init_data) {
#     pars <- names(init_data[[1]])
#     init_param_ests <- extract(stan_fit_simple_eitzen_initial, permuted=TRUE, par=pars)
#     # Convert arrays into numeric
#     last_iter_num <- nrow(init_param_ests[[1]])
#     for (n in 1:length(init_param_ests)) {
#         # Select the last estimate if the chain has more than one iteration:
#         if (length(dim(init_param_ests[[n]])) > 1) {
#             init_param_ests[[n]] <- as.numeric(init_param_ests[[n]][last_iter_num, ])
#         } else {
#             init_param_ests[[n]] <- as.numeric(init_param_ests[[n]])
#         }
#     }
#     # Stan wants inits to be a list of lists:
#     init_param_ests <- list(init_param_ests)
#     return(init_param_ests)
# }
#
# # Extract the final values from this chain and use them to initialize the next
# # chains.
# init_data_revised <- get_inits(stan_fit_simple_eitzen_initial, init_data)
# save(init_data_revised, file="init_data_revised.RData")

# Fit n_chains chains in parallel. Reuse same seed so that the chain_ids can be 
# used by stan to properly seed each chain differently.
cl <- makeCluster(n_chains)
registerDoParallel(cl)
sflist <- foreach(n=1:n_chains, .packages=c("rstan")) %dopar% {
    # Add 1 to n in order to ensure chain_id 1 is not reused
    stan(fit=stan_fit_simple_eitzen_initial, data=model_data, seed=seed, chains=1,
         iter=n_iter, chain_id=n+1, refresh=-1, init=init_data)
    # stan(fit=stan_fit_simple_eitzen_initial, data=model_data, seed=seed, chains=1,
    #      iter=n_iter, chain_id=n+1, refresh=-1, init=init_data_revised)
}
print("finished running stan models on cluster")
stopCluster(cl)

stan_fit_simple_eitzen_allchains <- sflist2stanfit(sflist)
save(stan_fit_simple_eitzen_allchains, file="stan_fit_simple_eitzen_allchains.RData")

# model_params <- c("log_dbh_latent",
#                   "inter",
#                   "slp_dbh",
#                   "sigma_obs",
#                   "sigma_proc",
#                   "sigma_ijk",
#                   "sigma_jk",
#                   "sigma_k",
#                   "b_ijk",
#                   "b_jk",
#                   "b_k")
