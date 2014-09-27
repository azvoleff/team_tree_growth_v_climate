library(rstan)
library(foreach)
library(doParallel)

model_file <- "full_model.stan" 
n_chains <- 12
n_iter <- 1000

load("stan_data.RData")
load("stan_init.RData")

seed <- 1638
fit <- stan(model_file, data=stan_data, iter=n_iter, chains=1,
            init=rep(stan_init, 1), refresh=1, seed=seed)

# Fit initial model, on a single CPU. Run only one iteration.
stan_fit_initial <- stan(model_file, data=stan_data, iter=0, chains=1,
                         init=stan_init, chain_id=1)
print("finished running initial stan model")
save(stan_fit_initial, file="stan_fit_initial.RData")

# # Function to extract the last parameter estimates from a chain in a stanfit 
# # object to use as initialization values for a new chain.
# get_inits <- function(stan_fit, stan_init) {
#     pars <- names(stan_init[[1]])
#     init_param_ests <- extract(stan_fit_initial, permuted=TRUE, par=pars)
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
# stan_init_revised <- get_inits(stan_fit_initial, stan_init)
# save(stan_init_revised, file="stan_init_revised.RData")

# Fit n_chains chains in parallel. Reuse same seed so that the chain_ids can be 
# used by stan to properly seed each chain differently.
cl <- makeCluster(n_chains)
registerDoParallel(cl)
sflist <- foreach(n=1:n_chains, .packages=c("rstan")) %dopar% {
    # Add 1 to n in order to ensure chain_id 1 is not reused
    stan(fit=stan_fit_initial, data=stan_data, seed=seed, chains=1,
         iter=n_iter, chain_id=n+1, refresh=-1, init=stan_init)
    # stan(fit=stan_fit_initial, data=stan_data, seed=seed, chains=1,
    #      iter=n_iter, chain_id=n+1, refresh=-1, init=stan_init_revised)
}
print("finished running stan models on cluster")
stopCluster(cl)

stan_fit_allchains <- sflist2stanfit(sflist)
save(stan_fit_allchains, file="stan_fit_allchains.RData")

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
