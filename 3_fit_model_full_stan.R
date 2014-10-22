library(rstan)
library(foreach)
library(doParallel)

model_file <- "full_model.stan" 

load("model_data_wide.RData")
load("init_data_with_ranefs.RData")

# n_B is number of fixed effects
model_data$n_B <- 7
# n_B_g is number of genus-level random effects
model_data$n_B_g <- 5
# W is prior scale for the inverse-Wishart
model_data$W <- diag(model_data$n_B_g)
model_data$WD_sq <- model_data$WD^2
model_data$mcwd_sq <- model_data$mcwd^2

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
model_data <- model_data[names(model_data) != "spi"]
# Use linear indexing to select observeds from dbh:
obs_linear_ind <- (model_data$obs_indices_period - 1) * nrow(model_data$dbh) + model_data$obs_indices_tree # From http://bit.ly/1rnKrC3
model_data$dbh_obs <- model_data$dbh[obs_linear_ind]
model_data <- model_data[names(model_data) != "dbh"]

# Use linear indexing to select latent_dbhs as initialization values for 
# missing dbhs:
miss_linear_ind <- (model_data$miss_indices_period - 1) * nrow(init_data$dbh_latent) + model_data$miss_indices_tree # From http://bit.ly/1rnKrC3
init_data$dbh_miss <- init_data$dbh_latent[miss_linear_ind]

# Stan doesn't allow NAs in input, so store 0 in the NAs of dbh_latent (these 
# cells will never be accessed anyways)
init_data$dbh_latent[is.na(init_data$dbh_latent)] <- 0
model_data$mcwd[is.na(model_data$mcwd)] <- 0
model_data$mcwd_sq[is.na(model_data$mcwd_sq)] <- 0

#init_data$xi <- rep(1, model_data$n_B_g)
#init_data$mu_B_g_raw <- rep(0, model_data$n_B_g)
init_data$mu_B_g_raw <- rnorm(model_data$n_B_g, 0, 1)
init_data$xi <- rnorm(model_data$n_B_g, 0, 1)
init_data$B_g_raw_temp <- rnorm(model_data$n_B_g, 0, 1)
init_data$int_ijk_std <- init_data$int_ijk / init_data$sigma_int_ijk
init_data$int_jk_std <- init_data$int_jk / init_data$sigma_int_jk
init_data$int_k_std <- init_data$int_k / init_data$sigma_int_k
init_data$int_t_std <- init_data$int_t / init_data$sigma_int_t

init_data <- init_data[names(init_data) != "int_ijk"]
init_data <- init_data[names(init_data) != "int_jk"]
init_data <- init_data[names(init_data) != "int_k"]
init_data <- init_data[names(init_data) != "int_t"]
init_data <- init_data[names(init_data) != "B_g_raw"]
#init_data <- init_data[names(init_data) != "Tau_B_g_raw"]

get_inits <- function() {
    c(init_data, list(
        # Fixed effects
        B=c(rnorm(model_data$n_B, 0, 1)),
        # Sigmas
        sigma_obs=runif(1, .05, .2), 
        sigma_proc=abs(rnorm(1, 0, 1))
        # Tau_B_g_raw=matrix(rWishart(1, model_data$n_B_g+1, 
        #                             diag(model_data$n_B_g)), 
        #                    nrow=model_data$n_B_g)
    ))
}

seed <- 1638

stan_fit <- stan(model_file, data=model_data, iter=20, chains=2)

stan_fit <- stan(model_file, data=model_data, iter=20, chains=2, inits=get_inits)

stan_fit <- stan(model_file, data=model_data, iter=2000, chains=4, 
                 inits=get_inits)
print("finished running test stan model")
save(stan_fit, file="stan_fit_full.RData")

# # Fit n_chains chains in parallel. Reuse same seed so that the chain_ids can be 
# # used by stan to properly seed each chain differently.
# # Fit initial model, on a single CPU. Run only one iteration.
# stan_fit_initial <- stan(model_file, data=model_data, chains=0, init=get_inits)
# print("finished running initial stan model")
#
# id_string <- paste0(Sys.info()[4], format(Sys.time(), "_%Y%m%d-%H%M%S"))
# n_chains <- 4
# n_iter <- 2000
# #cl <- makeCluster(min(n_chains, detectCores()))
# cl <- makeCluster(2)
# registerDoParallel(cl)
# sflist <- foreach(n=1:n_chains, .packages=c("rstan")) %dopar% {
#     # Add 1 to n in order to ensure chain_id 1 is not reused
#     sink(paste0("stan_", id_string, "_chain", n, ".txt"), append=TRUE)
#     stanfit <- stan(fit=stan_fit_initial, data=model_data, seed=seed, chains=1,
#          iter=n_iter, chain_id=n, init=get_inits)
#     save(stanfit, file=paste0("stan_", id_string, "_chain", n, ".RData"))
#     sink()
#     stanfit
# }
# print("finished running stan models on cluster")
# stopCluster(cl)
#
# stan_fit_p <- sflist2stanfit(sflist)
# save(stan_fit_p, file="stan_fit_full_parallel.RData")
