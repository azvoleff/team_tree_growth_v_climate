library(rstan)
library(foreach)
library(doParallel)

model_file <- "full_model_blocked.stan" 

temp_var <- "tmp_meanannual"
precip_var <- "mcwd_run12"
model_type <- "full"
#model_type <- "testing"
in_folder <- 'Data'
out_folder <- 'MCMC_Chains'

suffix <- paste0('_', model_type, '-', temp_var, '-', precip_var)

load(file.path(in_folder, paste0("model_data_wide_blocked", suffix, ".RData")))
load(file.path(in_folder, paste0("init_data_with_ranefs", suffix, ".RData")))

# n_B is number of fixed effects
model_data_blocked$n_B <- 2
# n_B_g is number of genus-level random effects
model_data_blocked$n_B_g <- 7
# n_B_T is number of terms in the temperature model
model_data_blocked$n_B_T <- 3
# W is prior scale for the inverse-Wishart
model_data_blocked$W <- diag(model_data_blocked$n_B_g)
model_data_blocked$WD_sq <- model_data_blocked$WD^2
model_data_blocked$precip_sq <- model_data_blocked$precip^2

model_data_blocked$n_miss <- nrow(model_data_blocked$miss_indices)
model_data_blocked$n_obs <- nrow(model_data_blocked$obs_indices)
# Stan can't hold matrices of ints, but can do lists of ints
model_data_blocked$obs_indices_tree <- model_data_blocked$obs_indices[, 1]
model_data_blocked$obs_indices_period <- model_data_blocked$obs_indices[, 2]
model_data_blocked$miss_indices_tree <- model_data_blocked$miss_indices[, 1]
model_data_blocked$miss_indices_period <- model_data_blocked$miss_indices[, 2]
model_data_blocked <- model_data_blocked[names(model_data_blocked) != "obs_indices"]
model_data_blocked <- model_data_blocked[names(model_data_blocked) != "miss_indices"]
model_data_blocked <- model_data_blocked[names(model_data_blocked) != "spi"]
# Use linear indexing to select observeds from dbh:
obs_linear_ind <- (model_data_blocked$obs_indices_period - 1) * nrow(model_data_blocked$dbh) + model_data_blocked$obs_indices_tree # From http://bit.ly/1rnKrC3
model_data_blocked$dbh_obs <- model_data_blocked$dbh[obs_linear_ind]
model_data_blocked <- model_data_blocked[names(model_data_blocked) != "dbh"]

# Use linear indexing to select latent_dbhs as initialization values for 
# missing dbhs:
miss_linear_ind <- (model_data_blocked$miss_indices_period - 1) * nrow(init_data$dbh_latent) + model_data_blocked$miss_indices_tree # From http://bit.ly/1rnKrC3
init_data$dbh_miss <- init_data$dbh_latent[miss_linear_ind]

# Stan doesn't allow NAs in input, so store 0 in the NAs of dbh_latent (these 
# cells will never be accessed anyways)
init_data$dbh_latent[is.na(init_data$dbh_latent)] <- 0
model_data_blocked$precip[is.na(model_data_blocked$precip)] <- 0
model_data_blocked$precip_sq[is.na(model_data_blocked$precip_sq)] <- 0
model_data_blocked$temp[is.na(model_data_blocked$temp)] <- 0

init_data$int_ijk_std <- init_data$int_ijk / init_data$sigma_int_ijk
init_data$int_jk_std <- init_data$int_jk / init_data$sigma_int_jk
init_data$int_k_std <- init_data$int_k / init_data$sigma_int_k
init_data$int_t_std <- init_data$int_t / init_data$sigma_int_t

# Convert vector of scalings (variances)
sigma_B_g_sigma <- sqrt(diag(init_data$sigma_B_g))
#sigma_B_g_sigma 

# Compute correlation matrix:
rho_B_g <- init_data$sigma_B_g / (sigma_B_g_sigma %*% t(sigma_B_g_sigma))

# Compute cholesky factor of correlation matrix
L_rho_B_g <- chol(rho_B_g)

##
## Below is just for checking my math
##
# Compute cholesky factor of final covariance matrix (equivalent to 
# diag_pre_multiply line in Stan code):
#L_B_g_sigma <- diag(sigma_B_g_sigma) %*% t(L_rho_B_g)
#
# Verify original covariance matrix is equal to L_B_g_sigma %*% t(L_B_g_sigma) 
# (with allowances for rounding error)
#L_B_g_sigma %*% t(L_B_g_sigma)
#table(abs(init_data$sigma_B_g - L_B_g_sigma %*% t(L_B_g_sigma)) < 1E-15)

init_data$sigma_B_g_sigma <- sigma_B_g_sigma
init_data$L_rho_B_g <- L_rho_B_g
# Means of the genus-level random effects
init_data$gamma_B_g <- apply(init_data$B_g_raw, 2, mean)
init_data$B_g_std <- solve(diag(sigma_B_g_sigma) %*% L_rho_B_g) %*% t(init_data$B_g_raw - init_data$gamma_B_g)

# Below is to verify this line in Stan code:
# B_g <- transpose(rep_matrix(gamma_B_g, n_B_g) + diag_pre_multiply(sigma_B_g_sigma, L_rho_B_g * B_g_std));
#
# First redefine as:
# B_g <- transpose(a + b);

# Note that in below code the matrix has to be transposed to match the behavior 
# of rep_mat in Stan
# a <- matrix(rep(init_data$gamma_B_g, model_data_blocked$n_genus), 
#             ncol=model_data_blocked$n_genus)
# b <- (diag(init_data$sigma_B_g_sigma) %*% L_rho_B_g) %*% init_data$B_g_std
# B_g_check <- t(a + b)
# head(B_g_check - init_data$B_g_raw)
# table(abs(B_g - init_data$B_g_raw) < 1E-12)

init_data <- init_data[names(init_data) != "int_ijk"]
init_data <- init_data[names(init_data) != "int_jk"]
init_data <- init_data[names(init_data) != "int_k"]
init_data <- init_data[names(init_data) != "int_t"]
init_data <- init_data[names(init_data) != "B_g_raw"]

get_inits <- function() {
    c(init_data, list(
        # Fixed effects
        B=c(rnorm(model_data_blocked$n_B, 0, 1)),
        # Temperature model
        B_T=c(rnorm(model_data_blocked$n_B_T, 0, 1)),
        # Sigmas
        sigma_obs=runif(1, .00026, .001), 
        sigma_proc=abs(rnorm(1, 0, 1))
    ))
}

# stan_fit <- stan(model_file, data=model_data_blocked, iter=20, chains=2, 
#                  inits=get_inits)
# print("finished running test stan model")
# run_id <- paste0(Sys.info()[4], format(Sys.time(), "_%Y%m%d%H%M%S"))
# out_name <- file.path(out_folder, paste0("stan_fit", suffix, '-', run_id, ".RData"))
# save(stan_fit, file=out_name)
# print(paste("Finished", out_name))

model_data_blocked <- model_data_blocked[names(model_data_blocked) != "first_obs_period"]
model_data_blocked <- model_data_blocked[names(model_data_blocked) != "last_obs_period"]

stan_fit_initial <- stan(model_file, data=model_data_blocked, chains=1, iter=10, 
                         sample_file="stan_test_samples.csv")
# Test run
stan_fit_initial <- stan(model_file, data=model_data_blocked, chains=1, iter=1000, 
                         init=get_inits, sample_file="stan_test_samples.csv")

# Fit n_chains chains in parallel. Reuse same seed so that the chain_ids can be 
# used by stan to properly seed each chain differently.
# Fit initial model, on a single CPU. Run only one iteration.
stan_fit_initial <- stan(model_file, data=model_data_blocked, chains=0, init=get_inits)
print("finished setting up stan model")

n_chains <- 3
n_cpu <- n_chains
n_iter <- 1000
cl <- makeCluster(n_cpu)
registerDoParallel(cl)
run_id <- paste0(Sys.info()[4], format(Sys.time(), "_%Y%m%d%H%M%S"))
sflist <- foreach(n=1:n_chains, .packages=c("rstan")) %dopar% {
    # Add 1 to n in order to ensure chain_id 1 is not reused
    sink(paste0("stan_", run_id, "_chain", n, ".txt"), append=TRUE)
    sample_file_csv <- paste0("stan_", run_id, "_chain", n, "_samples.csv")
    this_stanfit <- stan(fit=stan_fit_initial, data=model_data_blocked, seed=seed, 
                         chains=1, iter=n_iter, chain_id=n, init=get_inits, 
                         sample_file=sample_file_csv)
    save(this_stan_fit, file=file.path(out_folder, paste0("stan_fit", suffix, '-', run_id, "_chain", n,  ".RData")))
    sink()
    return(this_stanfit)
}
print("finished running stan models on cluster")
stopCluster(cl)

stan_fit <- sflist2stanfit(sflist)
out_name <- file.path(out_folder, paste0("stan_fit", suffix, '-', run_id, "_fullfit.RData"))
save(stan_fit, file=out_name)
print(paste("Finished", out_name))
