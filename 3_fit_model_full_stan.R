library(rstan)
library(foreach)
library(doParallel)

stan_path <- "C:/cmdstan"

source("0_settings.R")

model_file <- "full_model.stan" 
model_name <- "full_model" 

temp_var <- "tmp_meanannual"
precip_var <- "mcwd_run12"
#model_type <- "full"
model_type <- "testing"

suffix <- paste0('_', model_type, '-', temp_var, '-', precip_var)

data_folder <- file.path(prefix, "TEAM", "Tree_Growth", "Data")
init_folder <- file.path(prefix, "TEAM", "Tree_Growth", "Initialization")
mcmc_folder <- file.path(prefix, "TEAM", "Tree_Growth", "MCMC_Chains")

load(file.path(data_folder, paste0("model_data_wide_blocked", suffix, ".RData")))
load(file.path(init_folder, paste0("init_data_with_ranefs", suffix, ".RData")))
load(file.path(init_folder, paste0("init_data_blocked", suffix, ".RData")))

init_data$dbh_latent <- init_data_blocked$dbh_latent

# n_B is number of fixed effects
model_data_blocked$n_B <- 2
# n_B_g is number of genus-level random effects
model_data_blocked$n_B_g <- 7
# n_B_T is number of terms in the temperature model
model_data_blocked$n_B_T <- 2
# W is prior scale for the inverse-Wishart
model_data_blocked$W <- diag(model_data_blocked$n_B_g)

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
        B_T=matrix(rnorm(model_data_blocked$n_B_T*model_data_blocked$n_site, 0, 1),
                   ncol=model_data_blocked$n_B_T),
        # Sigmas
        sigma_obs=runif(1, model_data_blocked$sigma_obs_lower, .01), 
        sigma_proc=abs(rnorm(1, 0, 1))
    ))
}

# stan_fit <- stan(model_file, data=model_data_blocked, iter=20, chains=2, 
#                  inits=get_inits)
# print("finished running test stan model")
# run_id <- paste0(Sys.info()[4], format(Sys.time(), "_%Y%m%d%H%M%S"))
# out_name <- file.path(mcmc_folder, paste0("stan_fit", suffix, '-', run_id, ".RData"))
# save(stan_fit, file=out_name)
# print(paste("Finished", out_name))

model_data_blocked <- model_data_blocked[names(model_data_blocked) != "last_obs_period"]

# Test run
# stan_fit <- stan(model_file, data=model_data_blocked, chains=2, iter=20, 
#                  init=get_inits, control=list(refresh=1))

# Fit n_chains chains in parallel. Reuse same seed so that the chain_ids can be 
# used by stan to properly seed each chain differently.
# Fit initial model, on a single CPU. Run only one iteration.
# stan_fit_initial <- stan(model_file, data=model_data_blocked, chains=0, 
#                          init=get_inits)
print("finished setting up stan model")


# Run in CmdStan for memory reasons (so that output can be streamed directly to 
# disk.
n_chains <- 3
n_cpu <- n_chains
n_iter <- 1000
cl <- makeCluster(n_cpu)
registerDoParallel(cl)
run_id <- paste0(Sys.info()[4], format(Sys.time(), "_%Y%m%d%H%M%S"))

data_file <- file.path(mcmc_folder, paste0("model_data", suffix, ".R"))
attach(model_data_blocked)
stan_rdump(names(model_data_blocked), file=data_file)
detach(model_data_blocked)

sflist <- foreach(n=1:n_chains, .packages=c("rstan")) %dopar% {
    sink_file <- paste0("stan_", run_id, "_chain", n, "_sink.txt")
    init_file <- paste0("stan_", run_id, "_chain", n, "_inits.csv")
    sample_file <- paste0("stan_", run_id, "_chain", n, "_samples.csv")
    sink(sink_file, append=TRUE)
    inits <- get_inits()
    attach(inits)
    stan_rdump(names(inits), init_file)
    detach(inits)

                        " --o=", gsub(".stan", ".cpp", model_file),

    cpp_file <- file.path(stan_path, "mymodels", gsub(".stan", ".cpp", model_file))

    # Compile stan model (Windows)
    stanc_cmd <- paste0(file.path(stan_path, "bin", "stanc.exe"), " --name=", 
                        model_name, " --o=", cpp_file, " ", model_file)
    system(stanc_cmd)

    setwd(stan_path)
    make_cmd <- paste0("make ", gsub('.stan', '', model_file))
    system(make_cmd)


    # Sample
    sampl_cmd <- paste0("./", model_name, " sample random seed=12345 id=", n,
                        " data file=", data_file, " output file=", sample_file,
                        " &")


    # Windows
    sampl_cmd <- paste0("start /b ", model_name, " sample random seed=12345 id=", n,
                        " data file=", data_file, " output file=", sample_file)
    system(sampl_cmd)

    this_stanfit <- stan(fit=stan_fit_initial, data=model_data_blocked, seed=seed, 
                         chains=1, iter=n_iter, chain_id=n, init=get_inits, 
                         sample_file=sample_file_csv)

    save(this_stan_fit, file=file.path(mcmc_folder, paste0("stan_fit", suffix, '-', run_id, "_chain", n,  ".RData")))
    sink()
    return(this_stanfit)
}
print("finished running stan models on cluster")
stopCluster(cl)

stan_fit <- sflist2stanfit(sflist)
out_name <- file.path(mcmc_folder, paste0("stan_fit", suffix, '-', run_id, "_fullfit.RData"))
save(stan_fit, file=out_name)
print(paste("Finished", out_name))
