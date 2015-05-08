library(runjags)
library(rjags)

source("0_settings.R")

# Allow block-updating using glm module
load.module("glm")

#model_structure <- "simple"
#model_structure <- "correlated"
model_structure <- "interact"

temp_var <- 'tmn_meanannual'
#temp_var <- 'tmp_meanannual'
#temp_var <- 'tmx_meanannual'

precip_var <- 'mcwd_run12'

model_type <- "full"
#model_type <- "testing"

in_suffix <- paste0('_', model_type, '-', temp_var, '-', precip_var)
if (note != "") in_suffix <- paste0(in_suffix, '_', note)
out_suffix <- paste0(in_suffix, '_', model_structure)

monitored <- c("B",
               "sigma_obs",
               "sigma_proc",
               "sigma_int_jk",
               "sigma_int_k",
               "sigma_int_t",
               "int_jk",
               "int_k",
               "int_t",
               "B_k",
               "sigma_B_k",
               "B_g",
               "mu_B_g",
               "sigma_B_g",
               "rho_B_g")

data_folder <- file.path(prefix, "TEAM", "Tree_Growth", "Data")
init_folder <- file.path(prefix, "TEAM", "Tree_Growth", "Initialization")
mcmc_folder <- file.path(prefix, "TEAM", "Tree_Growth", "MCMC_Chains")

load(file.path(data_folder, paste0("model_data_wide", in_suffix, ".RData")))
load(file.path(data_folder, paste0("model_data_standardizing", in_suffix, ".RData")))

model_data <- model_data[names(model_data) != "WD"]

if (model_structure == "simple") {
    load(file.path(init_folder, paste0("init_data_with_ranefs", in_suffix, "_correlated.RData")))
    model_file <- "growth_model_simple.bug" 
    # n_B is number of fixed effects
    model_data$n_B <- 7
    init_data$B <- rnorm(model_data$n_B, 0, 1)
    monitored <- monitored[monitored != "B_g"]
    monitored <- monitored[monitored != "mu_B_g"]
    monitored <- monitored[monitored != "rho_B_g"]
    monitored <- monitored[monitored != "sigma_B_g"]
    model_data <- model_data[names(model_data) != "genus_ID"]
    model_data <- model_data[names(model_data) != "n_genus"]
    init_data <- init_data[names(init_data) != 'B_g_raw']
    init_data <- init_data[names(init_data) != 'sigma_B_g']
} else if (model_structure == "correlated") {
    load(file.path(init_folder, paste0("init_data_with_ranefs", in_suffix, "_correlated.RData")))
    model_file <- "growth_model_correlated.bug" 
    # n_B_g is number of genus-level random effects
    model_data$n_B_g <- 7
    monitored <- monitored[monitored != "B"]
    # W is the prior scale for the inverse-wishart
    model_data$W <- diag(model_data$n_B_g)
    # Initialize xi to the standard deviations of the genus-level effects, in
    # an effort to get the scale of mu_B_g_raw and Tau_B_g_raw in the right
    # ballpark
    init_data$xi <- apply(init_data$B_g_raw, 2, sd)
    init_data$mu_B_g_raw <- apply(init_data$B_g_raw, 2, mean) / init_data$xi
    # Initialize mean of group-level intercept to zero
    init_data$mu_B_g_raw[1] <- 0
    # Center the B_g_raw estimates
    init_data$B_g_raw <- init_data$B_g_raw - matrix(rep(init_data$mu_B_g_raw,
                                                        model_data$n_genus),
                                                    ncol=model_data$n_B_g,
                                                    byrow=TRUE)
    # Jags uses the inverse of the variance-covariance matrix to parameterize 
    # the wishart.
    #init_data$Tau_B_g_raw <- solve(diag(init_data$xi)) %*% init_data$sigma_B_g %*% solve(diag(init_data$xi))
    init_data <- init_data[!(names(init_data) %in% c("sigma_B_g"))]
    init_data <- init_data[!(names(init_data) %in% c("Tau_B_g_raw"))]
} else if (model_structure == "interact") {
    load(file.path(init_folder, paste0("init_data_with_ranefs", in_suffix, "_interact.RData")))
    model_file <- "growth_model_interact.bug"
    # n_B_g is number of genus-level random effects
    model_data$n_B_g <- 9
    monitored <- monitored[monitored != "B"]
    monitored <- monitored[monitored != "rho_B_g"]
    init_data <- init_data[names(init_data) != 'B_g_raw']
    init_data <- init_data[names(init_data) != 'sigma_B_g']
    model_type <- paste0(model_type, "_interact")
} else {
    stop(paste0('Unknown model_structure "', model_structure, '"'))
}

model_data$precip_sq <- model_data$precip^2

# Drop missing data indicators (not needed for JAGS)
model_data <- model_data[!(names(model_data) %in% c("miss_indices", "obs_indices"))]

# jags_fit <- run.jags(model=model_file, monitor=monitored, data=model_data, 
# inits=rep(list(init_data), 1),  n.chains=1, adapt=200, burnin=200, 
# sample=200)
# print("finished running single JAGS chain")
# run_id <- paste0(Sys.info()[4], format(Sys.time(), "_%Y%m%d%H%M%S"))
# out_name <- file.path(mcmc_folder, paste0("jags_fit", out_suffix, '-', run_id, ".RData"))
# save(jags_fit, file=out_name)
# print(paste("Finished", out_name))

n_chains <- 6
jags_fit <- run.jags(model=model_file, monitor=monitored,
                     data=model_data, inits=rep(list(init_data), n_chains),
                     n.chains=n_chains, method="parallel", adapt=1000,
                     burnin=10000, sample=200, thin=100, summarise=FALSE)
print("finished running JAGS chains in parallel")
run_id <- paste0(Sys.info()[4], format(Sys.time(), "_%Y%m%d%H%M%S"))
out_name <- file.path(mcmc_folder, paste0("jags_fit", out_suffix, '-', run_id, ".RData"))
save(jags_fit, file=out_name)
print(paste("Finished", out_name))
