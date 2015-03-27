library(runjags)
library(rjags)

source("0_settings.R")

# Allow block-updating using glm module
#load.module("glm")

# Include a random intercept by period?
model_structure <- "simple"
#model_structure <- "full_model"
#model_structure <- "full_model_no_t_effects"
#model_structure <- "full_model_no_t_effects_interact"

note <- 'genuslimits'

# temp_var <- "tmn_meanannual"
# temp_var <- "tmp_meanannual"
temp_var <- "tmx_meanannual"
precip_var <- "mcwd_run12"

#model_type <- "full"
model_type <- "testing"

in_suffix <- paste0('_', model_type, '-', temp_var, '-', precip_var)
out_suffix <- paste0(in_suffix, '_', note)

monitored <- c("B",
               "B_T_int",
               "B_T_lapse",
               "sigma_obs",
               "sigma_proc",
               "sigma_int_ijk",
               "sigma_int_jk",
               "sigma_int_k",
               "sigma_int_t",
               "int_jk",
               "int_k",
               "int_t",
               "B_g",
               "mu_B_g",
               "sigma_B_g",
               "rho_B_g")

data_folder <- file.path(prefix, "TEAM", "Tree_Growth", "Data")
init_folder <- file.path(prefix, "TEAM", "Tree_Growth", "Initialization")
mcmc_folder <- file.path(prefix, "TEAM", "Tree_Growth", "MCMC_Chains")

load(file.path(init_folder, paste0("init_data_with_ranefs", in_suffix,  ".RData")))
load(file.path(data_folder, paste0("model_data_wide", in_suffix, ".RData")))
load(file.path(data_folder, paste0("model_data_standardizing", in_suffix, ".RData")))

if (model_structure == "simple") {
    model_file <- "simple_model.bug" 
    # n_B is number of fixed effects
    model_data$n_B <- 9
    monitored <- monitored[monitored != "B_g"]
    monitored <- monitored[monitored != "mu_B_g"]
    monitored <- monitored[monitored != "rho_B_g"]
    monitored <- monitored[monitored != "sigma_B_g"]
    model_data <- model_data[names(model_data) != "genus_ID"]
    model_data <- model_data[names(model_data) != "n_genus"]
} else if (model_structure == "full_model") {
    model_file <- "full_model.bug" 
    # n_B_g is number of genus-level random effects
    n_B_g <- 7
    # W is the prior scale for the inverse-wishart
    model_data$W <- diag(n_B_g)
} else if (model_structure == "full_model_no_t_effects") {
    model_file <- "full_model_no_t_effects.bug"
    monitored <- monitored[monitored != "int_t"]
    monitored <- monitored[monitored != "sigma_int_t"]
    model_data <- model_data[names(model_data) != "n_period"]
    model_type <- paste0(model_type, "_no_t_effects")
    # n_B_g is number of genus-level random effects
    n_B_g <- 7
    # W is the prior scale for the inverse-wishart
    model_data$W <- diag(n_B_g)
} else if (model_structure == "full_model_no_t_effects_interact") {
    model_file <- "full_model_no_t_effects_interact.bug"
    monitored <- monitored[monitored != "int_t"]
    monitored <- monitored[monitored != "sigma_int_t"]
    monitored <- monitored[monitored != "rho_B_g"]
    init_data <- init_data[!grepl('int', names(init_data))]
    model_data <- model_data[names(model_data) != "n_period"]
    model_type <- paste0(model_type, "_no_t_effects_interact")
    n_B_g <- 11
} else {
    stop(paste0('Unknown model_structure "', model_structure, '"'))
}

if (model_structure != "simple") {
    # n_B_g is number of genus-level random effects
    model_data$n_B_g <- n_B_g
    # n_B is number of fixed effects
    model_data$n_B <- 2
}

# n_B_T is number of terms in the temperature model
n_B_T <- 2

model_data$WD_sq <- model_data$WD^2
model_data$precip_sq <- model_data$precip^2

# Drop missing data indicators (not needed for JAGS)
model_data <- model_data[!(names(model_data) %in% c("miss_indices", "obs_indices"))]
model_data <- model_data[!(names(model_data) %in% c("spi"))]

# Setup mean for lapse rate prior
model_data$lapse_mean <- -6.5 / temp_sd
# Recall the precision is 1 over the variance. Define the lapse rate prior to 
# have a standard deviation of 2.5 degrees
model_data$lapse_prec <- (2.5 / temp_sd)^-2

init_data$B <- rnorm(model_data$n_B, 0, 10)
init_data$B_T_int <- rnorm(model_data$n_site, 0, 10)
# Constrain lapse rate to be negative
init_data$B_T_lapse <- -abs(rnorm(model_data$n_site, model_data$lapse_mean, (model_data$lapse_prec)^-2))

# seq_n_chains <- 1
# jags_fit <- run.jags(model=model_file, monitor=monitored, data=model_data, 
#                      inits=rep(list(init_data), seq_n_chains), 
#                      n.chains=seq_n_chains, adapt=200, burnin=200, 
#                      sample=200)
# print("finished running single JAGS chain")
# run_id <- paste0(Sys.info()[4], format(Sys.time(), "_%Y%m%d%H%M%S"))
# out_name <- file.path(mcmc_folder, paste0("jags_fit", out_suffix, '-', run_id, ".RData"))
# save(jags_fit, file=out_name)
# print(paste("Finished", out_name))

jags_fit <- run.jags(model=model_file, monitor=monitored,
                     data=model_data, inits=rep(list(init_data), 3),
                     n.chains=3, method="parallel", adapt=1000,
                     burnin=2500, sample=2500, thin=4, summarise=FALSE)
print("finished running JAGS chains in parallel")
run_id <- paste0(Sys.info()[4], format(Sys.time(), "_%Y%m%d%H%M%S"))
out_name <- file.path(mcmc_folder, paste0("jags_fit", out_suffix, '-', run_id, ".RData"))
save(jags_fit, file=out_name)
print(paste("Finished", out_name))

# print(paste("Starting autorun", out_name))
# jags_fit <- autorun.jags(jags_fit, summarise=FALSE, max.time="10 days")
# autorun_out_name <- file.path(mcmc_folder, paste0("jags_fit", out_suffix, '-', run_id, "_autorun.RData"))
# save(jags_fit, file=autorun_out_name)
# print(paste("Finished", autorun_out_name))
