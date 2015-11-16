library(runjags)
library(rjags)

source("0_settings.R")

# Allow block-updating using glm module
load.module("glm")

#model_structure <- "simple"
#model_structure <- "correlated"
model_structure <- "interact"

precip_var <- 'mcwd_run12'

model_type <- "full"
#model_type <- "testing"

# Use the temp_min model just to get the base data
temp_var <- 'tmn_meanannual'
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

load(file.path(data_folder, paste0("model_data_wide_full-tmn_meanannual-mcwd_run12.RData")))
temp_min <- model_data$temp
load(file.path(data_folder, paste0("model_data_wide_full-tmp_meanannual-mcwd_run12.RData")))
temp_mean <- model_data$temp
load(file.path(data_folder, paste0("model_data_wide_full-tmx_meanannual-mcwd_run12.RData")))
temp_max <- model_data$temp

load(file.path(init_folder, paste0("init_data_with_ranefs_full-tmx_meanannual-mcwd_run12_interact.RData")))
model_file <- "growth_model_interact_allvars.bug"
# n_B_g is number of genus-level random effects
model_data$n_B_g <- 15
monitored <- monitored[monitored != "B"]
monitored <- monitored[monitored != "rho_B_g"]
init_data <- init_data[names(init_data) != 'B_g_raw']
init_data <- init_data[names(init_data) != 'sigma_B_g']
model_type <- paste0(model_type, "_interact")

model_data <- model_data[names(model_data) != "temp"]
model_data$temp_min <- temp_min
model_data$temp_mean <- temp_mean
model_data$temp_max <- temp_max
model_data$temp_min_sq <- temp_min^2
model_data$temp_mean_sq <- temp_mean^2
model_data$temp_max_sq <- temp_max^2
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
                     burnin=50000, sample=200, thin=100, summarise=FALSE)
print("finished running JAGS chains in parallel")
run_id <- paste0(Sys.info()[4], format(Sys.time(), "_%Y%m%d%H%M%S"))
out_name <- file.path(mcmc_folder, paste0("jags_fit", out_suffix, '-', run_id, ".RData"))
save(jags_fit, file=out_name)
print(paste("Finished", out_name))
