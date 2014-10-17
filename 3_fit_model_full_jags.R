library(runjags)
library(rjags)

# Allow block-updating using glm module
load.module("glm")

model_file <- "full_model.bug" 

load("model_data_wide.RData")
load("init_data.RData")

monitored <- c("int",
               "slp_dbh",
               "slp_dbh_sq",
               "slp_WD",
               "slp_WD_sq",
               "slp_mcwd",
               "slp_mcwd_sq",
               "sigma_obs",
               "sigma_proc",
               "sigma_int_ijk",
               "sigma_int_jk",
               "sigma_int_k",
               "sigma_int_t",
               "mu",
               "sigma_B",
               "rho_B")

init_data <- init_data[names(init_data) %in% c("dbh_latent")]

model_data$K <- 5
model_data$W <- diag(model_data$K)
model_data$WD_sq <- model_data$WD^2
model_data$mcwd_sq <- model_data$mcwd^2
# Drop missing data indicators (not needed for JAGS)
model_data <- model_data[!(names(model_data) %in% c("miss_indices", "obs_indices"))]

# seq_n_chains <- 1
# jags_fit <- run.jags(model=model_file, monitor=monitored, data=model_data, 
#                      inits=rep(list(init_data), seq_n_chains), 
#                      n.chains=seq_n_chains, adapt=100, burnin=100, 
#                      sample=100, thin=2)
# print("finished running single JAGS chain")
# run_id <- paste0(Sys.info()[4], format(Sys.time(), "_%Y%m%d-%H%M%S"))
# save(jags_fit, file=paste0("full_model_fit_", run_id, ".RData"))

init_data$int <- 0
init_data$slp_dbh <- 0
init_data$slp_dbh_sq <- 0
init_data$slp_WD <- 0
init_data$slp_WD_sq <- 0
init_data$slp_mcwd <- 0

jags_fit <- run.jags(model=model_file, monitor=monitored,
                     data=model_data, inits=rep(list(init_data), 4),
                     n.chains=4, method="parallel", adapt=2000,
                     burnin=5000, sample=10000)
print("finished running JAGS chains in parallel")
run_id <- paste0(Sys.info()[4], format(Sys.time(), "_%Y%m%d-%H%M%S"))
save(jags_fit, file=paste0("full_model_fit_parallel_", run_id, ".RData"))
