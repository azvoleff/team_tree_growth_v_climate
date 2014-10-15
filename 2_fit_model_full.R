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
               "slp_spi",
               "sigma_obs",
               "sigma_proc",
               "sigma_int_ijk",
               "sigma_int_jk",
               "sigma_int_k",
               "sigma_int_t",
               "mu_int_g",
               "mu_slp_g_spi",
               "mu_slp_g_dbh",
               "mu_slp_g_dbh_sq",
               "sigma_int_g",
               "sigma_slp_g_spi",
               "sigma_slp_g_dbh",
               "sigma_slp_g_dbh_sq",
               "rho_int_g_slp_g_spi",
               "rho_int_g_slp_g_dbh",
               "rho_int_g_slp_g_dbh_sq",
               "rho_slp_g_spi_slp_g_dbh",
               "rho_slp_g_spi_slp_g_dbh_sq",
               "rho_slp_g_dbh_slp_g_dbh_sq")

init_data <- init_data[names(init_data) %in% c("dbh_latent")]

model_data$W <- diag(4)
model_data$WD_sq <- model_data$WD^2
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
init_data$slp_spi <- 0
init_data$mu_int_g <- 0
init_data$mu_slp_g_spi <- 0
init_data$mu_slp_g_dbh <- 0
init_data$mu_slp_g_dbh_sq <- 0
init_data$rho_int_g_slp_g_spi <- 0
init_data$rho_int_g_slp_g_dbh <- 0
init_data$rho_int_g_slp_g_dbh_sq <- 0
init_data$rho_slp_g_spi_slp_g_dbh <- 0
init_data$rho_slp_g_spi_slp_g_dbh_sq <- 0

jags_fit <- run.jags(model=model_file, monitor=monitored, data=model_data, 
                     inits=rep(list(init_data), 4), n.chains=4, 
                     method="parallel", adapt=2000, burnin=10000, sample=5000, 
                     thin=2)
print("finished running JAGS chains in parallel")
run_id <- paste0(Sys.info()[4], format(Sys.time(), "_%Y%m%d-%H%M%S"))
save(jags_fit, file=paste0("full_model_fit_parallel_", run_id, ".RData"))
