library(runjags)
library(rjags)

# Allow block-updating using glm module
#load.module("glm")

model_file <- "full_model.bug" 

temp_var <- "tmn_meanannual"
precip_var <- "mcwd_run12"
model_type <- "full"
#model_type <- "testing"
in_folder <- 'Data'
out_folder <- 'MCMC_Chains'

suffix <- paste0('_', model_type, '-', temp_var, '-', precip_var)

load(file.path(in_folder, paste0("model_data_wide", suffix, ".RData")))
load(file.path(in_folder, paste0("init_data_with_ranefs", suffix, ".RData")))

monitored <- c("B",
               "B_T",
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

# n_B is number of fixed effects
model_data$n_B <- 2
# n_B_g is number of genus-level random effects
model_data$n_B_g <- 7
# n_B_T is number of terms in the temperature model
model_data$n_B_T <- 3
# W is prior scale for the inverse-Wishart
model_data$W <- diag(model_data$n_B_g)
model_data$WD_sq <- model_data$WD^2
model_data$precip_sq <- model_data$precip^2
# Drop missing data indicators (not needed for JAGS)
model_data <- model_data[!(names(model_data) %in% c("miss_indices", "obs_indices"))]
model_data <- model_data[!(names(model_data) %in% c("spi"))]

init_data$B <- rep(0, model_data$n_B)
init_data$B_T <- matrix(0, nrow=model_data$n_site, ncol=model_data$n_B_T)
# Initialize xi to the standard deviations of the genus-level effects, in an 
# effort to get the scale of mu_B_g_raw and Tau_B_g_raw in the right ballpark
init_data$xi <- apply(init_data$B_g_raw, 2, sd)
init_data$mu_B_g_raw <- apply(init_data$B_g_raw, 2, mean) / init_data$xi
# Center the B_g_raw estimates
init_data$B_g_raw <- init_data$B_g_raw - matrix(rep(init_data$mu_B_g_raw, 
                                                    model_data$n_genus), 
                                                ncol=model_data$n_B_g, 
                                                byrow=TRUE)

# Jags uses the inverse of the variance-covariance matrix to parameterize the 
# wishart.
init_data$Tau_B_g_raw <- solve(diag(init_data$xi)) %*% init_data$sigma_B_g %*% solve(diag(init_data$xi))
init_data <- init_data[!(names(init_data) %in% c("sigma_B_g"))]

# seq_n_chains <- 1
# jags_fit <- run.jags(model=model_file, monitor=monitored, data=model_data, 
#                      inits=rep(list(init_data), seq_n_chains), 
#                      n.chains=seq_n_chains, adapt=100, burnin=100, 
#                      sample=100)
# print("finished running single JAGS chain")
# run_id <- paste0(Sys.info()[4], format(Sys.time(), "_%Y%m%d%H%M%S"))
# save(jags_fit, file=file.path(out_folder, paste0("jags_fit", suffix, '-', run_id, ".RData")))

jags_fit <- run.jags(model=model_file, monitor=monitored,
                     data=model_data, inits=rep(list(init_data), 3),
                     n.chains=3, method="parallel", adapt=500,
                     burnin=2500, sample=2500, thin=4)
print("finished running JAGS chains in parallel")
run_id <- paste0(Sys.info()[4], format(Sys.time(), "_%Y%m%d%H%M%S"))
save(jags_fit, file=file.path(out_folder, paste0("jags_fit", suffix, '-', run_id, ".RData")))
