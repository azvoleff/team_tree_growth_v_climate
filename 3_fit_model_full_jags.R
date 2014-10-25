library(runjags)
library(rjags)

# Allow block-updating using glm module
#load.module("glm")

model_file <- "full_model.bug" 

load("model_data_wide.RData")
load("init_data_with_ranefs.RData")

monitored <- c("B",
               "sigma_obs",
               "sigma_proc",
               "sigma_int_ijk",
               "sigma_int_jk",
               "sigma_int_k",
               "sigma_int_t",
               "mu_B_g",
               "sigma_B_g",
               "rho_B_g")

# n_B is number of fixed effects
model_data$n_B <- 2
# n_B_g is number of genus-level random effects
model_data$n_B_g <- 5
# W is prior scale for the inverse-Wishart
model_data$W <- diag(model_data$n_B_g)
model_data$WD_sq <- model_data$WD^2
model_data$mcwd_sq <- model_data$mcwd^2
# Drop missing data indicators (not needed for JAGS)
model_data <- model_data[!(names(model_data) %in% c("miss_indices", "obs_indices"))]
model_data <- model_data[!(names(model_data) %in% c("spi"))]

init_data$B <- rep(0, model_data$n_B)
init_data$xi <- rep(1, model_data$n_B_g)
init_data$mu_B_g_raw <- apply(init_data$B_g_raw, 2, mean) / init_data$xi
# Jags uses the inverse of the variance-covariance matrix to parameterize the 
# wishart.
init_data$Tau_B_g_raw <- solve(init_data$sigma_B_g)
init_data <- init_data[!(names(init_data) %in% c("sigma_B_g"))]

# seq_n_chains <- 1
# jags_fit <- run.jags(model=model_file, monitor=monitored, data=model_data, 
#                      inits=rep(list(init_data), seq_n_chains), 
#                      n.chains=seq_n_chains, adapt=200, burnin=200, 
#                      sample=200)
# print("finished running single JAGS chain")
# run_id <- paste0(Sys.info()[4], format(Sys.time(), "_%Y%m%d-%H%M%S"))
# save(jags_fit, file=paste0("full_model_fit_", run_id, ".RData"))

jags_fit <- run.jags(model=model_file, monitor=monitored,
                     data=model_data, inits=rep(list(init_data), 3),
                     n.chains=3, method="parallel", adapt=500,
                     burnin=1250, sample=2000, thin=3)
print("finished running JAGS chains in parallel")
run_id <- paste0(Sys.info()[4], format(Sys.time(), "_%Y%m%d-%H%M%S"))
save(jags_fit, file=paste0("full_model_fit_parallel_", run_id, ".RData"))
