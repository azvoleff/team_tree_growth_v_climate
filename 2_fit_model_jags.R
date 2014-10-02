library(R2jags)

model_file <- "full_model.bug" 

load("model_data.RData")
load("init_data.RData")

seed <- 70276

jags_params <- c("intercept",
                 "slp_dbh",
                 "slp_dbh_sq",
                 "slp_WD",
                 "slp_WD_sq",
                 "slp_spi",
                 "inter_spi_dbh",
                 "inter_spi_WD",
                 "obs_sigma",
                 "proc_sigma",
                 "sigma_ijk",
                 "sigma_jk",
                 "sigma_k")

init_data[[1]] <- init_data[[1]][names(init_data[[1]]) %in% c("dbh_latent")]

model_data$WD_sq <- model_data$WD^2
# Drop missing data indicators (not needed for JAGS)
model_data <- model_data[!(names(model_data) %in% c("miss_indices", "obs_indices"))]

# set.seed(seed)
# seq_n_chains <- 10
# jags_fit <- jags(data=model_data, inits=rep(init_data, seq_n_chains),
#                  parameters.to.save=jags_params, n.chains=seq_n_chains, 
#                  n.iter=10000, model.file=model_file)
# print("finished running single JAGS chain")
# save(jags_fit, file="jags_fit_full.RData")

source("jagsparallel.R")
jags_fit_p <- jagsparallel(model_data, init_data,
                           parameters.to.save=jags_params, 
                           model.file="full_model.bug", n.chains=8, 
                           n.iter=100000, jags.seed=70276, n.cluster=8)
print("finished running JAGS chains in parallel")
save(jags_fit_p, file="jags_fit_full_parallel.RData")
