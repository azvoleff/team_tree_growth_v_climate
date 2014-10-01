library(R2jags)

model_file <- "simple_model_eitzen.bug" 
n_chains <- 8
n_burnin <- 10000
n_iter <- 10000
n_thin <- 10

load("model_data.RData")
load("init_data.RData")

seed <- 70276

jags_params <- c("dbh_latent",
                 "inter",
                 "slp_dbh",
                 "slp_dbh_sq",
                 "obs_sigma",
                 "proc_sigma")

# Drop missing data indicators (not needed for JAGS)
init_data[[1]] <- init_data[[1]][names(init_data[[1]]) %in% c("dbh_latent")]

model_data$dbh <- exp(model_data$log_dbh)
model_data <- model_data[names(model_data) %in% c("n_tree", 
                                                  "first_obs_period", 
                                                  "last_obs_period", "dbh")]
# misses <- calc_missings(model_data$dbh)
# model_data$dbh <- model_data$dbh[-unique(misses$miss[, 1]), ]
# model_data$first_obs_period <- model_data$first_obs_period[-unique(misses$miss[, 1])]
# model_data$last_obs_period <- model_data$last_obs_period[-unique(misses$miss[, 1])]
# model_data$n_tree <- nrow(model_data$dbh)

set.seed(seed)
jags_fit <- jags(data=model_data, inits=rep(init_data, 1),
                 parameters.to.save=jags_params, n.chains=1, 
                 n.burnin=1000, n.iter=1000, n.thin=n_thin, 
                 model.file=model_file)
print("finished running single JAGS chain")
save(jags_fit, file="jags_fit.RData")

jags_fit_p <- jags.parallel(data=model_data, inits=rep(init_data, n_chains), 
                            parameters.to.save=jags_params, 
                            model.file=model_file, n.chains=n_chains, 
                            n.burnin=n_burnin, n.iter=n_iter, n.thin=n_thin, 
                            jags.seed=seed)
print("finished running JAGS chains in parallel")
save(jags_fit_p, file="jags_fit_parallel.RData")
