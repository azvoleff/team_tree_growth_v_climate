library(R2jags)

model_file <- "full_model.bug" 
n_chains <- 12
n_burnin <- 12
n_iter <- 1000
n_thin <- 10

load("model_data.RData")
load("init_data.RData")

# Vector of zeros (see http://bit.ly/1uFTjlx)
zeros <- rep(0, n_obs)

seed <- 1638

jags_params <- c("log_dbh_latent",
                 "inter",
                 "slp_dbh",
                 "sigma_obs",
                 "sigma_proc",
                 "sigma_ijk",
                 "sigma_jk",
                 "sigma_k",
                 "b_ijk",
                 "b_jk",
                 "b_k")

jags_fit <- jags(data=jags_vars, inits=jags_inits, 
                 parameters.to.save=jags_params, n.chains=n_chains, 
                 n.burnin=n_burnin, n.iter=n_iter, n.thin=n_thin, 
                 model.file=model_file)
print("finished running single JAGS chain")
save(jags_fit_parallel, file="jags_fit_parallel.RData")


jags_fit_parallel <- jags.parallel(data=model_data, inits=model_inits, 
                                   parameters.to.save=jags_params, 
                                   model.file=model_file, n.chains=n_chains, 
                                   n.burnin=n_burnin, n.iter=n_iter, 
                                   n.thin=n_thin)
print("finished running JAGS chains in parallel")
save(jags_fit_parallel, file="jags_fit_parallel.RData")
