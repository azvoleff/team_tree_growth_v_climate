library(R2jags)

model_file <- "simple_model.bug" 
n_chains <- 8
n_burnin <- 10000
n_iter <- 10000
n_thin <- 10

load("model_data.RData")
load("init_data.RData")

# Add a vector of zeros to support the "zeros trick" (see 
# http://bit.ly/1uFTjlx)
model_data$zeros <- matrix(0, nrow=model_data$n_tree, 
                           ncol=ncol(model_data$log_dbh))

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

# Drop missing data indicators (not needed for JAGS)
model_data <- model_data[!(names(model_data) %in% c("obs", "miss"))]

init_data[[1]] <- init_data[[1]][names(init_data[[1]]) %in% c("log_dbh_latent")]
model_data <- model_data[!(names(model_data) %in% c("n_plot", "n_site", "plot_ID", "site_ID", "WD", "spi"))]

set.seed(seed)
jags_fit_simple <- jags(data=model_data, inits=rep(init_data, 1),
                 parameters.to.save=jags_params, n.chains=1, 
                 n.burnin=1000, n.iter=1000, n.thin=n_thin, 
                 model.file=model_file)
print("finished running single JAGS chain")
save(jags_fit_simple, file="jags_fit_simple.RData")

jags_fit_simple_p <- jags.parallel(data=model_data, inits=rep(init_data, n_chains), 
                            parameters.to.save=jags_params, 
                            model.file=model_file, n.chains=n_chains, 
                            n.burnin=n_burnin, n.iter=n_iter, n.thin=n_thin, 
                            jags.seed=seed)
print("finished running JAGS chains in parallel")
save(jags_fit_simple_p, file="jags_fit_simple_parallel.RData")
