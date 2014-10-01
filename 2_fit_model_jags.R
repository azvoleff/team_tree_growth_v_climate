library(R2jags)

model_file <- "full_model.bug" 
n_chains <- 8
n_iter <- 20000

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

# Drop missing data indicators (not needed for JAGS)
init_data[[1]] <- init_data[[1]][names(init_data[[1]]) %in% c("dbh_latent")]

model_data$WD_sq <- model_data$WD^2
model_data <- model_data[!(names(model_data) %in% c("miss_indices", "obs_indices"))]

set.seed(seed)
jags_fit <- jags(data=model_data, inits=rep(init_data, 1),
                 parameters.to.save=jags_params, n.chains=1, 
                 n.iter=1000, model.file=model_file)
print("finished running single JAGS chain")
save(jags_fit, file="jags_fit.RData")

jags_fit_p <- jags.parallel(data=model_data, inits=rep(init_data, n_chains), 
                            parameters.to.save=jags_params, 
                            model.file=model_file, n.chains=n_chains, 
                            n.iter=n_iter, jags.seed=seed)
print("finished running JAGS chains in parallel")
save(jags_fit_p, file="jags_fit_parallel.RData")
