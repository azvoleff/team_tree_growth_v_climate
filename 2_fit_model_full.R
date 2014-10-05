library(runjags)

model_file <- "full_model.bug" 

load("model_data.RData")
load("init_data.RData")

monitored <- c("intercept",
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
               "sigma_k",
               "sigma_g",
               "sigma_t")

init_data[[1]] <- init_data[[1]][names(init_data[[1]]) %in% c("dbh_latent")]

model_data$WD_sq <- model_data$WD^2
# Drop missing data indicators (not needed for JAGS)
model_data <- model_data[!(names(model_data) %in% c("miss_indices", "obs_indices"))]

# seq_n_chains <- 2
# jags_fit <- run.jags(model=model_file, monitor=monitored, data=model_data, 
#                      inits=rep(init_data, seq_n_chains), n.chains=seq_n_chains, 
#                      sample=100)
# print("finished running single JAGS chain")
# save(jags_fit, file="jags_fit_full.RData")

jags_fit_p <- run.jags(model=model_file, monitor=monitored, data=model_data, 
                       inits=rep(init_data, 4), n.chains=4, method="parallel")
print("finished running JAGS chains in parallel")
save(jags_fit_p, file="jags_fit_full_parallel.RData")
