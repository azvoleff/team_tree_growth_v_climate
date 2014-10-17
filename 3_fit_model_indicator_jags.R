library(runjags)

model_file <- "full_model_indicator.bug" 

load("model_data_wide.RData")
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
               "b_k",
               "sigma_t",
               "b_t",
               "b_g",
               "mu_b_g",
               "sigma_b_g",
               "slp_spi_g",
               "mu_slp_spi_g",
               "sigma_slp_spi_g",
               "rho_g",
               "z_slp_dbh",
               "z_slp_dbh_sq",
               "z_slp_WD",
               "z_slp_WD_sq",
               "z_slp_spi",
               "z_inter_spi_dbh",
               "z_inter_spi_WD")

init_data[[1]] <- init_data[[1]][names(init_data[[1]]) %in% c("dbh_latent")]

model_data$WD_sq <- model_data$WD^2
# Drop missing data indicators (not needed for JAGS)
model_data <- model_data[!(names(model_data) %in% c("miss_indices", "obs_indices"))]

# seq_n_chains <- 1
# jags_fit <- run.jags(model=model_file, monitor=monitored, data=model_data, 
#                      inits=rep(init_data, seq_n_chains), n.chains=seq_n_chains, 
#                      sample=100)
# print("finished running single JAGS chain")
# save(jags_fit,
#      file=paste0("jags_fit_indicator_",
#                  format(Sys.time(), "%Y%m%d-%H%M%S"), ".RData"))

jags_fit_p <- run.jags(model=model_file, monitor=monitored, data=model_data, 
                       inits=rep(init_data, 4), n.chains=4, method="parallel",
                       burnin=10000, sample=10000, thin=3)
print("finished running JAGS chains in parallel")
save(jags_fit_p,
     file=paste0("jags_fit_indicator_parallel_",
                 format(Sys.time(), "%Y%m%d-%H%M%S"), ".RData"))
