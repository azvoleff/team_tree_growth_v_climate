library(R2jags)

model_file <- "simple_model.bug" 
n_chains <- 8
n_iter <- 20000

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
               "proc_sigma")

# Drop missing data indicators (not needed for JAGS)
init_data[[1]] <- init_data[[1]][names(init_data[[1]]) %in% c("dbh_latent")]

model_data$WD_sq <- model_data$WD^2
model_data <- model_data[!(names(model_data) %in% c("n_plot", "n_site", 
                                                    "plot_ID", "site_ID",
                                                    "miss_indices", "obs_indices"))]

seq_n_chains <- 1
jags_fit <- run.jags(model=model_file, monitor=monitored, data=model_data, 
                     inits=rep(init_data, seq_n_chains), n.chains=seq_n_chains, 
                     sample=100)
print("finished running single JAGS chain")
save(jags_fit,
     file=paste0("jags_fit_simple_",
                 format(Sys.time(), "%Y%m%d-%H%M%S"), ".RData"))

jags_fit_p <- run.jags(model=model_file, monitor=monitored, data=model_data, 
                       inits=rep(init_data, 4), n.chains=4, method="parallel")
print("finished running JAGS chains in parallel")
print("finished running JAGS chains in parallel")
save(jags_fit_p,
     file=paste0("jags_fit_simple_parallel_",
                 format(Sys.time(), "%Y%m%d-%H%M%S"), ".RData"))
