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

# jags.parallel only works when it is pulling the data from the R environment
n_tree <- model_data$n_tree
n_plot <- model_data$n_plot
n_site <- model_data$n_site
first_obs_period <- model_data$first_obs_period
last_obs_period <- model_data$last_obs_period
plot_ID <- model_data$plot_ID
site_ID <- model_data$site_ID
dbh <- model_data$dbh
WD <- model_data$WD
spi <- model_data$spi
WD_sq <- model_data$WD_sq
jags_fit_p <- jags.parallel(names(model_data), init_data,
                            parameters.to.save=jags_params, 
                            model.file="full_model.bug", n.chains=12, 
                            n.iter=100000, jags.seed=70276)
print("finished running JAGS chains in parallel")
save(jags_fit_p, file="jags_fit_full_parallel.RData")
