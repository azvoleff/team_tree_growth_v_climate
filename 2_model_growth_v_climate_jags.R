library(dplyr)
library(ggplot2)
library(ggmcmc)
library(R2jags)

load("growth_ctfsflagged_merged_detrended.RData")

img_height <- 4
img_width <- 3
img_dpi <- 300

growth <- tbl_df(growth)

growth <- filter(growth, ctfs_accept)
growth <- filter(growth, n_days > 200)
growth <- filter(growth, n_days < 550)

# Exclude NAK, CSN, and YAN - too little data at these sites
#growth <- filter(growth, !(sitecode %in% c("NAK", "CSN", "YAN")))

growth$SamplingPeriodID <- with(growth, factor(factor(sitecode):factor(SamplingPeriodNumber)))


###############################################################################
### TESTING ONLY
growth <- filter(growth, sitecode %in% c("VB", "CAX"))
###############################################################################

# n_burnin <- 20000
# n_chains <- 8
# n_iter <- 100000
# n_thin <- 10
n_burnin <- 200
n_chains <- 1
n_iter <- 1000
n_thin <- 10

###############################################################################
#  JAGS model based on Clark, 2007
# See http://bit.ly/1qFRrJb for an example
#
# See http://bit.ly/1pQc5B4 for more on fitting the random slopes and 
# intercepts variance-covariance matrix in JAGS

model_file <- "2_model_growth_v_climate_jags_model_clark.bug"

# Setup variables in JAGS model
ID_site <- as.integer(factor(growth$sitecode))
ID_plot <- as.integer(factor(growth$plot_ID))
ID_tree <- as.integer(factor(growth$SamplingUnitName))
ID_period <- as.integer(factor(growth$SamplingPeriodID))
dbh_st <- growth$diameter_start
ln_dbh_st <- log(growth$diameter_start)
ln_dbh_end <- log(growth$diameter_end)
WD <- growth$WD - mean(growth$WD)
spi_24 <- growth$spi_24
spi_24_sq <- growth$spi_24^2
n_obs <- length(ID_site)
n_site <- length(unique(ID_site))
n_plot <- length(unique(ID_plot))
n_tree <- length(unique(ID_tree))
n_period <- length(unique(ID_period))
# Vector of zeros (see http://bit.ly/1uFTjlx)
zeros <- rep(0, n_obs)

# Setup latent diameter and true growth rate inits
growth <- group_by(growth, SamplingUnitName) %>%
    arrange(SamplingPeriodNumber) %>%
    mutate(ln_dbh_st_latent=log(diameter_start),  
           ln_dbh_end_latent=log(ifelse(diameter_end > diameter_start, diameter_end, diameter_start*1.001)))

jags_inits <- list(list(ln_dbh_st_latent=growth$ln_dbh_st_latent,
                        ln_dbh_end_latent=growth$ln_dbh_end_latent))
#jags_inits <- NULL

jags_inits <- rep(jags_inits, n_chains)

jags_vars <- list("ID_site", "ID_plot", "ID_tree", "ID_period", "ln_dbh_st", 
                  "ln_dbh_end", "WD", "spi_24", "spi_24_sq", "n_obs", "n_site", 
                  "n_plot", "n_tree", "n_period", "zeros")

jags_params <- c("int", "slp_dbh", "slp_WD", "slp_spi", "slp_spi_sq", 
                 "sigma.dia_inc", "sigma.p", "sigma.u", "sigma.v1", "sigma.v2", 
                 "rho", "sigma.w", "sigma.z")

jags_m <- jags(data=jags_vars, inits=jags_inits, parameters.to.save=jags_params, 
               n.chains=n_chains, n.burnin=n_burnin, n.iter=n_iter, 
               n.thin=n_thin, model.file=model_file)

log(exp(growth$ln_dbh_end_latent[25]) - exp(growth$ln_dbh_st_latent[25]))

jags_m <- jags.parallel(data=jags_vars, inits=NULL, 
                        parameters.to.save=jags_params, model.file=model_file, 
                        n.chains=n_chains, n.burnin=n_burnin, n.iter=n_iter, 
                        n.thin=n_thin)

# TODO: Include continent factors?
# TODO: see coda.samples

summary(jags_m)
plot(jags_m)
traceplot(jags_m)

jags_m_mcmc <- as.mcmc(jags_m)
summary(jags_m_mcmc)
xyplot(jags_m_mcmc)
densityplot(jags_m_mcmc)

s <- coda.samples(jags_m, jags_vars, n.iter=n_iter-n_burnin)
s <- ggmcmc(jags_m)
