library(R2jags)
library(coda)

load("jags_fit_full.RData")

trace_vars <- c("intercept",
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
                "sigma_t",
                "sigma_g",
                "sigma_slp_spi_g",
                "rho")

traceplot(jags_fit, trace_vars)

jags_fit_mcmc <- as.mcmc(jags_fit)

# What Gelman-Rubin to be about 1. Above 1.5 indicates lack of convergence.
gelman.diag(jags_fit_mcmc)
