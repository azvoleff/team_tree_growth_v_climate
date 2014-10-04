library(R2jags)
library(coda)

load("jags_fit_full.RData")

traceplot(jags_fit)

jags_fit_mcmc <- as.mcmc(jags_fit)

# What Gelman-Rubin to be about 1. Above 1.5 indicates lack of convergence.
gelman.diag(jags_fit_mcmc)
