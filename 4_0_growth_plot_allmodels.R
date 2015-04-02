# Forest growth model MCMC results

prefixes <- c('D:/azvoleff/Data', # CI-TEAM
              'H:/Data', # Buffalo drive
              'O:/Data', # Blue drive
              '/localdisk/home/azvoleff/Data') # vertica1
prefix <- prefixes[match(TRUE, unlist(lapply(prefixes, function(x) file_test('-d', x))))]

library(mcmcplots)
library(runjags)
library(coda)

start_val <- 10000
thin_val <- 8

plot_mcmc <- function(fit, params, title) {
    mcmc_results <- as.mcmc.list(fit, c(params))
    mcmc_results <- window(mcmc_results, start=start_val, thin=thin_val)
    mcmcplot(mcmc_results, title=paste0(title, ": ", paste(params, collapse=", ")))
    #caterplot(mcmc_results)
}

###############################################################################
## Simple model (no random effects)
##


    mcmc_results <- as.mcmc.list(jags_fit, c('sigma_'))
    mcmc_results <- window(mcmc_results, start=start_val, thin=thin_val)
    plot(mcmc_results, ask=TRUE)




plot_simple_model <- function(title) {
    plot_mcmc(jags_fit, "B[", title)
    plot_mcmc(jags_fit, "sigma_", title)
    #plot_mcmc(jags_fit, "^int_k", title)
}

### Mean temperature
load(file.path(prefix, "TEAM", "Tree_Growth", "MCMC_Chains", "simple", 
               "simple_tmp.RData"))
plot_simple_model('simple mean temp')

### Minimum temperature
load(file.path(prefix, "TEAM", "Tree_Growth", "MCMC_Chains", "simple", 
               "simple_tmn.RData"))
plot_simple_model('simple min temp')

### Maximum temperature
load(file.path(prefix, "TEAM", "Tree_Growth", "MCMC_Chains", "simple", 
               "simple_tmx.RData"))
plot_simple_model('simple max temp')

###############################################################################
## Smaller model (UDZ, BIF, VB only)

plot_small_model <- function(title) {
    plot_mcmc(jags_fit, c("B[1]", "B[2]"), title)
    plot_mcmc(jags_fit, "mu_B_g", title)
    plot_mcmc(jags_fit, "sigma_", title)
}

### Mean temperature
load(file.path(prefix, "TEAM", "Tree_Growth", "MCMC_Chains", "highelev", 
               "highelev_tmp.RData"))
plot_small_model("small mean temp")

### Minimum temperature
load(file.path(prefix, "TEAM", "Tree_Growth", "MCMC_Chains", "highelev", 
               "highelev_tmn.RData"))
plot_small_model("small min temp")

### Maximum temperature
load(file.path(prefix, "TEAM", "Tree_Growth", "MCMC_Chains", "highelev", 
               "highelev_tmx.RData"))
plot_small_model("small max temp")

###############################################################################
## Full model (interactions, uncorrelated random effects)

plot_full_interact_model <- function() {
    plot_mcmc(jags_fit, c("B[1]", "B[2]"), title)
    plot_mcmc(jags_fit, "mu_B_g", title)
    plot_mcmc(jags_fit, "sigma_", title)
}

### Mean temperature
load(file.path(prefix, "TEAM", "Tree_Growth", "MCMC_Chains", "interact", 
               "interact_tmx.RData"))
plot_full_interact_model("interact mean temp")

### Minimum temperature
load(file.path(prefix, "TEAM", "Tree_Growth", "MCMC_Chains", "interact", 
               "interact_tmx.RData"))
plot_full_interact_model("interact min temp")

### Maximum temperature
load(file.path(prefix, "TEAM", "Tree_Growth", "MCMC_Chains", "interact", 
               "interact_tmx.RData"))
plot_full_interact_model("interact max temp")

###############################################################################
## Full model (correlated random effects)

plot_full_correlated_model <- function() {
}

### Minimum temperature

### Mean temperature

### Maximum temperature

