library(mcmcplots)
library(runjags)
library(foreach)
library(doParallel)

prefixes <- c('D:/azvoleff/Data', # CI-TEAM
              'H:/Data', # Buffalo drive
              'O:/Data', # Blue drive
              '/localdisk/home/azvoleff/Data') # vertica1
prefix <- prefixes[match(TRUE, unlist(lapply(prefixes, function(x) file_test('-d', x))))]

data_folder <- file.path(prefix, "TEAM", "Tree_Growth", "Data")

cl <- makeCluster(3)
registerDoParallel(cl)

plot_mcmc <- function(fit, params, model_type, temp_var) {
    mcmc_results <- as.mcmc.list(fit, c(params))
    mcmc_results <- window(mcmc_results, start=start_val, thin=thin_val)
    if (!file_test('-d', model_type)) dir.create(model_type)
    mcmcplot(mcmc_results, title=paste0(model_type, ' - ', temp_var), 
             dir=model_type, filename=temp_var)
}

###############################################################################
## Simple model (no random effects)
start_val <- 40000
thin_val <- 50
foreach(temp_var=c('tmn', 'tmp', 'tmx'), .inorder=FALSE,
        .packages=c('runjags', 'mcmcplots')) %dopar% {
    load(file.path(prefix, "TEAM", "Tree_Growth", "MCMC_Chains", "simple", 
                   paste0('simple_', temp_var, '.RData')))
    plot_mcmc(jags_fit, c("B", "sigma_", "int_", "B_k"), "simple", temp_var)
}

###############################################################################
## Full model (interactions, uncorrelated random effects)
start_val <- 40000
thin_val <- 50
foreach(temp_var=c('tmn', 'tmp', 'tmx'), .inorder=FALSE,
        .packages=c('runjags', 'mcmcplots')) %dopar% {
    load(file.path(prefix, "TEAM", "Tree_Growth", "MCMC_Chains", "interact", 
                   paste0('interact_', temp_var, '.RData')))
    plot_mcmc(jags_fit, c("mu_B_g", "sigma_", "int_", "B_k"), 'interact', 
              temp_var)
}

###############################################################################
## Full model (correlated random effects)
start_val <- 20000
thin_val <- 50
foreach(temp_var=c('tmn', 'tmp', 'tmx'), .inorder=FALSE,
        .packages=c('runjags', 'mcmcplots')) %dopar% {
    load(file.path(prefix, "TEAM", "Tree_Growth", "MCMC_Chains", "correlated", 
                   paste0('correlated_', temp_var, '.RData')))
    plot_mcmc(jags_fit, c("mu_B_g", "sigma_", "int_", "B_k", "rho_B_g"), 
              "correlated", temp_var)
}
