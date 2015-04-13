# Forest growth model MCMC results

prefixes <- c('D:/azvoleff/Data', # CI-TEAM
              'H:/Data', # Buffalo drive
              'O:/Data', # Blue drive
              '/localdisk/home/azvoleff/Data') # vertica1
prefix <- prefixes[match(TRUE, unlist(lapply(prefixes, function(x) file_test('-d', x))))]

library(mcmcplots)
library(runjags)
library(ggmcmc)
library(ggthemes)
library(coda)

start_val <- 30000
thin_val <- 10

plot_mcmc <- function(fit, params, title) {
    mcmc_results <- as.mcmc.list(fit, c(params))
    mcmc_results <- window(mcmc_results, start=start_val, thin=thin_val)
    mcmcplot(mcmc_results, title=paste0(title))
    #caterplot(mcmc_results)
}

###############################################################################
## Simple model (no random effects)
##

plot_simple_model <- function(title) {
    plot_mcmc(jags_fit, c("B", "mu_B_g", "sigma_", "int_jk", "ink_k"), title)
}

### Minimum temperature
load(file.path(prefix, "TEAM", "Tree_Growth", "MCMC_Chains", "simple", 
               "simple_tmn.RData"))
plot_simple_model('simple min temp')

### Mean temperature
load(file.path(prefix, "TEAM", "Tree_Growth", "MCMC_Chains", "simple", 
               "simple_tmp.RData"))
plot_simple_model('simple mean temp')

### Maximum temperature
load(file.path(prefix, "TEAM", "Tree_Growth", "MCMC_Chains", "simple", 
               "simple_tmx.RData"))
plot_simple_model('simple max temp')

###############################################################################
## Model with only sites with large elevation gradients (UDZ, BIF, VB only)

# plot_small_model <- function(title) {
#     plot_mcmc(jags_fit, c("B[1]", "B[2]"), title)
#     plot_mcmc(jags_fit, "mu_B_g", title)
#     plot_mcmc(jags_fit, "sigma_", title)
# }
#
# ### Mean temperature
# load(file.path(prefix, "TEAM", "Tree_Growth", "MCMC_Chains", "highelev", 
#                "highelev_tmp.RData"))
# plot_small_model("small mean temp")
#
# ### Minimum temperature
# load(file.path(prefix, "TEAM", "Tree_Growth", "MCMC_Chains", "highelev", 
#                "highelev_tmn.RData"))
# plot_small_model("small min temp")
#
# ### Maximum temperature
# load(file.path(prefix, "TEAM", "Tree_Growth", "MCMC_Chains", "highelev", 
#                "highelev_tmx.RData"))
# plot_small_model("small max temp")

###############################################################################
## Full model (interactions, uncorrelated random effects)

plot_full_interact_model <- function(title) {
    plot_mcmc(jags_fit, c("B[1]", "B[2]", "mu_B_g", "sigma_", "B_T_"), title)
}

### Min temperature
load(file.path(prefix, "TEAM", "Tree_Growth", "MCMC_Chains", "interact", 
               "interact_tmn.RData"))
plot_full_interact_model("interact min temp")

### Mean temperature
load(file.path(prefix, "TEAM", "Tree_Growth", "MCMC_Chains", "interact", 
               "interact_tmp.RData"))
plot_full_interact_model("interact mean temp")

### Maximum temperature
load(file.path(prefix, "TEAM", "Tree_Growth", "MCMC_Chains", "interact", 
               "interact_tmx.RData"))
plot_full_interact_model("interact max temp")

# interact_labels_mu_B_g <- data.frame(Parameter=paste0('mu_B_g[', 1:12, ']'),
#     Label=c("int",
#             "P",
#             "P^2",
#             "T",
#             "T^2",
#             "D",
#             "D^2",
#             "W*P",
#             "D*P",
#             "W*T",
#             "D*T",
#             "E"))
#
# mu_B_g_vals <- ggs(as.mcmc.list(jags_fit, "mu_B_g"), par_labels=interact_labels_mu_B_g)
# ggs_caterpillar(mu_B_g_vals)
# ggs_caterpillar(mu_B_g_vals) + theme_tufte()
# ggs_caterpillar(mu_B_g_vals) + theme_solarized_2()

###############################################################################
## Full model (correlated random effects)
start_val <- 20000
thin_val <- 10
plot_full_correlated_model <- function(title) {
    plot_mcmc(jags_fit, c("B[1]", "B[2]", "mu_B_g", "sigma_", "B_T_"), title)
}

### Min temperature
load(file.path(prefix, "TEAM", "Tree_Growth", "MCMC_Chains", "correlated", 
               "correlated_tmn.RData"))
plot_full_correlated_model("correlated min temp")

### Mean temperature
load(file.path(prefix, "TEAM", "Tree_Growth", "MCMC_Chains", "correlated", 
               "correlated_tmp.RData"))
plot_full_correlated_model("correlated mean temp")

### Maximum temperature
load(file.path(prefix, "TEAM", "Tree_Growth", "MCMC_Chains", "correlated", 
               "correlated_tmx.RData"))
plot_full_correlated_model("correlated max temp")
