library(runjags)
library(coda)
library(mcgibbsit)

# load("MCMC_Chains/full_model_fit_parallel_vertica1.team.sdsc.edu_20141025-173906_extended_1.RData")
# load("MCMC_Chains/full_model_fit_parallel_vertica1.team.sdsc.edu_20141025-173906_extended_2_with_bg.RData")
load("MCMC_Chains/full_model_fit_parallel_vertica1.team.sdsc.edu_20141025-173906_extended_3_with_bg_and_ints.RData")

# mcgibbsit(as.mcmc.list(jags_fit))

fixefs <- as.mcmc.list(jags_fit, c("B[1]", "B[2]"))
# plot(fixefs, ask=TRUE)
# gelman.diag(fixefs)
save(fixefs, file="jags_fit_full_model_fixefs.RData")

ranefs <- as.mcmc.list(jags_fit, c("^int"))
# plot(ranefs, ask=TRUE)
# gelman.diag(ranefs)
save(ranefs , file="jags_fit_full_model_ranefs.RData")

ranefs_sigmas <- as.mcmc.list(jags_fit, c("sigma_obs", "sigma_proc", "sigma_int"))
# plot(ranefs_sigmas, ask=TRUE)
# gelman.diag(ranefs_sigmas)
save(ranefs_sigmas , file="jags_fit_full_model_ranefs_sigmas.RData")

ranefs <- as.mcmc.list(jags_fit, c("sigma_obs", "sigma_proc", "sigma_int"))
# plot(ranefs, ask=TRUE)
# gelman.diag(ranefs_sigmas)
save(ranefs, file="jags_fit_full_model_ranefs_sigmas.RData")

ranefs_B_g <- as.mcmc.list(jags_fit, c("^B_g"))
# plot(ranefs_B_g, ask=TRUE)
# gelman.diag(ranefs_B_g)
save(ranefs_B_g, file="jags_fit_full_model_ranefs_B_g.RData")

ranefs_mu_B_g <- as.mcmc.list(jags_fit, c("mu_B_g"))
# plot(ranefs_mu_B_g, ask=TRUE)
# gelman.diag(ranefs_mu_B_g)
save(ranefs_mu_B_g, file="jags_fit_full_model_ranefs_mu_B_g.RData")

ranefs_sigma_B_g <- as.mcmc.list(jags_fit, c("sigma_B_g"))
# plot(ranefs_sigma_B_g, ask=TRUE)
# gelman.diag(ranefs_sigma_B_g)
save(ranefs_sigma_B_g, file="jags_fit_full_model_ranefs_sigma_B_g.RData")

# Below species the subset of the correlation matrix that is stochastic
raneg_rho_B_g_names <- c("rho_B_g[1,2]",
                         "rho_B_g[1,3]",
                         "rho_B_g[1,4]",
                         "rho_B_g[1,5]",
                         "rho_B_g[2,3]",
                         "rho_B_g[2,4]",
                         "rho_B_g[2,5]",
                         "rho_B_g[3,4]",
                         "rho_B_g[3,5]",
                         "rho_B_g[4,5]")
ranefs_rho_B_g <- as.mcmc.list(jags_fit, c("rho_B_g"))
# plot(as.mcmc.list(jags_fit, raneg_rho_B_g_names), ask=TRUE)
# gelman.diag(as.mcmc.list(jags_fit, raneg_rho_B_g_names))
save(ranefs_rho_B_g, file="jags_fit_full_model_ranefs_rho_B_g.RData")
