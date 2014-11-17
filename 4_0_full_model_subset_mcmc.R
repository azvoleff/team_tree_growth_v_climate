library(runjags)
library(coda)
library(mcgibbsit)

model_type <- "full"
#model_type <- "testing"
precip_var <- "mcwd_run12"

temp_var <- "tmx_meanannual"
run_ID <- "vertica1.team.sdsc.edu_20141110224426_extend3" 

# temp_var <- "tmp_meanannual"
# run_ID <- "vertica1.team.sdsc.edu_20141110152032_extend1"

# temp_var <- "tmn_meanannual"
# run_ID <- "vertica1.team.sdsc.edu_20141115021513"

in_folder <- 'MCMC_Chains'
suffix <- paste0(model_type, '-', temp_var, '-', precip_var)
load(file.path(in_folder, paste0("jags_fit_", suffix, "-", run_ID, ".RData")))

#mcgibbsit(as.mcmc.list(jags_fit))

fixefs <- as.mcmc.list(jags_fit, c("B[1]", "B[2]"))
# plot(fixefs, ask=TRUE)
# gelman.diag(fixefs)
# mcgibbsit(fixefs)
save(fixefs, file=file.path(in_folder, paste0(suffix, "jags_fit_full_model_fixefs.RData")))

ranefs <- as.mcmc.list(jags_fit, c("^int"))
# plot(ranefs, ask=TRUE)
# gelman.diag(ranefs)
save(ranefs, file=file.path(in_folder, paste0(suffix, "jags_fit_full_model_ranefs.RData")))

ranefs_sigmas <- as.mcmc.list(jags_fit, c("sigma_obs", "sigma_proc", "sigma_int"))
# plot(ranefs_sigmas, ask=TRUE)
# gelman.diag(ranefs_sigmas)
save(ranefs_sigmas , file=file.path(in_folder, paste0(suffix, "jags_fit_full_model_ranefs_sigmas.RData")))

ranefs_B_T <- as.mcmc.list(jags_fit, c("^B_T"))
# plot(ranefs_B_T, ask=TRUE)
# gelman.diag(ranefs_B_T)
# mcgibbsit(ranefs_B_T)
save(ranefs_B_T, file=file.path(in_folder, paste0(suffix, "jags_fit_full_model_ranefs_B_T.RData")))

ranefs_mu_B_g <- as.mcmc.list(jags_fit, c("mu_B_g"))
# plot(ranefs_mu_B_g, ask=TRUE)
# gelman.diag(ranefs_mu_B_g)
save(ranefs_mu_B_g, file=file.path(in_folder, paste0(suffix, "jags_fit_full_model_ranefs_mu_B_g.RData")))

ranefs_B_g <- as.mcmc.list(jags_fit, c("^B_g"))
# plot(ranefs_B_g, ask=TRUE)
# gelman.diag(ranefs_B_g)
save(ranefs_B_g, file=file.path(in_folder, paste0(suffix, "jags_fit_full_model_ranefs_B_g.RData")))

ranefs_sigma_B_g <- as.mcmc.list(jags_fit, c("sigma_B_g"))
# plot(ranefs_sigma_B_g, ask=TRUE)
# gelman.diag(ranefs_sigma_B_g)
save(ranefs_sigma_B_g, file=file.path(in_folder, paste0(suffix, "jags_fit_full_model_ranefs_sigma_B_g.RData")))

ranefs_rho_B_g <- as.mcmc.list(jags_fit, c("rho_B_g"))
# # Below list is of the subset of the correlation matrix that is stochastic
# ranef_rho_B_g_names <- c("rho_B_g[1,2]",
#                          "rho_B_g[1,3]",
#                          "rho_B_g[1,4]",
#                          "rho_B_g[1,5]",
#                          "rho_B_g[1,6]",
#                          "rho_B_g[1,7]",
#                          "rho_B_g[2,3]",
#                          "rho_B_g[2,4]",
#                          "rho_B_g[2,5]",
#                          "rho_B_g[2,6]",
#                          "rho_B_g[2,5]",
#                          "rho_B_g[3,4]",
#                          "rho_B_g[3,5]",
#                          "rho_B_g[3,6]",
#                          "rho_B_g[3,7]",
#                          "rho_B_g[4,5]",
#                          "rho_B_g[4,6]",
#                          "rho_B_g[4,7]",
#                          "rho_B_g[5,6]",
#                          "rho_B_g[5,7]",
#                          "rho_B_g[6,7]")
# plot(as.mcmc.list(jags_fit, ranef_rho_B_g_names), ask=TRUE)
# gelman.diag(as.mcmc.list(jags_fit, ranef_rho_B_g_names))
save(ranefs_rho_B_g, file=file.path(in_folder, paste0(suffix, "jags_fit_full_model_ranefs_rho_B_g.RData")))
