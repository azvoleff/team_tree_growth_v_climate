library(runjags)
library(coda)

#load("full_model_fit_CC-4B-51-8F-CD_20141019-205136.RData")
load("full_model_fit_parallel_vertica1.team.sdsc.edu_20141021-014706.RData")
load("full_model_fit_parallel_vertica1.team.sdsc.edu_20141022-jags_extended.RData")

fixefs <- as.mcmc.list(jags_fit, c("^int", "slp_"))
save(fixefs, file="jags_fit_full_model_fixefs.RData")
plot(fixefs, ask=TRUE)

ranefs_g_sigma <- as.mcmc.list(jags_fit, c("sigma_B_g"))
save(ranefs_g_sigma, file="jags_fit_full_model_ranefs_g_sigma.RData")
plot(ranefs_g_sigma, ask=TRUE)

ranefs_g_rho <- as.mcmc.list(jags_fit, c("rho_B_g"))
save(ranefs_g_rho, file="jags_fit_full_model_ranefs_g_rho.RData")
plot(ranefs_g_rho, ask=TRUE)

ranefs <- as.mcmc.list(jags_fit, c("sigma_obs", "sigma_proc", "sigma_int"))
save(ranefs , file="jags_fit_full_model_ranefs.RData")
plot(ranefs, ask=TRUE)

# # Calculate growth increment
# dbh_preds <- extend.jags(jags_fit_p,
#                          drop.monitor=c("b_k", "b_t", "b_h", "slp_mcwd_g"), 
#                          add.monitor="dbh_latent", sample=1000, thin=5)
#
