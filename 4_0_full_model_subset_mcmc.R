library(runjags)
library(coda)
library(mcgibbsit)

load("full_model_fit_parallel_vertica1.team.sdsc.edu_20141025-173906_extended.RData")

# mcgibbsit(as.mcmc.list(jags_fit))

fixefs <- as.mcmc.list(jags_fit, c("^B"))
save(fixefs, file="jags_fit_full_model_fixefs.RData")
#plot(fixefs, ask=TRUE)

ranefs <- as.mcmc.list(jags_fit, c("sigma_obs", "sigma_proc", "sigma_int"))
save(ranefs , file="jags_fit_full_model_ranefs.RData")
#plot(ranefs, ask=TRUE)

ranefs_mu_B_g <- as.mcmc.list(jags_fit, c("mu_B_g"))
save(ranefs_mu_B_g, file="jags_fit_full_model_ranefs_mu_B_g.RData")
#plot(ranefs_mu_B_g, ask=TRUE)

ranefs_sigma_B_g <- as.mcmc.list(jags_fit, c("sigma_B_g"))
save(ranefs_sigma_B_g, file="jags_fit_full_model_ranefs_sigma_B_g.RData")
#plot(ranefs_sigma_B_g, ask=TRUE)

ranefs_rho_B_g <- as.mcmc.list(jags_fit, c("rho_B_g"))
save(ranefs_rho_B_g, file="jags_fit_full_model_ranefs_rho_B_g.RData")
#plot(ranefs_rho_B_g, ask=TRUE)

# # Calculate growth increment
# dbh_preds <- extend.jags(jags_fit_p,
#                          drop.monitor=c("b_k", "b_t", "b_h", "slp_mcwd_g"), 
#                          add.monitor="dbh_latent", sample=1000, thin=5)
