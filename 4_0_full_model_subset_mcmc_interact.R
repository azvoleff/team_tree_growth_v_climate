library(runjags)
library(coda)
library(mcgibbsit)

source("0_settings.R")

model_type <- "full_no_t_effects_interact"
#model_type <- "testing"
precip_var <- "mcwd_run12"

temp_var <- "tmn_meanannual"
run_ID <- "vertica1.team.sdsc.edu_20150104150226"

# temp_var <- "tmp_meanannual"
# run_ID <- "vertica1.team.sdsc.edu_20150104123029"

suffix <- paste0(model_type, '-', temp_var, '-', precip_var)

mcmc_folder <- file.path(prefix, "TEAM", "Tree_Growth", "MCMC_Chains")
params_folder <- file.path(prefix, "TEAM", "Tree_Growth", "Extracted_Parameters")

load(file.path(mcmc_folder, paste0("jags_fit_", suffix, "-", run_ID, ".RData")))

#mcgibbsit(as.mcmc.list(jags_fit))

fixefs <- as.mcmc.list(jags_fit, c("B[1]", "B[2]"))
save(fixefs, file=file.path(params_folder, paste0(suffix, "jags_fit_full_model_fixefs.RData")))

ranefs <- as.mcmc.list(jags_fit, c("^int"))
save(ranefs, file=file.path(params_folder, paste0(suffix, "jags_fit_full_model_ranefs.RData")))

ranefs_sigmas <- as.mcmc.list(jags_fit, c("sigma_obs", "sigma_proc", "sigma_int"))
save(ranefs_sigmas , file=file.path(params_folder, paste0(suffix, "jags_fit_full_model_ranefs_sigmas.RData")))

ranefs_B_T <- as.mcmc.list(jags_fit, c("^B_T"))
save(ranefs_B_T, file=file.path(params_folder, paste0(suffix, "jags_fit_full_model_ranefs_B_T.RData")))

all_efs <- as.mcmc.list(jags_fit, c("B[1]", "B[2]", "mu_B_g"))
save(all_efs, file=file.path(params_folder, paste0(suffix, "jags_fit_full_model_all_efs.RData")))

ranefs_mu_B_g <- as.mcmc.list(jags_fit, c("mu_B_g"))
save(ranefs_mu_B_g, file=file.path(params_folder, paste0(suffix, "jags_fit_full_model_ranefs_mu_B_g.RData")))

ranefs_B_g <- as.mcmc.list(jags_fit, c("^B_g"))
save(ranefs_B_g, file=file.path(params_folder, paste0(suffix, "jags_fit_full_model_ranefs_B_g.RData")))

ranefs_sigma_B_g <- as.mcmc.list(jags_fit, c("sigma_B_g"))
save(ranefs_sigma_B_g, file=file.path(params_folder, paste0(suffix, "jags_fit_full_model_ranefs_sigma_B_g.RData")))

plot(all_efs, ask=TRUE)
gelman.diag(ranefs)
mcgibbsit(ranefs)

plot(ranefs_B_T, ask=TRUE)
gelman.diag(ranefs_B_T)
mcgibbsit(ranefs_B_T)

plot(ranefs_sigmas, ask=TRUE)
gelman.diag(ranefs_sigmas)
mcgibbsit(ranefs_sigmas)

plot(ranefs_sigma_B_g, ask=TRUE)
gelman.diag(ranefs_sigma_B_g)
mcgibbsit(ranefs_sigma_B_g)

## Warning - below takes a long time to plot
# plot(ranefs_B_g, ask=TRUE)
# gelman.diag(ranefs_B_g)
# mcgibbsit(ranefs_B_g)
