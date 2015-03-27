library(runjags)
library(coda)
library(mcgibbsit)

source("0_settings.R")

model_type <- "full"
#model_type <- "testing"

# model_structure <- "full_model"
#model_structure <- "full_model_no_t_effects"
model_structure <- "full_model_no_t_effects_interact"

precip_var <- "mcwd_run12"

note <- 'genuslimits'

temp_var <- "tmn_meanannual"
run_id <- "vertica1.team.sdsc.edu_20150327024228"

# temp_var <- "tmp_meanannual"
# run_id <- "vertica1.team.sdsc.edu_20150327044629"

# temp_var <- "tmx_meanannual"
# run_id <- "vertica1.team.sdsc.edu_20150326223806"

suffix <- paste0(model_type, '-', temp_var, '-', precip_var, '_', note, '-', run_id)

mcmc_folder <- file.path(prefix, "TEAM", "Tree_Growth", "MCMC_Chains")
params_folder <- file.path(prefix, "TEAM", "Tree_Growth", "Extracted_Parameters")

load(file.path(mcmc_folder, paste0("jags_fit_", suffix, ".RData")))

#mcgibbsit(as.mcmc.list(jags_fit))

fixefs <- as.mcmc.list(jags_fit, c("B[1]", "B[2]"))
save(fixefs, file=file.path(params_folder, paste0(suffix, "_fixefs.RData")))

ranefs <- as.mcmc.list(jags_fit, c("^int"))
save(ranefs, file=file.path(params_folder, paste0(suffix, "_ranefs.RData")))

ranefs_sigmas <- as.mcmc.list(jags_fit, c("sigma_obs", "sigma_proc", "sigma_int"))
save(ranefs_sigmas , file=file.path(params_folder, paste0(suffix, "_ranefs_sigmas.RData")))

ranefs_B_T <- as.mcmc.list(jags_fit, c("^B_T"))
save(ranefs_B_T, file=file.path(params_folder, paste0(suffix, "_ranefs_B_T.RData")))

ranefs_mu_B_g <- as.mcmc.list(jags_fit, c("mu_B_g"))
save(ranefs_mu_B_g, file=file.path(params_folder, paste0(suffix, "_ranefs_mu_B_g.RData")))

ranefs_all_vars <- as.mcmc.list(jags_fit, c("mu_B_g", "B[1]", "B[2]"))
save(ranefs_all_vars, file=file.path(params_folder, paste0(suffix, "_ranefs_all_vars.RData")))

ranefs_B_g <- as.mcmc.list(jags_fit, c("^B_g"))
save(ranefs_B_g, file=file.path(params_folder, paste0(suffix, "_ranefs_B_g.RData")))

ranefs_sigma_B_g <- as.mcmc.list(jags_fit, c("sigma_B_g"))
save(ranefs_sigma_B_g, file=file.path(params_folder, paste0(suffix, "_ranefs_sigma_B_g.RData")))

if (model_structure != "full_model_no_t_effects_interact") {
    ranefs_rho_B_g <- as.mcmc.list(jags_fit, c("rho_B_g"))
    save(ranefs_rho_B_g, file=file.path(params_folder, paste0(suffix, "_ranefs_rho_B_g.RData")))
}

plot(fixefs, ask=TRUE)
gelman.diag(fixefs)
mcgibbsit(fixefs)

plot(ranefs, ask=TRUE)
gelman.diag(ranefs)
mcgibbsit(ranefs)

plot(ranefs_mu_B_g, ask=TRUE)
gelman.diag(ranefs_mu_B_g)
mcgibbsit(ranefs_mu_B_g)

plot(ranefs_B_T, ask=TRUE)
gelman.diag(ranefs_B_T)
mcgibbsit(ranefs_B_T)

plot(ranefs_sigmas, ask=TRUE)
gelman.diag(ranefs_sigmas)
mcgibbsit(ranefs_sigmas)

plot(ranefs_sigma_B_g, ask=TRUE)
gelman.diag(ranefs_sigma_B_g)
mcgibbsit(ranefs_sigma_B_g)

plot(ranefs_B_g, ask=TRUE)
gelman.diag(ranefs_B_g)

# Subset of the correlation matrix that is stochastic:
ranef_rho_B_g_names <- c("rho_B_g[1,2]",
                         "rho_B_g[1,3]",
                         "rho_B_g[1,4]",
                         "rho_B_g[1,5]",
                         "rho_B_g[1,6]",
                         "rho_B_g[1,7]",
                         "rho_B_g[2,3]",
                         "rho_B_g[2,4]",
                         "rho_B_g[2,5]",
                         "rho_B_g[2,6]",
                         "rho_B_g[2,5]",
                         "rho_B_g[3,4]",
                         "rho_B_g[3,5]",
                         "rho_B_g[3,6]",
                         "rho_B_g[3,7]",
                         "rho_B_g[4,5]",
                         "rho_B_g[4,6]",
                         "rho_B_g[4,7]",
                         "rho_B_g[5,6]",
                         "rho_B_g[5,7]",
                         "rho_B_g[6,7]")
plot(as.mcmc.list(jags_fit, ranef_rho_B_g_names), ask=TRUE)
gelman.diag(as.mcmc.list(jags_fit, ranef_rho_B_g_names))
mcgibbsit(as.mcmc.list(jags_fit, ranef_rho_B_g_names))
