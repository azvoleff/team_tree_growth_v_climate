library(runjags)
library(coda)

load("jags_fit_ranslope_parallel.RData")

fixef_names <- c("intercept",
                 "slp_WD",
                 "slp_WD_sq",
                 "slp_dbh",
                 "slp_dbh_sq",
                 "'slp_spi'",
                 "inter_spi_WD",
                 "inter_spi_dbh")
fixefs <- as.mcmc.list(jags_fit_p, fixef_names)
save(fixefs, file="jags_fit_ranslope_parallel_fixefs.RData")

ranef_names <- c("obs_sigma",
                 "proc_sigma",
                 "sigma_ijk",
                 "sigma_jk",
                 "sigma_k",
                 "sigma_t",
                 "sigma_b_g",
                 "sigma_slp_spi_g",
                 "rho_g")
ranefs <- as.mcmc.list(jags_fit_p, ranef_names)
save(ranefs, file="jags_fit_ranslope_parallel_ranefs.RData")

genus_ranef_names <- c("slp_spi_g")
genus_ranefs <- as.mcmc.list(jags_fit_p, genus_ranef_names)
save(genus_ranefs, file="jags_fit_ranslope_parallel_genus_ranefs.RData")

# Calculate growth increment
dbh_preds <- extend.jags(jags_fit_p,
                         drop.monitor=c("b_k", "b_t", "b_h", "slp_spi_g"), 
                         add.monitor="dbh_latent", sample=1000, thin=5)

