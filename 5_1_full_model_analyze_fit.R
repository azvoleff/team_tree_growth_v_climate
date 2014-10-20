library(runjags)
library(coda)
library(dplyr)
library(reshape2)
library(ggplot2)
library(corrplot)

img_dpi <- 300

load("model_data_standardizing.RData")

plot_estimates <- function(mcmc_ests, pars, pars_pretty=NULL, xmin=NULL, xmax=NULL, 
                           scaling=NULL) {
    stopifnot(length(pars) == length(pars_pretty))
    suppressMessages(mcmc_ests <- melt(mcmc_ests))
    # Melt will convert the variable names to be syntactically valid R variable 
    # names. So do the same to the names supplied in pars
    pars <- make.names(pars)
    # Filter to only include in plot variables that are in pars
    stopifnot(all(pars %in% mcmc_ests$variable))
    mcmc_ests <- mcmc_ests[mcmc_ests$variable %in% pars, ]
    if (!is.null(scaling)) {
        stopifnot(length(pars) == length(scaling))
        for (i in 1:length(pars)) {
            these_rows <- mcmc_ests$variable == pars[i]
            mcmc_ests[these_rows, ]$value <- mcmc_ests[these_rows, ]$value * scaling[i]
        }
    }
    mcmc_ests <- group_by(mcmc_ests, variable) %>%
        summarize(median=median(value),
                  q2pt5=quantile(value, .025),
                  q97pt5=quantile(value, .975))

    if (is.null(pars_pretty)) {
        mcmc_ests$variable <- ordered(mcmc_ests$variable, 
                                      levels=rev(pars))
    } else {
        mcmc_ests$variable <- ordered(mcmc_ests$variable, 
                                      levels=rev(pars),
                                      labels=rev(pars_pretty))
    }

    if (is.null(xmin)) xmin <- min(mcmc_ests$q2pt5)
    if (is.null(xmax)) xmax <- max(mcmc_ests$q97pt5)
    p <- ggplot(mcmc_ests, aes(median, variable)) +
        theme_bw(base_size=10) +
        geom_point() + 
        geom_errorbarh(aes(xmin=q2pt5, xmax=q97pt5), height=0.5) +
        geom_vline(aes(xintercept=0), color="grey") +
        xlim(c(xmin, xmax)) +
        xlab("") + ylab("") +
        scale_y_discrete(breaks=pars_pretty, labels=pars_pretty)
    return(p)
}

###############################################################################
### Plot fixed effects
###############################################################################

load("jags_fit_full_model_fixefs.RData")

fixefs_names <- c("int",
                 "slp_mcwd",
                 "slp_mcwd_sq",
                 "slp_dbh",
                 "slp_dbh_sq",
                 "slp_WD",
                 "slp_WD_sq")
fixefs_names_pretty <- c("Intercept",
                        "MCWD",
                        expression(MCWD^2),
                        "Diameter",
                        expression(Diameter^2),
                        "Density",
                        expression(Density^2))
fixefs_scaling <- c(dbh_sd,
                    (dbh_sd/mcwd_sd) * 100, # Convert from mm to 10s of cm
                    (dbh_sd/mcwd_sd) * 100, # Convert from mm to 10s of cm
                    dbh_sd/dbh_sd,
                    dbh_sd/dbh_sd,
                    dbh_sd/WD_sd,
                    dbh_sd/WD_sd)

# gelman.diag(fixefs)
# fixefs <- window(fixefs, start=9001, end=15000, thin=50)
plot(fixefs, ask=TRUE)
#autocorr.plot(fixefs)
fixefs_comb <- data.frame(combine.mcmc(fixefs))

plot_estimates(fixefs_comb, fixefs_names, fixefs_names_pretty,
               scaling=fixefs_scaling)
ggsave("growth_model_fixefs.png", width=3, height=4, dpi=img_dpi)

###############################################################################
### Plot random intercepts and process and observation error
###############################################################################

load("jags_fit_full_model_ranefs.RData")

ranef_names <- c("sigma_obs",
                 "sigma_proc",
                 "sigma_int_ijk",
                 "sigma_int_jk",
                 "sigma_int_k",
                 "sigma_int_t")
ranef_names_pretty <- c(expression(sigma[obs]),
                        expression(sigma[proc]),
                        expression(sigma[tree]),
                        expression(sigma[plot]),
                        expression(sigma[site]),
                        expression(sigma[period]))

# gelman.diag(ranefs)
# ranefs <- window(ranefs, start=9001, end=15000, thin=50)
#plot(ranefs, ask=TRUE)
#autocorr.plot(ranefs)
ranefs_comb <- data.frame(combine.mcmc(ranefs))

plot_estimates(ranefs_comb, ranef_names, ranef_names_pretty, xmin=0, xmax=.03)
ggsave("growth_model_ranefs.png", width=3, height=5, dpi=img_dpi)

###############################################################################
### Plot genus-level random effects variances
###############################################################################

load("jags_fit_full_model_ranefs_g_sigma.RData")

ranefs_g_sigma_names <- c("sigma_B_g[1]",
                          "sigma_B_g[2]",
                          "sigma_B_g[3]",
                          "sigma_B_g[4]",
                          "sigma_B_g[5]")
ranefs_g_sigma_names_pretty <- c(expression(sigma[int,g]),
                                 expression(sigma[MCWD,g]),
                                 expression(sigma[MCWD^2,g]),
                                 expression(sigma[DBH,g]),
                                 expression(sigma[DBH^2,g]))

# gelman.diag(ranefs_g_sigma)
# ranefs_g_sigma <- window(ranefs_g_sigma, start=9001, end=15000, thin=50)
#plot(ranefs_g_sigma, ask=TRUE)
#autocorr.plot(ranefs_g_sigma)
ranefs_g_sigma_comb <- data.frame(combine.mcmc(ranefs_g_sigma))

plot_estimates(ranefs_g_sigma_comb, ranefs_g_sigma_names,
               ranefs_g_sigma_names_pretty, xmin=0, xmax=.05)
ggsave("growth_model_ranefs_g_sigma.png", width=3, height=5, dpi=img_dpi)

###############################################################################
### Plot correlation matrix for genus-level random effects
###############################################################################

load("jags_fit_full_model_ranefs_g_rho.RData")

ranefs_g_rho_names <- c("rho_B_g[1,2]",
                        "rho_B_g[1,3]",
                        "rho_B_g[1,4]",
                        "rho_B_g[1,5]",
                        "rho_B_g[2,3]",
                        "rho_B_g[2,4]",
                        "rho_B_g[2,5]",
                        "rho_B_g[3,4]",
                        "rho_B_g[3,5]",
                        "rho_B_g[4,5]")
ranefs_g_rho_names_pretty <- c(expression(rho[list(int,MCWD)]),
                        expression(rho[list(int,MCWD^2)]),
                        expression(rho[list(int,DBH[t-1])]),
                        expression(rho[list(int,DBH[t-1]^2)]),
                        expression(rho[list(MCWD,MCWD^2)]),
                        expression(rho[list(MCWD,DBH[t-1])]),
                        expression(rho[list(MCWD,DBH[t-1]^2)]),
                        expression(rho[list(MCWD^2,DBH[t-1])]),
                        expression(rho[list(MCWD^2,DBH[t-1]^2)]),
                        expression(rho[list(DBH[t-1],DBH[t-1]^2)]))

# gelman.diag(ranefs_g_rho)
# ranefs_g_rho <- window(ranefs_g_rho, start=9001, end=15000, thin=50)
#plot(ranefs_g_rho, ask=TRUE)
#autocorr.plot(ranefs_g_rho)
ranefs_g_rho_comb <- data.frame(combine.mcmc(ranefs_g_rho))

plot_estimates(ranefs_g_rho_comb, ranefs_g_rho_names,
               ranefs_g_rho_names_pretty, xmin=-1, xmax=1)
ggsave("growth_model_ranefs_g_rho.png", width=3, height=6, dpi=img_dpi)

corrplot
