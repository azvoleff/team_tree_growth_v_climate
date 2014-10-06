library(runjags)
library(coda)
library(dplyr)
library(reshape2)
library(ggplot2)

img_dpi <- 300

load("jags_fit_full_parallel.RData")
#load("jags_fit_ranslope.RData")
plot_estimates <- function(mcmc_ests, pars, pars_pretty=NULL, xmin=-1, xmax=1) {
    mcmc_ests <- melt(mcmc_ests)
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

load("jags_fit_ranslope_parallel_fixefs.RData")

fixef_names <- c("intercept",
                 "slp_WD",
                 "slp_WD_sq",
                 "slp_dbh",
                 "slp_dbh_sq",
                 "slp_spi",
                 "inter_spi_WD",
                 "inter_spi_dbh")
fixef_names_pretty <- c("Intercept",
                        "Density",
                        expression(Density^2),
                        "Diameter",
                        expression(Diameter^2),
                        "SPI",
                        "Density:SPI",
                        "Diameter:SPI")
gelman.diag(fixefs)
fixefs <- window(fixefs, start=9001, end=15000, thin=50)
#plot(fixefs, ask=TRUE)
#autocorr.plot(fixefs)
fixefs_comb <- data.frame(combine.mcmc(fixefs))

plot_estimates(fixefs_comb, fixef_names, fixef_names_pretty)
ggsave("growth_model_ranslope_fixefs.png",
       width=3, height=4, dpi=img_dpi)

###############################################################################
### Plot random effects
###############################################################################

load("jags_fit_ranslope_parallel_ranefs.RData")

ranef_names <- c("obs_sigma",
                 "proc_sigma",
                 "sigma_ijk",
                 "sigma_jk",
                 "sigma_k",
                 "sigma_t",
                 "sigma_b_g",
                 "sigma_slp_spi_g",
                 "rho_g")
ranef_names_pretty <- c(expression(sigma[obs]),
                    expression(sigma[proc]),
                    expression(sigma[tree]),
                    expression(sigma[plot]),
                    expression(sigma[site]),
                    expression(sigma[period]),
                    expression(sigma[genus]),
                    expression(sigma[list(SPI,genus)]),
                    expression(rho[list(SPI,genus)]))

gelman.diag(ranefs)
ranefs <- window(ranefs, start=9001, end=15000, thin=50)
#plot(ranefs, ask=TRUE)
#autocorr.plot(ranefs)
ranefs_comb <- data.frame(combine.mcmc(ranefs))

plot_estimates(ranefs_comb, ranef_names, ranef_names_pretty, xmin=0, xmax=3)
ggsave("growth_model_ranslope_ranefs.png",
       width=3, height=5, dpi=img_dpi)
