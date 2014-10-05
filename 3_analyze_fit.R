library(runjags)
library(coda)
library(dplyr)
library(reshape2)
library(ggplot2)

img_height <- 4
img_width <- 3
img_dpi <- 300

load("jags_fit_full_parallel.RData")
#load("jags_fit_ranslope.RData")

###############################################################################
### Plot fixed effects
###############################################################################

fixed_names <- c("intercept",
                   "slp_WD",
                   "slp_WD_sq",
                   "slp_dbh",
                   "slp_dbh_sq",
                   "'slp_spi'",
                   "inter_spi_WD",
                   "inter_spi_dbh")
fixed_names_pretty <-c("Intercept",
              "Density",
              "Density^2",
              "Diameter",
              "Diameter^2",
              "SPI",
              "Density:SPI",
              "Diameter:SPI")

# print(jags_fit_p, vars=fixed_names)
# plot(jags_fit_p, vars=fixed_names, layout=c(4, 4))
fixefs <- as.mcmc.list(jags_fit_p, fixed_names)
#gelman.diag(fixefs)
fixefs <- window(fixefs, start=9001, end=15000, thin=50)
#autocorr.plot(fixefs)
fixefs <- data.frame(combine.mcmc(fixefs))
fixefs <- melt(fixefs)
fixefs <- group_by(fixefs, variable) %>%
    summarize(median=median(value),
              q2pt5=quantile(value, .025),
              q97pt5=quantile(value, .975))
fixefs$variable <- ordered(fixefs$variable, 
                           levels=rev(fixed_names),
                           labels=rev(fixed_names_pretty))

ggplot(fixefs, aes(median, variable)) +
    theme_bw(base_size=10) +
    geom_point() + 
    geom_errorbarh(aes(xmin=q2pt5, xmax=q97pt5), height=0.5) +
    geom_vline(aes(xintercept=0), color="grey") +
    #xlim(c(-1, 1)) +
    xlab("") + ylab("")
ggsave("growth_model_ranslope_fixedefs.png",
       width=img_width, height=img_height, dpi=img_dpi)

###############################################################################
### Plot random effects
###############################################################################
random_names <- c("obs_sigma",
                    "proc_sigma",
                    "sigma_ijk",
                    "sigma_jk",
                    "sigma_k",
                    "sigma_t",
                    "sigma_g",
                    "sigma_slp_spi_g",
                    "rho")
random_names_pretty <- c("sigma_obs",
                    "sigma_proc",
                    "sigma_ijk",
                    "sigma_jk",
                    "sigma_k",
                    "sigma_t",
                    "sigma_g",
                    "sigma_(SPI g)",
                    "rho")

# print(jags_fit_p, vars=random_names)
# plot(jags_fit_p, vars=random_names, layout=c(4, 4))
ranefs <- as.mcmc.list(jags_fit_p, random_names)
#gelman.diag(ranefs)
ranefs <- window(ranefs, start=9001, end=15000, thin=50)
#autocorr.plot(ranefs)
ranefs <- data.frame(combine.mcmc(ranefs))
ranefs <- melt(ranefs)
ranefs <- group_by(ranefs, variable) %>%
    summarize(median=median(value),
              q2pt5=quantile(value, .025),
              q97pt5=quantile(value, .975))
ranefs$variable <- ordered(ranefs$variable, 
                           levels=rev(random_names),
                           labels=rev(random_names_pretty))

ggplot(ranefs, aes(median, variable)) +
    theme_bw(base_size=10) +
    geom_point() + 
    geom_errorbarh(aes(xmin=q2pt5, xmax=q97pt5), height=0.5) +
    geom_vline(aes(xintercept=0), color="grey") +
    #xlim(c(-1, 1)) +
    xlab("") + ylab("")
ggsave("growth_model_ranslope_randomefs.png",
       width=img_width, height=img_height, dpi=img_dpi)
