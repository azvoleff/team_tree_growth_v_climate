library(runjags)
library(coda)
library(dplyr)
library(reshape2)
library(ggplot2)
library(grid) # for unit

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
        geom_errorbarh(aes(xmin=q2pt5, xmax=q97pt5), height=0.25) +
        geom_vline(aes(xintercept=0), color="grey") +
        xlim(c(xmin, xmax)) +
        xlab("") + ylab("") +
        scale_y_discrete(breaks=pars_pretty, labels=pars_pretty)
    return(p)
}

load("jags_fit_full_model_fixefs.RData")
load("jags_fit_full_model_ranefs.RData")
load("jags_fit_full_model_ranefs_g_sigma.RData")
load("jags_fit_full_model_ranefs_g_rho.RData")

start_val <- 20001
thin_val <- 20

fixefs <- window(fixefs, start=start_val, thin=thin_val)
ranefs <- window(ranefs, start=start_val, thin=thin_val)
ranefs_g_sigma <- window(ranefs_g_sigma, start=start_val, thin=thin_val)
ranefs_g_rho <- window(ranefs_g_rho, start=start_val, thin=thin_val)

# gelman.diag(fixefs)
# gelman.diag(ranefs)
# gelman.diag(ranefs_g_sigma)
# gelman.diag(ranefs_g_rho)

# plot(fixefs, ask=TRUE)
# plot(ranefs, ask=TRUE)
# plot(ranefs_g_sigma, ask=TRUE)
# plot(ranefs_g_rho, ask=TRUE)

# autocorr.plot(fixefs)
# autocorr.plot(ranefs)
# autocorr.plot(ranefs_g_sigma)
# autocorr.plot(ranefs_g_rho)

###############################################################################
### Plot fixed effects
###############################################################################


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

fixefs_comb <- data.frame(combine.mcmc(fixefs))

plot_estimates(fixefs_comb, fixefs_names, fixefs_names_pretty,
               scaling=fixefs_scaling)
ggsave("growth_model_fixefs.png", width=3, height=3.75, dpi=img_dpi)

###############################################################################
### Plot random intercepts and process and observation error
###############################################################################


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

ranefs_comb <- data.frame(combine.mcmc(ranefs))

plot_estimates(ranefs_comb, ranef_names, ranef_names_pretty)
ggsave("growth_model_ranefs.png", width=3, height=3.5, dpi=img_dpi)

###############################################################################
### Plot genus-level random effects variances
###############################################################################

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


ranefs_g_sigma_comb <- data.frame(combine.mcmc(ranefs_g_sigma))

plot_estimates(ranefs_g_sigma_comb, ranefs_g_sigma_names,
               ranefs_g_sigma_names_pretty)
ggsave("growth_model_ranefs_g_sigma.png", width=3, height=3, dpi=img_dpi)

###############################################################################
### Plot correlation matrix for genus-level random effects
###############################################################################

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

ranefs_g_rho_comb <- data.frame(combine.mcmc(ranefs_g_rho))

plot_estimates(ranefs_g_rho_comb, ranefs_g_rho_names,
               ranefs_g_rho_names_pretty, xmin=-1, xmax=1)
ggsave("growth_model_ranefs_g_rho.png", width=4, height=5, dpi=img_dpi)

corr_matrix_names <- c("sigma[int,g]", "sigma[MCWD,g]", 
                 "sigma[MCWD^2,g]", "sigma[DBH,g]", 
                 "sigma[DBH^2,g]")

corr_matrix <- matrix(apply(ranefs_g_rho_comb, 2, mean), nrow=5)

corr_matrix_q2pt5 <- matrix(apply(ranefs_g_rho_comb, 2, quantile, .025), nrow=5)
corr_matrix_q97pt5 <- matrix(apply(ranefs_g_rho_comb, 2, quantile, .975), nrow=5)
signif_mat <- (corr_matrix_q2pt5 > 0) | (corr_matrix_q97pt5 < 0)

corrplot <- function(corr_mat, labels=NULL, signif=NULL, base_size=10, 
                     colours=c('dodgerblue4', 'steelblue2', 'white', 
                               'orangered', 'red3')) {
    if (is.null(labels)) labels <- seq(1:nrow(corr_mat))
    colnames(corr_mat) <- labels
    rownames(corr_mat) <- labels
    corr_mat <- corr_mat*upper.tri(corr_mat)
    corr_mat_long <- melt(corr_mat)
    corr_mat_long$Var1 <- ordered(corr_mat_long$Var1, 
                                   levels=(unique(corr_mat_long$Var1)))
    corr_mat_long$Var2 <- ordered(corr_mat_long$Var2, 
                                   levels=rev(unique(corr_mat_long$Var2)))
    corr_mat_long$value[is.na(corr_mat_long$value)] <- 0
    corr_labels <- data.frame(label=labels)
    corr_labels$xpos=seq(1, length.out=nrow(corr_labels))
    corr_labels$ypos=seq(nrow(corr_labels), by=-1, 
                         length.out=nrow(corr_labels))
    p <- ggplot(corr_mat_long, aes(Var1, Var2)) +
        theme_bw(base_size=base_size) +
        geom_tile(fill="white") +
        geom_point(aes(x=Var1, y=Var2, colour=value, size=abs(value)), 
                   data=corr_mat_long[corr_mat_long$value != 0, ]) +
        scale_size_area(max_size=20/12*base_size, guide=FALSE) +
        scale_colour_gradientn(expression(rho), colours=colours, limits=c(-1, 
                                                                          1)) +
        scale_fill_gradientn(expression(rho), colours=colours,
                             limits=c(-1, 1)) +
        theme(axis.text=element_blank(), axis.ticks=element_blank(),
              axis.title=element_blank(), panel.background=element_blank(), 
              panel.border=element_blank(), panel.grid=element_blank(),
              plot.background=element_blank(), 
              legend.title.align=.2, legend.position=c(.7, .8),
              plot.margin=unit(c(0, 0, 0, 0), 'cm'),
              panel.margin=unit(c(0, 0, 0, 0), 'cm'),
              axis.ticks.length=unit(0,"null"),
              axis.ticks.margin=unit(0,"null"),
              legend.direction="horizontal") +
        geom_text(aes(label=label, x=xpos, y=ypos, hjust=.5, vjust=.5), parse=TRUE, 
                  data=corr_labels)
    if (is.null(signif_mat)) {
        geom_text(aes(label=round(value, 2), x=Var1, y=Var2, hjust=.5, vjust=.5), 
                  parse=TRUE, data=corr_matrix_long[corr_matrix_long$value != 0, ])
    } else {
        signif_mat <- signif_mat*upper.tri(signif_mat)
        signif_mat_long <- melt(signif_mat)
        unbold_labels <- corr_mat_long[signif_mat_long$value == FALSE & 
                                       corr_mat_long$value != 0, ]
        p <- p + geom_text(aes(label=round(value, 2), x=Var1, y=Var2, hjust=.5, 
                               vjust=.5), size=4/12*base_size,, alpha=.4, 
                           parse=TRUE, data=unbold_labels)
        bold_labels <- corr_mat_long[signif_mat_long$value == TRUE & 
                                     corr_mat_long$value != 0, ]
        p <- p + geom_text(aes(label=round(value, 2), x=Var1, y=Var2, hjust=.5, 
                               vjust=.5, face="bold"), size=4/12*base_size, 
                           parse=TRUE, data=bold_labels)
    }
    return(p)
}

corrplot(corr_matrix, corr_matrix_names, base_size=10)
ggsave("growth_model_ranefs_g_rho_corrplot.png", width=3, height=3, dpi=img_dpi)
