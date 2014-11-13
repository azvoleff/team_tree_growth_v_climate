library(runjags)
library(coda)
library(dplyr)
library(reshape2)
library(ggplot2)
library(grid) # for unit

img_dpi <- 300

n_B_g <- 7

model_type <- "full"
#model_type <- "testing"
temp_var <- "tmp_meanannual"
precip_var <- "mcwd_run12"
run_ID <- "vertica1.team.sdsc.edu_20141110152032_extend1"

in_folder <- 'MCMC_Chains'
suffix <- paste0(model_type, '-', temp_var, '-', precip_var)

load(file.path('Data', paste0("model_data_standardizing_", suffix, ".RData")))

plot_estimates <- function(mcmc_ests, pars=NULL, pars_pretty=NULL, xmin=NULL, 
                           xmax=NULL, scaling=NULL) {
    if (is.null(pars)) pars <- names(mcmc_ests)
    if (is.null(pars_pretty)) pars_pretty <- pars
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

    mcmc_ests$variable <- ordered(mcmc_ests$variable, 
                                  levels=rev(pars),
                                  labels=rev(pars_pretty))

    if (is.null(xmin)) xmin <- min(mcmc_ests$q2pt5)
    if (is.null(xmax)) xmax <- max(mcmc_ests$q97pt5)
    p <- ggplot(mcmc_ests, aes(median, variable)) +
        theme_bw(base_size=10) +
        geom_point(size=1.5) + 
        geom_errorbarh(aes(xmin=q2pt5, xmax=q97pt5), height=0.25) +
        geom_vline(aes(xintercept=0), color="grey") +
        xlim(c(xmin, xmax)) +
        xlab("") + ylab("") +
        scale_y_discrete(breaks=pars_pretty, labels=pars_pretty)
    return(p)
}

load(file.path(in_folder, paste0(suffix, "jags_fit_full_model_fixefs.RData")))
load(file.path(in_folder, paste0(suffix, "jags_fit_full_model_ranefs_sigmas.RData")))
load(file.path(in_folder, paste0(suffix, "jags_fit_full_model_ranefs_B_T.RData")))
load(file.path(in_folder, paste0(suffix, "jags_fit_full_model_ranefs_mu_B_g.RData")))
load(file.path(in_folder, paste0(suffix, "jags_fit_full_model_ranefs_sigma_B_g.RData")))
load(file.path(in_folder, paste0(suffix, "jags_fit_full_model_ranefs_rho_B_g.RData")))

start_val <- 18000
thin_val <- 4

fixefs <- window(fixefs, start=start_val, thin=thin_val)
ranefs_sigmas <- window(ranefs_sigmas, start=start_val, thin=thin_val)
ranefs_mu_B_g <- window(ranefs_mu_B_g, start=start_val, thin=thin_val)
ranefs_B_T <- window(ranefs_B_T, start=start_val, thin=thin_val)
ranefs_sigma_B_g <- window(ranefs_sigma_B_g, start=start_val, thin=thin_val)
ranefs_rho_B_g <- window(ranefs_rho_B_g, start=start_val, thin=thin_val)

# gelman.diag(fixefs)
# gelman.diag(ranefs_sigmas)
# gelman.diag(ranefs_sigma_B_g)
# gelman.diag(ranefs_mu_B_g)
# gelman.diag(ranefs_rho_B_g)

# plot(fixefs, ask=TRUE)
# plot(ranefs_sigmas, ask=TRUE)
# plot(ranefs_sigma_B_g, ask=TRUE)
# plot(ranefs_mu_B_g, ask=TRUE)
# plot(ranefs_rho_B_g, ask=TRUE)

# gelman.plot(fixefs)
# gelman.plot(ranefs_sigmas)
# gelman.plot(ranefs_sigma_B_g)
# gelman.plot(ranefs_mu_B_g)

# HPDinterval(fixefs)
# HPDinterval(ranefs_sigmas)

# autocorr.plot(fixefs)
# autocorr.plot(ranefs)
# autocorr.plot(ranefs_mu_B_g)
# autocorr.plot(ranefs_sigma_B_g)
# autocorr.plot(ranefs_rho_B_g)

###############################################################################
### Plot fixed effects
###############################################################################

fixefs_names <- c("B[1]",
                  "B[2]")
fixefs_names_pretty <- c("Density",
                         expression(Density^2))
fixefs_scaling <- c(dbh_sd/WD_sd,
                    dbh_sd/WD_sd)

fixefs_comb <- data.frame(combine.mcmc(fixefs))

plot_estimates(fixefs_comb, fixefs_names, fixefs_names_pretty,
               scaling=fixefs_scaling)
ggsave("growth_model_fixefs.png", width=3, height=1.5, dpi=img_dpi)

###############################################################################
### Plot temperature model
###############################################################################

ranefs_B_T_names <- c("B[1]",
                      "B[2]")
ranefs_B_T_names_pretty <- c("Density",
                         expression(Density^2))
ranefs_B_T_scaling <- c(dbh_sd/WD_sd,
                        dbh_sd/WD_sd)

ranefs_B_T_comb <- data.frame(combine.mcmc(ranefs_B_T))

plot_estimates(ranefs_B_T_comb)
ggsave("growth_model_ranefs_B_T.png", width=3, height=1.5, dpi=img_dpi)

###############################################################################
### Plot random intercepts and process and observation error
###############################################################################

ranef_sigmas_names <- c("sigma_obs",
                        "sigma_proc",
                        "sigma_int_ijk",
                        "sigma_int_jk",
                        "sigma_int_k",
                        "sigma_int_t")
ranef_sigmas_names_pretty <- c(expression(sigma[obs]),
                               expression(sigma[proc]),
                               expression(sigma[tree]),
                               expression(sigma[plot]),
                               expression(sigma[site]),
                               expression(sigma[period]))
ranefs_sigmas_scaling <- c(dbh_sd,
                           dbh_sd,
                           dbh_sd,
                           dbh_sd,
                           dbh_sd,
                           dbh_sd)

ranefs_sigmas_comb <- data.frame(combine.mcmc(ranefs_sigmas))

plot_estimates(ranefs_sigmas_comb, ranef_sigmas_names, ranef_sigmas_names_pretty, 
               scaling=ranefs_sigmas_scaling)
ggsave("growth_model_ranefs_sigmas.png", width=3, height=3.5, dpi=img_dpi)

###############################################################################
### Plot genus-level random effects means
###############################################################################

ranefs_mu_B_g_names <- c("mu_B_g[1]",
                         "mu_B_g[2]",
                         "mu_B_g[3]",
                         "mu_B_g[4]",
                         "mu_B_g[5]",
                         "mu_B_g[6]",
                         "mu_B_g[7]")
ranefs_mu_B_g_names_pretty <- c(expression(mu[int,g]),
                                expression(mu[P,g]),
                                expression(mu[P^2,g]),
                                expression(mu[T,g]),
                                expression(mu[T^2,g]),
                                expression(mu[D,g]),
                                expression(mu[D^2,g]))
ranef_mu_b_g_scaling <- c(dbh_sd,
                          dbh_sd/temp_sd,
                          dbh_sd/temp_sd,
                          (dbh_sd/precip_sd) * 100, # Convert from mm to 10s of cm
                          (dbh_sd/precip_sd) * 100, # Convert from mm to 10s of cm
                          dbh_sd/dbh_sd,
                          dbh_sd/dbh_sd)

ranefs_mu_B_g_comb <- data.frame(combine.mcmc(ranefs_mu_B_g))

plot_estimates(ranefs_mu_B_g_comb, ranefs_mu_B_g_names,
               ranefs_mu_B_g_names_pretty, scaling=ranef_mu_b_g_scaling)
ggsave("growth_model_ranefs_mu_B_g.png", width=3, height=3, dpi=img_dpi)

###############################################################################
### Plot genus-level random effects variances
###############################################################################

ranefs_sigma_B_g_names <- c("sigma_B_g[1]",
                            "sigma_B_g[2]",
                            "sigma_B_g[3]",
                            "sigma_B_g[4]",
                            "sigma_B_g[5]",
                            "sigma_B_g[6]",
                            "sigma_B_g[7]")
ranefs_sigma_B_g_names_pretty <- c(expression(sigma[int,g]),
                                   expression(sigma[P,g]),
                                   expression(sigma[P^2,g]),
                                   expression(sigma[T,g]),
                                   expression(sigma[T^2,g]),
                                   expression(sigma[D,g]),
                                   expression(sigma[D^2,g]))
ranefs_sigma_B_g_scaling <- c(dbh_sd,
                              dbh_sd/temp_sd,
                              dbh_sd/temp_sd,
                              (dbh_sd/precip_sd) * 100, # Convert from mm to 10s of cm
                              (dbh_sd/precip_sd) * 100, # Convert from mm to 10s of cm
                              dbh_sd/dbh_sd,
                              dbh_sd/dbh_sd)

ranefs_sigma_B_g_comb <- data.frame(combine.mcmc(ranefs_sigma_B_g))

plot_estimates(ranefs_sigma_B_g_comb, ranefs_sigma_B_g_names,
               ranefs_sigma_B_g_names_pretty, scaling=ranefs_sigma_B_g_scaling)
ggsave("growth_model_ranefs_sigma_B_g.png", width=3, height=3, dpi=img_dpi)

###############################################################################
### Plot correlation matrix for genus-level random effects
###############################################################################

ranefs_rho_B_g_names <- c("rho_B_g[1,2]",
                          "rho_B_g[1,3]",
                          "rho_B_g[1,4]",
                          "rho_B_g[1,5]",
                          "rho_B_g[1,6]",
                          "rho_B_g[1,7]",
                          "rho_B_g[2,3]",
                          "rho_B_g[2,4]",
                          "rho_B_g[2,5]",
                          "rho_B_g[2,6]",
                          "rho_B_g[2,7]",
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
ranefs_rho_B_g_names_pretty <- c(expression(rho[list(int,P)]),
                                 expression(rho[list(int,P^2)]),
                                 expression(rho[list(int,T)]),
                                 expression(rho[list(int,T^2)]),
                                 expression(rho[list(int,D)]),
                                 expression(rho[list(int,D^2)]),
                                 expression(rho[list(P,P^2)]),
                                 expression(rho[list(P,T)]),
                                 expression(rho[list(P,T^2)]),
                                 expression(rho[list(P,D)]),
                                 expression(rho[list(P,D^2)]),
                                 expression(rho[list(P^2,T)]),
                                 expression(rho[list(P^2,T^2)]),
                                 expression(rho[list(P^2,D)]),
                                 expression(rho[list(P^2,D^2)]),
                                 expression(rho[list(T,T^2)]),
                                 expression(rho[list(T,D)]),
                                 expression(rho[list(T,D^2)]),
                                 expression(rho[list(T^2,D)]),
                                 expression(rho[list(T^2,D^2)]),
                                 expression(rho[list(D,D^2)]))
ranefs_rho_B_g_comb <- data.frame(combine.mcmc(ranefs_rho_B_g))

plot_estimates(ranefs_rho_B_g_comb, ranefs_rho_B_g_names,
               ranefs_rho_B_g_names_pretty, xmin=-1, xmax=1)
ggsave("growth_model_ranefs_rho_B_g.png", width=4, height=5, dpi=img_dpi)

# Make a color coded plot of the correlation matrix
corr_matrix_names <- c("Intercept", "P", "P^2", "T", "T^2", "D", "D^2")
corr_matrix <- matrix(apply(ranefs_rho_B_g_comb, 2, mean), nrow=n_B_g)
corr_matrix_q2pt5 <- matrix(apply(ranefs_rho_B_g_comb, 2, quantile, .025), nrow=5)
corr_matrix_q97pt5 <- matrix(apply(ranefs_rho_B_g_comb, 2, quantile, .975), nrow=5)
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
        geom_text(aes(label=label, x=xpos, y=ypos, hjust=.5, vjust=.3), parse=TRUE, 
                  data=corr_labels, size=4/12*base_size)
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
ggsave("growth_model_ranefs_rho_B_g_corrplot.png", width=4, height=4, dpi=img_dpi)
