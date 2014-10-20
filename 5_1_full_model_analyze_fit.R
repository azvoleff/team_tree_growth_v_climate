library(runjags)
library(coda)
library(dplyr)
library(reshape2)
library(ggplot2)
#library(GGally)
library(corrplot)
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

ranef_names <- c("sigma[int,g]", "sigma[MCWD,g]", 
                 "sigma[MCWD^2,g]", "sigma[DBH,g]", 
                 "sigma[DBH^2,g]")
corr_matrix <- matrix(apply(ranefs_g_rho_comb, 2, mean), nrow=5)
colnames(corr_matrix) <- ranef_names
rownames(corr_matrix) <- ranef_names
corr_matrix <- t(corr_matrix*lower.tri(corr_matrix))
corr_matrix_long <- melt(corr_matrix)
corr_matrix_long$Var1 <- ordered(corr_matrix_long$Var1, 
                                 levels=(unique(corr_matrix_long$Var1)))
corr_matrix_long$Var2 <- ordered(corr_matrix_long$Var2, 
                                 levels=rev(unique(corr_matrix_long$Var2)))
corr_matrix_long$value[is.na(corr_matrix_long$value)] <- 0
corr_labels <- data.frame(label=c("sigma[int,g]", "sigma[MCWD,g]", 
                                  "sigma[MCWD^2,g]", "sigma[DBH,g]", 
                                  "sigma[DBH^2,g]"))
corr_labels$xpos=seq(1, length.out=nrow(corr_labels))
corr_labels$ypos=seq(nrow(corr_labels), by=-1, length.out=nrow(corr_labels))
ggplot(corr_matrix_long, aes(Var1, Var2)) +
    geom_tile(fill="white") +
    geom_point(aes(x=Var1, y=Var2, colour=value, size=abs(value)), data=corr_matrix_long[corr_matrix_long$value != 0, ]) +
    scale_size_area(max_size=20, guide=FALSE) +
    theme_bw(base_size=10) +
    scale_colour_gradientn(expression(rho),
                         colours=c('firebrick4', 'darkorange3', 'white', 'steelblue2', 'dodgerblue4'),
                         limits=c(-1, 1)) +
    scale_fill_gradientn(expression(rho),
                         colours=c('firebrick4', 'darkorange3', 'white', 'steelblue2', 'dodgerblue4'),
                         limits=c(-1, 1)) +
    theme(axis.text=element_blank(), axis.ticks=element_blank(),
          axis.title=element_blank(), 
          panel.background=element_blank(), panel.border=element_blank(),
          panel.grid=element_blank(),
          plot.background=element_blank(), 
          legend.title.align=.2, legend.position=c(.75, .75),
          plot.margin=unit(c(.1, .1, .1, .1), 'cm'),
          legend.direction="horizontal") +
    geom_text(aes(label=label, x=xpos, y=ypos, hjust=.5, vjust=.5), parse=TRUE, 
              data=corr_labels) +
    geom_text(aes(label=round(value, 2), x=Var1, y=Var2, hjust=.5, vjust=.5), 
              parse=TRUE, data=corr_matrix_long[corr_matrix_long$value != 0, ])
ggsave("growth_model_ranefs_g_rho_corrplot.png", width=4, height=4, dpi=img_dpi)

corrplot.mixed(corr_matrix)

# gelman.diag(ranefs_g_sigma)
# ranefs_g_sigma <- window(ranefs_g_sigma, start=9001, end=15000, thin=50)
#plot(ranefs_g_sigma, ask=TRUE)
#autocorr.plot(ranefs_g_sigma)

# D <- round(corr_matrix, 1)
# D <- D * lower.tri(D)
# D <- as.data.frame(D)
# D <- data.frame(row = colnames(corr_matrix), D)
# D <- melt(D, id.vars = "row")
# M <- corr_matrix * lower.tri(corr_matrix)
# M <- as.data.frame(M)
# M <- data.frame(row = colnames(corr_matrix), M)
# names(M) <- c("row", colnames(corr_matrix))
# M <- melt(M, id.vars = "row")
# M$value[M$value == 0] <- NA
# s <- seq(-1, 1, by = 0.25)
# M$value <- cut(M$value, breaks = s, include.lowest = TRUE, 
#                label = cut(s, breaks = s)[-1])
# M$row <- factor(M$row, levels=corr_mat_names)
# M$num <- as.numeric(M$value)
# diag <- subset(M, row == variable)
# po.nopanel <- list(theme(panel.background = element_blank(), 
#                          panel.grid.minor = element_blank(), panel.grid.major = element_blank(), 
#                          axis.text.x = element_text(angle = -90)))
# p = ggplot(M, aes(row, variable)) 
# p = p + scale_colour_brewer(expression(rho), palette = "RdYlGn") + 
#     scale_size_area(expression(rho), max_size = 6,
#                     labels = 
#                     levels(M$value)[table(M$value) > 0]) +
#     geom_point(aes(size = num, colour = value))
# p = p + geom_text(data = diag, aes(label = variable)) +
#     scale_x_discrete(breaks = NULL) +
#     scale_y_discrete(breaks = NULL) + 
#     labs(x = NULL, y = NULL) + coord_equal() + po.nopanel
