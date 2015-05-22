# Forest growth model MCMC results
library(RColorBrewer)
library(ggmcmc)
library(gridExtra)
library(foreach)
library(dplyr)

prefixes <- c('D:/azvoleff/Data', # CI-TEAM
              'H:/Data', # Buffalo drive
              'O:/Data', # Blue drive
              '/localdisk/ci_share/azvoleff/Data', # vertica1
              '/localdisk/home/azvoleff/Data') # vertica1
prefix <- prefixes[match(TRUE, unlist(lapply(prefixes, function(x) file_test('-d', x))))]

base_folder <- file.path(prefix, "TEAM", "Tree_Growth")

plot_width <- 3.5
plot_height <- 3
plot_dpi <- 300

# How many mm of precip per 1 unit change? Convert precip from mm to cm in 
# order to make interpretation of results easier
mm_per_unit <- 10

# What thinning was used? This does not change the thinning - it only is used 
# to properly align the plots.
nThin <- 100

# What burn-in was used? This does not change the thinning - it only is used to 
# properly align the plots.
nBurnin  <- 100000

caterpillar <- function(mods, labels=NULL) {
    cis <- group_by(mods, model) %>%
        rename(Parameter=param) %>%
        do(ci(.)) %>%
        rename(param=Parameter)
    # cis$model <- ordered(cis$model, levels=rev(unique(cis$model)))
    # # Flip parameters so first parameters show up on top
    # cis$param <- ordered(cis$param, levels=rev(unique(cis$param)))
    p <- ggplot(cis, aes(x=reorder(param, rev(1:length(param))), 
                         y=median, colour=model, fill=model)) +
        theme_bw(base_size=8) +
        geom_point(position=position_dodge(width=.4), size=1.25) +
        geom_linerange(aes(ymin=Low, ymax=High), size=.75, position=position_dodge(width=.4)) +
        geom_linerange(aes(ymin=low, ymax=high), size=.25, position=position_dodge(width=.4)) +
        xlab('Parameter') +
        ylab('Effect (cm)') +
        scale_colour_manual(values=c('#882255', '#117733', '#88CCEE'))
    if (!is.null(labels)) {
        p <- p + scale_x_discrete(labels=labels) +
            theme(axis.text.x=element_text(vjust=0))
    }
    p <- p + coord_flip()
    p
}

load('B_g_betas_weighted_overall.RData')

B_g_labels <- c('mu_B_g[1]'='int.',
                'mu_B_g[2]'='P',
                'mu_B_g[3]'=expression(P^2),
                'mu_B_g[4]'='T',
                'mu_B_g[5]'=expression(T^2),
                'mu_B_g[6]'='D',
                'mu_B_g[7]'=expression(D^2),
                'mu_B_g[8]'=expression(D%*%P),
                'mu_B_g[9]'=expression(D%*%T))

# Change MCWD units from mm to cm
B_g_betas$value[B_g_betas$param == "B_g_2_mean"] <- B_g_betas$value[B_g_betas$param == "B_g_2_mean"] * mm_per_unit
B_g_betas$value[B_g_betas$param == "B_g_3_mean"] <- B_g_betas$value[B_g_betas$param == "B_g_3_mean"] * mm_per_unit^2
B_g_betas$value[B_g_betas$param == "B_g_8_mean"] <- B_g_betas$value[B_g_betas$param == "B_g_8_mean"] * mm_per_unit

p <- caterpillar(filter(B_g_betas, param %in% paste0('B_g_', c(2:9), '_median')),
                 labels=c('B_g_2_median'='P',
                          'B_g_3_median'=expression(P^2),
                          'B_g_4_median'='T',
                          'B_g_5_median'=expression(T^2),
                          'B_g_6_median'='D',
                          'B_g_7_median'=expression(D^2),
                          'B_g_8_median'=expression(D%*%P),
                          'B_g_9_median'=expression(D%*%T)))
ggsave('caterpillar_climate_interact_weighted.png', p, width=plot_width*1.5, 
       height=plot_height, dpi=plot_dpi)

p <- caterpillar(filter(B_g_betas, param %in% paste0('B_g_', c(2:5, 8:9), '_median')),
                 labels=c('B_g_2_median'='P',
                          'B_g_3_median'=expression(P^2),
                          'B_g_4_median'='T',
                          'B_g_5_median'=expression(T^2),
                          'B_g_8_median'=expression(D%*%P),
                          'B_g_9_median'=expression(D%*%T)))
ggsave('caterpillar_climate_interact_weighted_noD.png', p, width=plot_width*1.5, 
       height=plot_height, dpi=plot_dpi)

# Temp plots
mins <- c(11:24)

# Plot Rhats for each model
foreach(this_model=unique(params$model)) %do% {
    pars <- filter(params, model == this_model)
    attributes(pars)$nChains <- max(params$chain)
    attributes(pars)$nIterations <- max(params$iteration)
    attributes(pars)$nThin <- nThin
    attributes(pars)$nBurnin <- nBurnin

    ggs_Rhat(pars, 'B_k')
    ggsave(paste0('Rhat_interact_B_k_', this_model, '.png'), width=plot_width*2, 
           height=plot_height*2, dpi=plot_dpi)

    ggs_traceplot(pars, 'mu_B_g')
    ggsave(paste0('traceplot_interact_mu_B_g_', this_model, '.png'), width=plot_width*2, 
           height=plot_height*4, dpi=plot_dpi)
    ggs_Rhat(pars, 'mu_B_g')
    ggsave(paste0('Rhat_interact_mu_B_g_', this_model, '.png'), width=plot_width*2, 
           height=plot_height*2, dpi=plot_dpi)

    ggs_traceplot(pars, 'sigma_')
    ggsave(paste0('traceplot_interact_sigma_', this_model, '.png'), width=plot_width*2, 
           height=plot_height*6, dpi=plot_dpi)
    ggs_Rhat(pars, 'sigma_')
    ggsave(paste0('Rhat_interact_sigma_', this_model, '.png'), width=plot_width*2, 
           height=plot_height*2, dpi=plot_dpi)

    # ggs_traceplot(pars, 'int_')
    # ggsave(paste0('traceplot_interact_int_', this_model, '.png'), width=plot_width*2, 
    #        height=plot_height*3, dpi=plot_dpi)
    ggs_Rhat(pars, 'int_')
    ggsave(paste0('Rhat_interact_int_', this_model, '.png'), width=plot_width*2, 
           height=plot_height*3, dpi=plot_dpi)
}

# Function to make a growth prediction as a function of a particular climate 
# variable
make_preds <- function(B, temp, precip, dbhs) {
    X <- matrix(rep(1, length(dbhs)), ncol=1)
    X <- cbind(X, precip)
    X <- cbind(X, precip^2)
    X <- cbind(X, temp)
    X <- cbind(X, temp^2)
    X <- cbind(X, dbhs)
    X <- cbind(X, dbhs^2)
    X <- cbind(X, precip*dbhs)
    X <- cbind(X, temp*dbhs)
    # dbhs is subtracted so that output is growth increment
    these_preds <- X %*% B - rep(dbhs, ncol(B))
    medians <- apply(these_preds, 1, median)
    q2pt5 <- apply(these_preds, 1, quantile, .025)
    q97pt5 <- apply(these_preds, 1, quantile, .975)
    data.frame(median=medians, q2pt5, q97pt5)
}

dbhs <- seq(10, 120, 1)
preds <- foreach(this_model=unique(B_g_betas$model), .combine=rbind) %do% {
    # Need to center the predictor variables
    load(file.path(base_folder, 'Data',
                   paste0("model_data_standardizing_full-", this_model, 
                          "_meanannual-mcwd_run12.RData")))

    # Model effects of temp variation:
    B <- filter(B_g_betas, model == this_model) %>%
        ungroup() %>%
        arrange(param) %>%
        select(param, value)
    # Convert to matrix for linear algebra
    B <- matrix(B$value, nrow=length(unique(B$param)), byrow=TRUE)

    dbhs_centered <- dbhs - dbh_mean

    # Remember temp is standardized, so round it then center it just so the 
    # units look pretty
    temps <- round(temp_mean) - temp_mean
    temps <- c(temps - 1, temps, temps + 1)
    temp_preds <- foreach(temp=temps, .combine=rbind) %do% {
        temp_preds <- make_preds(B, temp, (150 - precip_mean)/mm_per_unit, dbhs_centered)
        temp_preds <- cbind(Panel="Temperature", clim=temp + temp_mean, 
                            dbh=dbhs,
                            temp_preds)
        return(temp_preds)
    }

    # Center, and set units
    precips <- (c(75, 150, 150 + 75) - precip_mean)/mm_per_unit
    precip_preds <- foreach(precip=precips, .combine=rbind) %do% {
        precip_preds <- make_preds(B, 0, precip, dbhs_centered)
        precip_preds <- cbind(Panel="MCWD",
                              clim=precip + precip_mean/mm_per_unit,
                              dbh=dbhs,
                              precip_preds)
        return(precip_preds)
    }

    preds <- rbind(temp_preds, precip_preds)
    preds <- cbind(model=this_model, preds)

    return(preds)
}

preds$clim <- factor(sprintf('%.02f', preds$clim))
preds$model <- factor(preds$model, levels=c('tmn', 'tmp', 'tmx'), 
                      labels=c('Min. Temp. Model',
                               'Mean Temp. Model',
                               'Max. Temp. Model'))
ann_text <- group_by(preds, Panel, model) %>%
    summarise(lab=paste0('Lower=', unique(clim)[1], '\nUpper=', unique(clim)[2]))
ann_text$dbh <- 90
ann_text$mean <- 5

# ggplot(preds, aes(x=dbh, y=mean)) +
#     geom_line(aes(colour=clim)) +
#     geom_ribbon(aes(ymin=q2pt5, ymax=q97pt5, fill=clim), alpha=.2) +
#     facet_grid(Panel~model) +
#     ylab('Increment (cm)') +
#     geom_text(data=ann_text, aes(label=lab))
#
# head(preds)

ps <- foreach(this_panel=unique(preds$Panel), .combine=c) %:%
    foreach(this_model=unique(preds$model)) %do% {
    p <- ggplot(filter(preds, this_model == model, this_panel == Panel), aes(x=dbh, y=median)) +
        geom_line(aes(colour=clim)) +
        geom_ribbon(aes(ymin=q2pt5, ymax=q97pt5, fill=clim), alpha=.2) +
        facet_wrap(~model) +
        xlab('Initial size (cm)') +
        ylab('Growth increment (cm)') +
        coord_cartesian(ylim=c(0, .7)) +
        theme(legend.position=c(.20, .2)) +
        guides(fill=guide_legend(this_panel),
               colour=guide_legend(this_panel))
    return(p)
}

p <- arrangeGrob(ps[[1]], ps[[2]], ps[[3]],
                 ps[[4]], ps[[5]], ps[[6]], ncol=3, nrow=2)
ggsave(paste0('predicted_growth_increments.png'), width=plot_width*2.5,
       height=plot_height*2.5, dpi=plot_dpi, plot=p)
ggsave(paste0('predicted_growth_increments.svg'), width=plot_width*2.5,
       height=plot_height*2.5, dpi=plot_dpi, plot=p)

###############################################################################
## Full model (correlated random effects)

# load(file.path(base_folder, "Extracted_Parameters", "parameter_estimates_correlated.RData"))
