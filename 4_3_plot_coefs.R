# Forest growth model MCMC results

prefixes <- c('D:/azvoleff/Data', # CI-TEAM
              'H:/Data', # Buffalo drive
              'O:/Data', # Blue drive
              '/localdisk/home/azvoleff/Data') # vertica1
prefix <- prefixes[match(TRUE, unlist(lapply(prefixes, function(x) file_test('-d', x))))]

base_folder <- file.path(prefix, "TEAM", "Tree_Growth")

library(ggmcmc)
library(gridExtra)
library(foreach)

plot_width <- 3.5
plot_height <- 3
plot_dpi <- 300

# Calculate weights for each genus ID (doesn't matter which temp_var is used 
# since the genus IDs and frequencies are the same across all simulations).
load(file.path(base_folder, 'Data', 
               paste0("model_data_wide_full-tmn_meanannual-mcwd_run12.RData")))
merged <- tbl_df(data.frame(site_ID=model_data$site_ID, 
                            genus_ID=as.integer(model_data$genus_ID)))
genus_weights <- group_by(merged, genus_ID) %>%
    summarize(n=n()) %>%
    ungroup() %>%
    mutate(weight=n/sum(n)) %>%
    select(-n) %>%
    arrange(desc(weight))

multimodel_caterpillar <- function(mods, labels=NULL) {
    cis <- group_by(mods, Model) %>%
        do(ci(.))
    p <- ggplot(cis, aes(x=Parameter, y=median, colour=Model)) +
        theme_bw(base_size=8) +
        geom_point(position=position_dodge(width=.4), size=1.25) +
        geom_linerange(aes(ymin=Low, ymax=High), size=.75, position=position_dodge(width=.4)) +
        geom_linerange(aes(ymin=low, ymax=high), size=.25, position=position_dodge(width=.4)) +
        ylab('Effect (cm)')

    if (!is.null(labels)) {
        p <- p + scale_x_discrete(labels=labels) +
            theme(axis.text.x=element_text(vjust=0))
    }
    p
}

# Function to calculated weighted coefficients from an array of MCMC results.  
# Weights should be a 2 column data.frame with IDs in the first column, and 
# weights in the second. D should be a ggs object with parameter IDs added 
# using the jags_param_ids function
weight_coef <- function(d, w) {
    d <- left_join(d, w)
    d <- group_by(d, Model, Chain, Iteration, param_ID) %>%
        summarise(Parameter=paste(Parameter_Base[1], param_ID[1], 'mean', sep='_'),
                  value=sum(weight * value))
    return(d)
}

###############################################################################
## Simple model (no random effects)
##

# load(file.path(base_folder, "Extracted_Parameters", "parameter_estimates_simple.RData"))
# params$Model <- factor(params$Model, levels=c('tmn', 'tmp', 'tmx'),
#                        labels=c('Min. Temp.', 'Mean Temp.', 'Max. Temp.'))
#
# multimodel_caterpillar(filter(params, Parameter %in% paste0('B[', 1:7, ']')))
# ggsave('caterpillar_climate_simple.png', p, width=plot_width, 
#        height=plot_height, dpi=plot_dpi)
#
# p <- multimodel_caterpillar(filter(params, Parameter %in% c('B[2]', 'B[3]', 
#                                                             'B[4]', 'B[5]')),
#                             labels=c('B[2]'='P',
#                                      'B[3]'=expression(P^2), 
#                                      'B[4]'='T',
#                                      'B[5]'=expression(T^2)))
# ggsave('caterpillar_climate_simple.png', p, width=plot_width, 
#        height=plot_height, dpi=plot_dpi)

###############################################################################
## Full model (interactions, uncorrelated random effects)

load(file.path(base_folder, "Extracted_Parameters", "parameter_estimates_interact.RData"))

B_g_betas <- weight_coef(filter(params, Parameter_Base == 'B_g') %>% 
                         rename(genus_ID=row_ID, param_ID=col_ID), 
                         genus_weights)


multimodel_caterpillar(filter(params, Parameter %in% paste0('mu_B_g[', 1:9, ']')))

p <- multimodel_caterpillar(filter(params, Parameter %in% paste0('mu_B_g[', c(2:9), ']')),
                            labels=c('mu_B_g[2]'='P',
                                     'mu_B_g[3]'=expression(P^2),
                                     'mu_B_g[4]'='T',
                                     'mu_B_g[5]'=expression(T^2),
                                     'mu_B_g[6]'='D',
                                     'mu_B_g[7]'=expression(D^2),
                                     'mu_B_g[8]'=expression(D%*%P),
                                     'mu_B_g[9]'=expression(D%*%T)))
ggsave('caterpillar_climate_interact_unweighted.png', p, width=plot_width*1.5, 
       height=plot_height, dpi=plot_dpi)

p <- multimodel_caterpillar(filter(B_g_betas, Parameter %in% paste0('B_g_', c(2:9), '_mean')),
                            labels=c('B_g_2_mean'='P',
                                     'B_g_3_mean'=expression(P^2),
                                     'B_g_4_mean'='T',
                                     'B_g_5_mean'=expression(T^2),
                                     'B_g_6_mean'='D',
                                     'B_g_7_mean'=expression(D^2),
                                     'B_g_8_mean'=expression(D%*%P),
                                     'B_g_9_mean'=expression(D%*%T)))
ggsave('caterpillar_climate_interact_weighted.png', p, width=plot_width*1.5, 
       height=plot_height, dpi=plot_dpi)

# Plot Rhats for each model
foreach(Model=unique(params$Model) %do% {
    pars <- filter(params, Model == Model)
    attributes(pars)$nChains <- 3
    attributes(pars)$nIterations <- max(params$Iteration)
    attributes(pars)$nThin <- 100
    attributes(pars)$nBurnin <- 60000

    ggs_Rhat(pars, 'B_k')
    ggsave(paste0('Rhat_interact_B_k_', model, '.png'), width=plot_width*2, 
           height=plot_height*2, dpi=plot_dpi)

    ggs_traceplot(pars, 'mu_B_g')
    ggsave(paste0('traceplot_interact_mu_B_g_', model, '.png'), width=plot_width*2, 
           height=plot_height*4, dpi=plot_dpi)
    ggs_Rhat(pars, 'mu_B_g')
    ggsave(paste0('Rhat_interact_mu_B_g_', model, '.png'), width=plot_width*2, 
           height=plot_height*2, dpi=plot_dpi)

    ggs_traceplot(pars, 'sigma_')
    ggsave(paste0('traceplot_interact_sigma_', model, '.png'), width=plot_width*2, 
           height=plot_height*6, dpi=plot_dpi)
    ggs_Rhat(pars, 'sigma_')
    ggsave(paste0('Rhat_interact_sigma_', model, '.png'), width=plot_width*2, 
           height=plot_height*2, dpi=plot_dpi)

    # ggs_traceplot(pars, 'int_')
    # ggsave(paste0('traceplot_interact_int_', model, '.png'), width=plot_width*2, 
    #        height=plot_height*3, dpi=plot_dpi)
    ggs_Rhat(pars, 'int_')
    ggsave(paste0('Rhat_interact_int_', model, '.png'), width=plot_width*2, 
           height=plot_height*3, dpi=plot_dpi)
}

# Function to make a growth prediction as a function of a particular climate 
# variable
make_preds <- function(B, clim_low, clim_high, clim_mean, dbhs, dbh_mean, panel) {
    dbhs_centered <- dbhs - dbh_mean
    clim_low_centered <- ceiling(clim_low) - clim_mean
    clim_high_centered <- floor(clim_high) - clim_mean
    preds <- foreach(clim=c(clim_low_centered, clim_high_centered), .combine=rbind) %do% {
        # X <- matrix(rep(clim, length(dbhs_centered)), ncol=1)
        X <- matrix(rep(1, length(dbhs_centered)), ncol=1)
        X <- cbind(X, clim)
        X <- cbind(X, clim^2)
        X <- cbind(X, dbhs_centered)
        X <- cbind(X, dbhs_centered^2)
        X <- cbind(X, clim*dbhs_centered)
        # dbhs_centered is subtracted so that output is growth increment
        these_preds <- X %*% B - rep(dbhs_centered, ncol(B))
        means <- apply(these_preds, 1, mean)
        q2pt5 <- apply(these_preds, 1, quantile, .025)
        q97pt5 <- apply(these_preds, 1, quantile, .975)
        data.frame(panel=panel, dbh=dbhs, clim=clim + clim_mean, mean=means, 
                   q2pt5, q97pt5)
    }
    return(preds)
}

dbhs <- seq(10, 125, 1)
preds <- foreach(Model=unique(B_g_betas$Model), .combine=rbind) %do% {
    # Need to center the predictor variables
    load(file.path(base_folder, 'Data',
                   paste0("model_data_standardizing_full-", Model, 
                          "_meanannual-mcwd_run12.RData")))

    # Model effects of temp variation:
    B_temp <- filter(B_g_betas, Parameter %in% paste0('B_g_', c(1,4:7,9), '_mean'),
                     Model == Model) %>%
        ungroup() %>%
        arrange(Parameter, Iteration) %>%
        select(Parameter, value)
    # Convert to matrix for linear algebra
    B_temp <- matrix(B_temp$value, nrow=length(unique(B_temp$Parameter)), byrow=TRUE)
    temp_preds <- make_preds(B_temp, temp_min, temp_max, temp_mean, dbhs, 
                             dbh_mean, panel="Temperature")

    # Model effects of MCWD variation:
    B_mcwd <- filter(B_g_betas, Parameter %in% paste0('B_g_', c(1:3,6:7,8), '_mean'),
                     Model == Model) %>%
        ungroup() %>%
        arrange(Parameter, Iteration) %>%
        select(Parameter, value)
    # Convert to matrix for linear algebra
    B_mcwd <- matrix(B_mcwd$value, nrow=length(unique(B_mcwd$Parameter)), byrow=TRUE)
    precip_preds <- make_preds(B_mcwd, precip_min, precip_max, precip_mean, 
                               dbhs, dbh_mean, panel="MCWD")

    preds <- rbind(temp_preds, precip_preds)
    preds <- cbind(Model=Model, preds)

    return(preds)
}

preds$clim <- factor(sprintf('%.02f', preds$clim))
preds$Model <- factor(preds$Model, levels=c('tmn', 'tmp', 'tmx'), 
                      labels=c('Min. Temp. Model',
                               'Mean Temp. Model',
                               'Max. Temp. Model'))
ann_text <- group_by(preds, panel, Model) %>%
    summarise(lab=paste0('Lower=', unique(clim)[1], '\nUpper=', unique(clim)[2]))
ann_text$dbh <- 110
ann_text$mean <- 25

ggplot(preds, aes(x=dbh, y=mean)) +
    geom_line(aes(colour=clim)) +
    geom_ribbon(aes(ymin=q2pt5, ymax=q97pt5, fill=clim), alpha=.2) +
    facet_grid(panel~Model) +
    ylab('Increment (cm)') +
    geom_text(data=ann_text, aes(label=lab))

head(preds)

ps <- foreach(this_panel=unique(preds$panel), .combine=c) %:%
    foreach(this_model=unique(preds$Model)) %do% {
    p <- ggplot(filter(preds, this_model == Model, this_panel == panel), aes(x=dbh, y=mean)) +
        geom_line(aes(colour=clim)) +
        geom_ribbon(aes(ymin=q2pt5, ymax=q97pt5, fill=clim), alpha=.2) +
        facet_wrap(~Model) +
        xlab('Initial size (cm)') +
        ylab('Growth increment (cm)') +
        coord_cartesian(ylim=c(-20, 25)) +
        theme(legend.position=c(.8, .85)) +
        guides(fill=guide_legend(this_panel),
               colour=guide_legend(this_panel))
    return(p)
}

grid.arrange(ps[[1]], ps[[2]], ps[[3]],
             ps[[4]], ps[[5]], ps[[6]], ncol=3, nrow=2)

###############################################################################
## Full model (correlated random effects)

# load(file.path(base_folder, "Extracted_Parameters", "parameter_estimates_correlated.RData"))
