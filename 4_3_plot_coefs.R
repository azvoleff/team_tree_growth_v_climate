# Forest growth model MCMC results

prefixes <- c('D:/azvoleff/Data', # CI-TEAM
              'H:/Data', # Buffalo drive
              'O:/Data', # Blue drive
              '/localdisk/home/azvoleff/Data') # vertica1
prefix <- prefixes[match(TRUE, unlist(lapply(prefixes, function(x) file_test('-d', x))))]

base_folder <- file.path(prefix, "TEAM", "Tree_Growth")

library(RColorBrewer)
library(ggmcmc)
library(gridExtra)
library(foreach)

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
nBurnin  <- 20000

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

caterpillar <- function(mods, labels=NULL) {
    cis <- group_by(mods, Model) %>%
        do(ci(.))
    # cis$Model <- ordered(cis$Model, levels=rev(unique(cis$Model)))
    # # Flip parameters so first parameters show up on top
    # cis$Parameter <- ordered(cis$Parameter, levels=rev(unique(cis$Parameter)))
    p <- ggplot(cis, aes(x=reorder(Parameter, rev(1:length(Parameter))), 
                         y=median, colour=Model, fill=Model)) +
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
# caterpillar(filter(params, Parameter %in% paste0('B[', 2:7, ']')))
# ggsave('caterpillar_climate_simple.png', p, width=plot_width, 
#        height=plot_height, dpi=plot_dpi)
#
# p <- caterpillar(filter(params, Parameter %in% c('B[2]', 'B[3]', 'B[4]', 'B[5]')),
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
params$value[params$Parameter == "mu_B_g[2]"] <- params$value[params$Parameter == "mu_B_g[2]"] * mm_per_unit
params$value[params$Parameter == "mu_B_g[3]"] <- params$value[params$Parameter == "mu_B_g[3]"] * mm_per_unit^2
params$value[params$Parameter == "mu_B_g[8]"] <- params$value[params$Parameter == "mu_B_g[8]"] * mm_per_unit
B_g_betas$value[B_g_betas$Parameter == "B_g_2_mean"] <- B_g_betas$value[B_g_betas$Parameter == "B_g_2_mean"] * mm_per_unit
B_g_betas$value[B_g_betas$Parameter == "B_g_3_mean"] <- B_g_betas$value[B_g_betas$Parameter == "B_g_3_mean"] * mm_per_unit^2
B_g_betas$value[B_g_betas$Parameter == "B_g_8_mean"] <- B_g_betas$value[B_g_betas$Parameter == "B_g_8_mean"] * mm_per_unit

caterpillar(filter(params, Parameter %in% paste0('mu_B_g[', 1:9, ']')), 
            labels=B_g_labels)

p <- caterpillar(filter(params, Parameter %in% paste0('mu_B_g[', c(2:9), ']')), 
                 labels=B_g_labels[2:9])

ggsave('caterpillar_climate_interact_unweighted.png', p, width=plot_width*1.5, 
       height=plot_height, dpi=plot_dpi)

p <- caterpillar(filter(B_g_betas, Parameter %in% paste0('B_g_', c(2:9), '_mean')),
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
foreach(this_model=unique(params$Model)) %do% {
    pars <- filter(params, Model == this_model)
    attributes(pars)$nChains <- max(params$Chain)
    attributes(pars)$nIterations <- max(params$Iteration)
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
    means <- apply(these_preds, 1, mean)
    q2pt5 <- apply(these_preds, 1, quantile, .025)
    q97pt5 <- apply(these_preds, 1, quantile, .975)
    data.frame(mean=means, q2pt5, q97pt5)
}

dbhs <- seq(10, 120, 1)
preds <- foreach(this_model=unique(B_g_betas$Model), .combine=rbind) %do% {
    # Need to center the predictor variables
    load(file.path(base_folder, 'Data',
                   paste0("model_data_standardizing_full-", this_model, 
                          "_meanannual-mcwd_run12.RData")))

    # Model effects of temp variation:
    B <- filter(B_g_betas, Model == this_model) %>%
        ungroup() %>%
        arrange(Parameter, Iteration) %>%
        select(Parameter, value)
    # Convert to matrix for linear algebra
    B <- matrix(B$value, nrow=length(unique(B$Parameter)), byrow=TRUE)

    dbhs_centered <- dbhs - dbh_mean

    # Remember temp is standardized, so round it then center it just so the 
    # units look pretty
    temps <- round(temp_mean) - temp_mean
    temps <- c(temps - 2, temps, temps + 2)
    temp_preds <- foreach(temp=temps, .combine=rbind) %do% {
        temp_preds <- make_preds(B, temp, (0 - precip_mean)/mm_per_unit, dbhs_centered)
        temp_preds <- cbind(Panel="Temperature", clim=temp + temp_mean, 
                            dbh=dbhs,
                            temp_preds)
        return(temp_preds)
    }

    # Center, and set units
    precips <- (c(0, 150, 150 + 100) - precip_mean)/mm_per_unit
    precip_preds <- foreach(precip=precips, .combine=rbind) %do% {
        precip_preds <- make_preds(B, 0, precip, dbhs_centered)
        precip_preds <- cbind(Panel="MCWD",
                              clim=precip + precip_mean/mm_per_unit,
                              dbh=dbhs,
                              precip_preds)
        return(precip_preds)
    }

    preds <- rbind(temp_preds, precip_preds)
    preds <- cbind(Model=this_model, preds)

    return(preds)
}

preds$clim <- factor(sprintf('%.02f', preds$clim))
preds$Model <- factor(preds$Model, levels=c('tmn', 'tmp', 'tmx'), 
                      labels=c('Min. Temp. Model',
                               'Mean Temp. Model',
                               'Max. Temp. Model'))
ann_text <- group_by(preds, Panel, Model) %>%
    summarise(lab=paste0('Lower=', unique(clim)[1], '\nUpper=', unique(clim)[2]))
ann_text$dbh <- 90
ann_text$mean <- 5

# ggplot(preds, aes(x=dbh, y=mean)) +
#     geom_line(aes(colour=clim)) +
#     geom_ribbon(aes(ymin=q2pt5, ymax=q97pt5, fill=clim), alpha=.2) +
#     facet_grid(Panel~Model) +
#     ylab('Increment (cm)') +
#     geom_text(data=ann_text, aes(label=lab))
#
# head(preds)

ps <- foreach(this_panel=unique(preds$Panel), .combine=c) %:%
    foreach(this_model=unique(preds$Model)) %do% {
    p <- ggplot(filter(preds, this_model == Model, this_panel == Panel), aes(x=dbh, y=mean)) +
        geom_line(aes(colour=clim)) +
        geom_ribbon(aes(ymin=q2pt5, ymax=q97pt5, fill=clim), alpha=.2) +
        facet_wrap(~Model) +
        xlab('Initial size (cm)') +
        ylab('Growth increment (cm)') +
        coord_cartesian(ylim=c(0, 6)) +
        theme(legend.position=c(.8, .75)) +
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
