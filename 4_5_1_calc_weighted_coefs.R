# Forest growth model MCMC results
library(RColorBrewer)
library(ggmcmc)
library(gridExtra)
library(foreach)
library(dplyr)
library(matrixStats)
library(RPostgreSQL)
library(doParallel)

cl <- makeCluster(3)
registerDoParallel(cl)

pgsqlpwd <- as.character(read.table('~/pgsqlpwd')[[1]])
pgsqluser <- as.character(read.table('~/pgsqluser')[[1]])

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
        do(ci(.))
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

###############################################################################
## Simple model (no random effects)
##

# load(file.path(base_folder, "Extracted_Parameters", "parameter_estimates_simple.RData"))
# params$model <- factor(params$model, levels=c('tmn', 'tmp', 'tmx'),
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

pg_src <- src_postgres('tree_growth', user=pgsqluser, password=pgsqlpwd)

stems <- tbl(pg_src, paste0('stems'))

# Calculate site-level genus weights, weighting all plots equally within sites
stems %>%
    group_by(model, site_id, plot_id, genus_id) %>%
    # Count num indiv of each genus within plot
    summarize(n=n()) %>%
    ungroup() %>%
    group_by(model, site_id, plot_id) %>%
    # Convert counts to plot-level weights
    mutate(weight=n/sum(n)) %>%
    collect() %>%
    group_by(model, site_id, genus_id) %>%
    # Convert plot-level weights to site-level weights for each genus
    summarise(weight=sum(weight)) %>%
    group_by(model, site_id) %>%
    # Normalize weights by site
    mutate(weight=weight/sum(weight)) -> genus_weights_bysite

# Calculate overall genus weights, weighting all sites equally
genus_weights <- group_by(genus_weights_bysite, genus_id) %>%
    # Sum genus weights across sites
    summarise(weight=sum(weight)) %>%
    # Normalize weights overall
    mutate(weight=weight/sum(weight))

# Check that weighting worked correctly (totals should be three as there are 
# three models:
#group_by(genus_weights, site_id) %>% summarise(sum(weight))

B_g_betas <- foreach(this_model=c('tmn', 'tmp', 'tmx'), .combine=rbind,
                     .packages=c('dplyr', 'RPostgreSQL', 'reshape2', 
                                 'matrixStats'), .inorder=FALSE) %dopar% {
    params <- tbl(src_postgres('tree_growth', user=pgsqluser, 
                               password=pgsqlpwd), 'interact')
    Bs <- filter(params, parameter_base == 'B_g', model == this_model) %>% 
        rename(genus_id=row_id, param_id=col_id) %>%
        collapse() %>%
        mutate(param=sql("parameter_base || '_' || param_id"),
               estimate_id=sql("chain || '_' || iteration")) %>%
        collapse() %>%
        select(model, param, estimate_id, genus_id, value) %>%
        arrange(model, param, estimate_id, genus_id) %>%
        collect() %>%
        dcast(model + param + estimate_id ~ genus_id, value.var="value") %>%
        select(-estimate_id)
    # Remember columns 1 and 2 are the model and parameter IDs
    stopifnot(names(Bs)[c(-1, -2)] == genus_weights$genus_id)
    data.frame(model=Bs$model,
               param=paste0(Bs$param, '_median'),
               value=rowWeightedMedians(as.matrix(Bs[c(-1, -2)]), 
                                        w=genus_weights$weight))
}

save(B_g_betas, 'B_g_betas_weighted_overall.RData')

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

p <- caterpillar(filter(params, param %in% paste0('int_k[', c(1:13), ']')))
ggsave('caterpillar_int_k_interact.png', p, width=plot_width*1.5, 
       height=plot_height*1.5, dpi=plot_dpi)

p <- caterpillar(filter(params, param %in% paste0('int_jk[', c(1:82), ']')))
ggsave('caterpillar_int_jk_interact.png', p, width=plot_width*1.5, 
       height=plot_height*2, dpi=plot_dpi)

p <- caterpillar(filter(params, param %in% paste0('B_k[', c(1:13), ']')))
ggsave('caterpillar_B_k_interact.png', p, width=plot_width*1.5, 
       height=plot_height*1.5, dpi=plot_dpi)

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
        arrange(param, iteration) %>%
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
#        coord_cartesian(ylim=c(0, 3)) +
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
