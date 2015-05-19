# Forest growth model MCMC results

prefixes <- c('D:/azvoleff/Data', # CI-TEAM
              'H:/Data', # Buffalo drive
              'O:/Data', # Blue drive
              '/localdisk/home/azvoleff/Data') # vertica1
prefix <- prefixes[match(TRUE, unlist(lapply(prefixes, function(x) file_test('-d', x))))]

base_folder <- file.path(prefix, "TEAM", "Tree_Growth")
data_folder <- file.path(base_folder, "Data")

library(RColorBrewer)
library(ggmcmc)
library(gridExtra)
library(scales) # for 'alpha'
library(foreach)
library(doParallel)

cl <- makeCluster(3)
registerDoParallel(cl)

plot_width <- 3.5
plot_height <- 3
plot_dpi <- 300

# How many mm of precip per 1 unit change? Convert precip from mm to cm in 
# order to make interpretation of results easier
mm_per_unit <- 10

# Calculate weights for each genus ID (doesn't matter which temp_var is used 
# since the genus IDs and frequencies are the same across all simulations).
load(file.path(data_folder,
               paste0("model_data_wide_full-tmn_meanannual-mcwd_run12.RData")))
merged <- tbl_df(data.frame(site_ID=model_data$site_ID, 
                            plot_ID=model_data$plot_ID, 
                            genus_ID=model_data$genus_ID))

# Function to calculated weighted coefficients from an array of MCMC results.  
# Weights should be a 2 column data.frame with IDs in the first column, and 
# weights in the second. D should be a ggs object with parameter IDs added 
# using the jags_param_ids function
weight_coef <- function(d, w) {
    left_join(d, w) %>%
        group_by(Model, Chain, Iteration, param_ID) %>%
        summarise(Parameter=paste(Parameter_Base[1], param_ID[1], 'median', sep='_'),
                  value=sum(weight * value))
}

###############################################################################
## Full model (interactions, uncorrelated random effects)

load(file.path(base_folder, "Extracted_Parameters", "parameter_estimates_interact.RData"))

# Function to make a growth prediction as a function of a particular climate 
# variable
#
# TODO: Test the intercept addition to make sure I did the rep correctly
make_dbh_preds <- function(B, temp, precip, dbhs, intercept) {
    X <- matrix(rep(1, length(dbhs)), ncol=1) # genus-level intercept
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
    these_preds <- these_preds + matrix(rep(intercept, length(dbhs)), 
                                        ncol=ncol(these_preds), byrow=TRUE)
    medians <- apply(these_preds, 1, median)
    q2pt5 <- apply(these_preds, 1, quantile, .025)
    q97pt5 <- apply(these_preds, 1, quantile, .975)
    data.frame(median=medians, q2pt5, q97pt5)
}

make_temp_preds <- function(B, temps, precip, dbh, intercept) {
    X <- matrix(rep(1, length(temps)), ncol=1) # genus-level intercept
    X <- cbind(X, precip)
    X <- cbind(X, precip^2)
    X <- cbind(X, temps)
    X <- cbind(X, temps^2)
    X <- cbind(X, dbh)
    X <- cbind(X, dbh^2)
    X <- cbind(X, precip*dbh)
    X <- cbind(X, temps*dbh)
    # dbh is subtracted so that output is growth increment
    these_preds <- X %*% B - rep(dbh, ncol(B))
    these_preds <- these_preds + matrix(rep(intercept, length(temps)), 
                                        ncol=ncol(these_preds), byrow=TRUE)
    medians <- apply(these_preds, 1, median)
    q2pt5 <- apply(these_preds, 1, quantile, .025)
    q97pt5 <- apply(these_preds, 1, quantile, .975)
    data.frame(median=medians, q2pt5, q97pt5)
}

# Genus-level random effects
B_g_params <- filter(params, Parameter_Base == 'B_g') %>% 
    rename(genus_ID=row_ID, param_ID=col_ID)

# Site-level intercepts
int_k_params <- filter(params, Parameter_Base == 'int_k') %>% 
    rename(site_ID=row_ID) %>% select(-col_ID)


# Genus-level random effects
B_g_params <- filter(params, Parameter_Base == 'B_g') %>% 
    rename(genus_ID=row_ID, param_ID=col_ID)

# Site-level intercepts
int_k_params <- filter(params, Parameter_Base == 'int_k') %>% 
    rename(site_ID=row_ID) %>% select(-col_ID)

# Plot-level intercepts
int_jk_params <- filter(params, Parameter_Base == 'int_jk') %>% 
    rename(plot_ID=row_ID) %>% select(-col_ID)

# Elevation difference between plot and CRU cell center
B_k_params <- filter(params, Parameter_Base == 'B_k') %>% 
    rename(site_ID=row_ID) %>% select(-col_ID)

dbhs <- seq(10, 120, 1)

# this_model <- unique(params$Model)[1]
# this_plot <- unique(merged$plot_ID)[1]

preds <- foreach(this_model=unique(params$Model), .combine=rbind, 
                 .final=tbl_df) %:%
    foreach(this_plot=unique(merged$plot_ID), .combine=rbind,
            .packages=c('foreach', 'dplyr')) %dopar% {

    suffix <- paste0("_full-", this_model, "_meanannual-mcwd_run12")
    # Need to key plot_ID to site_ID
    plot_ID_factor_key <- read.csv(file=file.path(data_folder,
                                                  paste0("plot_ID_factor_key", 
                                                         suffix, ".csv")))
    # Need to key plot_ID to site_ID
    site_ID_factor_key <- read.csv(file=file.path(data_folder, 
                                                  paste0("site_ID_factor_key", 
                                                         suffix, ".csv")))
    # Need elevation data
    load(file.path(data_folder, paste0("model_data_elev_key", suffix, ".RData")))

    # Need to center the predictor variables
    load(file.path(data_folder, paste0("model_data_standardizing", suffix, 
                                       ".RData")))

    plot_ID_factor_key$site_ID_char <- gsub('(VG)|([0-9]*)', '', plot_ID_factor_key$plot_ID_char)

    site_key <- left_join(plot_ID_factor_key, rename(elev_key, plot_ID_char=plot_ID)) %>%
        left_join(site_ID_factor_key, by="site_ID_char") %>%
        rename(plot_ID=plot_ID_numeric, site_ID=site_ID_numeric)

    this_site <- site_key$site_ID[match(this_plot, site_key$plot_ID)]
    this_site_char <- site_key$site_ID_char[match(this_plot, site_key$plot_ID)]
    this_plot_char <- site_key$plot_ID_char[match(this_plot, site_key$plot_ID)]

    this_elev_diff <- site_key$elev_diff[match(this_plot, site_key$plot_ID)]

    genus_weights <- filter(merged, plot_ID == this_plot) %>%
        group_by(genus_ID) %>%
        summarize(n=n()) %>%
        ungroup() %>%
        mutate(weight=n/sum(n)) %>%
        select(-n) %>%
        arrange(desc(weight))

    B_g_betas <- weight_coef(filter(B_g_params, genus_ID %in% 
                                    genus_weights$genus_ID,
                                    Model == this_model), genus_weights)

    # Change MCWD units from mm to cm
    B_g_betas$value[B_g_betas$Parameter == "B_g_2_median"] <- B_g_betas$value[B_g_betas$Parameter == "B_g_2_median"] * mm_per_unit
    B_g_betas$value[B_g_betas$Parameter == "B_g_3_median"] <- B_g_betas$value[B_g_betas$Parameter == "B_g_3_median"] * mm_per_unit^2
    B_g_betas$value[B_g_betas$Parameter == "B_g_8_median"] <- B_g_betas$value[B_g_betas$Parameter == "B_g_8_median"] * mm_per_unit

    # Calculate intercept
    intercept <- filter(int_k_params, site_ID == this_site, Model == this_model)$value +
        filter(int_jk_params, plot_ID == this_plot, Model == this_model)$value

    # Calculate elevation adjustment for this plot
    elev_adjust <- filter(B_k_params, site_ID == this_site, Model == this_model)$value * (this_elev_diff - elev_diff_mean)

    intercept <- intercept + elev_adjust

    # Model effects of temp variation:
    B <- filter(B_g_betas, Model == this_model) %>%
        ungroup() %>%
        arrange(Parameter, Iteration) %>%
        select(Parameter, value)
    # Convert to matrix for linear algebra
    B <- matrix(B$value, nrow=length(unique(B$Parameter)), byrow=TRUE)

    dbhs_centered <- dbhs - dbh_mean

    temps <- round(temp_mean) - temp_mean
    temps <- c(temps - 2, temps, temps + 2)
    temp_preds <- foreach(temp=temps, .combine=rbind) %do% {
        temp_preds <- make_dbh_preds(B, temp, (150 - precip_mean)/mm_per_unit, 
                                 dbhs_centered, intercept)
        temp_preds <- cbind(Panel="Temperature", clim=temp + temp_mean, 
                            dbh=dbhs,
                            temp_preds)
        return(temp_preds)
    }

    precips <- (c(75, 150, 150 + 75) - precip_mean)/mm_per_unit
    precip_preds <- foreach(precip=precips, .combine=rbind) %do% {
        precip_preds <- make_dbh_preds(B, 0, precip, dbhs_centered, intercept)
        precip_preds <- cbind(Panel="MCWD",
                              clim=precip + precip_mean/mm_per_unit,
                              dbh=dbhs,
                              precip_preds)
        return(precip_preds)
    }

    preds <- rbind(temp_preds, precip_preds)
    preds <- cbind(Model=this_model, plot_ID=this_plot_char, 
                   site_ID=this_site_char, preds)

    return(preds)
}

preds$clim <- factor(sprintf('%.02f', preds$clim))
preds$clim <- relevel(preds$clim, '7.50')
preds$Model <- factor(preds$Model, levels=c('tmn', 'tmp', 'tmx'), 
                      labels=c('Min. Temp. Model',
                               'Mean Temp. Model',
                               'Max. Temp. Model'))

save(preds, file='growth_predictions_byplot.RData')

preds_bysite <- group_by(preds, Model, Panel, clim, dbh, site_ID) %>%
    summarise(median=median(median),
              q2pt5=median(q2pt5),
              q97pt5=median(q97pt5))

preds_overall <- group_by(preds, Model, Panel, clim, dbh) %>%
    summarise(median=median(median),
              q2pt5=median(q2pt5),
              q97pt5=median(q97pt5))

# Make site-by-site plots
foreach(this_site_ID=unique(preds_bysite$site_ID)) %do% {
    ggplot(filter(preds_bysite, Model == "Mean Temp. Model",
                  Panel == "Temperature", site_ID == this_site_ID),
           aes(x=dbh, y=median)) +
        theme_bw(base_size=8) +
        geom_line(aes(colour=clim)) +
        geom_ribbon(aes(ymin=q2pt5, ymax=q97pt5, fill=clim), alpha=.2) +
        coord_cartesian(ylim=c(-1, 3)) +
        scale_colour_brewer(palette = "Dark2", guide=FALSE) + 
        scale_fill_brewer(palette = "Dark2", guide=FALSE) + 
        ggtitle(this_site_ID) +
        theme(legend.position=c(.8, .85),
              plot.title=element_text(size=8),
              panel.background=element_rect(fill='transparent', colour=NA),
              plot.background=element_rect(fill='transparent', colour=NA),
              axis.title=element_blank(),
              plot.margin=unit(c(.3, .4, -.5, -.5), 'lines'),
              panel.grid.minor=element_blank(),
              panel.border=element_blank())
    ggsave(paste0('ArcMap_plot_predgrowth_', this_site_ID, '.png'), 
           width=1, height=.75, dpi=300, bg='transparent')
}

ps <- foreach(this_panel=unique(preds_overall$Panel), .combine=c) %:%
    foreach(this_model=unique(preds_overall$Model)) %do% {
    p <- ggplot(filter(preds_overall, this_model == Model, this_panel == Panel), aes(x=dbh, y=median)) +
        geom_line(aes(colour=clim)) +
        geom_ribbon(aes(ymin=q2pt5, ymax=q97pt5, fill=clim), alpha=.2) +
        facet_wrap(~Model) +
        xlab('Initial size (cm)') +
        ylab('Growth increment (cm)') +
        coord_cartesian(ylim=c(0, 1)) +
        theme(legend.position=c(.8, .7)) +
        guides(fill=guide_legend(this_panel),
               colour=guide_legend(this_panel))
    return(p)
}

p <- arrangeGrob(ps[[1]], ps[[2]], ps[[3]],
                 ps[[4]], ps[[5]], ps[[6]], ncol=3, nrow=2)
ggsave('predicted_growth_increments_frombyplot.png', width=plot_width*2,
       height=plot_height*2, dpi=plot_dpi, plot=p)
ggsave('predicted_growth_increments_frombyplot.svg', width=plot_width*2,
       height=plot_height*2, dpi=plot_dpi, plot=p)

panels <- c(rep('Temperature', 3), 'MCWD')
models <- c('Min. Temp. Model', 'Mean Temp. Model',
            'Max. Temp. Model', 'Min. Temp. Model')
legend_titles <- c('Min. Temp.', 'Mean Temp.', 'Max. Temp', 'MCWD')
ps <- foreach(this_panel=panels, this_model=models,
              legend_title=legend_titles) %do% {
    p <- ggplot(filter(preds_overall, this_model == Model, this_panel == Panel), aes(x=dbh, y=median)) +
        theme_bw(base_size=10) +
        geom_line(aes(colour=clim)) +
        geom_ribbon(aes(ymin=q2pt5, ymax=q97pt5, fill=clim), alpha=.2) +
        geom_line(aes(y=q2pt5, colour=clim), alpha=.2) +
        geom_line(aes(y=q97pt5, colour=clim), alpha=.2) +
        xlab('Initial size (cm)') +
        ylab('Growth increment (cm)') +
        coord_cartesian(ylim=c(0, 1)) +
        theme(legend.position=c(.8, .7)) +
        guides(fill=guide_legend(legend_title),
               colour=guide_legend(legend_title)) +
        scale_colour_brewer(palette = "Dark2", guide=FALSE) + 
        scale_fill_brewer(palette = "Dark2", guide=FALSE)
    return(p)
}
p <- arrangeGrob(ps[[1]], ps[[2]], ps[[3]], ps[[4]], ncol=2, nrow=2)
ggsave('predicted_growth_increments_frombyplot_2x2.png', width=plot_width*2,
       height=plot_height*2, dpi=plot_dpi, plot=p)
ggsave('predicted_growth_increments_frombyplot_2x2.svg', width=plot_width*2,
       height=plot_height*2, dpi=plot_dpi, plot=p)

p <- ggplot(filter(preds_bysite, Model == "Mean Temp. Model", Panel == "Temperature"),
           aes(x=dbh, y=median)) +
        theme_bw(base_size=10) +
        geom_line(aes(colour=clim)) +
        geom_ribbon(aes(ymin=q2pt5, ymax=q97pt5, fill=clim), alpha=.2) +
        facet_wrap(~site_ID) +
        xlab('Initial size (cm)') +
        ylab('Growth increment (cm)') +
        scale_colour_brewer("Mean Temp.", palette = "Dark2") + 
        scale_fill_brewer("Mean Temp.", palette = "Dark2")
ggsave('predicted_growth_increments_sites.png', width=plot_width*2,
       height=plot_height*2, dpi=plot_dpi, plot=p)

p <- ggplot(filter(preds_bysite, Model == "Mean Temp. Model", Panel == 
                   "Temperature", site_ID != 'KRP'),
           aes(x=dbh, y=median)) +
        theme_bw(base_size=10) +
        geom_line(aes(colour=clim)) +
        geom_ribbon(aes(ymin=q2pt5, ymax=q97pt5, fill=clim), alpha=.2) +
        facet_wrap(~site_ID) +
        xlab('Initial size (cm)') +
        ylab('Growth increment (cm)') +
        coord_cartesian(ylim=c(0, 1.5)) +
        scale_colour_brewer("Mean Temp.", palette = "Dark2") + 
        scale_fill_brewer("Mean Temp.", palette = "Dark2")
ggsave('predicted_growth_increments_sites_mean_temp_noKRP.png', width=plot_width*2,
       height=plot_height*2, dpi=plot_dpi, plot=p)

p <- ggplot(filter(preds_bysite, Model == "Min. Temp. Model", Panel == 
                   "Temperature", site_ID != 'KRP'),
           aes(x=dbh, y=median)) +
        theme_bw(base_size=10) +
        geom_line(aes(colour=clim)) +
        geom_ribbon(aes(ymin=q2pt5, ymax=q97pt5, fill=clim), alpha=.2) +
        facet_wrap(~site_ID) +
        xlab('Initial size (cm)') +
        ylab('Growth increment (cm)') +
        coord_cartesian(ylim=c(0, 1.5)) +
        scale_colour_brewer("Min. Temp.", palette = "Dark2") + 
        scale_fill_brewer("Min. Temp.", palette = "Dark2")
ggsave('predicted_growth_increments_sites_min_temp_noKRP.png', width=plot_width*2,
       height=plot_height*2, dpi=plot_dpi, plot=p)

p <- ggplot(filter(preds_bysite, Model == "Max. Temp. Model", Panel == 
                   "Temperature", site_ID != 'KRP'),
           aes(x=dbh, y=median)) +
        theme_bw(base_size=10) +
        geom_line(aes(colour=clim)) +
        geom_ribbon(aes(ymin=q2pt5, ymax=q97pt5, fill=clim), alpha=.2) +
        facet_wrap(~site_ID) +
        xlab('Initial size (cm)') +
        ylab('Growth increment (cm)') +
        coord_cartesian(ylim=c(0, 1.5)) +
        scale_colour_brewer("Max. Temp.", palette = "Dark2") + 
        scale_fill_brewer("Max. Temp.", palette = "Dark2")
ggsave('predicted_growth_increments_sites_max_temp_noKRP.png', width=plot_width*2,
       height=plot_height*2, dpi=plot_dpi, plot=p)

p <- ggplot(filter(preds_bysite, Model == "Min. Temp. Model", Panel == 
                   "MCWD", site_ID != 'KRP'),
           aes(x=dbh, y=median)) +
        theme_bw(base_size=10) +
        geom_line(aes(colour=clim)) +
        geom_ribbon(aes(ymin=q2pt5, ymax=q97pt5, fill=clim), alpha=.2) +
        facet_wrap(~site_ID) +
        xlab('Initial size (cm)') +
        ylab('Growth increment (cm)') +
        coord_cartesian(ylim=c(0, 1.5)) +
        scale_colour_brewer("MCWD", palette = "Dark2") + 
        scale_fill_brewer("MCWD", palette = "Dark2")
ggsave('predicted_growth_increments_sites_MCWD_noKRP.png', width=plot_width*2,
       height=plot_height*2, dpi=plot_dpi, plot=p)
