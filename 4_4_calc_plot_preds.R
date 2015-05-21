# Forest growth model MCMC results

library(dplyr)
library(RColorBrewer)
library(ggmcmc)
library(gridExtra)
library(scales) # for 'alpha'
library(foreach)
library(doParallel)
library(RPostgreSQL)

cl <- makeCluster(20)
registerDoParallel(cl)

plot_width <- 3.5
plot_height <- 3
plot_dpi <- 300

# How many mm of precip per 1 unit change? Convert precip from mm to cm in 
# order to make interpretation of results easier
mm_per_unit <- 10

prefixes <- c('D:/azvoleff/Data', # CI-TEAM
              'H:/Data', # Buffalo drive
              'O:/Data', # Blue drive
              '/localdisk/home/azvoleff/Data') # vertica1
prefix <- prefixes[match(TRUE, unlist(lapply(prefixes, function(x) file_test('-d', x))))]

base_folder <- file.path(prefix, "TEAM", "Tree_Growth")
data_folder <- file.path(base_folder, "Data")

pgsqlpwd <- as.character(read.table('~/pgsqlpwd')[[1]])
pgsqluser <- as.character(read.table('~/pgsqluser')[[1]])

# con <- dbConnect(PostgreSQL(), dbname='tree_growth', user=pgsqluser, 
# password=pgsqlpwd)
# dbSendQuery(con, paste0("DROP TABLE IF EXISTS stems"))
# dbSendQuery(con, paste0("CREATE TABLE stems (site_id integer, plot_id integer, genus_id integer)"))

###############################################################################
## Full model (interactions, uncorrelated random effects)

dbhs <- seq(10, 120, 1)
precips <- c(7.5, 15, 22.5)
temps <- c(10:32)

this_model <- 'tmn'
this_plot <- 1

preds <- foreach(this_model=c('tmn', 'tmp', 'tmx'), .combine=rbind) %:%
    foreach(this_plot=c(1:82), .combine=rbind,
            .packages=c('foreach', 'dplyr', 'RPostgreSQL'),
            .inorder=FALSE) %dopar% {
        pg_src <- src_postgres('tree_growth', user=pgsqluser, 
                               password=pgsqlpwd)
        params <- tbl(pg_src, 'interact')
        params <- filter(params, model == this_model)
        stems <- tbl(pg_src, paste0('stems_', this_model))

        # Genus-level random effects
        B_g_params <- filter(params, parameter_base == 'B_g') %>% 
            rename(genus_id=row_id, param_id=col_id)

        # Site-level intercepts
        int_k_params <- filter(params, parameter_base == 'int_k') %>% 
            rename(site_id=row_id) %>% select(-col_id) %>% collect()

        # Plot-level intercepts
        int_jk_params <- filter(params, parameter_base == 'int_jk') %>% 
            rename(plot_id=row_id) %>% select(-col_id) %>% collect()

        # Elevation difference between plot and CRU cell center
        B_k_params <- filter(params, parameter_base == 'B_k') %>% 
            rename(site_id=row_id) %>% select(-col_id) %>% collect()


        suffix <- paste0("_full-", this_model, "_meanannual-mcwd_run12")
        # Need to key plot_id to site_id
        plot_id_factor_key <- read.csv(file=file.path(data_folder,
                                                      paste0("plot_ID_factor_key", 
                                                             suffix, ".csv")))
        # Need to key plot_id to site_id
        site_id_factor_key <- read.csv(file=file.path(data_folder, 
                                                      paste0("site_ID_factor_key", 
                                                             suffix, ".csv")))
        # Need elevation data
        load(file.path(data_folder, paste0("model_data_elev_key", suffix, ".RData")))

        # Need to center the predictor variables
        load(file.path(data_folder, paste0("model_data_standardizing", suffix, 
                                           ".RData")))

        plot_id_factor_key$site_ID_char <- gsub('(VG)|([0-9]*)', '', plot_id_factor_key$plot_ID_char)

        site_key <- left_join(plot_id_factor_key, rename(elev_key, plot_ID_char=plot_ID)) %>%
            left_join(site_id_factor_key, by="site_ID_char") %>%
            rename(plot_id=plot_ID_numeric, site_id=site_ID_numeric)

        this_site <- site_key$site_id[match(this_plot, site_key$plot_id)]
        this_site_char <- site_key$site_ID_char[match(this_plot, site_key$plot_id)]
        this_plot_char <- site_key$plot_ID_char[match(this_plot, site_key$plot_id)]

        this_elev_diff <- site_key$elev_diff[match(this_plot, site_key$plot_id)]

        genus_weights <- filter(stems, plot_id == this_plot) %>%
            group_by(genus_id) %>%
            summarize(n=n()) %>%
            ungroup() %>%
            mutate(weight=n/sum(n)) %>%
            select(-n)

        B_g_betas <- filter(left_join(genus_weights, B_g_params)) %>%
            mutate(param=sql("parameter_base || '_' || param_id || '_median'")) %>%
            group_by(chain, iteration, param) %>%
            summarise(value=sum(weight * value)) %>%
            ungroup() %>%
            arrange(iteration, param) %>%
            collect()

        # Calculate intercept
        intercept <- filter(int_k_params, site_id == this_site, model == this_model)$value +
            filter(int_jk_params, plot_id == this_plot, model == this_model)$value

        # Calculate elevation adjustment for this plot
        elev_adjust <- filter(B_k_params, site_id == this_site)$value * (this_elev_diff - elev_diff_mean)

        intercept <- intercept + elev_adjust

        # Model effects of temp variation:
        # Convert to matrix for linear algebra
        B <- matrix(B_g_betas$value, nrow=length(unique(B_g_betas$param)), byrow=TRUE)

        dbhs_centered <- dbhs - dbh_mean
        temps_centered <- temps - temp_mean
        # Change MCWD units from cm to mm
        precips_centered <- (precips * mm_per_unit) - precip_mean

        X <- foreach(dbh=dbhs_centered, .combine=rbind) %:%
            foreach(temp=temps_centered, .combine=rbind) %:%
            foreach(precip=precips_centered, .combine=rbind) %do% {
                X <- matrix(1, ncol=1) # genus-level intercept
                X <- cbind(X, precip)
                X <- cbind(X, precip^2)
                X <- cbind(X, temp)
                X <- cbind(X, temp^2)
                X <- cbind(X, dbh)
                X <- cbind(X, dbh^2)
                X <- cbind(X, precip*dbh)
                X <- cbind(X, temp*dbh)
        }
        these_preds <- X %*% B 
        these_preds <- these_preds + matrix(rep(intercept, nrow(X)), 
                                            nrow=nrow(X), byrow=TRUE)
        # X[, 6] (dbh) is subtracted so that output is growth increment
        these_preds <- these_preds - matrix(rep(X[, 6], each=ncol(B)), 
                                            ncol=ncol(B), byrow=TRUE)
        medians <- apply(these_preds, 1, median)
        q2pt5 <- apply(these_preds, 1, quantile, .025)
        q97pt5 <- apply(these_preds, 1, quantile, .975)
        these_preds <- data.frame(model=this_model,
                                  plot_id=this_plot_char, 
                                  site_id=this_site_char,
                                  precip=(X[, 2] + precip_mean)/mm_per_unit,
                                  temp=X[, 4] + temp_mean,
                                  dbh=X[, 6] + dbh_mean,
                                  median=medians,
                                  q2pt5,
                                  q97pt5)

        head(filter(these_preds, precip == 15, temp == 20))

        return(these_preds)
}

save(preds, file='growth_predictions_byplot.RData')

con <- dbConnect(PostgreSQL(), dbname='tree_growth', user=pgsqluser, 
                 password=pgsqlpwd)
if (dbExistsTable(con, "preds_interact")) dbRemoveTable(con, "preds_interact")
dbWriteTable(con, "preds_interact", preds)
dbSendQuery(con, paste0("VACUUM ANALYZE preds_interact;"))
idx_qry <- paste0("CREATE INDEX ON preds_interact (model);",
    "CREATE INDEX ON preds_interact (plot_id);",
    "CREATE INDEX ON preds_interact (site_id);",
    "CREATE INDEX ON preds_interact (precip);",
    "CREATE INDEX ON preds_interact (temp);",
    "CREATE INDEX ON preds_interact (dbh);")
dbSendQuery(con, idx_qry)

###############################################################################
# Now plot results

preds$clim <- factor(sprintf('%.02f', preds$clim))
preds$clim <- relevel(preds$clim, '7.50')
preds$model <- factor(preds$model, levels=c('tmn', 'tmp', 'tmx'), 
                      labels=c('Min. Temp. Model',
                               'Mean Temp. Model',
                               'Max. Temp. Model'))

preds_bysite <- group_by(preds, model, Panel, clim, dbh, site_id) %>%
    summarise(median=median(median),
              q2pt5=median(q2pt5),
              q97pt5=median(q97pt5))

preds_overall <- group_by(preds, model, Panel, clim, dbh) %>%
    summarise(median=median(median),
              q2pt5=median(q2pt5),
              q97pt5=median(q97pt5))

# Make site-by-site plots
foreach(this_site_id=unique(preds_bysite$site_id)) %do% {
    ggplot(filter(preds_bysite, model == "Mean Temp. Model",
                  Panel == "Temperature", site_id == this_site_id),
           aes(x=dbh, y=median)) +
        theme_bw(base_size=8) +
        geom_line(aes(colour=clim)) +
        geom_ribbon(aes(ymin=q2pt5, ymax=q97pt5, fill=clim), alpha=.2) +
        coord_cartesian(ylim=c(-1, 3)) +
        scale_colour_brewer(palette = "Dark2", guide=FALSE) + 
        scale_fill_brewer(palette = "Dark2", guide=FALSE) + 
        ggtitle(this_site_id) +
        theme(legend.position=c(.8, .85),
              plot.title=element_text(size=8),
              panel.background=element_rect(fill='transparent', colour=NA),
              plot.background=element_rect(fill='transparent', colour=NA),
              axis.title=element_blank(),
              plot.margin=unit(c(.3, .4, -.5, -.5), 'lines'),
              panel.grid.minor=element_blank(),
              panel.border=element_blank())
    ggsave(paste0('ArcMap_plot_predgrowth_', this_site_id, '.png'), 
           width=1, height=.75, dpi=300, bg='transparent')
}

ps <- foreach(this_panel=unique(preds_overall$Panel), .combine=c) %:%
    foreach(this_model=unique(preds_overall$model)) %do% {
    p <- ggplot(filter(preds_overall, this_model == model, this_panel == Panel), aes(x=dbh, y=median)) +
        geom_line(aes(colour=clim)) +
        geom_ribbon(aes(ymin=q2pt5, ymax=q97pt5, fill=clim), alpha=.2) +
        facet_wrap(~model) +
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
    p <- ggplot(filter(preds_overall, this_model == model, this_panel == Panel), aes(x=dbh, y=median)) +
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

p <- ggplot(filter(preds_bysite, model == "Mean Temp. Model", Panel == "Temperature"),
           aes(x=dbh, y=median)) +
        theme_bw(base_size=10) +
        geom_line(aes(colour=clim)) +
        geom_ribbon(aes(ymin=q2pt5, ymax=q97pt5, fill=clim), alpha=.2) +
        facet_wrap(~site_id) +
        xlab('Initial size (cm)') +
        ylab('Growth increment (cm)') +
        scale_colour_brewer("Mean Temp.", palette = "Dark2") + 
        scale_fill_brewer("Mean Temp.", palette = "Dark2")
ggsave('predicted_growth_increments_sites.png', width=plot_width*2,
       height=plot_height*2, dpi=plot_dpi, plot=p)

p <- ggplot(filter(preds_bysite, model == "Mean Temp. Model", Panel == 
                   "Temperature", site_id != 'KRP'),
           aes(x=dbh, y=median)) +
        theme_bw(base_size=10) +
        geom_line(aes(colour=clim)) +
        geom_ribbon(aes(ymin=q2pt5, ymax=q97pt5, fill=clim), alpha=.2) +
        facet_wrap(~site_id) +
        xlab('Initial size (cm)') +
        ylab('Growth increment (cm)') +
        coord_cartesian(ylim=c(0, 1.5)) +
        scale_colour_brewer("Mean Temp.", palette = "Dark2") + 
        scale_fill_brewer("Mean Temp.", palette = "Dark2")
ggsave('predicted_growth_increments_sites_mean_temp_noKRP.png', width=plot_width*2,
       height=plot_height*2, dpi=plot_dpi, plot=p)

p <- ggplot(filter(preds_bysite, model == "Min. Temp. Model", Panel == 
                   "Temperature", site_id != 'KRP'),
           aes(x=dbh, y=median)) +
        theme_bw(base_size=10) +
        geom_line(aes(colour=clim)) +
        geom_ribbon(aes(ymin=q2pt5, ymax=q97pt5, fill=clim), alpha=.2) +
        facet_wrap(~site_id) +
        xlab('Initial size (cm)') +
        ylab('Growth increment (cm)') +
        coord_cartesian(ylim=c(0, 1.5)) +
        scale_colour_brewer("Min. Temp.", palette = "Dark2") + 
        scale_fill_brewer("Min. Temp.", palette = "Dark2")
ggsave('predicted_growth_increments_sites_min_temp_noKRP.png', width=plot_width*2,
       height=plot_height*2, dpi=plot_dpi, plot=p)

p <- ggplot(filter(preds_bysite, model == "Max. Temp. Model", Panel == 
                   "Temperature", site_id != 'KRP'),
           aes(x=dbh, y=median)) +
        theme_bw(base_size=10) +
        geom_line(aes(colour=clim)) +
        geom_ribbon(aes(ymin=q2pt5, ymax=q97pt5, fill=clim), alpha=.2) +
        facet_wrap(~site_id) +
        xlab('Initial size (cm)') +
        ylab('Growth increment (cm)') +
        coord_cartesian(ylim=c(0, 1.5)) +
        scale_colour_brewer("Max. Temp.", palette = "Dark2") + 
        scale_fill_brewer("Max. Temp.", palette = "Dark2")
ggsave('predicted_growth_increments_sites_max_temp_noKRP.png', width=plot_width*2,
       height=plot_height*2, dpi=plot_dpi, plot=p)

p <- ggplot(filter(preds_bysite, model == "Min. Temp. Model", Panel == 
                   "MCWD", site_id != 'KRP'),
           aes(x=dbh, y=median)) +
        theme_bw(base_size=10) +
        geom_line(aes(colour=clim)) +
        geom_ribbon(aes(ymin=q2pt5, ymax=q97pt5, fill=clim), alpha=.2) +
        facet_wrap(~site_id) +
        xlab('Initial size (cm)') +
        ylab('Growth increment (cm)') +
        coord_cartesian(ylim=c(0, 1.5)) +
        scale_colour_brewer("MCWD", palette = "Dark2") + 
        scale_fill_brewer("MCWD", palette = "Dark2")
ggsave('predicted_growth_increments_sites_MCWD_noKRP.png', width=plot_width*2,
       height=plot_height*2, dpi=plot_dpi, plot=p)
