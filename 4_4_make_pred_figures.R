library(dplyr)
library(foreach)
library(RColorBrewer)
library(scales) # for 'alpha'
library(ggplot2)
library(gridExtra)
library(RPostgreSQL)

preds <- tbl(src_postgres('tree_growth', user=pgsqluser, password=pgsqlpwd), 
             'preds')

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
