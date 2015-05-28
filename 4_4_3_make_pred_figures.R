library(dplyr)
library(foreach)
library(RColorBrewer)
library(scales) # for 'alpha'
library(ggplot2)
library(gridExtra)
library(RPostgreSQL)

plot_width <- 3
plot_height <- 1.5
plot_dpi <- 300

pgsqlpwd <- as.character(read.table('~/pgsqlpwd')[[1]])
pgsqluser <- as.character(read.table('~/pgsqluser')[[1]])

preds <- tbl(src_postgres('tree_growth', user=pgsqluser, password=pgsqlpwd), 
             'preds_interact') %>% collect()

preds$model <- factor(preds$model, levels=c('tmn', 'tmp', 'tmx'), 
                      labels=c('Min. Temp.',
                               'Mean Temp.',
                               'Max. Temp.'))

# Fix inaccuracies due to roounding
preds <- mutate(preds, temp_diff=round(temp_diff))

preds_bysite <- group_by(preds, model, site_id, precip, temp_diff, dbh) %>%
    summarise(median=median(median),
              q2pt5=median(q2pt5),
              q97pt5=median(q97pt5))

preds_overall <- group_by(preds_bysite, model, precip, temp_diff, dbh) %>%
    summarise(median=median(median),
              q2pt5=median(q2pt5),
              q97pt5=median(q97pt5))

################################################################################
# Plot growth increment vs dbh with lines for precip
preds_overall_pchg <- group_by(preds, model, precip, dbh) %>%
    # Filter out precips that aren't the plot-level means
    filter(precip %in% c(7.5, 15, 22.5)) %>%
    # Filter out a set of temp diffs
    filter(temp_diff == 0) %>%
    summarise(median=median(median),
              q2pt5=median(q2pt5),
              q97pt5=median(q97pt5))
# preds_overall_pchg$temp_diff <- factor(sprintf('%.02f', preds_overall_pchg$temp_diff))
preds_overall_pchg$precip <- factor(preds_overall_pchg$precip)

p <- ggplot(preds_overall_pchg, aes(x=dbh, y=median)) +
    theme_bw(base_size=8) +
    geom_line(aes(colour=precip), size=.25) +
    geom_ribbon(aes(ymin=q2pt5, ymax=q97pt5, fill=precip), alpha=.2) +
    facet_wrap(~model) +
    xlab('Initial size (cm)') +
    ylab('Growth increment (cm)') +
    coord_cartesian(ylim=c(-.1, 1)) +
    scale_colour_brewer("MCWD", palette = "Dark2") + 
    scale_fill_brewer("MCWD", palette = "Dark2") +
    theme(legend.key.size=unit(.3, 'cm'))
          #legend.title=element_text(size=6, face='plain'))
          # legend.direction='horizontal',
          # legend.position=c(.1, .7))
ggsave('predicted_growth_increments_pchg.png', width=plot_width*1.5,
       height=plot_height, dpi=plot_dpi, plot=p)

################################################################################
# Plot growth increment vs dbh with lines for temperature increase at site 
# level
preds_overall_tchg_site <- group_by(preds, site_id, model, temp_diff, dbh) %>%
    # Filter out precips that aren't the plot-level means
    filter(!(precip %in% c(7.5, 15, 22.5))) %>%
    summarise(median=median(median),
              q2pt5=median(q2pt5),
              q97pt5=median(q97pt5))
# preds_overall_tchg_site$temp_diff <- factor(sprintf('%.02f', preds_overall_tchg_site$temp_diff))
preds_overall_tchg_site$temp_diff_factor <- factor(paste0('+', preds_overall_tchg_site$temp_diff, '* degree * C'))
preds_overall_tchg_site$dbh_factor <- factor(preds_overall_tchg_site$dbh)

p <- ggplot(filter(preds_overall_tchg_site, temp_diff %in% c(0, 2, 4)),
            aes(x=dbh, y=median)) +
    theme_bw(base_size=8) +
    geom_line(aes(colour=temp_diff_factor), size=.25) +
    geom_ribbon(aes(ymin=q2pt5, ymax=q97pt5, fill=temp_diff_factor), alpha=.2) +
    facet_grid(site_id~model) +
    xlab('Initial size (cm)') +
    ylab('Growth increment (cm)') +
    coord_cartesian(ylim=c(-.1, 1)) +
    theme(legend.key.size=unit(.3, 'cm')) +
    scale_colour_brewer("Temp.\nChange", labels=parse_format(), palette = "Dark2") + 
    scale_fill_brewer("Temp.\nChange", labels=parse_format(), palette = "Dark2")
ggsave('predicted_growth_increments_tchg_bydbh_bysite.png', width=plot_width*3,
       height=plot_height*5, dpi=plot_dpi, plot=p)

################################################################################
# Plot growth increment vs dbh with lines for temperature increase
preds_overall_tchg <- group_by(preds, model, temp_diff, dbh) %>%
    # Filter out precips that aren't the plot-level means
    filter(!(precip %in% c(7.5, 15, 22.5))) %>%
    summarise(median=median(median),
              q2pt5=median(q2pt5),
              q97pt5=median(q97pt5))
# preds_overall_tchg$temp_diff <- factor(sprintf('%.02f', preds_overall_tchg$temp_diff))
preds_overall_tchg$temp_diff_factor <- factor(paste0('+', preds_overall_tchg$temp_diff, '* degree * C'))
preds_overall_tchg$dbh_factor <- factor(preds_overall_tchg$dbh)

p <- ggplot(filter(preds_overall_tchg, temp_diff %in% c(0, 2, 4)),
            aes(x=dbh, y=median)) +
    theme_bw(base_size=8) +
    geom_line(aes(colour=temp_diff_factor), size=.25) +
    geom_ribbon(aes(ymin=q2pt5, ymax=q97pt5, fill=temp_diff_factor), alpha=.2) +
    facet_wrap(~model) +
    xlab('Initial size (cm)') +
    ylab('Growth increment (cm)') +
    coord_cartesian(ylim=c(-.1, 1)) +
    theme(legend.key.size=unit(.3, 'cm')) +
    scale_colour_brewer("Temp.\nChange", labels=parse_format(), palette = "Dark2") + 
    scale_fill_brewer("Temp.\nChange", labels=parse_format(), palette = "Dark2")
ggsave('predicted_growth_increments_tchg_bydbh.png', width=plot_width*1.5,
       height=plot_height, dpi=plot_dpi, plot=p)

################################################################################
# Plot growth increment vs temperature increase with lines for dbh

p <- ggplot(filter(preds_overall_tchg, dbh_factor %in% c(12.5, 22.5, 55)), 
            aes(x=temp_diff, y=median)) +
    theme_bw(base_size=8) +
    geom_line(aes(colour=dbh_factor), size=.25) +
    geom_ribbon(aes(ymin=q2pt5, ymax=q97pt5, fill=dbh_factor), alpha=.2) +
    facet_wrap(~model) +
    xlab(expression('Temp. Change (' * degree * 'C)')) +
    ylab('Growth increment (cm)') +
    coord_cartesian(ylim=c(-.1, 1)) +
    theme(legend.key.size=unit(.3, 'cm')) +
    scale_colour_brewer('Initial size (cm)', palette = "Dark2") + 
    scale_fill_brewer('Initial size (cm)', palette = "Dark2")
ggsave('predicted_growth_increments_tchg_byT.png', width=plot_width*1.5,
       height=plot_height, dpi=plot_dpi, plot=p)

preds_overall_tchg_overall <- group_by(preds, model, temp_diff) %>%
    # Filter out precips that aren't the plot-level means
    filter(!(precip %in% c(7.5, 15, 22.5))) %>%
    summarise(median=median(median),
              q2pt5=median(q2pt5),
              q97pt5=median(q97pt5))
# preds_overall_tchg$temp_diff <- factor(sprintf('%.02f', preds_overall_tchg$temp_diff))
preds_overall_tchg_overall$temp_diff_factor <- factor(paste0('+', preds_overall_tchg_overall$temp_diff, '* degree * C'))

p <- ggplot(preds_overall_tchg_overall, aes(x=temp_diff, y=median)) +
    theme_bw(base_size=8) +
    geom_line(aes(colour=model), size=.25) +
    geom_ribbon(aes(ymin=q2pt5, ymax=q97pt5, fill=model), alpha=.2) +
    xlab(expression('Temp. Change (' * degree * 'C)')) +
    ylab('Growth increment (cm)') +
    coord_cartesian(ylim=c(-.1, 1)) +
    theme(legend.key.size=unit(.3, 'cm')) +
    scale_colour_brewer('Model', palette = "Dark2") + 
    scale_fill_brewer('Model', palette = "Dark2")
ggsave('predicted_growth_increments_tchg_byT_overall.png', width=plot_width,
       height=plot_height, dpi=plot_dpi, plot=p)
