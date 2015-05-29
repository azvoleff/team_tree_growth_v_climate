library(dplyr)
library(foreach)
library(RColorBrewer)
library(scales) # for 'alpha'
library(ggplot2)
library(gridExtra)
library(RPostgreSQL)

plot_width <- 3
plot_height <- 1.5
plot_dpi <- 600

plot_theme <- theme_bw(base_size=8) +
    theme(legend.key.size=unit(.3, 'cm'),
          axis.line=element_line(size=.3, color='black'),
          plot.margin=unit(c(0, 0, 0, 0), "mm"))

pgsqlpwd <- as.character(read.table('~/pgsqlpwd')[[1]])
pgsqluser <- as.character(read.table('~/pgsqluser')[[1]])

preds <- tbl(src_postgres('tree_growth', user=pgsqluser, password=pgsqlpwd), 
             'preds_interact') %>% collect()

wide_dbh_classes <-  c('[10,20]', '(20,30]', '(30,50]', '(50,80]', '(80,120]')

preds$dbh_type <- 'narrow'
preds$dbh_type[preds$dbh_class %in% wide_dbh_classes] <- 'wide'

preds$model <- factor(preds$model, levels=c('tmn', 'tmp', 'tmx'), 
                      labels=c('Min. Temp.',
                               'Mean Temp.',
                               'Max. Temp.'))

# Fix inaccuracies arising due to rounding
preds <- mutate(preds, temp_diff=round(temp_diff))

# These are the comparison MCWDs that are used in the precip plot
precip_test_levels <- c(7.5, 15, 22.5)
# These are the comparison temps that are used in the temp plot
temp_test_levels <- c(0, 2, 4)

################################################################################
# Plot growth increment vs dbh with lines for precip
preds_overall_pchg <- filter(preds, dbh_type == 'narrow') %>%
    group_by(model, precip, dbh) %>%
    # Filter out the test precips
    filter(precip %in% precip_test_levels) %>%
    # Filter out a set of temp diffs
    filter(temp_diff == 0) %>%
    summarise(median=median(median),
              q2pt5=median(q2pt5),
              q97pt5=median(q97pt5))
preds_overall_pchg$precip <- factor(preds_overall_pchg$precip, 
                                    levels=precip_test_levels, 
                                    labels=paste(precip_test_levels, 'cm'))

p <- ggplot(preds_overall_pchg, aes(x=dbh, y=median)) +
    geom_line(aes(colour=precip, linetype=precip), size=.25) +
    geom_ribbon(aes(ymin=q2pt5, ymax=q97pt5, fill=precip), alpha=.2) +
    facet_wrap(~model) +
    xlab('Initial size (cm)') +
    ylab('Growth (cm)') +
    coord_cartesian(ylim=c(-.1, 1)) +
    scale_linetype_manual("MCWD", values=c(3, 2, 4)) +
    scale_colour_brewer("MCWD", palette = "Dark2") + 
    scale_fill_brewer("MCWD", palette = "Dark2") +
    plot_theme
ggsave('predicted_growth_increments_pchg.png', width=plot_width*1.5,
       height=plot_height, dpi=plot_dpi, plot=p)
ggsave('predicted_growth_increments_pchg.pdf', width=plot_width*1.5,
       height=plot_height, dpi=plot_dpi, plot=p)

################################################################################
# Plot growth increment vs dbh with lines for temperature increase at site 
# level
preds_overall_tchg_site <-  filter(preds, dbh_type == 'narrow') %>%
    group_by(site_id, model, temp_diff, dbh) %>%
    # Filter out precips that aren't the plot-level means
    filter(!(precip %in% precip_test_levels)) %>%
    summarise(median=median(median),
              q2pt5=median(q2pt5),
              q97pt5=median(q97pt5))
preds_overall_tchg_site$temp_diff_factor <- factor(paste0('+', preds_overall_tchg_site$temp_diff, '* degree * C'))
preds_overall_tchg_site$dbh_factor <- factor(preds_overall_tchg_site$dbh)

p <- ggplot(filter(preds_overall_tchg_site, temp_diff %in% temp_test_levels),
            aes(x=dbh, y=median)) +
    geom_line(aes(colour=temp_diff_factor), size=.25) +
    geom_ribbon(aes(ymin=q2pt5, ymax=q97pt5, fill=temp_diff_factor), alpha=.2) +
    facet_grid(site_id~model) +
    xlab('Initial size (cm)') +
    ylab('Growth (cm)') +
    coord_cartesian(ylim=c(-.1, 1)) +
    scale_colour_brewer("Temp.\nchange", labels=parse_format(), palette = "Dark2") + 
    scale_fill_brewer("Temp.\nchange", labels=parse_format(), palette = "Dark2") +
    plot_theme
ggsave('predicted_growth_increments_tchg_bydbh_bysite.png', width=plot_width*3,
       height=plot_height*5, dpi=plot_dpi, plot=p)
ggsave('predicted_growth_increments_tchg_bydbh_bysite.pdf', width=plot_width*3,
       height=plot_height*5, dpi=plot_dpi, plot=p)

################################################################################
# Plot growth increment vs dbh with lines for temperature increase
preds_overall_tchg <- group_by(preds, model, temp_diff, dbh_type, dbh_class, dbh) %>%
    # Filter out precips that aren't the plot-level means
    filter(!(precip %in% precip_test_levels)) %>%
    summarise(median=median(median),
              q2pt5=median(q2pt5),
              q97pt5=median(q97pt5))
preds_overall_tchg$temp_diff_factor <- factor(paste0('+', preds_overall_tchg$temp_diff, '* degree * C'))
preds_overall_tchg$dbh_factor <- factor(preds_overall_tchg$dbh)
preds_overall_tchg_narrow <- filter(preds_overall_tchg, dbh_type == 'narrow')
preds_overall_tchg_wide <- filter(preds_overall_tchg, dbh_type == 'wide')

p <- ggplot(filter(preds_overall_tchg_narrow, temp_diff %in% temp_test_levels),
            aes(x=dbh, y=median)) +
    geom_line(aes(colour=temp_diff_factor, linetype=temp_diff_factor), size=.25) +
    geom_ribbon(aes(ymin=q2pt5, ymax=q97pt5, fill=temp_diff_factor), alpha=.2) +
    facet_wrap(~model) +
    xlab('Initial size (cm)') +
    ylab('Growth (cm)') +
    coord_cartesian(ylim=c(-.1, 1)) +
    scale_linetype_manual("Temp.\nChange", values=c(3, 2, 4), labels=parse_format()) + 
    scale_colour_brewer("Temp.\nChange", labels=parse_format(), palette = "Dark2") + 
    scale_fill_brewer("Temp.\nChange", labels=parse_format(), palette = "Dark2") +
    plot_theme
ggsave('predicted_growth_increments_tchg_bydbh.png', width=plot_width*1.5,
       height=plot_height, dpi=plot_dpi, plot=p)
ggsave('predicted_growth_increments_tchg_bydbh.pdf', width=plot_width*1.5,
       height=plot_height, dpi=plot_dpi, plot=p)

################################################################################
# Plot growth increment vs temperature increase with lines for dbh
preds_overall_tchg_wide$dbh_class <- factor(preds_overall_tchg_wide$dbh_class, levels=wide_dbh_classes)
p <- ggplot(preds_overall_tchg_wide, aes(x=temp_diff, y=median)) +
    geom_line(aes(colour=dbh_class), size=.25) +
    geom_ribbon(aes(ymin=q2pt5, ymax=q97pt5, fill=dbh_class), alpha=.2) +
    facet_wrap(~model) +
    xlab(expression('Temperature change (' * degree * 'C)')) +
    ylab('Growth (cm)') +
    coord_cartesian(ylim=c(-.1, 1)) +
    scale_colour_brewer('Initial\nsize (cm)', palette = "Dark2") + 
    scale_fill_brewer('Initial\nsize (cm)', palette = "Dark2") +
    plot_theme
ggsave('predicted_growth_increments_tchg_byT.png', width=plot_width*1.5,
       height=plot_height, dpi=plot_dpi, plot=p)
ggsave('predicted_growth_increments_tchg_byT.pdf', width=plot_width*1.5,
       height=plot_height, dpi=plot_dpi, plot=p)

################################################################################
# Plot growth increment vs temperature increase with lines for models
preds_overall_tchg_overall <-  filter(preds, dbh_type == 'narrow') %>%
    group_by(model, temp_diff) %>%
    # Filter out precips that aren't the plot-level means
    filter(!(precip %in% precip_test_levels)) %>%
    summarise(median=median(median),
              q2pt5=median(q2pt5),
              q97pt5=median(q97pt5))
preds_overall_tchg_overall$temp_diff_factor <- factor(paste0('+', preds_overall_tchg_overall$temp_diff, '* degree * C'))

t_overall <- ggplot(preds_overall_tchg_overall, aes(x=temp_diff, y=median)) +
    geom_line(aes(colour=model, linetype=model), size=.25) +
    geom_ribbon(aes(ymin=q2pt5, ymax=q97pt5, fill=model), alpha=.2) +
    xlab(expression('Temperature change (' * degree * 'C)')) +
    ylab('Growth (cm)') +
    coord_cartesian(ylim=c(-.1, 1)) +
    scale_linetype_manual("Model", values=c(3, 2, 4)) +
    scale_colour_brewer('Model', palette = "Dark2") + 
    scale_fill_brewer('Model', palette = "Dark2") +
    plot_theme +
    theme(plot.margin=unit(c(2, 0, 0, 0), "mm"))
ggsave('predicted_growth_increments_tchg_byT_overall.png', width=plot_width,
       height=plot_height, dpi=plot_dpi, plot=t_overall)
ggsave('predicted_growth_increments_tchg_byT_overall.pdf', width=plot_width,
       height=plot_height, dpi=plot_dpi, plot=t_overall)

################################################################################
# Plot growth increment vs MCWD with lines for models
preds_overall_pchg_overall <- filter(preds, dbh_type == 'narrow') %>%
    group_by(model, precip) %>%
    # Filter out the test precips
    filter(precip %in% c(0, precip_test_levels)) %>%
    # Filter to include only plot-level mean temps
    filter(temp_diff == 0) %>%
    summarise(median=median(median),
              q2pt5=median(q2pt5),
              q97pt5=median(q97pt5))

p_overall <- ggplot(preds_overall_pchg_overall, aes(x=precip, y=median)) +
    geom_line(aes(colour=model, linetype=model), size=.25) +
    geom_ribbon(aes(ymin=q2pt5, ymax=q97pt5, fill=model), alpha=.2) +
    xlab('MCWD (cm)') +
    ylab('Growth (cm)') +
    coord_cartesian(ylim=c(-.1, 1)) +
    scale_linetype_manual("Model", values=c(3, 2, 4)) +
    scale_colour_brewer("Model", palette = "Dark2") + 
    scale_fill_brewer("Model", palette = "Dark2") +
    plot_theme +
    theme(plot.margin=unit(c(2, 0, 0, 0), "mm"))
ggsave('predicted_growth_increments_pchg_byP_overall.png', width=plot_width,
       height=plot_height, dpi=plot_dpi, plot=p_overall)
ggsave('predicted_growth_increments_pchg_byP_overall.pdf', width=plot_width,
       height=plot_height, dpi=plot_dpi, plot=p_overall)

p <- arrangeGrob(t_overall, p_overall, nrow=1)
ggsave('predicted_growth_increments_overall_combined.png', width=plot_width*2,
       height=plot_height, dpi=plot_dpi, plot=p)
ggsave('predicted_growth_increments_overall_combined.pdf', width=plot_width*2,
       height=plot_height, dpi=plot_dpi, plot=p)
ggsave('predicted_growth_increments_overall_combined.svg', width=plot_width*2,
       height=plot_height, dpi=plot_dpi, plot=p)
