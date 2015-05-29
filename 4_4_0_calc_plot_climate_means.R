# Calculates site-level climate normals for use in predictive models

library(dplyr)
library(lubridate)
library(reshape2)
library(RPostgreSQL)
library(ggplot2)
library(gridExtra)
library(RColorBrewer)
library(foreach)

plot_width <- 3
plot_height <- 1.5
plot_dpi <- 600

plot_theme <- theme_bw(base_size=8) +
    theme(legend.key.size=unit(.3, 'cm'),
          axis.line=element_line(size=.3, color='black'),
          plot.margin=unit(c(2.5, 2, 0, 0), "mm"))

pgsqlpwd <- as.character(read.table('~/pgsqlpwd')[[1]])
pgsqluser <- as.character(read.table('~/pgsqluser')[[1]])

prefixes <- c('D:/azvoleff/Data', # CI-TEAM
              'H:/Data', # Buffalo drive
              'O:/Data', # Blue drive
              '/localdisk/home/azvoleff/Data') # vertica1
prefix <- prefixes[match(TRUE, unlist(lapply(prefixes, function(x) file_test('-d', x))))]
base_folder <- file.path(prefix, "TEAM", "Tree_Growth")

load('vg_plot_cwds.RData')
cwds <- melt(cwds, id.vars=c('sitecode', 'plot_ID', 'date'))

load('vg_plot_cru.RData')
cru <- select(cru, -plot_num) %>%
    rename(variable=dataset)

clim <- full_join(cru, cwds)

clim$year <- year(clim$date)

clim <- group_by(clim, sitecode, plot_ID, year, variable) %>%
    summarise(value=mean(value))

clim <- filter(clim, year >= 1985, year <= 2014)

clim_stats <- group_by(clim, variable, sitecode, plot_ID, year) %>%
    summarise(plot_annual_mean=mean(value)) %>%
    group_by(variable, sitecode, plot_ID) %>%
    summarise(mean=mean(plot_annual_mean),
              sd=sd(plot_annual_mean)) 

# Add in the site IDs used in the database for the tree_growth paper
tmn_site_key <- read.csv(file.path(base_folder, 'Data', 
                              'site_ID_factor_key_full-tmn_meanannual-mcwd_run12.csv'))
tmp_site_key <- read.csv(file.path(base_folder, 'Data', 
                              'site_ID_factor_key_full-tmp_meanannual-mcwd_run12.csv'))
tmx_site_key <- read.csv(file.path(base_folder, 'Data', 
                              'site_ID_factor_key_full-tmx_meanannual-mcwd_run12.csv'))
stopifnot(tmn_site_key == tmp_site_key)
stopifnot(tmx_site_key == tmp_site_key)
clim_stats <- left_join(rename(tmn_site_key, site_id=site_ID_numeric),
                        rename(clim_stats, site_ID_char=sitecode))

# Add in the plot IDs used in the database for the tree_growth paper
tmn_plot_key <- read.csv(file.path(base_folder, 'Data', 
                              'plot_ID_factor_key_full-tmn_meanannual-mcwd_run12.csv'))
tmp_plot_key <- read.csv(file.path(base_folder, 'Data', 
                              'plot_ID_factor_key_full-tmp_meanannual-mcwd_run12.csv'))
tmx_plot_key <- read.csv(file.path(base_folder, 'Data', 
                              'plot_ID_factor_key_full-tmx_meanannual-mcwd_run12.csv'))
stopifnot(tmn_plot_key == tmp_plot_key)
stopifnot(tmx_plot_key == tmp_plot_key)
clim_stats <- left_join(rename(tmn_plot_key, plot_id=plot_ID_numeric),
                        rename(clim_stats, plot_ID_char=plot_ID))

clim_stats <- rename(clim_stats, site_id_char=site_ID_char, plot_id_char=plot_ID_char)

con <- dbConnect(PostgreSQL(), dbname='tree_growth', user=pgsqluser, password=pgsqlpwd)
if (dbExistsTable(con, 'climate')) dbRemoveTable(con, 'climate')
dbWriteTable(con, 'climate', clim_stats)
dbSendQuery(con, paste0("VACUUM ANALYZE climate;"))
idx_qry <- paste0("CREATE INDEX ON climate (variable);",
    "CREATE INDEX ON climate (site_id);",
    "CREATE INDEX ON climate (plot_id);",
    "CREATE INDEX ON climate (site_id_char);",
    "CREATE INDEX ON climate (plot_id_char);")
dbSendQuery(con, idx_qry)

sites <- read.csv(file.path(prefix, 'TEAM/Sitecode_Key/sitecode_key.csv'))
site_means <- group_by(clim, sitecode, year, variable) %>%
    summarise(value=mean(value)) %>%
    filter(variable %in% c('mcwd_run12', 'tmn', 'tmp', 'tmx'))
site_means <- left_join(site_means, sites)
site_means$year <- as.Date(paste0(site_means$year, '/1/1'))

site_means$variable <- factor(site_means$variable,
                              levels=c('mcwd_run12', 'tmn', 'tmp', 'tmx'),
                              labels=c('MCWD', 'Minimum Temperature', 'Mean Temperature', 'Maximum Temperature'))

mcwd_plot <- ggplot(filter(site_means, variable == 'MCWD', year >= 
                   as.Date('2003/1/1'),
                   !(sitecode %in% c('CSN', 'NAK', 'YAN'))),
       aes(year, value/10, colour=sitename_short, linetype=sitename_short)) +
    geom_line(size=.25) + geom_point(aes(shape=sitename_short), size=1) +
    scale_linetype_manual("Site", values=rep(c(3, 4, 2), each=5, length.out=13)) +
    scale_colour_manual("Site", values=rep(brewer.pal(3, "Dark2"), each=5, length.out=13)) +
    scale_shape_manual("Site", values=rep(c(3, 4, 16, 2, 8), length.out=13)) +
    xlab('Year') +
    ylab('MCWD (cm)') +
    guides(colour=guide_legend(ncol=3),
           linetype=guide_legend(ncol=3), 
           shape=guide_legend(ncol=3)) +
    plot_theme + 
    theme(legend.key.height=unit(.3, 'cm'),
          legend.key.width=unit(.7, 'cm'))
ggsave('climate_mcwd_run12.png', width=plot_width*2,
       height=plot_height, dpi=plot_dpi, plot=mcwd_plot)
ggsave('climate_mcwd_run12.pdf', width=plot_width*2,
       height=plot_height, dpi=plot_dpi, plot=mcwd_plot)
ggsave('climate_mcwd_run12.svg', width=plot_width*2,
       height=plot_height, dpi=plot_dpi, plot=mcwd_plot)

mcwd_plot_nolegend <- mcwd_plot + guides(colour=FALSE, linetype=FALSE, shape=FALSE)
ggsave('climate_mcwd_run12_nolegend.png', width=plot_width*2/3,
       height=plot_height, dpi=plot_dpi, plot=mcwd_plot_nolegend)
ggsave('climate_mcwd_run12_nolegend.pdf', width=plot_width*2/3,
       height=plot_height, dpi=plot_dpi, plot=mcwd_plot_nolegend)
ggsave('climate_mcwd_run12_nolegend.svg', width=plot_width*2/3,
       height=plot_height, dpi=plot_dpi, plot=mcwd_plot_nolegend)

temp_plot <- ggplot(filter(site_means,
                           variable %in% c('Minimum Temperature', 'Mean Temperature', 'Maximum Temperature'),
                           year >= as.Date('2003/1/1'),
                   !(sitecode %in% c('CSN', 'NAK', 'YAN'))),
       aes(year, value, colour=sitename_short, linetype=sitename_short)) +
    geom_line(size=.25) + geom_point(aes(shape=sitename_short), size=1) +
    scale_linetype_manual("Site", values=rep(c(3, 4, 2), each=5, length.out=13)) +
    scale_colour_manual("Site", values=rep(brewer.pal(3, "Dark2"), each=5, length.out=13)) +
    scale_shape_manual("Site", values=rep(c(3, 4, 16, 2, 8), length.out=13)) +
    facet_wrap(~variable) +
    xlab('Year') +
    ylab(expression('Temperature (' * degree * 'C)')) +
    guides(colour=guide_legend(ncol=3),
           linetype=guide_legend(ncol=3), 
           shape=guide_legend(ncol=3)) +
    plot_theme + 
    theme(legend.key.height=unit(.3, 'cm'),
          legend.key.width=unit(.7, 'cm'),
          legend.position='bottom',
          legend.title=element_blank())
ggsave('climate_temp.png', width=plot_width*2,
       height=plot_height*1.75, dpi=plot_dpi, plot=temp_plot)
ggsave('climate_temp.pdf', width=plot_width*2,
       height=plot_height*1.75, dpi=plot_dpi, plot=temp_plot)
ggsave('climate_temp.svg', width=plot_width*2,
       height=plot_height*1.75, dpi=plot_dpi, plot=temp_plot)

temp_plot_nolegend <- temp_plot + guides(colour=FALSE, linetype=FALSE, shape=FALSE)
ggsave('climate_temp_nolegend.png', width=plot_width*2,
       height=plot_height, dpi=plot_dpi, plot=temp_plot_nolegend)
ggsave('climate_temp_nolegend.pdf', width=plot_width*2,
       height=plot_height, dpi=plot_dpi, plot=temp_plot_nolegend)
ggsave('climate_temp_nolegend.svg', width=plot_width*2,
       height=plot_height, dpi=plot_dpi, plot=temp_plot_nolegend)

# Make a combined Temp and MCWD plot
temp_plots <- foreach(this_variable=c('Minimum Temperature', 'Mean Temperature', 'Maximum Temperature')) %do% {
    ggplot(filter(site_means, variable == this_variable,
                  year >= as.Date('2003/1/1'),
                       !(sitecode %in% c('CSN', 'NAK', 'YAN'))),
           aes(year, value, colour=sitename_short, linetype=sitename_short)) +
        geom_line(size=.25) + geom_point(aes(shape=sitename_short), size=1) +
        scale_linetype_manual("Site", values=rep(c(3, 4, 2), each=5, length.out=13)) +
        scale_colour_manual("Site", values=rep(brewer.pal(3, "Dark2"), each=5, length.out=13)) +
        scale_shape_manual("Site", values=rep(c(3, 4, 16, 2, 8), length.out=13)) +
        facet_wrap(~variable) +
        xlab('Year') +
        ylab(expression('Temperature (' * degree * 'C)')) +
        guides(colour=FALSE,
               linetype=FALSE,
               shape=FALSE) +
        plot_theme
}

mcwd_plot <- ggplot(filter(site_means, variable == 'MCWD', year >= 
                   as.Date('2003/1/1'),
                   !(sitecode %in% c('CSN', 'NAK', 'YAN'))),
       aes(year, value/10, colour=sitename_short, linetype=sitename_short)) +
    geom_line(size=.25) + geom_point(aes(shape=sitename_short), size=1) +
    scale_linetype_manual("Site", values=rep(c(3, 4, 2), each=5, length.out=13)) +
    scale_colour_manual("Site", values=rep(brewer.pal(3, "Dark2"), each=5, length.out=13)) +
    scale_shape_manual("Site", values=rep(c(3, 4, 16, 2, 8), length.out=13)) +
    facet_wrap(~variable) +
    xlab('Year') +
    ylab('MCWD (cm)') +
    guides(colour=FALSE,
           linetype=FALSE,
           shape=FALSE) +
    plot_theme

p <- arrangeGrob(temp_plots[[1]], temp_plots[[2]], temp_plots[[3]], 
                 mcwd_plot, nrow=2, ncol=3)
ggsave('climate_combined.png', width=plot_width*2,
       height=plot_height*2, dpi=plot_dpi, plot=p)
ggsave('climate_combined.pdf', width=plot_width*2,
       height=plot_height*2, dpi=plot_dpi, plot=p)
ggsave('climate_combined.svg', width=plot_width*2,
       height=plot_height*2, dpi=plot_dpi, plot=p)
