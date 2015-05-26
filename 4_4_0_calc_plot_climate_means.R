# Calculates site-level climate normals for use in predictive models

library(dplyr)
library(lubridate)
library(reshape2)
library(RPostgreSQL)

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
