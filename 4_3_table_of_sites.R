library(dplyr)
library(reshape2)
library(ggplot2)

img_dpi <- 300

sites <- read.csv('H:/Data/TEAM/Sitecode_Key/sitecode_key.csv')

load("init_data.RData")
load("model_data_wide.RData")
load("model_data_standardizing.RData")

site_ID_key <- read.csv("site_ID_factor_key.csv")
period_ID_key <- read.csv("period_ID_factor_key.csv")

dbh_latent_long <- melt(init_data$dbh_latent, varnames=c("tree_ID", "period_ID"), 
                        value.name="dbh_latent_end")
dbh_latent_long <- dbh_latent_long[complete.cases(dbh_latent_long), ]
mcwd_long <- melt(model_data$mcwd, varnames=c("tree_ID", "period_ID"), 
                  value.name="mcwd")
merged_data <- merge(dbh_latent_long, mcwd_long, all.x=TRUE)
WD <- data.frame(tree_ID=seq(1, model_data$n_tree), WD=model_data$WD)
merged_data <- merge(merged_data, WD, by="tree_ID")
site_ID <- data.frame(tree_ID=seq(1, model_data$n_tree), site_ID=model_data$site_ID)
merged_data <- merge(site_ID, merged_data)
plot_ID <- data.frame(tree_ID=seq(1, model_data$n_tree), plot_ID=model_data$plot_ID)
merged_data <- merge(plot_ID, merged_data)
genus_ID <- data.frame(tree_ID=seq(1, model_data$n_tree), genus_ID=model_data$genus_ID)
merged_data <- merge(genus_ID, merged_data)
merged_data <- merge(site_ID_key, merged_data, by.x="site_ID_numeric", 
                     by.y="site_ID")
merged_data <- merge(period_ID_key, merged_data, by.x="period_ID_numeric", 
                     by.y="period_ID")
merged_data <- arrange(merged_data, site_ID_char, plot_ID, genus_ID, tree_ID, 
                       period_ID_numeric)

merged_data <- merge(merged_data, sites, by.x="site_ID_char", by.y="sitecode")

summary_table <- group_by(merged_data, continent, country, sitename_abbrev) %>%
    summarize(period_first=min(period_ID_char),
              period_last=max(period_ID_char),
              n_plot=length(unique(plot_ID)),
              n_stem=length(unique(tree_ID)),
              n_genera=length(unique(genus_ID)))
write.csv(summary_table, file="growth_data_1ha_plot_summary.csv", row.names=FALSE)

var_stats <- data.frame(variable=c("dbh", "mcwd", "WD"),
           mean=c(dbh_mean, mcwd_mean, WD_mean),
           min=c(dbh_min, mcwd_min, WD_min),
           max=c(dbh_max, mcwd_max, WD_max),
           sd=c(dbh_sd, mcwd_sd, WD_sd))
write.csv(var_stats, file="growth_data_covariate_summary.csv", row.names=FALSE)
