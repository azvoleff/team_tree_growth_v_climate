library(dplyr)
library(foreach)
library(reshape2)
library(ggplot2)

source("0_settings.R")

img_dpi <- 300

data_folder <- file.path(prefix, "TEAM", "Tree_Growth", "Data")
init_folder <- file.path(prefix, "TEAM", "Tree_Growth", "Initialization")

# model_type <- model_types[1]
# temp_var <- temp_vars[1]
# precip_var <- precip_vars[1]

ret <- foreach (model_type=model_types) %:%
    foreach (temp_var=temp_vars) %:%
        foreach (precip_var=precip_vars,
                 .packages=c("reshape2", "dplyr", "ggplot2"),
                 .inorder=FALSE) %do% {
    suffix <- paste0('_', model_type, '-', temp_var, '-', precip_var)

    sites <- read.csv(file.path(prefix, "TEAM", "Sitecode_Key", "sitecode_key.csv"))

    load(file.path(init_folder, paste0("init_data", suffix, ".RData")))
    load(file.path(data_folder, paste0("model_data_wide", suffix, ".RData")))
    load(file.path(data_folder, paste0("model_data_standardizing", suffix, ".RData")))

    site_ID_factor_key <- read.csv(file.path(data_folder,
                                             paste0("site_ID_factor_key", suffix, 
                                                    ".csv")))
    period_num_factor_key <- read.csv(file.path(data_folder,
                                             paste0("period_num_factor_key", suffix, 
                                                    ".csv")))

    dbh_latent_long <- melt(init_data$dbh_latent, varnames=c("tree_ID", "period_num"), 
                            value.name="dbh_latent_end")
    dbh_latent_long <- dbh_latent_long[complete.cases(dbh_latent_long), ]
    precip_long <- melt(model_data$precip, varnames=c("tree_ID", "period_num"), 
                      value.name="precip")
    merged_data <- merge(dbh_latent_long, precip_long, all.x=TRUE)
    temp_long <- melt(model_data$temp, varnames=c("tree_ID", "period_num"), 
                      value.name="temp")
    merged_data <- merge(dbh_latent_long, temp_long, all.x=TRUE)
    WD <- data.frame(tree_ID=seq(1, model_data$n_tree), WD=model_data$WD)
    merged_data <- merge(merged_data, WD, by="tree_ID")
    site_ID <- data.frame(tree_ID=seq(1, model_data$n_tree), site_ID=model_data$site_ID)
    merged_data <- merge(site_ID, merged_data)
    plot_ID <- data.frame(tree_ID=seq(1, model_data$n_tree), plot_ID=model_data$plot_ID)
    merged_data <- merge(plot_ID, merged_data)
    genus_ID <- data.frame(tree_ID=seq(1, model_data$n_tree), genus_ID=model_data$genus_ID)
    merged_data <- merge(genus_ID, merged_data)
    merged_data <- merge(site_ID_factor_key, merged_data, by.x="site_ID_numeric", 
                         by.y="site_ID")
    merged_data <- merge(period_num_factor_key, merged_data, by.x="period_num_numeric", 
                         by.y="period_num")
    merged_data <- arrange(merged_data, site_ID_char, plot_ID, genus_ID, tree_ID, 
                           period_num_numeric)

    merged_data <- merge(merged_data, sites, by.x="site_ID_char", by.y="sitecode")

    # # Check the number of stems per plot
    # summary_check <- group_by(merged_data, continent, country, sitename_abbrev, 
    #                           plot_ID) %>%
    #     summarise(n=n())

    summary_table <- group_by(merged_data, continent, country, sitename_abbrev) %>%
        summarize(period_first=min(period_num_char),
                  period_last=max(period_num_char),
                  n_plot=length(unique(plot_ID)),
                  n_stem=length(unique(tree_ID)),
                  n_genera=length(unique(genus_ID)))
    write.csv(summary_table, file=paste0("growth_data_summary_sites", suffix, 
                                         ".csv"), row.names=FALSE)

    var_stats <- data.frame(variable=c("dbh", "precip", "temp", "WD"),
               mean=c(dbh_mean, precip_mean, temp_mean, WD_mean),
               min=c(dbh_min, precip_min, temp_min, WD_min),
               max=c(dbh_max, precip_max, temp_max, WD_max),
               sd=c(dbh_sd, precip_sd, temp_sd, WD_sd))
    write.csv(var_stats, file=paste0("growth_data_summary_covariates", suffix, 
                                     ".csv"), row.names=FALSE)
}
