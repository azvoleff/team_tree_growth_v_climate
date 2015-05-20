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


plot_ID_factor_key <- read.csv(file=file.path(data_folder, 
    paste0("plot_ID_factor_key_full-tmn_meanannual-mcwd_run12.csv")))
# Need to key plot_ID to site_ID
site_ID_factor_key <- read.csv(file=file.path(data_folder, 
    paste0("site_ID_factor_key_full-tmn_meanannual-mcwd_run12.csv")))
plot_ID_factor_key$site_ID_char <- gsub('(VG)|([0-9]*)', '', plot_ID_factor_key$plot_ID_char)

site_key <- left_join(plot_ID_factor_key, site_ID_factor_key, by="site_ID_char") %>%
    rename(plot_ID=plot_ID_numeric, site_ID=site_ID_numeric)

n_sites <- length(unique(merged$site_ID))
genus_weights <- group_by(merged, site_ID, genus_ID) %>%
    summarize(n=n()) %>%
    group_by(site_ID) %>%
    mutate(weight=n/sum(n)) %>%
    group_by(genus_ID) %>%
    summarise(weight=sum(weight) / n_sites) %>%
    arrange(desc(weight))


genus_weights$weight_cum <- cumsum(genus_weights$weight)
filter(genus_weights, weight_cum <= .95)

###############################################################################
## Full model (interactions, uncorrelated random effects)

load(file.path(base_folder, "Extracted_Parameters", "parameter_estimates_interact.RData"))

B_g_params <- filter(params, Parameter_Base == 'B_g') %>% 
    rename(genus_ID=row_ID, param_ID=col_ID) %>%
    group_by(Model, genus_ID, param_ID) %>%
    summarise(est=mean(value),
              q2pt5=quantile(value, .025),
              q97pt5=quantile(value, .975))

ggplot(B_g_params) +
    geom_pointrange(aes(genus_ID, est, ymin=q2pt5, ymax=q97pt5)) +
    facet_grid(Model ~ param_ID, scales="free")

B_g_params$Model <- factor(B_g_params$Model, levels=c('tmn', 'tmp', 'tmx'), 
                           labels=c('Minimum~~Temperature~~Model',
                                    'Mean~~Temperature~~Model',
                                    'Maximum~~Temperature~~Model'))
B_g_param_labels <- c('Intercept',
                      'MCWD', 'MCWD^2',
                      'Temp', 'Temp^2',
                      'DBH', 'DBH^2',
                      'DBH%*%MCWD', 'DBH%*%Temp')
# Figure out a plotting order for the data so that it is always organized by 
# increasing value:
B_g_params <- group_by(B_g_params, Model, param_ID) %>%
    arrange(est) %>%
    mutate(plot_order=1:n())
B_g_params$label <- factor(B_g_params$param_ID, levels=1:9, labels=B_g_param_labels)

B_g_params <- left_join(B_g_params, genus_weights)

p <- ggplot(filter(B_g_params, weight_cum <= .95)) +
    theme_bw(base_size=8) +
    geom_segment(aes(x=plot_order, xend=plot_order, y=q2pt5, yend=q97pt5, 
                     colour=weight), alpha=.5, size=.3) +
    geom_point(aes(x=plot_order, est, colour=weight), size=1) +
    facet_grid(label ~ Model, scales="free", labeller=label_parsed) +
    scale_colour_gradient(low="grey", high="blue") +
    theme(plot.background=element_rect(fill='transparent', colour=NA),
          panel.grid.minor.x=element_blank(),
          panel.grid.major.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.x=element_blank(),
          axis.title=element_blank())
ggsave(paste0('genus_level_estimates.png'), width=9,
       height=6.5, dpi=plot_dpi, plot=p)
ggsave(paste0('genus_level_estimates.svg'), width=9,
       height=6.5, dpi=plot_dpi, plot=p)

# Replot but filter to only show genera with over 5% abundance globally
