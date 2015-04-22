# Forest growth model MCMC results

prefixes <- c('D:/azvoleff/Data', # CI-TEAM
              'H:/Data', # Buffalo drive
              'O:/Data', # Blue drive
              '/localdisk/home/azvoleff/Data') # vertica1
prefix <- prefixes[match(TRUE, unlist(lapply(prefixes, function(x) file_test('-d', x))))]

base_folder <- file.path(prefix, "TEAM", "Tree_Growth")

library(ggmcmc)
library(foreach)

plot_width <- 3.5
plot_height <- 3
plot_dpi <- 300

# Calculate weights for each genus ID (doesn't matter which temp_var is used 
# since the genus IDs and frequencies are the same across all simulations).
load(file.path(data_folder, paste0("model_data_wide_full-tmn_meanannual-mcwd_run12.RData")))
merged <- tbl_df(data.frame(site_ID=model_data$site_ID, 
                            genus_ID=as.integer(model_data$genus_ID)))
genus_weights <- group_by(merged, genus_ID) %>%
    summarize(n=n()) %>%
    ungroup() %>%
    mutate(weight=n/sum(n)) %>%
    select(-n) %>%
    arrange(desc(weight))

multimodel_caterpillar <- function(mods, labels=NULL) {
    cis <- group_by(mods, Model) %>%
        do(ci(.))
    p <- ggplot(cis, aes(x=Parameter, y=median, colour=Model)) +
        theme_bw(base_size=8) +
        geom_point(position=position_dodge(width=.4), size=1.25) +
        geom_linerange(aes(ymin=Low, ymax=High), size=.75, position=position_dodge(width=.4)) +
        geom_linerange(aes(ymin=low, ymax=high), size=.25, position=position_dodge(width=.4)) +
        ylab('Effect (cm)')

    if (!is.null(labels)) {
        p <- p + scale_x_discrete(labels=labels) +
            theme(axis.text.x=element_text(vjust=0))
    }
    p
}

# Function to calculated weighted coefficients from an array of MCMC results.  
# Weights should be a 2 column data.frame with IDs in the first column, and 
# weights in the second. D should be a ggs object with parameter IDs added 
# using the jags_param_ids function
weight_coef <- function(d, w) {
    d <- left_join(d, w)
    d <- group_by(d, Model, Chain, Iteration, param_ID) %>%
        summarise(Parameter=paste(Parameter_Base[1], param_ID[1], 'mean', sep='_'),
                  value=sum(weight * value))
    return(d)
}

###############################################################################
## Simple model (no random effects)
##

load(file.path(base_folder, "Extracted_Parameters", "parameter_estimates_simple.RData"))
params$Model <- factor(params$Model, levels=c('tmn', 'tmp', 'tmx'),
                       labels=c('Min. Temp.', 'Mean Temp.', 'Max. Temp.'))

multimodel_caterpillar(filter(params, Parameter %in% paste0('B[', 1:7, ']')))
ggsave('caterpillar_climate_simple.png', p, width=plot_width, 
       height=plot_height, dpi=plot_dpi)

p <- multimodel_caterpillar(filter(params, Parameter %in% c('B[2]', 'B[3]', 
                                                            'B[4]', 'B[5]')),
                            labels=c('B[2]'='P',
                                     'B[3]'=expression(P^2), 
                                     'B[4]'='T',
                                     'B[5]'=expression(T^2)))
ggsave('caterpillar_climate_simple.png', p, width=plot_width, 
       height=plot_height, dpi=plot_dpi)

###############################################################################
## Full model (interactions, uncorrelated random effects)

load(file.path(base_folder, "Extracted_Parameters", "parameter_estimates_interact.RData"))

multimodel_caterpillar(filter(params, Parameter %in% paste0('mu_B_g[', 1:9, ']')))

p <- multimodel_caterpillar(filter(params, Parameter %in% paste0('mu_B_g[', c(2:5, 8:9), ']')),
                            labels=c('mu_B_g[2]'='P',
                                     'mu_B_g[3]'=expression(P^2),
                                     'mu_B_g[4]'='T',
                                     'mu_B_g[5]'=expression(T^2),
                                     'mu_B_g[8]'=expression(D%*%P),
                                     'mu_B_g[9]'=expression(D%*%T)))
ggsave('caterpillar_climate_interact_unweighted.png', p, width=plot_width, 
       height=plot_height, dpi=plot_dpi)

B_g_betas <- weight_coef(filter(params, Parameter_Base == 'B_g') %>% 
                         rename(genus_ID=row_ID, param_ID=col_ID), 
                         genus_weights)
p <- multimodel_caterpillar(filter(B_g_betas, Parameter %in% paste0('B_g_', c(2:5, 8:9), '_mean')),
                            labels=c('B_g_2_mean'='P',
                                     'B_g_3_mean'=expression(P^2),
                                     'B_g_4_mean'='T',
                                     'B_g_5_mean'=expression(T^2),
                                     'B_g_8_mean'=expression(D%*%P),
                                     'B_g_9_mean'=expression(D%*%T)))
ggsave('caterpillar_climate_interact_weighted.png', p, width=plot_width, 
       height=plot_height, dpi=plot_dpi)

# Plot Rhats for each model
foreach(model=c('tmn', 'tmp', 'tmx')) %do% {
    pars <- filter(params, Model == model)
    attributes(pars)$nChains <- 3
    attributes(pars)$nIterations <- max(params$Iteration)
    ggs_Rhat(pars, 'mu_B_g')
    ggsave(paste0('Rhat_interact_mu_B_g_', model, '.png'), width=plot_width*2, 
           height=plot_height*2, dpi=plot_dpi)
    ggs_Rhat(pars, 'sigma_')
    ggsave(paste0('Rhat_interact_sigma_', model, '.png'), width=plot_width*2, 
           height=plot_height*2, dpi=plot_dpi)
    ggs_Rhat(pars, 'int_')
    ggsave(paste0('Rhat_interact_int_', model, '.png'), width=plot_width*2, 
           height=plot_height*3, dpi=plot_dpi)
}

###############################################################################
## Full model (correlated random effects)

load(file.path(base_folder, "Extracted_Parameters", "parameter_estimates_correlated.RData"))
