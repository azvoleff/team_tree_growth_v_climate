# Forest growth model MCMC results

prefixes <- c('D:/azvoleff/Data', # CI-TEAM
              'H:/Data', # Buffalo drive
              'O:/Data', # Blue drive
              '/localdisk/home/azvoleff/Data') # vertica1
prefix <- prefixes[match(TRUE, unlist(lapply(prefixes, function(x) file_test('-d', x))))]

data_folder <- file.path(prefix, "TEAM", "Tree_Growth", "Data")

library(mcmcplots)
library(runjags)
library(ggmcmc)
library(foreach)
library(stringr)
library(doParallel)
library(ggthemes)
library(coda)
library(Rcpp)

# Below is required for using regex as it must link (so can't work solely 
# through BH package. This is easiest to get working on a Linux machine.
Sys.setenv("PKG_LIBS"="-lboost_regex")
sourceCpp('jags_param_ids.cpp')
# # Below is test code for jags_param_ids function
# jags_param_ids(c('B_g', 'B_g[1]', 'B_g[2,3]', 'B_g[3,4]'))
# jags_param_ids(c('', 'B_g', 'B_g[1]', 'B_g[2,3]', 'B_g[3,4]'))

cl <- makeCluster(3)
registerDoParallel(cl)

plot_width <- 3.5
plot_height <- 3
plot_dpi <- 300

start_val <- 40000
thin_val <- 100

plot_mcmc <- function(fit, params, title) {
    mcmc_results <- as.mcmc.list(fit, c(params))
    mcmc_results <- window(mcmc_results, start=start_val, thin=thin_val)
    mcmcplot(mcmc_results, title=paste0(title))
    #caterplot(mcmc_results)
}

multimodel_caterpillar <- function(mods, labels=NULL) {
    cis <- group_by(mods, Model) %>%
        do(calc_cis(.))
    p <- ggplot(cis, aes(x=median, y=Parameter, colour=Model)) +
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

# Function to "destandardize" standardized coefficients in an array of MCMC 
# results
destd <- function(d, rows, multiplier) {
    d[rows, ]$value_destd <- d[rows, ]$value * multiplier
    return(d)
}

# Function to calculated weighted coefficients from an array of MCMC results.  
# Weights should be a 2 column data.frame with IDs in the first column, and 
# weights in the second. D should be a ggs object with parameter IDs added 
# using the jags_param_ids function
weight_coef <- function(d, w) {
    d <- left_join(d, w)
    d <- group_by(d, Model, Chain, Iteration, param_ID) %>%
        summarise(Parameter=paste(Parameter_Base[1], param_ID[1], 'mean', sep='_'),
                  value=mean(weight * value))
    return(d)
}

# ###############################################################################
# ## Simple model (no random effects)
# ##
#
# plot_simple_model <- function(title) {
#     plot_mcmc(jags_fit, c("B", "sigma_", "int_"), title)
# }
#
# simple_labels_B <- data.frame(Parameter=paste0('B[', 1:10, ']'),
#     Label=c("int",
#             "P" ,
#             "P^2",
#             "T",
#             "T^2",
#             "D",
#             "D^2"))
#
# load(file.path(data_folder, paste0("model_data_standardizing_full-tmn_meanannual-mcwd_run12.RData")))
#
# ### Minimum temperature
# load(file.path(prefix, "TEAM", "Tree_Growth", "MCMC_Chains", "simple", 
#                "simple_tmn.Rdata"))
# plot_simple_model('simple min temp')
# simple_tmn_B_vals <- ggs(as.mcmc.list(jags_fit, "B"), par_labels=simple_labels_B)
# simple_tmn_B_vals$Model <- "Min temp"
#
# ### Mean temperature
# load(file.path(prefix, "TEAM", "Tree_Growth", "MCMC_Chains", "simple", 
#                "simple_tmp.Rdata"))
# plot_simple_model('simple mean temp')
# simple_tmp_B_vals <- ggs(as.mcmc.list(jags_fit, "B"), par_labels=simple_labels_B)
# simple_tmp_B_vals$Model <- "Mean temp"
#
# ### Maximum temperature
# load(file.path(prefix, "TEAM", "Tree_Growth", "MCMC_Chains", "simple", 
#                "simple_tmx.Rdata"))
# plot_simple_model('simple max temp')
# simple_tmx_B_vals <- ggs(as.mcmc.list(jags_fit, "B"), par_labels=simple_labels_B)
# simple_tmx_B_vals$Model <- "Max temp"
#
# simple_B <- rbind(simple_tmn_B_vals, simple_tmp_B_vals, simple_tmx_B_vals)
#
# load(file.path(data_folder, paste0("model_data_standardizing_full-tmn_meanannual-mcwd_run12.RData")))
# simple_B[simple_b$parameter == 'p'] <- simple_B[simple_B$Parameter == 'P'] * (dbh_sd/precip_sd) * 100 # Convert from mm to 10s of cm
# simple_B[simple_B$Parameter == 'P^2'] <- simple_B[simple_B$Parameter == 'P^2'] * (dbh_sd/precip_sd) * 100 # Convert from mm to 10s of cm
# simple_B[simple_B$Parameter == 'T'] <- simple_B[simple_B$Parameter == 'T'] * (dbh_sd/temp_sd)
# simple_B[simple_B$Parameter == 'T^2'] <- simple_B[simple_B$Parameter == 'T^2'] * (dbh_sd/temp_sd)
#
# p <- multimodel_caterpillar(filter(simple_B, Parameter %in% c('P', 'P^2', 'T', 'T^2')),
#                             labels=c('P'='P',
#                                      'P^2'=expression(P^2), 
#                                      'T'='T',
#                                      'T^2'=expression(T^2)))
# ggsave('simple_caterpillar_climate.png', p, width=plot_width, 
#        dpi=plot_dpi)

###############################################################################
## Full model (interactions, uncorrelated random effects)

interact_labels <- data.frame(Label=c("int",
            "P",
            "P^2",
            "T",
            "T^2",
            "D",
            "D^2",
            "D*P",
            "D*T"))
interact_labels$Parameter <- paste0('mu_B_g[', 1:nrow(interact_labels), ']')

int_mods <- foreach(temp_var=c('tmn', 'tmp', 'tmx'), .combine=rbind, .inorder=FALSE,
                    .packages=c('runjags', 'ggmcmc', 'mcmcplots')) %dopar% {
    load(file.path(prefix, "TEAM", "Tree_Growth", "MCMC_Chains", "interact", 
                   paste0('interact_', temp_var, '.RData')))
    jags_fit <- window(as.mcmc.list(jags_fit), start=start_val, thin=thin_val)
    plot_mcmc(jags_fit, c("mu_B_g", "sigma_", "int_", "B_k"), paste0("interact - ", temp_var))
}

# Note the below is not parallel since running it in parallel breaks the Rcpp 
# function.
int_mods <- foreach(temp_var=c('tmn', 'tmp', 'tmx'), .combine=rbind, .inorder=FALSE,
                    .packages=c('runjags', 'ggmcmc', 'mcmcplots', 'Rcpp',
                                'dplyr')) %do% {
    load(file.path(prefix, "TEAM", "Tree_Growth", "MCMC_Chains", "interact", 
                   paste0('interact_', temp_var, '.RData')))
    jags_fit <- window(as.mcmc.list(jags_fit), start=start_val, thin=thin_val)
    mod <- ggs(jags_fit)

    ids <- jags_param_ids(mod$Parameter)
    names(ids) <- c("Parameter_Base", "genus_ID", "param_ID")

    mod <- cbind(mod, ids)

    # Return variables to original units ("destandardize")
    mod$value_destd <- NA
    load(file.path(data_folder, paste0("model_data_standardizing_full-", 
                                       temp_var, 
                                       "_meanannual-mcwd_run12.RData")))

    mod <- destd(mod, mod$param_ID == 1, dbh_sd)
    mod <- destd(mod, mod$param_ID == 2, (dbh_sd/precip_sd) * 100) # Convert from mm to 10s of cm
    mod <- destd(mod, mod$param_ID == 3, (dbh_sd/precip_sd) * 100) # Convert from mm to 10s of cm
    mod <- destd(mod, mod$param_ID == 4, dbh_sd/temp_sd)
    mod <- destd(mod, mod$param_ID == 5, dbh_sd/temp_sd)
    mod <- destd(mod, mod$param_ID == 6, dbh_sd/dbh_sd)
    mod <- destd(mod, mod$param_ID == 7, dbh_sd/dbh_sd)
    # TODO: How to de-standardize interactions?
    #mod <- destd(mod, mod$param_ID == 8, dbh_sd/dbh_sd)
    #mod <- destd(mod, mod$param_ID == 9, dbh_sd/dbh_sd)

    mod$Model <- temp_var 
    return(tbl_df(mod))
}
#save(int_mods, file="int_mods.RData")

# Calculate weights for each genus ID (doesn't matter which temp_var is used 
# since the genus IDs and frequencies are the same across all simulations).
load(file.path(data_folder, paste0("model_data_wide_full-tmn_meanannual-mcwd_run12.RData")))
merged <- tbl_df(data.frame(site_ID=model_data$site_ID, genus_ID=as.integer(model_data$genus_ID)))
genus_weights <- group_by(merged, genus_ID) %>%
    summarize(n=n()) %>%
    ungroup() %>%
    mutate(weight=n/sum(n)) %>%
    select(-n) %>%
    arrange(desc(weight))

g_betas <- weight_coef(filter(int_mods, Parameter_Base == 'B_g'), genus_weights)

clim_betas <- filter(int_mods, Parameter %in% paste0('mu_B_g[', 2:5, ']')) %>%
    mutate(value=value_destd) %>% select(-value_destd)

p <- multimodel_caterpillar(clim_betas,
                            labels=c('mu_B_g[2]'='P',
                                     'mu_B_g[3]'=expression(P^2), 
                                     'mu_B_g[4]'='T',
                                     'mu_B_g[5]'=expression(T^2)))
ggsave('interact_caterpillar_climate.png', p, width=plot_width, height=plot_height, 
       dpi=plot_dpi)

# mu_B_g_vals <- ggs(as.mcmc.list(jags_fit, "mu_B_g"), par_labels=interact_labels_mu_B_g)
# ggs_caterpillar(mu_B_g_vals)
# ggs_caterpillar(mu_B_g_vals) + theme_tufte()
# ggs_caterpillar(mu_B_g_vals) + theme_solarized_2()

###############################################################################
## Full model (correlated random effects)
start_val <- 20000
thin_val <- 10
plot_full_correlated_model <- function(title) {
    plot_mcmc(jags_fit, c("B[1]", "B[2]", "mu_B_g", "sigma_", "int_"), title)
}

correlated_labels_mu_B_g <- data.frame(Parameter=paste0('mu_B_g[', 1:8, ']'),
    Label=c("int",
            "P",
            "P^2",
            "T",
            "T^2",
            "D",
            "D^2",
            "E"))

### Min temperature
load(file.path(prefix, "TEAM", "Tree_Growth", "MCMC_Chains", "correlated", 
               "correlated_tmn.Rdata"))
plot_full_correlated_model("correlated min temp")
correlated_tmn_B_g_vals <- ggs(as.mcmc.list(jags_fit, "mu_B_g"), par_labels=correlated_labels_mu_B_g)

### Mean temperature
load(file.path(prefix, "TEAM", "Tree_Growth", "MCMC_Chains", "correlated", 
               "correlated_tmp.Rdata"))
plot_full_correlated_model("correlated mean temp")
correlated_tmp_B_g_vals <- ggs(as.mcmc.list(jags_fit, "mu_B_g"), par_labels=correlated_labels_mu_B_g)

### Maximum temperature
load(file.path(prefix, "TEAM", "Tree_Growth", "MCMC_Chains", "correlated", 
               "correlated_tmx.Rdata"))
plot_full_correlated_model("correlated max temp")
correlated_tmx_B_g_vals <- ggs(as.mcmc.list(jags_fit, "mu_B_g"), par_labels=correlated_labels_mu_B_g)

correlated_climate <- list("Min temp"=filter(correlated_tmn_B_g_vals, Parameter %in% c('P', 'P^2', 'T', 'T^2')),
                           "Mean temp"=filter(correlated_tmp_B_g_vals, Parameter %in% c('P', 'P^2', 'T', 'T^2')),
                           "Max temp"=filter(correlated_tmx_B_g_vals, Parameter %in% c('P', 'P^2', 'T', 'T^2')))
p <- multimodel_caterpillar(correlated_climate,
                            labels=c('P'='P',
                                     'P^2'=expression(P^2), 
                                     'T'='T',
                                     'T^2'=expression(T^2)))
ggsave('correlated_caterpillar_climate.png', p, width=plot_width, height=plot_height, 
       dpi=plot_dpi)

