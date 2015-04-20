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
# through BH package.
#Sys.setenv("PKG_LIBS"="-lboost_regex")
# Use version of C++ with regex support:
Sys.setenv("PKG_CXXFLAGS"="-std=c++0x")

sourceCpp('regex_test.cpp')
jags_param_ids(c('asdf', 'fda', 'werq', 'awee'))

sourceCpp('jags_param_ids.cpp')

cl <- makeCluster(3)
registerDoParallel(cl)

plot_width <- 3.5
plot_height <- 3
plot_dpi <- 300

start_val <- 30000
thin_val <- 50

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
destd <- function(d, regex, multiplier) {
    rows <- grepl(regex, d$Parameter)
    d[rows, ]$value_destd <- d[rows, ]$value * multiplier
    return(d)
}

# Function to calculated weighted coefficients from an array of MCMC results.  
# Weights should be a 2 column data.frame with IDs in the first column, and 
# weights in the second
weight_coef <- function(d, w) {
    par_grp_name <- unique(str_extract(d$Parameter, '[a-zA-Z_]*'))
    stopifnot(length(par_grp_name) == 1)
    d$grp_ID <- gsub('[\\[,]*', '', str_extract(d$Parameter, '\\[[0-9]*,'))
    d$par_ID <- gsub('[\\],]*', '', str_extract(d$Parameter, ',[0-9]*\\]'))
    stopifnot(ncol(w) == 2)
    stopifnot(nrow(w) == length(unique(d$grp_ID)))
    names(w) <- c('grp_ID', 'weight')
    d <- left_join(d, w)
    d <- group_by(d, Chain, par_ID) %>%
        summarise(Parameter=paste0(value=mean(value * weight))
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

# Calculate weights for each genus ID (doesn't matter which temp_var is used).  
load(file.path(data_folder, paste0("model_data_wide_full-tmn_meanannual-mcwd_run12.RData")))
merged <- tbl_df(data.frame(site_ID=model_data$site_ID, genus_ID=as.integer(model_data$genus_ID)))
genus_weights <- group_by(merged, genus_ID) %>%
    summarize(n=n()) %>%
    ungroup() %>%
    mutate(weight=n/sum(n)) %>%
    select(-n) %>%
    arrange(desc(weight))

int_mods <- foreach(temp_var=c('tmn', 'tmp', 'tmx'), .combine=rbind, .inorder=FALSE,
                    .packages=c('runjags', 'ggmcmc', 'mcmcplots')) %dopar% {
    load(file.path(prefix, "TEAM", "Tree_Growth", "MCMC_Chains", "interact", 
                   paste0('interact_', temp_var, '.Rdata')))
    plot_mcmc(jags_fit, c("mu_B_g", "sigma_", "int_", "B_k"), paste0("interact - ", temp_var))
    jags_fit <- as.mcmc.list(jags_fit)
    mod <- ggs(window(jags_fit, start=start_val, thin=thin_val))

    # Return variables to original units ("destandardize")
    mod$value_destd <- NA
    load(file.path(data_folder, paste0("model_data_standardizing_full-", 
                                       temp_var, 
                                       "_meanannual-mcwd_run12.RData")))

    mod <- destd(mod, 'B_g\\[([0-9]*,)?1\\]', dbh_sd)
    mod <- destd(mod, 'B_g\\[([0-9]*,)?2\\]', (dbh_sd/precip_sd) * 100) # Convert from mm to 10s of cm
    mod <- destd(mod, 'B_g\\[([0-9]*,)?3\\]', (dbh_sd/precip_sd) * 100) # Convert from mm to 10s of cm
    mod <- destd(mod, 'B_g\\[([0-9]*,)?4\\]', dbh_sd/temp_sd)
    mod <- destd(mod, 'B_g\\[([0-9]*,)?5\\]', dbh_sd/temp_sd)
    mod <- destd(mod, 'B_g\\[([0-9]*,)?6\\]', dbh_sd/dbh_sd)
    mod <- destd(mod, 'B_g\\[([0-9]*,)?7\\]', dbh_sd/dbh_sd)

    mod$Model <- temp_var 
    return(mod)
}
#save(int_mods, file="int_mods.RData")

g_betas <- filter(int_mods, grepl('^B_g\\[[0-9]*,[0-9]*\\]$', Parameter))

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

