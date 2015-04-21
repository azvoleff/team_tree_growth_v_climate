# This script needs to be run a machine with the C++ Boost libraries installed

prefixes <- c('D:/azvoleff/Data', # CI-TEAM
              'H:/Data', # Buffalo drive
              'O:/Data', # Blue drive
              '/localdisk/home/azvoleff/Data') # vertica1
prefix <- prefixes[match(TRUE, unlist(lapply(prefixes, function(x) file_test('-d', x))))]

data_folder <- file.path(prefix, "TEAM", "Tree_Growth", "Data")

library(runjags)
library(ggmcmc)
library(foreach)
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

###############################################################################
## Simple model (no random effects)
##

# Note the below is not parallel since running it in parallel breaks the Rcpp 
# function.
smp_mods <- foreach(temp_var=c('tmn', 'tmp', 'tmx'), .combine=rbind, .inorder=FALSE,
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

    mod <- destd(mod, mod$Parameter_Base == "B_g" & mod$param_ID == 1, dbh_sd)
    mod <- destd(mod, mod$Parameter_Base == "B_g" & mod$param_ID == 2, (dbh_sd/precip_sd) * 100) # Convert from mm to 10s of cm
    mod <- destd(mod, mod$Parameter_Base == "B_g" & mod$param_ID == 3, (dbh_sd/precip_sd) * 100) # Convert from mm to 10s of cm
    mod <- destd(mod, mod$Parameter_Base == "B_g" & mod$param_ID == 4, dbh_sd/temp_sd)
    mod <- destd(mod, mod$Parameter_Base == "B_g" & mod$param_ID == 5, dbh_sd/temp_sd)
    mod <- destd(mod, mod$Parameter_Base == "B_g" & mod$param_ID == 6, dbh_sd/dbh_sd)
    mod <- destd(mod, mod$Parameter_Base == "B_g" & mod$param_ID == 7, dbh_sd/dbh_sd)
    # TODO: How to de-standardize interactions?
    # mod <- destd(mod, mod$Parameter_Base == "B_g" & mod$param_ID == 8, dbh_sd/dbh_sd)
    # mod <- destd(mod, mod$Parameter_Base == "B_g" & mod$param_ID == 9, dbh_sd/dbh_sd)

    mod$Model <- temp_var 
    return(tbl_df(mod))
}
save(smp_mods, file="smp_mods.RData")

simple_B[simple_b$parameter == 'p'] <- simple_B[simple_B$Parameter == 'P'] * (dbh_sd/precip_sd) * 100 # Convert from mm to 10s of cm
simple_B[simple_B$Parameter == 'P^2'] <- simple_B[simple_B$Parameter == 'P^2'] * (dbh_sd/precip_sd) * 100 # Convert from mm to 10s of cm
simple_B[simple_B$Parameter == 'T'] <- simple_B[simple_B$Parameter == 'T'] * (dbh_sd/temp_sd)
simple_B[simple_B$Parameter == 'T^2'] <- simple_B[simple_B$Parameter == 'T^2'] * (dbh_sd/temp_sd)

###############################################################################
## Full model (interactions, uncorrelated random effects)

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

    mod <- destd(mod, mod$Parameter_Base == "B_g" & mod$param_ID == 1, dbh_sd)
    mod <- destd(mod, mod$Parameter_Base == "B_g" & mod$param_ID == 2, (dbh_sd/precip_sd) * 100) # Convert from mm to 10s of cm
    mod <- destd(mod, mod$Parameter_Base == "B_g" & mod$param_ID == 3, (dbh_sd/precip_sd) * 100) # Convert from mm to 10s of cm
    mod <- destd(mod, mod$Parameter_Base == "B_g" & mod$param_ID == 4, dbh_sd/temp_sd)
    mod <- destd(mod, mod$Parameter_Base == "B_g" & mod$param_ID == 5, dbh_sd/temp_sd)
    mod <- destd(mod, mod$Parameter_Base == "B_g" & mod$param_ID == 6, dbh_sd/dbh_sd)
    mod <- destd(mod, mod$Parameter_Base == "B_g" & mod$param_ID == 7, dbh_sd/dbh_sd)
    # TODO: How to de-standardize interactions?
    # mod <- destd(mod, mod$Parameter_Base == "B_g" & mod$param_ID == 8, dbh_sd/dbh_sd)
    # mod <- destd(mod, mod$Parameter_Base == "B_g" & mod$param_ID == 9, dbh_sd/dbh_sd)

    mod$Model <- temp_var 
    return(tbl_df(mod))
}
save(int_mods, file="int_mods.RData")

###############################################################################
## Full model (correlated random effects)
start_val <- 20000
thin_val <- 10

# Note the below is not parallel since running it in parallel breaks the Rcpp 
# function.
cor_mods <- foreach(temp_var=c('tmn', 'tmp', 'tmx'), .combine=rbind, .inorder=FALSE,
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

    mod <- destd(mod, mod$Parameter_Base == "B_g" & mod$param_ID == 1, dbh_sd)
    mod <- destd(mod, mod$Parameter_Base == "B_g" & mod$param_ID == 2, (dbh_sd/precip_sd) * 100) # Convert from mm to 10s of cm
    mod <- destd(mod, mod$Parameter_Base == "B_g" & mod$param_ID == 3, (dbh_sd/precip_sd) * 100) # Convert from mm to 10s of cm
    mod <- destd(mod, mod$Parameter_Base == "B_g" & mod$param_ID == 4, dbh_sd/temp_sd)
    mod <- destd(mod, mod$Parameter_Base == "B_g" & mod$param_ID == 5, dbh_sd/temp_sd)
    mod <- destd(mod, mod$Parameter_Base == "B_g" & mod$param_ID == 6, dbh_sd/dbh_sd)
    mod <- destd(mod, mod$Parameter_Base == "B_g" & mod$param_ID == 7, dbh_sd/dbh_sd)
    # TODO: How to de-standardize interactions?
    # mod <- destd(mod, mod$Parameter_Base == "B_g" & mod$param_ID == 8, dbh_sd/dbh_sd)
    # mod <- destd(mod, mod$Parameter_Base == "B_g" & mod$param_ID == 9, dbh_sd/dbh_sd)

    mod$Model <- temp_var 
    return(tbl_df(mod))
}
save(cor_mods, file="cor_mods.RData")
