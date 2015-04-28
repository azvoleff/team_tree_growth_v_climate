# This script needs to be run a machine with the C++ Boost libraries installed

prefixes <- c('D:/azvoleff/Data', # CI-TEAM
              'H:/Data', # Buffalo drive
              'O:/Data', # Blue drive
              '/localdisk/home/azvoleff/Data') # vertica1
prefix <- prefixes[match(TRUE, unlist(lapply(prefixes, function(x) file_test('-d', x))))]

base_folder <- file.path(prefix, "TEAM", "Tree_Growth")

library(runjags)
library(ggmcmc)
library(foreach)
library(doParallel)
library(coda)
library(Rcpp)
library(dplyr)

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

# Functions to "destandardize" standardized coefficients in an array of MCMC 
# results
destd <- function(d, rows, multiplier) {
    stopifnot(sum(rows) > 0)
    d[rows, ]$value <- d[rows, ]$value * multiplier
    return(d)
}

destandardize <- function(d, model_type, temp_var) {
    load(file.path(base_folder, 'Data',
                   paste0("model_data_standardizing_full-", temp_var, 
                          "_meanannual-mcwd_run12.RData")))
    if (model_type == "simple") {
        d <- destd(d, d$Parameter_Base == "B" & d$row_ID == 1, dbh_sd)
        d <- destd(d, d$Parameter_Base == "B" & d$row_ID == 2, dbh_sd/precip_sd)
        d <- destd(d, d$Parameter_Base == "B" & d$row_ID == 3, dbh_sd/(precip_sd^2))
        d <- destd(d, d$Parameter_Base == "B" & d$row_ID == 4, dbh_sd/temp_sd)
        d <- destd(d, d$Parameter_Base == "B" & d$row_ID == 5, dbh_sd/(temp_sd^2))
        d <- destd(d, d$Parameter_Base == "B" & d$row_ID == 6, 1)
        d <- destd(d, d$Parameter_Base == "B" & d$row_ID == 7, 1/dbh_sd)
        d <- destd(d, d$Parameter_Base == "int_jk", dbh_sd)
        d <- destd(d, d$Parameter_Base == "int_k", dbh_sd)
        d <- destd(d, d$Parameter_Base == "int_t", dbh_sd)
        d <- destd(d, d$Parameter_Base == "sigma_obs", dbh_sd)
        d <- destd(d, d$Parameter_Base == "sigma_proc", dbh_sd)
        d <- destd(d, d$Parameter_Base == "sigma_B_k", dbh_sd)
        d <- destd(d, d$Parameter_Base == "sigma_int_jk", dbh_sd)
        d <- destd(d, d$Parameter_Base == "sigma_int_k", dbh_sd)
        d <- destd(d, d$Parameter_Base == "sigma_int_t", dbh_sd)
        d <- destd(d, d$Parameter_Base == "B_k", dbh_sd)
    } else if (model_type == "interact") {
        B_g_rows <- d$Parameter_Base == "B_g"
        B_g_sigma_mu_rows <- d$Parameter_Base %in% c("mu_B_g", "sigma_B_g")
        d <- destd(d, (B_g_sigma_mu_rows & d$row_ID == 1) | (B_g_rows & d$col_ID == 1), dbh_sd)
        d <- destd(d, (B_g_sigma_mu_rows & d$row_ID == 2) | (B_g_rows & d$col_ID == 2), dbh_sd/precip_sd)
        d <- destd(d, (B_g_sigma_mu_rows & d$row_ID == 3) | (B_g_rows & d$col_ID == 3), dbh_sd/(precip_sd^2))
        d <- destd(d, (B_g_sigma_mu_rows & d$row_ID == 4) | (B_g_rows & d$col_ID == 4), dbh_sd/temp_sd)
        d <- destd(d, (B_g_sigma_mu_rows & d$row_ID == 5) | (B_g_rows & d$col_ID == 5), dbh_sd/(temp_sd^2))
        d <- destd(d, (B_g_sigma_mu_rows & d$row_ID == 6) | (B_g_rows & d$col_ID == 6), 1)
        d <- destd(d, (B_g_sigma_mu_rows & d$row_ID == 7) | (B_g_rows & d$col_ID == 7), 1/dbh_sd)
        d <- destd(d, (B_g_sigma_mu_rows & d$row_ID == 8) | (B_g_rows & d$col_ID == 8), 1/precip_sd)
        d <- destd(d, (B_g_sigma_mu_rows & d$row_ID == 9) | (B_g_rows & d$col_ID == 9), 1/temp_sd)
        d <- destd(d, d$Parameter_Base == "int_jk", dbh_sd)
        d <- destd(d, d$Parameter_Base == "int_k", dbh_sd)
        d <- destd(d, d$Parameter_Base == "int_t", dbh_sd)
        d <- destd(d, d$Parameter_Base == "sigma_obs", dbh_sd)
        d <- destd(d, d$Parameter_Base == "sigma_proc", dbh_sd)
        d <- destd(d, d$Parameter_Base == "sigma_B_k", dbh_sd)
        d <- destd(d, d$Parameter_Base == "sigma_int_jk", dbh_sd)
        d <- destd(d, d$Parameter_Base == "sigma_int_k", dbh_sd)
        d <- destd(d, d$Parameter_Base == "sigma_int_t", dbh_sd)
        d <- destd(d, d$Parameter_Base == "B_k", dbh_sd)
    } else if (model_type == "correlated") {
        B_g_rows <- d$Parameter_Base == "B_g"
        B_g_sigma_mu_rows <- d$Parameter_Base %in% c("mu_B_g", "sigma_B_g")
        d <- destd(d, (B_g_sigma_mu_rows & d$row_ID == 1) | (B_g_rows & d$col_ID == 1), dbh_sd)
        d <- destd(d, (B_g_sigma_mu_rows & d$row_ID == 2) | (B_g_rows & d$col_ID == 2), dbh_sd/precip_sd)
        d <- destd(d, (B_g_sigma_mu_rows & d$row_ID == 3) | (B_g_rows & d$col_ID == 3), dbh_sd/(precip_sd^2))
        d <- destd(d, (B_g_sigma_mu_rows & d$row_ID == 4) | (B_g_rows & d$col_ID == 4), dbh_sd/temp_sd)
        d <- destd(d, (B_g_sigma_mu_rows & d$row_ID == 5) | (B_g_rows & d$col_ID == 5), dbh_sd/(temp_sd^2))
        d <- destd(d, (B_g_sigma_mu_rows & d$row_ID == 6) | (B_g_rows & d$col_ID == 6), 1)
        d <- destd(d, (B_g_sigma_mu_rows & d$row_ID == 7) | (B_g_rows & d$col_ID == 7), 1/dbh_sd)
        d <- destd(d, d$Parameter_Base == "int_jk", dbh_sd)
        d <- destd(d, d$Parameter_Base == "int_k", dbh_sd)
        d <- destd(d, d$Parameter_Base == "int_t", dbh_sd)
        d <- destd(d, d$Parameter_Base == "sigma_obs", dbh_sd)
        d <- destd(d, d$Parameter_Base == "sigma_proc", dbh_sd)
        d <- destd(d, d$Parameter_Base == "sigma_B_k", dbh_sd)
        d <- destd(d, d$Parameter_Base == "sigma_int_jk", dbh_sd)
        d <- destd(d, d$Parameter_Base == "sigma_int_k", dbh_sd)
        d <- destd(d, d$Parameter_Base == "sigma_int_t", dbh_sd)
        d <- destd(d, d$Parameter_Base == "B_k", dbh_sd)
    } else {
        stop('unrecognized model type')
    }
    return(d)
}

# model_type <- "interact"
# temp_var <- "tmn"

start_val <- 100000
thin_val <- 100
foreach(model_type=c('simple', 'interact', 'correlated')) %do% {
    # Note the below is not parallel since running it in parallel breaks the Rcpp 
    # function.
    params <- foreach(temp_var=c('tmn', 'tmp', 'tmx'), .combine=rbind, 
                      .inorder=FALSE,
                      .packages=c('runjags', 'ggmcmc', 'mcmcplots', 'Rcpp',
                                  'dplyr')) %do% {
        message(paste0("Processing ", model_type, "/", temp_var))
        load(file.path(base_folder, "MCMC_Chains", model_type,
                       paste0(model_type, '_', temp_var, '.RData')))
        jags_fit <- window(as.mcmc.list(jags_fit), start=start_val, thin=thin_val)
        these_params <- ggs(jags_fit)

        ids <- jags_param_ids(these_params$Parameter)
        these_params <- cbind(these_params, ids)

        these_params <- destandardize(these_params, model_type, temp_var)

        these_params$Model <- temp_var 
        return(tbl_df(these_params))
    }
    save(params, file=file.path(base_folder, "Extracted_Parameters", 
                                paste0("parameter_estimates_", model_type, 
                                       ".RData")))
}
