library(dplyr)
library(reshape2)
library(lme4)

library(foreach)
library(doParallel)

source("0_settings.R")

cl <- makeCluster(8)
registerDoParallel(cl)

runmodels <- TRUE

data_folder <- file.path(prefix, "TEAM", "Tree_Growth", "Data")
init_folder <- file.path(prefix, "TEAM", "Tree_Growth", "Initialization")

model_type <- model_types[1]
precip_var <- precip_vars[1]
temp_var <- temp_vars[3]

#note <- ""
note <- "highelev"

ret <- foreach (model_type=model_types) %:%
    foreach (temp_var=temp_vars) %:%
        foreach (precip_var=precip_vars,
                 .packages=c("reshape2", "dplyr", "lme4"),
                 .inorder=FALSE) %dopar% {

    in_suffix <- paste0('_', model_type, '-', temp_var, '-', precip_var)
    if (note != "") out_suffix <- paste0(in_suffix, '_', note)

    load(file.path(init_folder, paste0("init_data", in_suffix, ".RData")))
    load(file.path(data_folder, paste0("model_data_wide", in_suffix, ".RData")))
    load(file.path(data_folder, paste0("model_data_long", in_suffix, ".RData")))

    precip_long <- melt(model_data$precip, varnames=c("tree_ID", "period_num"), 
                      value.name="precip")
    dbh_latent_long <- melt(init_data$dbh_latent, varnames=c("tree_ID", "period_num"), 
                            value.name="dbh_latent_end")
    temp_long <- melt(model_data$temp, varnames=c("tree_ID", "period_num"), 
                      value.name="temp")
    dbh_latent_long <- group_by(dbh_latent_long, tree_ID) %>%
        arrange(period_num) %>%
        mutate(dbh_latent_start=c(NA, dbh_latent_end[1:length(dbh_latent_end) - 1]))
    calib_data <- merge(dbh_latent_long, precip_long)
    calib_data <- merge(calib_data, temp_long)
    calib_data <- calib_data[complete.cases(calib_data), ]

    WD <- data.frame(tree_ID=seq(1, model_data$n_tree), WD=model_data$WD)
    calib_data <- merge(calib_data, WD)
    site_ID <- data.frame(tree_ID=seq(1, model_data$n_tree), site_ID=model_data$site_ID)
    calib_data <- merge(calib_data, site_ID)
    plot_ID <- data.frame(tree_ID=seq(1, model_data$n_tree), plot_ID=model_data$plot_ID)
    calib_data <- merge(calib_data, plot_ID)

    get_variance <- function(model, grp, var1, var2) {
        vc <- as.data.frame(VarCorr(model))
        if (missing(var1) & missing(var2)) {
            if (grp != "Residual") {
                stop('must specify a var1 if grp is not "Residual"')
            } else {
                return(vc[vc$grp == "Residual", ]$vcov)
            }
        } else if (missing(var1) & !missing(var2)) {
            stop('must specify a var1 if var2 is specified')
        } else {
            vc <- vc[vc$grp == grp & vc$var1== var1, ]
            if (missing(var2)) {
                return(vc[is.na(vc$var2), ]$vcov)
            } else {
                return(vc$vcov[match(var2, vc$var2)])
            }
        }
    }

    ###########################################################################
    # Inits with period random intercepts
    if (runmodels) {
        calib_model  <- lmer(dbh_latent_end ~ (1|site_ID) + (1|plot_ID) + 
                             (1|tree_ID) + (1|period_num), data=calib_data,
                             control=lmerControl(optCtrl=list(maxfun=10*35^2)))
        save(calib_model, file=file.path(init_folder, paste0("calib_model", out_suffix, ".RData")))
    } else {
        load(file.path(init_folder, paste0("calib_model", out_suffix, ".RData")))
    }

    init_data$int_ijk <- as.numeric(unlist(ranef(calib_model)$tree_ID))
    init_data$int_jk <- as.numeric(unlist(ranef(calib_model)$plot_ID))
    init_data$int_k <- as.numeric(unlist(ranef(calib_model)$site_ID))
    init_data$int_t <- as.numeric(unlist(ranef(calib_model)$period))
    init_data$sigma_int_ijk <- sqrt(get_variance(calib_model, "tree_ID", "(Intercept)"))
    init_data$sigma_int_jk <- sqrt(get_variance(calib_model, "plot_ID", "(Intercept)"))
    init_data$sigma_int_k <- sqrt(get_variance(calib_model, "site_ID", "(Intercept)"))
    init_data$sigma_int_t <- sqrt(get_variance(calib_model, "period_num", "(Intercept)"))
    save(init_data, file=file.path(init_folder, paste0("init_data_with_ranefs", out_suffix, ".RData")))

    return(TRUE)
}

stopCluster(cl)
