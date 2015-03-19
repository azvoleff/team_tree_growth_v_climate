library(dplyr)
library(reshape2)
library(zoo) # for na.approx
library(Rcpp)
library(inline)
library(foreach)
library(doParallel)

source("0_settings.R")

cl <- makeCluster(3)
registerDoParallel(cl)

data_folder <- file.path(prefix, "TEAM", "Tree_Growth", "Data")
init_folder <- file.path(prefix, "TEAM", "Tree_Growth", "Initialization")

# temp_var <- 'tmn_meanannual'
# precip_var <- 'mcwd_run12'
# model_type <- 'full'

foreach (model_type=model_types) %:%
    foreach (temp_var=temp_vars) %:%
        foreach (precip_var=precip_vars,
                 .packages=c("dplyr", "Rcpp", "inline", "zoo", "reshape2"),
                 .inorder=FALSE) %dopar% {
    sourceCpp('calc_missings.cpp')
    sourceCpp('left_align.cpp')

    load("growth_ctfsflagged_merged_detrended.RData")

    suffix <- paste0('_', model_type, '-', temp_var, '-', precip_var)

    growth <- tbl_df(growth)

    # table(growth$ctfs_accept)
    growth <- filter(growth, ctfs_accept)

    # Exclude CSN, YAN, and NAK as genera are poorly identified
    growth <- filter(growth, !(sitecode %in% c('NAK', 'CSN', 'YAN')))

    # Remove very rare genera (those with fewer than 20 individuals) to avoid
    # convergence problems.
    #
    # n_per_genus <- group_by(growth, Genus) %>%
    #     summarize(n=length(unique(SamplingUnitName)))
    # sum(n_per_genus$n > 0)
    # sum(n_per_genus$n > 5)
    # sum(n_per_genus$n > 10)
    # sum(n_per_genus$n > 20)
    growth <- group_by(growth, Genus) %>%
        mutate(n_indiv_per_genus=length(unique(SamplingUnitName))) %>%
        filter(n_indiv_per_genus >= 20)

    setup_factors <- function(growth) {
        growth$tree_ID_char <- factor(growth$SamplingUnitName)
        growth$plot_ID_char <- factor(growth$plot_ID)
        growth$site_ID_char <- factor(growth$sitecode)
        growth$period_ID_char <- factor(growth$SamplingPeriodEnd)
        growth$genus_ID_char <- factor(growth$Genus)
        growth$tree_ID <- factor(growth$SamplingUnitName)
        growth$plot_ID <- factor(growth$plot_ID)
        growth$site_ID <- factor(growth$sitecode)
        growth$genus_ID <- factor(growth$Genus)
        growth$period_ID <- as.integer(ordered(growth$SamplingPeriodEnd))
        return(growth)
    }
    growth <- setup_factors(growth)

    # Subsample dataset for running testing models
    if (model_type == "testing") {
        set.seed(1)
        # Generate a subsample including only 5% of the trees included in the 
        # initial sampling period at each site
        initial_trees <- group_by(growth, site_ID) %>%
            filter(period_ID == min(period_ID)) %>%
            group_by(site_ID) %>%
            sample_frac(.1)
        growth <- filter(growth, tree_ID %in% initial_trees$tree_ID)
        growth <- group_by(growth, Genus) %>%
            mutate(n_indiv_per_genus=length(unique(SamplingUnitName))) %>%
            filter(n_indiv_per_genus >= 20)
        growth <- setup_factors(growth)
    }

    # Rename chosen precip var to be named "precip" and chosen temp var to be 
    # "temp" for easier handling in the later code.
    stopifnot(sum(names(growth) == temp_var) == 1)
    stopifnot(sum(names(growth) == precip_var) == 1)
    stopifnot(sum(names(growth) == 'temp') == 0)
    stopifnot(sum(names(growth) == 'precip') == 0)
    names(growth)[names(growth) ==  precip_var] <- "precip"
    names(growth)[names(growth) ==  temp_var] <- "temp"

    # Setup timeseries formatted dataframe (with NAs for covariates at time 1 since 
    # first growth measurement can only be calculated from time 1 to time 2)
    dbh_ts <- arrange(growth, tree_ID, period_ID) %>%
        select(tree_ID, plot_ID, site_ID, period_ID, genus_ID, WD, 
               dbh=diameter_end, diameter_start, precip, temp)

    # Setup an array of initial diameters for all the trees, by taking the 
    # dbh_start
    # value from the first period with an available observation
    dbh_time_0 <- arrange(dbh_ts, tree_ID, period_ID) %>%
        group_by(tree_ID) %>%
        filter(period_ID==min(period_ID)) %>%
        select(tree_ID, plot_ID, site_ID, period_ID, genus_ID, WD, diameter_start)
    dbh_time_0$dbh <- dbh_time_0$diameter_start
    dbh_time_0 <- select(dbh_time_0, -diameter_start)
    dbh_ts <- select(dbh_ts, -diameter_start)

    dbh_time_0$precip <- NA
    dbh_time_0$temp <- NA
    dbh_time_0$period_ID <- dbh_time_0$period_ID - 1

    dbh_ts <- merge(dbh_ts, dbh_time_0, all=TRUE)
    dbh_ts <- arrange(dbh_ts, tree_ID, period_ID)

    #dbh_ts$precip <- log(1 + dbh_ts$precip)

    # Standardize outcome and predictors.
    dbh_mean <- mean(dbh_ts$dbh)
    dbh_min <- min(dbh_ts$dbh)
    dbh_max <- max(dbh_ts$dbh)
    dbh_sd <- sd(dbh_ts$dbh)
    dbh_ts$dbh <- (dbh_ts$dbh - dbh_mean) / dbh_sd
    precip_mean <- mean(dbh_ts$precip, na.rm=TRUE)
    precip_min <- min(dbh_ts$precip, na.rm=TRUE)
    precip_max <- max(dbh_ts$precip, na.rm=TRUE)
    precip_sd <- sd(dbh_ts$precip, na.rm=TRUE)
    dbh_ts$precip <- (dbh_ts$precip - precip_mean) / precip_sd
    temp_mean <- mean(dbh_ts$temp, na.rm=TRUE)
    temp_min <- min(dbh_ts$temp, na.rm=TRUE)
    temp_max <- max(dbh_ts$temp, na.rm=TRUE)
    temp_sd <- sd(dbh_ts$temp, na.rm=TRUE)
    dbh_ts$temp <- (dbh_ts$temp - temp_mean) / temp_sd
    WD_mean <- mean(dbh_ts$WD)
    WD_min <- min(dbh_ts$WD)
    WD_max <- max(dbh_ts$WD)
    WD_sd <- sd(dbh_ts$WD)
    WD <- (dbh_time_0$WD - WD_mean) / WD_sd

    genus_ID <- dbh_time_0$genus_ID
    sum(genus_ID == "Unknown") / length(genus_ID)

    # Setup wide format dbh, precip, and temp dataframes
    dbh <- dcast(dbh_ts, tree_ID + plot_ID + site_ID ~ period_ID, value.var="dbh")
    precip <- dcast(dbh_ts, tree_ID + plot_ID + site_ID ~ period_ID, value.var="precip")
    temp <- dcast(dbh_ts, tree_ID + plot_ID + site_ID ~ period_ID, value.var="temp")
    tree_ID <- dbh$tree_ID
    plot_ID <- dbh$plot_ID
    site_ID <- dbh$site_ID
    # Eliminate the ID columns
    dbh <- dbh[!grepl('_ID$', names(dbh))]
    precip <- precip[!grepl('_ID$', names(precip))]
    temp <- temp[!grepl('_ID$', names(temp))]

    # Calculate per-plot elevation (not time-varying)
    elev <- group_by(growth, plot_ID) %>%
        summarize(elev=elev[1])
    stopifnot(levels(plot_ID) == elev$plot_ID)
    elev <- elev$elev
    elev_mean <- mean(elev)
    elev_min <- min(elev)
    elev_max <- max(elev)
    elev_sd <- sd(elev)
    elev <- (elev - elev_mean) / elev_sd

    # Save sd and means so the variables can be unstandardized later
    save(dbh_mean, dbh_sd, dbh_min, dbh_max,
         WD_mean, WD_sd, WD_min, WD_max,
         elev_mean, elev_sd, elev_min, elev_max,
         temp_mean, temp_sd, temp_min, temp_max,
         precip_mean, precip_min, precip_max, precip_sd, 
         file=file.path(data_folder, paste0("model_data_standardizing", suffix, 
                                     ".RData")))

    # Setup data
    n_tree <- length(unique(dbh_ts$tree_ID))
    n_site <- length(unique(dbh_ts$site_ID))
    n_plot <- length(unique(dbh_ts$plot_ID))
    # One less period than there are max number of observations (since n_period
    # relates to number of growth measurements).
    n_period <- length(unique(dbh_ts$period_ID)) - 1
    n_genus <- length(unique(dbh_ts$genus_ID))

    # Calculate the first and last observation for each tree. The +1's below are 
    # because the period_ID variable starts at zero.
    obs_per_tree <- group_by(dbh_ts, tree_ID) %>%
        summarize(first_obs_period=(min(period_ID) + 1),
                  last_obs_period=(max(period_ID) + 1))

    # precip observations are missing for periods when dbh observations are 
    # missing.  Fill these observations using the mean precip for the appropriate 
    # period.
    #
    # Note that the same centering and scaling needs to be done here (and log) that 
    # was done for the other precip data.
    precip_means <- group_by(growth, plot_ID, period_ID) %>%
         summarize(precip_means=mean((precip - precip_mean) / precip_sd, na.rm=TRUE))
    precip_means <- data.frame(precip_means) # Fix for indexing bug in dplyr 0.3.2
    precip <- as.matrix(precip)
    precip_missings <- calc_missings(precip)$miss
    # Use linear indexing to replace NAs in precips
    precip_miss_linear_ind <- (precip_missings[, 2] - 1) * nrow(precip) + precip_missings[, 1] # From http://bit.ly/1rnKrC3
    precip[precip_miss_linear_ind] <- precip_means[match(paste(plot_ID[precip_missings[, 1]], precip_missings[, 2]),
                                                paste(precip_means$plot_ID, precip_means$period_ID)), 3]
    stopifnot(is.null(calc_missings(precip)$miss))

    # Fill in temp (same issue as above)
    temp_means <- group_by(growth, plot_ID, period_ID) %>%
         summarize(temp_means=mean((temp - temp_mean) / temp_sd, na.rm=TRUE))
    temp_means <- data.frame(temp_means) # Fix for indexing bug in dplyr 0.3.2
    temp <- as.matrix(temp)
    temp_missings <- calc_missings(temp)$miss
    # Use linear indexing to replace NAs in temps
    temp_miss_linear_ind <- (temp_missings[, 2] - 1) * nrow(temp) + temp_missings[, 1] # From http://bit.ly/1rnKrC3
    temp[temp_miss_linear_ind] <- temp_means[match(paste(plot_ID[temp_missings[, 1]], temp_missings[, 2]),
                                                paste(temp_means$plot_ID, temp_means$period_ID)), 3]
    stopifnot(is.null(calc_missings(temp)$miss))

    # Calculate indices of missing and observed data (using the function 
    # defined in calc_missings.cpp), so that observed and missing data can be 
    # modeled separately in Stan.
    missings <- calc_missings(as.matrix(dbh))

    stopifnot(length(elev) == n_plot)
    stopifnot(length(WD) == n_tree)
    stopifnot(length(genus_ID) == n_tree)
    stopifnot(length(plot_ID) == n_tree)
    stopifnot(length(site_ID) == n_tree)

    # Calculate precision of diameter tape (1 mm) in standardized units:
    sigma_obs_lower <- .1 / dbh_sd / sqrt(12)

    ###############################################################################
    # Output wide format data
    model_data <- list(n_tree=n_tree,
                       n_plot=n_plot,
                       n_site=n_site,
                       n_period=n_period,
                       n_genus=n_genus,
                       first_obs_period=obs_per_tree$first_obs_period,
                       last_obs_period=obs_per_tree$last_obs_period,
                       plot_ID=as.integer(factor(plot_ID)),
                       site_ID=as.integer(factor(site_ID)),
                       genus_ID=as.integer(factor(genus_ID)),
                       sigma_obs_lower=sigma_obs_lower,
                       dbh=as.matrix(dbh),
                       WD=WD,
                       elev=elev,
                       temp=temp,
                       precip=precip,
                       obs_indices=missings$obs,
                       miss_indices=missings$miss)
    save(model_data, file=file.path(data_folder, paste0("model_data_wide", 
                                                       suffix, ".RData")))

    ###############################################################################
    # Output long format data
    model_data_long <- data.frame(tree_ID=factor(tree_ID),
                                  plot_ID=factor(plot_ID),
                                  site_ID=factor(site_ID),
                                  genus_ID=factor(genus_ID),
                                  sigma_obs_lower=sigma_obs_lower,
                                  WD=WD)
    # Standardize the diameter variables and precip variable the same way they are 
    # standardized in the wide dataset
    merge_data <- data.frame(tree_ID=growth$tree_ID,
                             growth_rgr=growth$growth_rgr,
                             diameter_start=(growth$diameter_start - dbh_mean) / dbh_sd,
                             diameter_end=(growth$diameter_end - dbh_mean) / dbh_sd,
                             n_days=growth$n_days,
                             elev=(growth$elev - elev_mean) / elev_sd,
                             temp=(growth$temp - temp_mean) / temp_sd,
                             precip=(growth$precip - precip_mean) / precip_sd)
    model_data_long <- merge(model_data_long, merge_data, by="tree_ID", all=TRUE)
    save(model_data_long, file=file.path(data_folder, paste0("model_data_long", 
                                                            suffix, ".RData")))

    genus_ID_factor_levels <- levels(factor(genus_ID))
    genus_ID_factor_key <- cbind(genus_ID_char=as.character(genus_ID_factor_levels), 
                                 genus_ID_numeric=seq(1:length(genus_ID_factor_levels)))
    write.csv(genus_ID_factor_key, file=file.path(data_folder, paste0("genus_ID_factor_key", 
                                                         suffix, ".csv")), 
                                                  row.names=FALSE)

    site_ID_factor_levels <- levels(factor(site_ID))
    site_ID_factor_key <- cbind(site_ID_char=as.character(site_ID_factor_levels), 
                                site_ID_numeric=seq(1:length(site_ID_factor_levels)))
    write.csv(site_ID_factor_key, file=file.path(data_folder, 
                                                 paste0("site_ID_factor_key", 
                                                        suffix, ".csv")), 
                                                 row.names=FALSE)

    period_ID_factor_levels <- levels(ordered(growth$SamplingPeriodEnd))
    period_ID_factor_key <- data.frame(period_ID_char=as.character(period_ID_factor_levels), 
                                       period_ID_numeric=seq(1:length(period_ID_factor_levels)),
                                       stringsAsFactors=FALSE)
    # Add in initial period (coded as zero, but not in the SamplingPeriodEnd 
    # vector)
    initial_period <- paste0(as.numeric(substr(period_ID_factor_key$period_ID_char[1], 1, 4)) - 1, ".01")
    period_ID_factor_key <- rbind(c(initial_period, 0), period_ID_factor_key)
    write.csv(period_ID_factor_key,
              file=file.path(data_folder, paste0("period_ID_factor_key", suffix, 
                                                ".csv")), row.names=FALSE)

    ###############################################################################
    # Setup inits

    # Setup wide format dbh_latent dataframes
    dbh_latent <- dcast(dbh_ts, tree_ID + plot_ID + site_ID ~ period_ID, value.var="dbh")
    dbh_latent <- dbh_latent[!grepl('_ID$', names(dbh_latent))]
    # Smooth dbh observations with spline to calculate latent dbh
    smooth_dbh_obs <- function(x) {
        first_obs <- min(which(!is.na(x)))
        last_obs <- max(which(!is.na(x)))
        if (sum(!is.na(x)) >= 4) {
            spline_smoother <- smooth.spline(x=c(1:length(x))[!is.na(x)], y=x[!is.na(x)])
            preds <- predict(spline_smoother, c(1:length(x)))$y
            x[first_obs:last_obs] <- preds[first_obs:last_obs]
        } else {
            x[first_obs:last_obs] <- na.approx(as.numeric(x[first_obs:last_obs]))
        }
        return(x)
    }
    dbh_latent <- t(apply(dbh_latent, 1, smooth_dbh_obs))

    init_data <- list(dbh_latent=as.matrix(dbh_latent))
    save(init_data, file=file.path(init_folder, paste0("init_data", suffix, 
                                                      ".RData")))

    ###############################################################################
    # Output data in a format that is clustered by the number of periods of 
    # observations there are per individual - this format allows Stan to run 
    # faster.
    head(model_data$temp)

    dbh_la <- left_align(model_data$dbh)
    dbh_latent_la <- left_align(dbh_latent)
    temp_la <- left_align(model_data$temp, indent=1)
    precip_la <- left_align(model_data$precip, indent=1)

    # Reorder covariates and independent variables by number of periods
    last_obs_dbh_latent <- data.frame(tree_ID=as.integer(factor(tree_ID)),
                                      last_obs=apply(dbh_latent_la, 1, function(x) match(NA, x) - 1))
    last_obs_dbh <- data.frame(tree_ID=as.integer(factor(tree_ID)),
                               last_obs=apply(dbh_latent_la, 1, function(x) match(NA, x) - 1))
    stopifnot(identical(last_obs_dbh, last_obs_dbh_latent))
    new_order <- order(last_obs_dbh$last_obs, last_obs_dbh$tree_ID)
    last_obs_dbh <- last_obs_dbh[new_order, ]

    dbh_la <- dbh_la[new_order, ]
    dbh_latent_la <- dbh_latent_la[new_order, ]
    temp_la <- temp_la[new_order, ]
    precip_la <- precip_la[new_order, ]
    tree_ID_la <- as.integer(factor(tree_ID))[new_order]

    bl_size <- unique(last_obs_dbh$last_obs)
    # Drop the NA (which is for rows without ANY missing observations - meaning 
    # they have no first NA)
    bl_size <- bl_size[!is.na(bl_size)]
    bl_st <- match(bl_size, last_obs_dbh$last_obs)
    bl_end <- c(bl_st[2:length(bl_st)] - 1, nrow(last_obs_dbh))

    missings_wide <- calc_missings(as.matrix(dbh_la))

    model_data_blocked <- list(n_tree=n_tree,
                               n_plot=n_plot,
                               n_site=n_site,
                               n_period=n_period,
                               n_genus=n_genus,
                               n_blocks=length(bl_size),
                               bl_size=bl_size,
                               bl_st=bl_st,
                               bl_end=bl_end,
                               tree_ID=tree_ID_la,
                               first_obs_period=obs_per_tree$first_obs_period, # doesn't need to be left aligned since it is indexed by tree_ID
                               last_obs_period=obs_per_tree$last_obs_period, # doesn't need to be left aligned since it is indexed by tree_ID
                               plot_ID=as.integer(as.factor(plot_ID)), # doesn't need to be left aligned since it is indexed by tree_ID
                               site_ID=as.integer(as.factor(site_ID)), # doesn't need to be left aligned since it is indexed by tree_ID
                               genus_ID=as.integer(as.factor(genus_ID)), # doesn't need to be left aligned since it is indexed by tree_ID
                               sigma_obs_lower=sigma_obs_lower,
                               dbh=as.matrix(dbh_la),
                               WD=WD, # doesn't need to be left aligned since it is indexed by tree_ID
                               elev=elev[plot_ID], # doesn't need to be left aligned since it is indexed by tree_ID
                               temp=temp_la,
                               precip=precip_la,
                               obs_indices=missings_wide$obs,
                               miss_indices=missings_wide$miss)
    save(model_data_blocked,
         file=file.path(data_folder, paste0("model_data_wide_blocked", suffix, ".RData")))

    init_data_blocked <- list(dbh_latent_la=as.matrix(dbh_latent_la))
    save(init_data_blocked,
         file=file.path(init_folder, paste0("init_data_blocked", suffix, ".RData")))

}

stopCluster(cl)
