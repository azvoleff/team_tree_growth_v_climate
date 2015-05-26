library(dplyr)
library(reshape2)
library(zoo) # for na.approx
library(Rcpp)
library(inline)
library(foreach)
library(doParallel)

source("0_settings.R")

cl <- makeCluster(6)
registerDoParallel(cl)

data_folder <- file.path(prefix, "TEAM", "Tree_Growth", "Data")
init_folder <- file.path(prefix, "TEAM", "Tree_Growth", "Initialization")

foreach (model_type=model_types) %:%
    foreach (temp_var=temp_vars) %:%
        foreach (precip_var=precip_vars,
                 .packages=c("dplyr", "Rcpp", "inline", "zoo", "reshape2"),
                 .inorder=FALSE) %dopar% {
    sourceCpp('calc_missings.cpp')

    load("growth_ctfsflagged_merged_detrended.RData")

    suffix <- paste0('_', model_type, '-', temp_var, '-', precip_var)
    if (note != "") suffix <- paste0(suffix, '_', note)

    growth <- tbl_df(growth)

    # table(growth$ctfs_accept)
    growth <- filter(growth, ctfs_accept)

    # Exclude CSN, YAN, and NAK as genera are poorly identified
    growth <- filter(growth, !(sitecode %in% c('NAK', 'CSN', 'YAN')))

    if (note == "highelev") growth <- filter(growth, sitecode %in% c('VB', 'UDZ', 'BIF', 'YAN', 'NAK'))

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
        growth$tree_ID <- factor(growth$SamplingUnitName)
        growth$plot_ID <- factor(growth$plot_ID)
        growth$site_ID <- factor(growth$sitecode)
        growth$genus_ID <- factor(growth$Genus)
        growth$period_num <- as.integer(ordered(growth$SamplingPeriodEnd))
        return(growth)
    }
    growth <- setup_factors(growth)

    # Subsample dataset for running testing models
    if (model_type == "testing") {
        set.seed(1)
        # Generate a subsample including only 5% of the trees included in the 
        # initial sampling period at each site
        initial_trees <- group_by(growth, site_ID) %>%
            filter(period_num == min(period_num)) %>%
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
    dbh_ts <- arrange(growth, tree_ID, period_num) %>%
        select(tree_ID, plot_ID, site_ID, period_num, genus_ID, WD, 
               dbh=diameter_end, diameter_start, precip, temp)

    # Setup an array of initial diameters for all the trees, by taking the 
    # dbh_start
    # value from the first period with an available observation
    dbh_time_0 <- arrange(dbh_ts, tree_ID, period_num) %>%
        group_by(tree_ID) %>%
        filter(period_num==min(period_num)) %>%
        select(tree_ID, plot_ID, site_ID, period_num, genus_ID, WD, diameter_start)
    dbh_time_0$dbh <- dbh_time_0$diameter_start
    dbh_time_0 <- select(dbh_time_0, -diameter_start)
    dbh_ts <- select(dbh_ts, -diameter_start)

    dbh_time_0$precip <- NA
    dbh_time_0$temp <- NA
    dbh_time_0$period_num <- dbh_time_0$period_num - 1

    dbh_ts <- merge(dbh_ts, dbh_time_0, all=TRUE)
    dbh_ts <- arrange(dbh_ts, tree_ID, period_num)

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

    # Setup wide format dbh, precip, and temp dataframes
    dbh <- dcast(dbh_ts, tree_ID + plot_ID + site_ID + genus_ID ~ period_num, value.var="dbh")
    precip <- dcast(dbh_ts, tree_ID + plot_ID + site_ID + genus_ID ~ period_num, value.var="precip")
    temp <- dcast(dbh_ts, tree_ID + plot_ID + site_ID + genus_ID ~ period_num, value.var="temp")
    tree_ID <- dbh$tree_ID
    plot_ID <- dbh$plot_ID
    site_ID <- dbh$site_ID
    genus_ID <- dbh$genus_ID
    # Eliminate the ID columns
    dbh <- dbh[!grepl('_ID$', names(dbh))]
    precip <- precip[!grepl('_ID$', names(precip))]
    temp <- temp[!grepl('_ID$', names(temp))]

    # Calculate per-plot elevation as difference mean elevation of CRU grid 
    # cell (not time-varying)
    growth$elev_diff <- growth$elev - growth$elev_cru
    elev_diff <- group_by(growth, plot_ID) %>%
        summarize(elev_diff=elev_diff[1])
    stopifnot(levels(plot_ID) == elev_diff$plot_ID)
    elev_diff <- elev_diff$elev_diff
    elev_diff_mean <- mean(elev_diff)
    elev_diff_min <- min(elev_diff)
    elev_diff_max <- max(elev_diff)
    elev_diff_sd <- sd(elev_diff)
    elev_diff <- (elev_diff - elev_diff_mean) / elev_diff_sd

    # Save sd and means so the variables can be unstandardized later
    save(dbh_mean, dbh_sd, dbh_min, dbh_max,
         WD_mean, WD_sd, WD_min, WD_max,
         elev_diff_mean, elev_diff_sd, elev_diff_min, elev_diff_max,
         temp_mean, temp_sd, temp_min, temp_max,
         precip_mean, precip_min, precip_max, precip_sd, 
         file=file.path(data_folder, paste0("model_data_standardizing", suffix, 
                                     ".RData")))

    # Save a key that gives elevation difference by plot ID and site - will be 
    # used later when making predictions.
    elev_key <- data.frame(plot_ID=unique(growth$plot_ID))
    elev_key$elev_diff <- elev_diff # Note elev_diff is sorted by plot_ID
    save(elev_key, file=file.path(data_folder, paste0("model_data_elev_key", 
                                                      suffix, ".RData")))

    # Setup data
    n_tree <- length(unique(dbh_ts$tree_ID))
    n_site <- length(unique(dbh_ts$site_ID))
    n_plot <- length(unique(dbh_ts$plot_ID))
    # One less period than there are max number of observations (since n_period
    # relates to number of growth measurements).
    n_period <- length(unique(dbh_ts$period_num)) - 1
    n_genus <- length(unique(dbh_ts$genus_ID))

    # Calculate the first and last observation for each tree. The +1's below are 
    # because the period_num variable starts at zero.
    obs_per_tree <- group_by(dbh_ts, tree_ID) %>%
        summarize(first_obs_period=(min(period_num) + 1),
                  last_obs_period=(max(period_num) + 1))

    # precip observations are missing for periods when dbh observations are 
    # missing.  Fill these observations using the mean precip for the appropriate 
    # period.
    #
    # Note that the same centering and scaling needs to be done here (and log) that 
    # was done for the other precip data.
    precip_means <- group_by(growth, plot_ID, period_num) %>%
         summarize(precip_means=mean((precip - precip_mean) / precip_sd, na.rm=TRUE))
    precip_means <- data.frame(precip_means) # Fix for indexing bug in dplyr 0.3.2
    precip <- as.matrix(precip)
    precip_missings <- calc_missings(precip)$miss
    # Use linear indexing to replace NAs in precips
    precip_miss_linear_ind <- (precip_missings[, 2] - 1) * nrow(precip) + precip_missings[, 1] # From http://bit.ly/1rnKrC3
    precip[precip_miss_linear_ind] <- precip_means[match(paste(plot_ID[precip_missings[, 1]], precip_missings[, 2]),
                                                paste(precip_means$plot_ID, precip_means$period_num)), 3]
    stopifnot(is.null(calc_missings(precip)$miss))

    # Fill in temp (same issue as above)
    temp_means <- group_by(growth, plot_ID, period_num) %>%
         summarize(temp_means=mean((temp - temp_mean) / temp_sd, na.rm=TRUE))
    temp_means <- data.frame(temp_means) # Fix for indexing bug in dplyr 0.3.2
    temp <- as.matrix(temp)
    temp_missings <- calc_missings(temp)$miss
    # Use linear indexing to replace NAs in temps
    temp_miss_linear_ind <- (temp_missings[, 2] - 1) * nrow(temp) + temp_missings[, 1] # From http://bit.ly/1rnKrC3
    temp[temp_miss_linear_ind] <- temp_means[match(paste(plot_ID[temp_missings[, 1]], temp_missings[, 2]),
                                                paste(temp_means$plot_ID, temp_means$period_num)), 3]
    stopifnot(is.null(calc_missings(temp)$miss))

    # Calculate indices of missing and observed data (using the function 
    # defined in calc_missings.cpp), so that observed and missing data can be 
    # modeled separately in Stan.
    missings <- calc_missings(as.matrix(dbh))

    stopifnot(length(elev_diff) == n_plot)
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
                       elev_diff=elev_diff,
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
                             period_num=growth$period_num,
                             growth_rgr=growth$growth_rgr,
                             diameter_start=(growth$diameter_start - dbh_mean) / dbh_sd,
                             diameter_end=(growth$diameter_end - dbh_mean) / dbh_sd,
                             n_days=growth$n_days,
                             elev_diff=(growth$elev_diff - elev_diff_mean) / elev_diff_sd,
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

    plot_ID_factor_levels <- levels(factor(plot_ID))
    plot_ID_factor_key <- cbind(plot_ID_char=as.character(plot_ID_factor_levels), 
                                plot_ID_numeric=seq(1:length(plot_ID_factor_levels)))
    write.csv(plot_ID_factor_key, file=file.path(data_folder, 
                                                 paste0("plot_ID_factor_key", 
                                                        suffix, ".csv")), 
                                                 row.names=FALSE)


    site_ID_factor_levels <- levels(factor(site_ID))
    site_ID_factor_key <- cbind(site_ID_char=as.character(site_ID_factor_levels), 
                                site_ID_numeric=seq(1:length(site_ID_factor_levels)))
    write.csv(site_ID_factor_key, file=file.path(data_folder, 
                                                 paste0("site_ID_factor_key", 
                                                        suffix, ".csv")), 
                                                 row.names=FALSE)

    period_num_factor_levels <- levels(ordered(growth$SamplingPeriodEnd))
    period_num_factor_key <- data.frame(period_num_char=as.character(period_num_factor_levels), 
                                       period_num_numeric=seq(1:length(period_num_factor_levels)),
                                       stringsAsFactors=FALSE)
    # Add in initial period (coded as zero, but not in the SamplingPeriodEnd 
    # vector)
    initial_period <- paste0(as.numeric(substr(period_num_factor_key$period_num_char[1], 1, 4)) - 1, ".01")
    period_num_factor_key <- rbind(c(initial_period, 0), period_num_factor_key)
    write.csv(period_num_factor_key,
              file=file.path(data_folder, paste0("period_num_factor_key", suffix, 
                                                ".csv")), row.names=FALSE)

    ###############################################################################
    # Setup inits

    # Setup wide format dbh_latent dataframes
    dbh_latent <- dcast(dbh_ts, tree_ID + plot_ID + site_ID ~ period_num, value.var="dbh")
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

}

stopCluster(cl)
