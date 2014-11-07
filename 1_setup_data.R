library(dplyr)
library(reshape2)
library(zoo) # for na.approx
library(Rcpp)
library(inline)

sourceCpp('calc_missings.cpp')

load("growth_ctfsflagged_merged_detrended.RData")

suffixes <- c("", "_testing")

for (suffix in suffixes) {
    growth <- tbl_df(growth)

    #TODO Fix WD data
    growth <- filter(growth, !is.na(WD))

    # table(growth$ctfs_accept)
    growth <- filter(growth, ctfs_accept)

    # Exclude Pasoh until Pasoh plot coordinates are confirmed
    growth <- filter(growth, !(sitecode %in% c("PSH")))

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

    ###############################################################################
    # Subsample dataset fro running testing models
    if (suffix == "_testing") {
        set.seed(1)
        # Generate a subsample including only 5% of the trees included in the 
        # initial sampling period at each site
        initial_trees <- group_by(growth, site_ID) %>%
            filter(period_ID == min(period_ID)) %>%
            group_by(site_ID) %>%
            sample_frac(.05)
        growth <- filter(growth, tree_ID %in% initial_trees$tree_ID)
        growth <- setup_factors(growth)
    }
    ###############################################################################

    # Define maximum cumulative water deficit as a positive number to make 
    # interpretation easier
    growth$cwd <- abs(growth$cwd)
    growth$mcwd12 <- abs(growth$mcwd12)
    growth$cwd_run12 <- abs(growth$cwd_run12)
    growth$mcwd_run12 <- abs(growth$mcwd_run12)

    # Setup timeseries formatted dataframe (with NAs for covariates at time 1 since 
    # first growth measurement can only be calculated from time 1 to time 2)
    dbh_ts <- arrange(growth, tree_ID, period_ID) %>%
        select(tree_ID, site_ID, plot_ID, period_ID, genus_ID, spi=spi_24, 
               mcwd=mcwd_run12, temp=tmn_meanannual, WD, dbh=diameter_end)

    # Setup an array of initial diameters for all the trees, by taking the dbh_st 
    # value from the first period with an available observation
    dbh_time_0 <- arrange(growth, tree_ID, period_ID) %>%
        group_by(tree_ID) %>%
        filter(period_ID==min(period_ID)) %>%
        select(tree_ID, site_ID, plot_ID, period_ID, genus_ID, spi=spi_24, 
               mcwd=mcwd_run12, temp=tmn_meanannual, WD, dbh=diameter_start)
    dbh_time_0$spi <- NA
    dbh_time_0$mcwd <- NA
    dbh_time_0$temp <- NA
    dbh_time_0$period_ID <- dbh_time_0$period_ID - 1

    dbh_ts <- rbind(dbh_ts, dbh_time_0)
    dbh_ts <- arrange(dbh_ts, tree_ID, period_ID)

    #dbh_ts$mcwd <- log(1 + dbh_ts$mcwd)

    # Standardize outcome and predictors.
    dbh_mean <- mean(dbh_ts$dbh)
    dbh_min <- min(dbh_ts$dbh)
    dbh_max <- max(dbh_ts$dbh)
    dbh_sd <- sd(dbh_ts$dbh)
    dbh_ts$dbh <- (dbh_ts$dbh - dbh_mean) / dbh_sd
    mcwd_mean <- mean(dbh_ts$mcwd, na.rm=TRUE)
    mcwd_min <- min(dbh_ts$mcwd, na.rm=TRUE)
    mcwd_max <- max(dbh_ts$mcwd, na.rm=TRUE)
    mcwd_sd <- sd(dbh_ts$mcwd, na.rm=TRUE)
    dbh_ts$mcwd <- (dbh_ts$mcwd - mcwd_mean) / mcwd_sd
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

    # Setup wide format dbh and spi dataframes
    dbh <- dcast(dbh_ts, tree_ID + plot_ID + site_ID ~ period_ID, value.var="dbh")
    spi <- dcast(dbh_ts, tree_ID + plot_ID + site_ID ~ period_ID, value.var="spi")
    mcwd <- dcast(dbh_ts, tree_ID + plot_ID + site_ID ~ period_ID, value.var="mcwd")
    temp <- dcast(dbh_ts, tree_ID + plot_ID + site_ID ~ period_ID, value.var="temp")
    tree_ID <- dbh$tree_ID
    plot_ID <- dbh$plot_ID
    site_ID <- dbh$site_ID
    # Eliminate the ID columns
    dbh <- dbh[!grepl('_ID$', names(dbh))]
    spi <- spi[!grepl('_ID$', names(spi))]
    mcwd <- mcwd[!grepl('_ID$', names(mcwd))]
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
         mcwd_mean, mcwd_min, mcwd_max, mcwd_sd, 
         file=paste0("model_data_standardizing", suffix, ".RData"))

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

    # mcwd observations are missing for periods when dbh observations are 
    # missing.  Fill these observations using the mean mcwd for the appropriate 
    # period.
    #
    # Note that the same centering and scaling needs to be done here (and log) that 
    # was done for the other mcwd data.
    # mcwd_means <- group_by(growth, plot_ID, period_ID) %>%
    #     summarize(mcwd_means=mean((log(1 + mcwd_run12) - mcwd_mean) / mcwd_sd, 
    #     na.rm=TRUE))
    mcwd_means <- group_by(growth, plot_ID, period_ID) %>%
         summarize(mcwd_means=mean((mcwd_run12 - mcwd_mean) / mcwd_sd, na.rm=TRUE))
    mcwd_means <- data.frame(mcwd_means) # Fix for indexing bug in dplyr 0.3.2
    mcwd <- as.matrix(mcwd)
    mcwd_missings <- calc_missings(mcwd)$miss
    # Use linear indexing to replace NAs in mcwds
    mcwd_miss_linear_ind <- (mcwd_missings[, 2] - 1) * nrow(mcwd) + mcwd_missings[, 1] # From http://bit.ly/1rnKrC3
    mcwd[mcwd_miss_linear_ind] <- mcwd_means[match(paste(plot_ID[mcwd_missings[, 1]], mcwd_missings[, 2]),
                                                paste(mcwd_means$plot_ID, mcwd_means$period_ID)), 3]
    stopifnot(is.null(calc_missings(mcwd)$miss))

    # Fill in temp (same issue as above)
    temp_means <- group_by(growth, plot_ID, period_ID) %>%
         summarize(temp_means=mean((tmn_meanannual - temp_mean) / temp_sd, na.rm=TRUE))
    temp_means <- data.frame(temp_means) # Fix for indexing bug in dplyr 0.3.2
    temp <- as.matrix(temp)
    temp_missings <- calc_missings(temp)$miss
    # Use linear indexing to replace NAs in temps
    temp_miss_linear_ind <- (temp_missings[, 2] - 1) * nrow(temp) + temp_missings[, 1] # From http://bit.ly/1rnKrC3
    temp[temp_miss_linear_ind] <- temp_means[match(paste(plot_ID[temp_missings[, 1]], temp_missings[, 2]),
                                                paste(temp_means$plot_ID, temp_means$period_ID)), 3]
    stopifnot(is.null(calc_missings(temp)$miss))

    # SPI observations are missing for periods when dbh observations are missing.  
    # Fill these observations using the mean SPI for the appropriate period.
    spi_means <- group_by(growth, plot_ID, period_ID) %>%
        summarize(spi_means=mean(spi_24, na.rm=TRUE))
    spi_means <- data.frame(spi_means) # Fix for indexing bug in dplyr 0.3.2
    spi <- as.matrix(spi)
    spi_missings <- calc_missings(spi)$miss
    # Use linear indexing to replace NAs in SPIs
    spi_miss_linear_ind <- (spi_missings[, 2] - 1) * nrow(spi) + spi_missings[, 1] # From http://bit.ly/1rnKrC3
    spi[spi_miss_linear_ind] <- spi_means[match(paste(plot_ID[spi_missings[, 1]], spi_missings[, 2]),
                                                paste(spi_means$plot_ID, spi_means$period_ID)), 3]
    stopifnot(is.null(calc_missings(spi)$miss))

    # Calculate indices of missing and observed data (using the function defined in 
    # calc_missings.cpp), so that observed and missing data can be modeled 
    # separately in Stan.
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
                       spi=spi,
                       elev=elev,
                       temp=temp,
                       mcwd=mcwd,
                       obs_indices=missings$obs,
                       miss_indices=missings$miss)
    save(model_data, file=paste0("model_data_wide", suffix, ".RData"))

    ###############################################################################
    # Output long format data
    model_data_long <- data.frame(tree_ID=factor(tree_ID),
                                  plot_ID=factor(plot_ID),
                                  site_ID=factor(site_ID),
                                  genus_ID=factor(genus_ID),
                                  sigma_obs_lower=sigma_obs_lower,
                                  WD=WD)
    # Standardize the diameter variables and mcwd variable the same way they are 
    # standardized in the wide dataset
    merge_data <- data.frame(tree_ID=growth$tree_ID,
                             growth_rgr=growth$growth_rgr,
                             diameter_start=(growth$diameter_start - dbh_mean) / dbh_sd,
                             diameter_end=(growth$diameter_end - dbh_mean) / dbh_sd,
                             n_days=growth$n_days,
                             elev=(growth$elev - elev_mean) / elev_sd,
                             temp=(growth$tmn_meanannual - temp_mean) / temp_sd,
                             mcwd=(growth$mcwd_run12 - mcwd_mean) / mcwd_sd)
    model_data_long <- merge(model_data_long, merge_data, by="tree_ID", all=TRUE)
    save(model_data_long, file=paste0("model_data_long", suffix, ".RData"))

    genus_ID_factor_levels <- levels(factor(genus_ID))
    genus_ID_factor_key <- cbind(genus_ID_char=as.character(genus_ID_factor_levels), 
                                 genus_ID_numeric=seq(1:length(genus_ID_factor_levels)))
    write.csv(genus_ID_factor_key, file=paste0("genus_ID_factor_key", suffix, 
                                              ".csv"), row.names=FALSE)

    site_ID_factor_levels <- levels(factor(site_ID))
    site_ID_factor_key <- cbind(site_ID_char=as.character(site_ID_factor_levels), 
                                site_ID_numeric=seq(1:length(site_ID_factor_levels)))
    write.csv(site_ID_factor_key, file=paste0("site_ID_factor_key", suffix, 
                                              ".csv"), row.names=FALSE)

    period_ID_factor_levels <- levels(ordered(growth$SamplingPeriodEnd))
    period_ID_factor_key <- data.frame(period_ID_char=as.character(period_ID_factor_levels), 
                                       period_ID_numeric=seq(1:length(period_ID_factor_levels)),
                                       stringsAsFactors=FALSE)
    # Add in initial period (coded as zero, but not in the SamplingPeriodEnd 
    # vector)
    initial_period <- paste0(as.numeric(substr(period_ID_factor_key$period_ID_char[1], 1, 4)) - 1, ".01")
    period_ID_factor_key <- rbind(c(initial_period, 0), period_ID_factor_key)
    write.csv(period_ID_factor_key,
              file=paste0("period_ID_factor_key", suffix, ".csv"), 
              row.names=FALSE)

    ###############################################################################
    # Setup inits

    # test_m <- lm(growth$diameter_end ~ growth$diameter_start + 
    # I(growth$diameter_start^2) + growth$WD + I(growth$WD^2) + growth$spi_24 +
    #    growth$diameter_start * growth$spi_24 +
    #    growth$WD * growth$spi_24)
    # #summary(test_m)

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
    save(init_data, file=paste0("init_data", suffix, ".RData"))
}
