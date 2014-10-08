library(dplyr)
library(reshape2)
library(zoo) # for na.approx
library(Rcpp)
library(inline)

sourceCpp('calc_missings.cpp')

load("growth_ctfsflagged_merged_detrended.RData")

growth <- tbl_df(growth)

#TODO Fix WD data
growth <- filter(growth, !is.na(WD))

# table(growth$ctfs_accept)
growth <- filter(growth, ctfs_accept)

###############################################################################
### TESTING ONLY
# growth <- filter(growth, sitecode %in% c("VB", "CAX"))
# growth <- filter(growth, SamplingPeriodEnd %in% c("2011.01", "2012.01", "2013.01"))

###############################################################################

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

# Setup timeseries formatted dataframe (with NAs for covariates at time 1 since 
# first growth measurement can only be calculated from time 1 to time 2)
dbh_ts <- arrange(growth, tree_ID, period_ID) %>%
    select(tree_ID, site_ID, plot_ID, period_ID, genus_ID, spi=spi_24, WD, 
           dbh=diameter_end)

# Setup an array of initial diameters for all the trees, by taking the dbh_st 
# value from the first period with an available observation
dbh_time_0 <- arrange(growth, tree_ID, period_ID) %>%
    group_by(tree_ID) %>%
    filter(period_ID==min(period_ID)) %>%
    select(tree_ID, site_ID, plot_ID, period_ID, genus_ID, spi=spi_24, WD, 
           dbh=diameter_start)
dbh_time_0$spi <- NA
dbh_time_0$period_ID <- dbh_time_0$period_ID - 1

dbh_ts <- rbind(dbh_ts, dbh_time_0)
dbh_ts <- arrange(dbh_ts, tree_ID, period_ID)

# Standardize outcome and predictors.
dbh_mean <- mean(dbh_ts$dbh)
dbh_sd <- sd(dbh_ts$dbh)
dbh_ts$dbh <- (dbh_ts$dbh - dbh_mean) / dbh_sd
WD_mean <- mean(dbh_ts$WD)
WD_sd <- sd(dbh_ts$WD)
WD <- (dbh_time_0$WD - WD_mean) / WD_sd
# Save sd and means so the variables can be unstandardized later
save(dbh_mean, dbh_sd, WD_mean, WD_sd, file="model_data_standardizing.RData")

# Calculate precision of diameter tape (1 mm) in standardized units:
.1 / dbh_sd

genus_ID <- dbh_time_0$genus_ID
sum(genus_ID == "Unknown") / length(genus_ID)

# Add latent growth inits
calc_latent_dbh <- function(dbh) {
    dbh_latent <- rep(NA, length(dbh))
    dbh_latent[1] <- dbh[1]
    for (i in 2:length(dbh)) {
        dbh_latent[i] <- ifelse(dbh[i] > dbh_latent[i - 1],
            dbh[i], dbh_latent[i - 1] + .005)
    }
    return(dbh_latent)
}
dbh_ts <- arrange(dbh_ts, period_ID) %>%
    group_by(tree_ID) %>%
    mutate(dbh_latent=calc_latent_dbh(dbh)) %>%
    arrange(tree_ID, period_ID)

# Setup wide format dbh and dbh_latent dataframes
dbh <- dcast(dbh_ts, tree_ID + plot_ID + site_ID ~ period_ID, value.var="dbh")
dbh_latent <- dcast(dbh_ts, tree_ID + plot_ID + site_ID ~ period_ID, value.var="dbh_latent")
spi <- dcast(dbh_ts, tree_ID + plot_ID + site_ID ~ period_ID, value.var="spi")
tree_ID <- dbh$tree_ID
plot_ID <- dbh$plot_ID
site_ID <- dbh$site_ID
# Eliminate the ID columns
dbh <- dbh[!grepl('_ID$', names(dbh))]
dbh_latent <- dbh_latent[!grepl('_ID$', names(dbh_latent))]
spi <- spi[!grepl('_ID$', names(spi))]

# Interpolate missing values in dbh_latent
interp_dbh_obs <- function(x) {
    first_obs <- min(which(!is.na(x)))
    last_obs <- max(which(!is.na(x)))
    x[first_obs:last_obs] <- na.approx(as.numeric(x[first_obs:last_obs]))
    return(x)
}
dbh_latent <- t(apply(dbh_latent, 1, interp_dbh_obs))

# Setup data
n_tree <- length(unique(dbh_ts$tree_ID))
n_site <- length(unique(dbh_ts$site_ID))
n_plot <- length(unique(dbh_ts$plot_ID))
n_period <- length(unique(dbh_ts$period_ID))
n_genus <- length(unique(dbh_ts$genus_ID))

# Calculate the first and last observation for each tree. The +1's below are 
# because the period_ID variable starts at zero.
obs_per_tree <- group_by(dbh_ts, tree_ID) %>%
    summarize(first_obs_period=(min(period_ID) + 1),
              last_obs_period=(max(period_ID) + 1))

# SPI observations are missing for periods when dbh observations are missing.  
# Fill these observations using the mean SPI for the appropriate period.
spi_means <- group_by(growth, plot_ID, period_ID) %>%
    summarize(spi_means=mean(spi_24, na.rm=TRUE))
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

stopifnot(length(WD) == n_tree)
stopifnot(length(genus_ID) == n_tree)
stopifnot(length(plot_ID) == n_tree)
stopifnot(length(site_ID) == n_tree)

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
                   dbh=as.matrix(dbh),
                   WD=WD,
                   spi=spi,
                   obs_indices=missings$obs,
                   miss_indices=missings$miss)
save(model_data, file="model_data.RData")

test_m <- lm(growth$diameter_end ~ growth$diameter_start + I(growth$diameter_start^2) + 
   growth$WD + I(growth$WD^2) + growth$spi_24 +
   growth$diameter_start * growth$spi_24 +
   growth$WD * growth$spi_24)
#summary(test_m)

# Setup inits
init_data <- list(dbh_latent=as.matrix(dbh_latent))
save(init_data, file="init_data.RData")
