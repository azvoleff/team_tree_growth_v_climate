library(dplyr)
library(reshape2)
library(zoo) # for na.approx
library(Rcpp)
library(inline)

sourceCpp('calc_missings.cpp')

load("growth_ctfsflagged_merged_detrended.RData")

lm(growth$diameter_end ~ growth$diameter_start + I(growth$diameter_start^2))

growth <- tbl_df(growth)

#TODO Fix WD data
growth <- filter(growth, !is.na(WD))

# table(growth$ctfs_accept)
# table(growth$n_days < 200)
# table(growth$n_days > 550)
# table(growth$n_days > 1000)

growth <- filter(growth, ctfs_accept)
# growth <- filter(growth, n_days > 200)
# growth <- filter(growth, n_days < 550)

# Exclude NAK, CSN, and YAN - too little data at these sites
#growth <- filter(growth, !(sitecode %in% c("NAK", "CSN", "YAN")))

###############################################################################
# Fix issue with having two periods in 2005 

# Check how many trees were observed in both 2005.01 and 2005.02
# table(growth$SamplingPeriodStart == "2005.01")
# table(growth$SamplingPeriodEnd == "2005.02")
# table(growth$SamplingPeriodStart == "2005.01" & growth$SamplingPeriodEnd == "2005.02")
# table(growth$SamplingPeriodStart == "2005.01" & growth$SamplingPeriodEnd == "2006.01")
# table(growth$SamplingPeriodStart == "2005.01" & growth$SamplingPeriodEnd == "2006.01")
# Eliminate the 2005.01 sampling period end dbhs into the 2005.02 period start 
# dbhs.
# dbh_2005_01 <- growth$SamplingPeriodStart == "2005.01"
# dbh_2005_02 <- growth$SamplingPeriodStart == "2005.02"
# growth[dbh_2005_01, ]$diameter_end <- NA
# growth[dbh_2005_01, ]$diameter_end <- growth[dbh_2005_02, ]$diameter_end[match(growth[dbh_2005_01, ]$SamplingUnitName,
#                                                                                growth[dbh_2005_02, ]$SamplingUnitName)]
# growth[dbh_2005_01, ]$SamplingPeriodEnd <- "2006.01"
# # Now remove the no longer needed 2005.02 rows
# growth <- filter(growth, SamplingPeriodStart != "2005.02")
# growth <- filter(growth, SamplingPeriodEnd != "2005.02")
#
# stopifnot(sum(growth$SamplingPeriodStart == "2005.02") == 0)
# stopifnot(sum(growth$SamplingPeriodEnd == "2005.02") == 0)

###############################################################################
### TESTING ONLY
#growth <- filter(growth, sitecode %in% c("VB", "CAX"))
###############################################################################

growth$ID_tree_char <- factor(growth$SamplingUnitName)
growth$ID_plot_char <- factor(growth$plot_ID)
growth$ID_site_char <- factor(growth$sitecode)
growth$ID_period_char <- factor(growth$SamplingPeriodEnd)
growth$ID_tree <- as.integer(factor(growth$SamplingUnitName))
growth$ID_plot <- as.integer(factor(growth$plot_ID))
growth$ID_site <- as.integer(factor(growth$sitecode))
# The plus one below is to ensure that the first period (rbinded on below) will 
# have period number equal to 1, which will allow period numbers to be used in 
# JAGS and R with 1 based indexing.
growth$ID_period <- as.integer(ordered(growth$SamplingPeriodEnd))

# Setup timeseries formatted dataframe (with NAs for covariates at time 1 since 
# first growth measurement can only be calculated from time 1 to time 2)
dbh_ts <- arrange(growth, ID_tree, ID_period) %>%
    select(ID_tree, ID_site, ID_plot, ID_period, spi=spi_24, WD, 
           dbh=diameter_end)

# Setup an array of initial diameters for all the trees, by taking the dbh_st 
# value from the first period with an available observation
dbh_time_0 <- arrange(growth, ID_tree, ID_period) %>%
    group_by(ID_tree) %>%
    filter(ID_period==min(ID_period)) %>%
    select(ID_tree, ID_site, ID_plot, ID_period, spi=spi_24, WD, 
           dbh=diameter_start)
dbh_time_0$spi <- NA
dbh_time_0$ID_period <- dbh_time_0$ID_period - 1

dbh_ts <- rbind(dbh_ts, dbh_time_0)
dbh_ts <- arrange(dbh_ts, ID_tree, ID_period)

# Add latent growth inits
calc_latent_dbh <- function(dbh) {
    dbh_latent <- rep(NA, length(dbh))
    dbh_latent[1] <- dbh[1]
    for (i in 2:length(dbh)) {
        dbh_latent[i] <- ifelse(dbh[i] > dbh_latent[i - 1],
            dbh[i], dbh_latent[i - 1] + .01)
    }
    return(dbh_latent)
}
dbh_ts <- arrange(dbh_ts, ID_period) %>%
    group_by(ID_tree) %>%
    mutate(dbh_latent=calc_latent_dbh(dbh)) %>%
    arrange(ID_tree, ID_period)

# Setup wide format dbh and dbh_latent dataframes
dbh <- dcast(dbh_ts, ID_tree + ID_plot + ID_site ~ ID_period, value.var="dbh")
dbh_latent <- dcast(dbh_ts, ID_tree + ID_plot + ID_site ~ ID_period, value.var="dbh_latent")
spi <- dcast(dbh_ts, ID_tree + ID_plot + ID_site ~ ID_period, value.var="spi")
ID_tree <- dbh$ID_tree
ID_plot <- dbh$ID_plot
ID_site <- dbh$ID_site
# Eliminate the ID columns
dbh <- dbh[!grepl('^ID_', names(dbh))]
dbh_latent <- dbh_latent[!grepl('^ID_', names(dbh_latent))]
spi <- spi[!grepl('^ID_', names(spi))]

# Interpolate missing values in dbh_latent
interp_dbh_obs <- function(x) {
    first_obs <- min(which(!is.na(x)))
    last_obs <- max(which(!is.na(x)))
    x[first_obs:last_obs] <- na.approx(as.numeric(x[first_obs:last_obs]))
    return(x)
}
dbh_latent <- t(apply(dbh_latent, 1, interp_dbh_obs))

WD <- dbh_time_0$WD
WD <- (WD - mean(WD)) / sd(WD)

# Setup data
n_tree <- length(unique(dbh_ts$ID_tree))
n_site <- length(unique(dbh_ts$ID_site))
n_plot <- length(unique(dbh_ts$ID_plot))

# Calculate the first and last observation for each tree. The +1's below are 
# because the ID_period variable starts at zero.
obs_per_tree <- group_by(dbh_ts, ID_tree) %>%
    summarize(first_obs_period=(min(ID_period) + 1),
              last_obs_period=(max(ID_period) + 1))

# SPI observations are missing for periods when dbh observations are missing.  
# Fill these observations using the mean SPI for the appropriate period.
spi_means <- group_by(growth, ID_plot, ID_period) %>%
    summarize(spi_means=mean(spi_24, na.rm=TRUE))
spi <- as.matrix(spi)
spi_missings <- calc_missings(spi)$miss
# Use linear indexing to replace NAs in SPIs
spi_miss_linear_ind <- (spi_missings[, 2] - 1) * nrow(spi) + spi_missings[, 1] # From http://bit.ly/1rnKrC3
spi[spi_miss_linear_ind] <- spi_means[match(paste(ID_plot[spi_missings[, 1]], spi_missings[, 2]),
                                            paste(spi_means$ID_plot, spi_means$ID_period)), 3]
stopifnot(is.null(calc_missings(spi)$miss))

# # Some SPI measurements are still missing for period 3 (2005.01) from site 4 
# # (CAX).  Fill these in with correct values.
# table(growth[growth$sitecode == "CAX", ]$period_end_month, growth[growth$sitecode == "CAX", ]$SamplingPeriodEnd)
# load('../CHIRPS/vg_plot_spis.RData')
# spis$plot_ID <- as.character(spis$plot_ID)
# spis$period_end_month <- spis$date
#
# spi_missings <- calc_missings(spi)$miss
# ID_plot_char <- levels(growth$ID_plot_char)[ID_plot]
# for (n in 1:nrow(spi_missings)) {
#     this_tree_num <- spi_missings[n, 1]
#     this_period <- spi_missings[n, 2]
#     this_plot_ID_char <- ID_plot_char[this_tree_num]
#     # Get most common period
#     tmp <- table(growth[growth$ID_plot_char == this_plot_ID_char, ]$period_end_month)
#     names(tmp)[tmp == max(tmp)]
#     growth[(growth$ID_plot_char == this_plot_ID_char), ]$period_end_month
# }
#
#
# SPI_period_3_site_4 <- filter(spis, spi_period == 24, period_end_month == "2005-12-1") %>%
#     select(plot_ID, spi)
# SPI_period_3_site_4 <- match(levels(growth$ID_plot_char)
#
# spi_24 <- filter(spis, spi_period == 24, period_end_month > "2002-1-1") %>%
#     select(plot_ID, period_end_month, spi)
#
# names(spi_24)[names(spi_24) == "spi"] <- "spi_24"
#
# levels(growth$sitecode)[ID_site[calc_missings(spi)$miss[1, 1]]]
# levels(ordered(growth$SamplingPeriodStart))[3]
# ID_plot[calc_missings(spi)$miss[1, 1]]
# calc_missings(spi)$miss[1, 2]
# # plot 21, period 3
# spi_means[spi_means$ID_plot == 21, ]
#
# table(ID_plot[calc_missings(spi)$miss[, 1]])
# table(calc_missings(spi)$miss[, 2])
#
# Calculate indices of missing and observed data (using the function defined in 
# calc_missings.cpp), so that observed and missing data can be modeled 
# separately.

missings <- calc_missings(as.matrix(dbh))

model_data <- list(n_tree=n_tree,
                   n_plot=n_plot,
                   n_site=n_site,
                   first_obs_period=obs_per_tree$first_obs_period,
                   last_obs_period=obs_per_tree$last_obs_period,
                   plot_ID=as.integer(factor(ID_plot)),
                   site_ID=as.integer(factor(ID_site)),
                   dbh=as.matrix(dbh),
                   WD=WD,
                   spi=spi,
                   obs_indices=missings$obs,
                   miss_indices=missings$miss)
save(model_data, file="model_data.RData")

# Setup inits
init_data <- list(list(dbh_latent=as.matrix(dbh_latent),
                       inter=-.27,
                       slp_dbh=1.05,
                       slp_dbh_sq=-.001,
                       slp_spi=.1,
                       slp_spi_sq=-.1,
                       slp_WD=.1,
                       slp_WD_sq=-.1,
                       sigma_obs=2, 
                       sigma_proc=10, 
                       sigma_ijk=2,
                       sigma_jk=5,
                       sigma_k=10))
save(init_data, file="init_data.RData")
