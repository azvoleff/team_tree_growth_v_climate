library(dplyr)
library(reshape2)
library(zoo) # for na.approx

load("growth_ctfsflagged_merged_detrended.RData")

growth <- tbl_df(growth)

table(growth$ctfs_accept)
table(growth$n_days < 200)
table(growth$n_days > 550)
table(growth$n_days > 1000)

growth <- filter(growth, ctfs_accept)
growth <- filter(growth, n_days > 200)
growth <- filter(growth, n_days < 550)

# Exclude NAK, CSN, and YAN - too little data at these sites
#growth <- filter(growth, !(sitecode %in% c("NAK", "CSN", "YAN")))

growth$SamplingPeriodID <- with(growth, factor(factor(sitecode):factor(SamplingPeriodNumber)))

###############################################################################
### TESTING ONLY
#growth <- filter(growth, sitecode %in% c("VB", "CAX"))
###############################################################################

growth$ID_tree <- as.integer(factor(growth$SamplingUnitName))
growth$ID_plot <- as.integer(factor(growth$plot_ID))
growth$ID_site <- as.integer(factor(growth$sitecode))
growth$ID_period <- as.integer(factor(growth$SamplingPeriodNumber))

# Setup timeseries formatted dataframe (with NAs for covariates at time 1 since 
# first growth measurement can only be calculated from time 1 to time 2)
dbh_ts <- arrange(growth, ID_period) %>%
    group_by(ID_tree) %>%
    do(data.frame(ID_tree=c(.$ID_tree[1], .$ID_tree),
                  ID_plot=c(.$ID_plot[1], .$ID_plot),
                  ID_site=c(.$ID_site[1], .$ID_site),
                  ID_period=c(1, .$ID_period + 1), # Start period IDs at 1
                  obs_num=seq(0, length(.$dbh_st)),
                  spi=c(NA, .$spi_24),
                  WD=c(NA, .$WD),
                  dbh=c(.$diameter_start[1], .$diameter_end))) %>%
    arrange(ID_tree, ID_period)

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

WD <- growth$WD[match(ID_tree, growth$ID_tree)]
WD <- (WD - mean(WD)) / sd(WD)

dbh_missing <- is.na(dbh)

# Check how many trees have missing observations in the middle of their 
# observation periods (with each tree's observation period defined by the 
# periods with the first and last known diameter measurements).
calc_n_missings <- function(x) {
    first_obs <- min(which(!is.na(x)))
    last_obs <- max(which(!is.na(x)))
    sum(is.na(x[first_obs:last_obs]))
}
n_missing_periods <- t(apply(dbh, 1, calc_n_missings))
table(n_missing_periods)

# TEMPORARY: interpolate missing values in dbh (TODO: model NAs in Stan/JAGS)
interp_dbh_obs <- function(x) {
    first_obs <- min(which(!is.na(x)))
    last_obs <- max(which(!is.na(x)))
    x[first_obs:last_obs] <- na.approx(as.numeric(x[first_obs:last_obs]))
    return(x)
}
dbh <- t(apply(dbh, 1, interp_dbh_obs))

# Interpolate missing values in dbh_latent
dbh_latent <- t(apply(dbh_latent, 1, interp_dbh_obs))

# Setup data
n_tree <- length(unique(dbh_ts$ID_tree))
n_site <- length(unique(dbh_ts$ID_site))
n_plot <- length(unique(dbh_ts$ID_plot))

obs_per_tree <- group_by(dbh_ts, ID_tree) %>%
    summarize(first_obs_period=min(ID_period),
              last_obs_period=max(ID_period))

model_data <- list(n_tree=n_tree,
                   n_plot=n_plot,
                   n_site=n_site,
                   first_obs_period=obs_per_tree$first_obs_period,
                   last_obs_period=obs_per_tree$last_obs_period,
                   n_periods=ncol(dbh),
                   plot_ID=as.integer(factor(ID_plot)),
                   site_ID=as.integer(factor(ID_site)),
                   log_dbh=as.matrix(log(dbh)),
                   WD=WD,
                   spi=spi)
save(model_data, file="model_data.RData")

#var(log(dbh_ts$dbh))
# Setup inits
init_data <- list(list(log_dbh_latent=as.matrix(log(dbh_latent)),
                       inter=-1.7,
                       slp_dbh=.008,
                       slp_spi=.1,
                       slp_spi_sq=-.1,
                       slp_WD=.1,
                       slp_WD_sq=-.1,
                       sigma_obs=.02, 
                       sigma_proc=.02, 
                       sigma_ijk=.06,
                       sigma_jk=.12,
                       sigma_k=.24,
                       b_ijk_std=rep(.5, n_tree),
                       b_jk_std=rep(.3, n_plot),
                       b_k_std=rep(.2, n_site)))
save(init_data, file="init_data.RData")
