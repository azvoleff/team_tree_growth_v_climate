library(ggplot2)
library(stringr)
library(scales) # for alpha()
library(grid)
library(runjags)
library(coda)
library(dplyr)
library(reshape2)
library(foreach)
library(doParallel)

source("settings.R")

cl <- makeCluster(3)
registerDoParallel(cl)

model_type <- "full"
#model_type <- "testing"
precip_var <- "mcwd_run12"

temp_var <- "tmx_meanannual"
run_ID <- "vertica1.team.sdsc.edu_20141110224426_extend3" 

suffix <- paste0(model_type, '-', temp_var, '-', precip_var)

params_folder <- file.path(prefix, "TEAM", "Tree_Growth", "Extracted_Parameters")
data_folder <- file.path(prefix, "TEAM", "Tree_Growth", "Data")

load(file.path(data_folder, paste0("model_data_wide", suffix, ".RData")))
load(file.path(data_folder, paste0("model_data_standardizing_", suffix, ".RData")))

load(file.path(params_folder, paste0(suffix, "jags_fit_full_model_fixefs.RData")))

n_site <- model_data$n_site
n_genus <- model_data$n_genus

merged <- with(model_data, data.frame(site_ID=site_ID, 
                                      genus_ID=as.integer(genus_ID),
                                      WD=WD,
                                      dbh=apply(dbh, 1, mean, na.rm=TRUE)))
dbhs <- c(10, 20, 30, 40, 50, 75, 100, 150)
dbhs_std <- (dbhs - dbh_mean)/dbh_sd
dbh_class_midpoints <- diff(dbhs_std)/2 + dbhs_std[-length(dbhs_std)]
merged$dbh_class_ID <- as.numeric(cut(merged$dbh, dbhs_std, include.lowest=TRUE))

d_dbh <- 10
dbhs_smalld <- seq(10, 150, d_dbh)
dbhs_smalld_std <- (dbhs_smalld - dbh_mean)/dbh_sd
dbh_smalld_class_midpoints <- diff(dbhs_smalld_std)/2 + dbhs_smalld_std[-length(dbhs_smalld_std)]

# Check overlap of genera between sites
genera_overlap <- group_by(merged, genus_ID) %>%
    summarize(n=n(),
              n_sites=length(unique(site_ID))) %>%
    arrange(desc(n_sites))
genera_overlap$genus_ID_ordered <- with(genera_overlap,
                                        ordered(genus_ID, levels=genus_ID[order(n_sites, decreasing=TRUE)]))
genera_overlap$genus_ID_ordered <- as.numeric(genera_overlap$genus_ID_ordered)
ggplot(genera_overlap) +
    geom_bar(aes(genus_ID_ordered, n_sites), colour="black", fill="black", stat="identity") +
    xlab("Genus ID number") +
    ylab("Number of sites with observations")

genus_stats_overall <- group_by(merged, genus_ID) %>%
    summarize(n=n()) %>%
    arrange(desc(n))

genus_stats <- group_by(merged, site_ID, genus_ID) %>%
    summarize(min_obs_dbh_std=min(dbh),
              min_obs_dbh_std=min(dbh),
              min_obs_dbh=min(dbh)*dbh_sd + dbh_mean,
              max_obs_dbh=max(dbh)*dbh_sd + dbh_mean,
              n=n()) %>%
    group_by(site_ID) %>%
    mutate(freq=n/sum(n)) %>%
    group_by() %>%
    arrange(desc(freq))

###############################################################################
# Load MCMC samples

# Fixed effects on density
# Combine the MCMC chains and convert to data.frame with coefficients in 
    # columns and MCMC samples in rows
fixefs <- as.data.frame(combine.mcmc(fixefs))
n_mcmc <- nrow(fixefs)
B <- fixefs[grepl('^B\\[', names(fixefs))]
B <- B[order(names(B))]
B <- as.matrix(B)
n_B <- ncol(B)

# Random intercepts at plot, site, and period levels
load(file.path(params_folder, paste0(suffix, "jags_fit_full_model_ranefs.RData")))
ranefs <- as.data.frame(combine.mcmc(ranefs))
# Combine the MCMC chains and convert to data.frames with coefficients in 
# columns and MCMC samples in rows
int_jk <- ranefs[grepl('^int_jk\\[', names(ranefs))]
int_jk <- int_jk[order(names(int_jk))]
int_k <- ranefs[grepl('^int_k\\[', names(ranefs))]
int_k <- int_k[order(names(int_k))]
int_t <- ranefs[grepl('^int_t\\[', names(ranefs))]
int_t <- int_t[order(names(int_t))]

# Random effects at genus-level
load(file.path(params_folder, paste0(suffix, "jags_fit_full_model_ranefs_B_g.RData")))
B_g <- combine.mcmc(ranefs_B_g)
n_B_g <- ncol(B_g)/n_genus

# Make a series of vectors to check that the array conversion is done 
# correctly.
# 
# First vector is the set of coefficients from the 1st mcmc sample for genus 1
B_g_set1 <- B_g[1, c(0:(n_B_g - 1)) * n_genus+ 1]
# Second vector is the set of coefficients from the 5th mcmc sample for genus 7
B_g_set2 <- B_g[5, c(0:(n_B_g - 1)) * n_genus + 7]
# Now convert to a 3D array, with genera in rows, coefficients in columns, and 
# separate MCMC samples in third dimension
dim(B_g)
B_g_array <- t(matrix(t(B_g), byrow=TRUE, ncol=n_genus))
dim(B_g_array)
B_g_array <- array(B_g_array, dim=c(n_genus, n_B_g, n_mcmc))
dim(B_g_array)
stopifnot(all(B_g_array[1, , 1] == B_g_set1))
stopifnot(all(B_g_array[7, , 5] == B_g_set2))

###############################################################################
# Make predictions of mean genus-level growth
WD <- group_by(merged, genus_ID) %>%
    summarize(WD=WD[1])
# Construct fixed effects design matrix
x <- cbind(WD$WD, WD$WD^2)
# Construct genus-level random effects design matrix
#
x_g <- cbind(1, # genus-level intercept
             rep(0, nrow(WD)), # MCWD
             rep(0^2, nrow(WD)), # MCWD
             rep(dbh_class_midpoints[1], nrow(WD)),
             rep(dbh_class_midpoints[1]^2, nrow(WD)))
# Combine the two matrices:
x_all <- cbind(x, x_g)

# Convert B into a n_genus x n_beta x n_mcmc array. (This repeats the same B 
# for each genus row since B is a fixed effect).
B_rep <- matrix(rep(B, each=n_genus), ncol=2)
B_rep <- array(B_rep, dim=c(n_genus, n_B, n_mcmc))

###############################################################################
# Now make growth curves for each genus
# B <- array(B, dim=c(1, ncol(B), nrow(B)))

# dim(aperm(array(B), c(1)))
# B <- matrix(rep(t(B), nrow(WD)), ncol=ncol(B), byrow=TRUE)

# genera_IDs <- unique(this_genus_set$genus_ID)
# dbhs <- this_dbh
# mcwd <- 0

# Below function predictions the growth by genus for each of a set of MCMC 
# samples
pred_growth <- function(genera_IDs=1:n_genus, dbhs=dbh_class_midpoints,
                        mcwd=0) {
    n_genus <- length(genera_IDs)
    # Make predictions into a genus x dbh x n_mcmc array.
    growth_preds <- array(NA, dim=c(n_genus, length(dbhs), n_mcmc))
    for (n in 1:n_genus) {
        genus_ID <- genera_IDs[n]
        # Collapse betas for this genus into flat matrices (temporarily eliminating 
        # third dimension for efficient matrix multiplication)
        B_flat <- matrix(B_rep[genus_ID, , ], nrow=dim(B_rep)[3])
        B_g_flat <- matrix(B_g_array[genus_ID, , ], nrow=dim(B_g_array)[3], byrow=TRUE)
        B_all <- cbind(B_flat, B_g_flat)
        this_x_all <- matrix(rep(x_all[genus_ID, ], length(dbhs)), ncol=n_B + n_B_g, byrow=TRUE)
        this_x_all[, 4] <- mcwd 
        this_x_all[, 5] <- mcwd^2
        this_x_all[, 6] <- dbhs 
        this_x_all[, 7] <- dbhs^2
        growth_preds_flat <- this_x_all %*% t(B_all) - rep(dbhs, nrow(B_all))
        growth_preds[n, , ] <- growth_preds_flat
    }
    return(growth_preds)
}

# Below function takes output of pred_growth to calculate a genus-level mean 
# growth with quantiles
genus_mean_growth <- function(genera_IDs=1:n_genus, dbhs=dbh_class_midpoints, 
                              mcwds=0) {
    n_genus <- length(genera_IDs)
    growths <- foreach (mcwd=mcwds, .combine=rbind,
                        .export=c("pred_growth", "n_genus", "n_mcmc", "B_rep", 
                                  "B_g_array", "x_all", "n_B", "n_B_g", 
                                  "dbh_class_midpoints", "dbh_sd"),
                        .packages=c("reshape2", "dplyr")) %dopar% {
        growth_preds <- pred_growth(dbhs=dbhs, mcwd=mcwd)
        # Return to original units
        growth_preds <- growth_preds * dbh_sd
        growth_means <- apply(growth_preds, c(1, 2), mean)
        growth_q2pt5 <- apply(growth_preds, c(1, 2), quantile, .025)
        growth_q97pt5 <- apply(growth_preds, c(1, 2), quantile, .975)
        growth_means <- melt(growth_means, varnames=c("genus_ID", "dbh_class_ID"), 
                             value.name="mean")
        growth_q2pt5 <- melt(growth_q2pt5, varnames=c("genus_ID", "dbh_class_ID"), 
                             value.name="q2pt5")
        growth_q97pt5 <- melt(growth_q97pt5, varnames=c("genus_ID", "dbh_class_ID"), 
                              value.name="q97pt5")
        data.frame(cbind(growth_means, mcwd=mcwd, q2pt5=growth_q2pt5$q2pt5, 
                         q97pt5=growth_q97pt5$q97pt5))
    }
    return(growths)
}

# Calculate growth for 0, the mean, and the mean plus 2 and plus 3 SD
growth <- genus_mean_growth(mcwds=(c(0, 185, 185+(150*2), 185+(150*3))- 
                                   mcwd_mean)/mcwd_sd, 
                            dbhs=dbh_smalld_class_midpoints)

genus_ID_factor_key <- read.csv(file.path(data_folder,
                                         paste0("genus_ID_factor_key_", suffix, 
                                                ".csv")))
growth <- merge(growth, genus_ID_factor_key, by.x="genus_ID", by.y="genus_ID_numeric")
growth$dbh <- dbh_smalld_class_midpoints[growth$dbh_class_ID]*dbh_sd + dbh_mean
growth$mcwd <- (growth$mcwd*mcwd_sd) + mcwd_mean
growth <- tbl_df(growth)
growth <- arrange(growth, genus_ID, dbh_class_ID, mcwd)

# Plot subset of growth curves
ggplot(filter(growth, genus_ID %in% c(833, 566, 784, 261), mcwd == 0)) +
    geom_line(aes(x=dbh, y=mean, colour=genus_ID_char))+
    geom_ribbon(aes(x=dbh, ymin=q2pt5, ymax=q97pt5, fill=genus_ID_char), alpha=.25)+
    coord_cartesian(xlim=c(10, 150), ylim=c(0, 3)) +
    scale_fill_discrete("Genus") +
    scale_colour_discrete("Genus") +
    xlab("Initial diameter (cm)") +
    ylab("Growth (cm/yr)") +
    theme(legend.text=element_text(face="italic"))
ggsave("growth_genus_level_subset.png", width=4, height=2, dpi=300)

# Plot ALL growth curves
ggplot(filter(growth, mcwd == 0)) +
    geom_line(aes(x=dbh, y=mean, group=genus_ID), alpha=.2, size=.25)+
    #geom_ribbon(aes(x=dbh, ymin=q2pt5, ymax=q97pt5, group=genus_ID), 
    #alpha=.05)+
    coord_cartesian(xlim=c(10, 150), ylim=c(0, 3)) +
    xlab("Initial diameter (cm)") +
    ylab("Growth (cm/yr)")
ggsave("growth_genus_level_all.png", width=6, height=4, dpi=300)

# Plot specific curve versus MCWD
ggplot(filter(growth, genus_ID %in% c(784))) +
    geom_line(aes(x=dbh, y=mean, colour=factor(mcwd)))+
    geom_ribbon(aes(x=dbh, ymin=q2pt5, ymax=q97pt5, fill=factor(mcwd)), alpha=.25)+
    coord_cartesian(xlim=c(10, 150), ylim=c(0, 3)) +
    scale_fill_discrete("MCWD") +
    scale_colour_discrete("MCWD") +
    xlab("Initial diameter (cm)") +
    ylab("Growth (cm/yr)") +
    theme(legend.text=element_text(face="italic"))
ggsave("growth_genus_level_Syzygium.png", width=4, height=2, dpi=300)

# Plot specific curve versus MCWD
ggplot(filter(growth, genus_ID %in% c(261))) +
    geom_line(aes(x=dbh, y=mean, colour=factor(mcwd)))+
    geom_ribbon(aes(x=dbh, ymin=q2pt5, ymax=q97pt5, fill=factor(mcwd)), alpha=.25)+
    coord_cartesian(xlim=c(10, 150), ylim=c(0, 3)) +
    scale_fill_discrete("MCWD") +
    scale_colour_discrete("MCWD") +
    xlab("Initial diameter (cm)") +
    ylab("Growth (cm/yr)") +
    theme(legend.text=element_text(face="italic"))
ggsave("growth_genus_level_dipterocarpus.png", width=4, height=2, dpi=300)

# Plot specific curve versus MCWD
ggplot(filter(growth, genus_ID %in% c(566))) +
    geom_line(aes(x=dbh, y=mean, colour=factor(mcwd)))+
    geom_ribbon(aes(x=dbh, ymin=q2pt5, ymax=q97pt5, fill=factor(mcwd)), alpha=.25)+
    coord_cartesian(xlim=c(10, 150), ylim=c(0, 3)) +
    scale_fill_discrete("MCWD") +
    scale_colour_discrete("MCWD") +
    xlab("Initial diameter (cm)") +
    ylab("Growth (cm/yr)") +
    theme(legend.text=element_text(face="italic"))
ggsave("growth_genus_level_oenocarpus.png", width=4, height=2, dpi=300)

# Plot specific curve versus MCWD
ggplot(filter(growth, genus_ID %in% c(833))) +
    geom_line(aes(x=dbh, y=mean, colour=factor(mcwd)))+
    geom_ribbon(aes(x=dbh, ymin=q2pt5, ymax=q97pt5, fill=factor(mcwd)), alpha=.25)+
    coord_cartesian(xlim=c(10, 150), ylim=c(0, 3)) +
    scale_fill_discrete("MCWD") +
    scale_colour_discrete("MCWD") +
    xlab("Initial diameter (cm)") +
    ylab("Growth (cm/yr)") +
    theme(legend.text=element_text(face="italic"))
ggsave("growth_genus_level_unknown.png", width=4, height=2, dpi=300)

###############################################################################
# Predict site-level mean growth based on observed frequencies of each genus at 
# each site

genus_freq_by_dbh <- group_by(merged, site_ID, dbh_class_ID, genus_ID) %>%
    summarize(dbh=dbh_class_midpoints[dbh_class_ID[1]],
              n=n()) %>%
    group_by(site_ID, dbh_class_ID) %>%
    mutate(freq=n/sum(n)) %>%
    group_by(site_ID, genus_ID, dbh_class_ID) %>%
    arrange(site_ID, genus_ID, dbh_class_ID)

this_site_ID <- 1
this_dbh <- dbh_class_midpoints[1]

site_means <- foreach (this_site_ID=c(1:n_site), .combine=rbind) %:% 
        foreach (mcwd=c(0, 185, 185+(150*2), 185+(150*3)),
                 .packages=c("reshape2", "dplyr"),
                 .export=c("n_genus"),
                 .combine=rbind) %do% {
    # Need to standardize mcwd before using it as a predictor
    mcwd_std <- (mcwd - mcwd_mean)/mcwd_sd
    # Figure out the genera present at this site, and their frequencies
    this_genus_set <- filter(genus_freq_by_dbh, site_ID == this_site_ID)
    # Calculate predictions for all genera present at this site for all MCMC 
    # samples
    this_growth_preds <- pred_growth(unique(this_genus_set$genus_ID), dbhs=dbh_class_midpoints, mcwd=mcwd_std)
    this_growth_preds <- this_growth_preds * dbh_sd

    # Apply site-level intercept
    this_int_k <- array(rep(int_k[, this_site_ID]*dbh_sd, 
                            each=dim(this_growth_preds)[1] *
                                 dim(this_growth_preds)[2]), 
                        dim=dim(this_growth_preds))
    this_growth_preds <- this_growth_preds + this_int_k

    # Now weight each growth prediction by the frequency of each genus within 
    # each dbh class
    weight_matrix <- matrix(0, ncol=ncol(this_growth_preds), nrow=nrow(this_growth_preds))
    this_genus_set$genus_row_num <- as.numeric(factor(this_genus_set$genus_ID)) 
    for (i in 1:nrow(this_genus_set)) {
        row_num <- this_genus_set$genus_row_num[i]
        col_num <- this_genus_set$dbh_class_ID[i]
        weight_matrix[row_num, col_num] <- this_genus_set$freq[i]
    }
    stopifnot(all((apply(weight_matrix, 2, sum) == 0) |
                  (apply(weight_matrix, 2, sum) - 1 < 1e-10)))

    # Repeat weight_matrix in z direction to match this_growth_preds
    weight_matrix <- array(rep(weight_matrix, dim(this_growth_preds)[3]), dim=dim(this_growth_preds))
    # Weight growth predictions by frequencies of occurrence of individuals of 
    # each genera
    this_growth_preds <- this_growth_preds*weight_matrix
    this_growth_preds <- t(apply(this_growth_preds, c(2, 3), sum))

    growth_means <- apply(this_growth_preds, 2, mean)
    growth_q2pt5 <- apply(this_growth_preds, 2, quantile, .025)
    growth_q97pt5 <- apply(this_growth_preds, 2, quantile, .975)
    data.frame(cbind(site_ID=this_site_ID, mcwd=mcwd, 
                     dbh=dbh_class_midpoints*dbh_sd+dbh_mean, 
                     growth=growth_means, q2pt5=growth_q2pt5, 
                     q97pt5=growth_q97pt5))
}

site_ID_factor_key <- read.csv(file.path(data_folder,
                                         paste0("site_ID_factor_key_", suffix, 
                                                ".csv")))
site_means <- merge(site_means, site_ID_factor_key, by.x="site_ID", by.y="site_ID_numeric")

sitecode_key <- read.csv(file.path(prefix, "TEAM", "Sitecode_Key", "sitecode_key.csv")
site_means <- merge(site_means, sitecode_key, by.x="site_ID_char", by.y="sitecode")

ggplot(filter(site_means, mcwd == 0)) +
    geom_line(aes(x=dbh, y=growth, colour=site_ID_char)) +
    geom_point(aes(x=dbh, y=growth, colour=site_ID_char, shape=site_ID_char), size=2) +
    geom_ribbon(aes(x=dbh, ymin=q2pt5, ymax=q97pt5, fill=site_ID_char), alpha=.25)+
    coord_cartesian(xlim=c(15, 125), ylim=c(0, 1.5)) +
    scale_colour_discrete("Site") +
    scale_fill_discrete("Site") +
    scale_shape_manual("Site", values=rep(1:4, 4)) +
    xlab("Initial diameter (cm)") +
    ylab("Growth (cm/yr)")
ggsave("growth_site_level_mcwd_0.png", width=8, height=6, dpi=300)

ggplot(filter(site_means, mcwd == 0)) +
    geom_line(aes(x=dbh, y=growth, colour=site_ID_char)) +
    geom_point(aes(x=dbh, y=growth, colour=site_ID_char, shape=site_ID_char), size=2) +
    geom_ribbon(aes(x=dbh, ymin=q2pt5, ymax=q97pt5, fill=site_ID_char), alpha=.25)+
    coord_cartesian(xlim=c(15, 125), ylim=c(0, 1.5)) +
    facet_grid(continent~.) +
    scale_colour_discrete("Site") +
    scale_fill_discrete("Site") +
    scale_shape_manual("Site", values=rep(1:4, 4)) +
    xlab("Initial diameter (cm)") +
    ylab("Growth (cm/yr)")
ggsave("growth_site_level_mcwd_0_continent.png", width=8, height=6, dpi=300)

color_values <- c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3")
ggplot(filter(site_means, mcwd == 0)) +
    geom_line(aes(x=dbh, y=growth, colour=site_ID_char)) +
    geom_point(aes(x=dbh, y=growth, colour=site_ID_char, shape=site_ID_char), size=2) +
    geom_ribbon(aes(x=dbh, ymin=q2pt5, ymax=q97pt5, fill=site_ID_char), alpha=.25)+
    coord_cartesian(xlim=c(15, 125), ylim=c(0, 1.5)) +
    # scale_colour_discrete("Site") +
    # scale_fill_discrete("Site") +
    # scale_shape_discrete("Site") +
    scale_colour_manual("Site", values=rep(color_values, each=4)) +
    scale_fill_manual("Site", values=rep(color_values, each=4)) +
    scale_shape_manual("Site", values=rep(1:4, 4)) +
    xlab("Initial diameter (cm)") +
    ylab("Growth (cm/yr)")
ggsave("growth_site_level_mcwd_0_v1.png", width=8, height=6, dpi=300)

ggplot(filter(site_means)) +
    geom_line(aes(x=dbh, y=growth, colour=factor(mcwd))) +
    geom_ribbon(aes(x=dbh, ymin=q2pt5, ymax=q97pt5, fill=factor(mcwd)), alpha=.25)+
    coord_cartesian(xlim=c(15, 125), ylim=c(0, 1.5)) +
    facet_wrap(~site_ID_char) +
    scale_fill_discrete("MCWD") +
    scale_colour_discrete("MCWD") +
    xlab("Initial diameter (cm)") +
    ylab("Growth (cm/yr)")
ggsave("growth_site_level_mcwd_all.png", width=8, height=6, dpi=300)

ggplot(filter(site_means, !(site_ID_char %in% c("KRP", "YAS", "COU", "YAN")))) +
    geom_line(aes(x=dbh, y=growth, colour=factor(mcwd))) +
    geom_ribbon(aes(x=dbh, ymin=q2pt5, ymax=q97pt5, fill=factor(mcwd)), alpha=.25)+
    coord_cartesian(xlim=c(15, 125), ylim=c(0, 1.5)) +
    facet_wrap(~site_ID_char) +
    scale_fill_discrete("MCWD") +
    scale_colour_discrete("MCWD") +
    xlab("Initial diameter (cm)") +
    ylab("Growth (cm/yr)")
ggsave("growth_site_level_mcwd_all_without_short_records.png", width=8, height=6, dpi=300)

ggplot(filter(site_means)) +
    geom_line(aes(x=dbh, y=growth, colour=site_ID_char)) +
    geom_point(aes(x=dbh, y=growth, colour=site_ID_char, shape=site_ID_char), size=2) +
    geom_ribbon(aes(x=dbh, ymin=q2pt5, ymax=q97pt5, fill=site_ID_char), alpha=.25)+
    coord_cartesian(xlim=c(15, 125), ylim=c(0, 1.5)) +
    facet_grid(continent~mcwd) +
    scale_colour_discrete("Site") +
    scale_fill_discrete("Site") +
    scale_shape_manual("Site", values=rep(1:4, 4)) +
    xlab("Initial diameter (cm)") +
    ylab("Growth (cm/yr)")
ggsave("growth_site_level_mcwd_all_continent.png", width=8, height=6, dpi=300)

ggplot(filter(site_means, !(site_ID_char %in% c("KRP", "YAS", "COU", "YAN")))) +
    geom_line(aes(x=dbh, y=growth, colour=site_ID_char)) +
    geom_point(aes(x=dbh, y=growth, colour=site_ID_char, shape=site_ID_char), size=2) +
    geom_ribbon(aes(x=dbh, ymin=q2pt5, ymax=q97pt5, fill=site_ID_char), alpha=.25)+
    coord_cartesian(xlim=c(15, 125), ylim=c(0, 1.5)) +
    facet_grid(continent~mcwd) +
    scale_colour_discrete("Site") +
    scale_fill_discrete("Site") +
    scale_shape_manual("Site", values=rep(1:4, 4)) +
    xlab("Initial diameter (cm)") +
    ylab("Growth (cm/yr)")
ggsave("growth_site_level_mcwd_all_continent_sites_without_shortrecords.png", width=8, height=6, dpi=300)

###############################################################################
# Plot variation in MCWD slope versus wood density

B_g_scaling <- c(dbh_sd,
                 (dbh_sd/mcwd_sd) * 100, # Convert from mm to 100s of mm
                 (dbh_sd/mcwd_sd) * 100, # Convert from mm to 100s of mm
                 dbh_sd/dbh_sd,
                 dbh_sd/dbh_sd)

B_g_scaling <- matrix(rep(B_g_scaling, each=nrow(B_g_array)), ncol=ncol(B_g_array))
B_g_scaling <- array(rep(B_g_scaling, dim(B_g_array)[3]), dim=dim(B_g_array))
B_g_rescaled <- B_g_array * B_g_scaling

B_g_means <- data.frame(apply(B_g_rescaled, c(1, 2), mean))
B_g_means <- cbind(1:nrow(B_g_means), B_g_means)
names(B_g_means) <- c("genus_ID", "int", "MCWD", "MCWD^2", "DBH", "DBH^2")
B_g_means <- melt(B_g_means, id.vars="genus_ID", value.name="mean")

B_g_q2pt5 <- data.frame(apply(B_g_rescaled, c(1, 2), quantile, .025))
B_g_q2pt5 <- cbind(1:nrow(B_g_q2pt5), B_g_q2pt5)
names(B_g_q2pt5) <- c("genus_ID", "int", "MCWD", "MCWD^2", "DBH", "DBH^2")
B_g_q2pt5 <- melt(B_g_q2pt5, id.vars="genus_ID", value.name="q2pt5")

B_g_q97pt5 <- data.frame(apply(B_g_rescaled, c(1, 2), quantile, .975))
B_g_q97pt5 <- cbind(1:nrow(B_g_q97pt5), B_g_q97pt5)
names(B_g_q97pt5) <- c("genus_ID", "int", "MCWD", "MCWD^2", "DBH", "DBH^2")
B_g_q97pt5 <- melt(B_g_q97pt5, id.vars="genus_ID", value.name="q97pt5")

B_g_means <- merge(B_g_means, B_g_q2pt5)
B_g_means <- merge(B_g_means, B_g_q97pt5)
B_g_means <- merge(B_g_means, WD)
B_g_means$WD <- B_g_means$WD*WD_sd + WD_mean

cor(B_g_means$mean, B_g_means$WD)

B_g_means$WD_class <- cut(B_g_means$WD, c(0, .4, .7, 1.5), include.lowest=TRUE)
min(B_g_means$WD)
max(B_g_means$WD)
table(B_g_means$WD_class)
table(is.na(B_g_means$WD_class))

ggplot(filter(B_g_means, variable %in% c("MCWD", "MCWD^2")),
       aes(x=WD, y=mean)) +
    geom_point() + facet_wrap(~variable, scales="free") +
    geom_smooth()
ggsave("growth_MCWD_slope_vs_WD_points.png", width=6, height=6, dpi=300)

ggplot(filter(B_g_means, variable %in% c("MCWD", "MCWD^2"))) +
    geom_density(aes(mean, fill=WD_class), alpha=.3) + facet_grid(~variable)
ggsave("growth_MCWD_slope_vs_WD_densities.png", width=6, height=6, dpi=300)

ggplot(filter(B_g_means, variable %in% c("int")),
       aes(x=WD, y=mean)) +
    geom_point() + facet_wrap(~variable, scales="free")

    geom_density(aes(mean, fill=WD), alpha=.3) + facet_grid(~variable)

ggsave("growth_WD_slope_vs_WD_densities.png", width=6, height=6, dpi=300)


# Plot five dominant genera in each site
# Plot by region

###############################################################################
# Plot mean global growth curve
load(file.path(params_folder, paste0(suffix, "jags_fit_full_model_ranefs_mu_B_g.RData")))
ranefs_mu_B_g <- as.data.frame(combine.mcmc(ranefs_mu_B_g))

# Build coefficient matrix
B_all <- cbind(B, ranefs_mu_B_g)

#overall_growth_preds <- foreach(mcwd=(c(0, 185, 185+(150*2), 185+(150*3))- 
#overall_growth_preds <- foreach(mcwd=(c(185, 185+(150*2))- 
overall_growth_preds <- foreach(mcwd=(c(0, 185, 185+(150*2))- 
                                      mcwd_mean)/mcwd_sd,
                                .combine=rbind) %do% {
    # Build design matrix
    xs <- cbind(rep(0, length(dbh_smalld_class_midpoints)), # WD
                rep(0, length(dbh_smalld_class_midpoints)), # WD^2
                rep(1, length(dbh_smalld_class_midpoints)), # genus-level intercept
                rep(mcwd, length(dbh_smalld_class_midpoints)), # genus-level MCWD
                rep(mcwd^2, length(dbh_smalld_class_midpoints)), # genus-level MCWD^2
                dbh_smalld_class_midpoints,          # genus-level DBH
                dbh_smalld_class_midpoints^2)        # genus-level DBH^2
    xs <- matrix(xs, nrow=nrow(xs))
    preds <- t(xs %*% t(B_all) - rep(dbh_smalld_class_midpoints, 
                                     nrow(B_all)))*dbh_sd
    means <- apply(preds, 2, mean)
    q2pt5 <- apply(preds, 2, quantile, .025)
    q97pt5 <- apply(preds, 2, quantile, .975)
    data.frame(cbind(mcwd=mcwd*mcwd_sd + mcwd_mean, 
                     dbh=dbh_smalld_class_midpoints*dbh_sd+dbh_mean, 
                     growth=means, q2pt5=q2pt5, q97pt5=q97pt5))
}

ggplot(overall_growth_preds, aes(dbh, growth)) +
    theme_bw(base_size=10) +
    geom_line(aes(colour=factor(mcwd), linetype=factor(mcwd)), size=.4) +
    geom_ribbon(aes(ymin=q2pt5, ymax=q97pt5, fill=factor(mcwd)), alpha=.25)+
    coord_cartesian(xlim=c(15, 125), ylim=c(0, .75)) +
    scale_fill_discrete("MCWD", guide=guide_legend(title.position="top")) +
    scale_colour_discrete("MCWD", guide=guide_legend(title.position="top")) +
    scale_shape_discrete("MCWD", guide=guide_legend(title.position="top")) +
    scale_linetype_discrete("MCWD", guide=guide_legend(title.position="top")) +
    xlab("Initial diameter (cm)") +
    ylab("Growth (cm/yr)") +
    theme(legend.key.width=unit(.5, "cm"),
          #legend.position=c(.29,.84),
          legend.position=c(.7,.8),
          legend.direction="horizontal",
          legend.background=element_rect(fill=alpha('white', .5)))

ggsave("growth_overall_mcwd_all.png", width=3, height=2, dpi=300)

###############################################################################
# Plot site-level mean growth rates (as deviations from the overall mean)
int_k_mean <- apply(int_k, 2, mean) * dbh_sd
int_k_q2pt5 <- apply(int_k, 2, quantile, .025) * dbh_sd
int_k_q97pt5 <- apply(int_k, 2, quantile, .975) * dbh_sd
int_k_mean <- data.frame(site_ID_factor_key, mean=int_k_mean)
int_k_mean <- cbind(int_k_mean, q2pt5=int_k_q2pt5)
int_k_mean <- cbind(int_k_mean, q97pt5=int_k_q97pt5)
int_k_mean <- merge(int_k_mean, sitecode_key, by.x="site_ID_char", by.y="sitecode")
int_k_mean$positive <- int_k_mean$mean > 0

# Reorder factors
#int_k_mean$sitename_abbrev <- with(int_k_mean, factor(sitename_abbrev, levels=sitename_abbrev[order(continent, mean)]))
int_k_mean$sitename_abbrev <- with(int_k_mean, factor(sitename_abbrev, levels=sitename_abbrev[order(mean)]))
int_k_mean$site_ID_char <- with(int_k_mean, factor(site_ID_char, levels=site_ID_char[order(mean)]))

#int_k_mean$sitename_pretty_detailed <- gsub('\\(', '\\\n\\(', int_k_mean$sitename_pretty_detailed )
ggplot(int_k_mean, aes(site_ID_char, mean)) +
    theme_bw(base_size=8) +
    geom_point(stat="identity", aes(colour=continent_long)) +
    geom_errorbar(aes(ymin=q2pt5, ymax=q97pt5, colour=continent_long), width=.3) +
    scale_colour_discrete("Continent", guide=guide_legend(title.position="top")) +
    xlab("Site") +
    ylab("Mean deviation (cm/yr)") +
    theme(legend.key.width=unit(.5, "cm"),
          #legend.position=c(.29,.84),
          legend.position=c(.3,.8),
          legend.direction="horizontal",
          legend.background=element_rect(fill=alpha('white', .5)))
ggsave("growth_site_level_mean_mcwd_all.png", width=4.5, height=3, dpi=300)

int_k_long <- melt(int_k, variable.name="site_ID")
int_k_long$site_ID <- gsub('int_k\\[', '', int_k_long$site_ID)
int_k_long$site_ID <- gsub('\\]', '', int_k_long$site_ID)
int_k_long$site_ID <- as.numeric(int_k_long$site_ID)
int_k_long <- merge(int_k_long, site_ID_factor_key, by.x="site_ID", by.y="site_ID_numeric")
int_k_long <- merge(int_k_long, sitecode_key, by.x="site_ID_char", by.y="sitecode")
ggplot(int_k_long, aes(site_ID_char, value)) +
    theme_bw(base_size=8) +
    geom_violin(aes(fill=continent_long, colour=continent_long), width=.3) +
    scale_colour_discrete("Continent", guide=guide_legend(title.position="top")) +
    scale_fill_discrete("Continent", guide=guide_legend(title.position="top")) +
    xlab("Site") +
    ylab("Mean deviation (cm/yr)") +
    theme(legend.key.width=unit(.5, "cm"),
          #legend.position=c(.29,.84),
          legend.position=c(.3,.8),
          legend.direction="horizontal",
          legend.background=element_rect(fill=alpha('white', .5)))
ggsave("growth_site_level_mean_mcwd_all_violin.png", width=4.5, height=3, dpi=300)

###############################################################################
# Plot percentage of genera with significant drought effects by site

# Calculate relative frequency of each genus at each site
genus_freq_by_site <- group_by(merged, site_ID, genus_ID) %>%
    summarize(n=n()) %>%
    group_by(site_ID) %>%
    mutate(freq=n/sum(n)) %>%
    group_by(site_ID, genus_ID) %>%
    arrange(site_ID, genus_ID)

# Calculate mean MCWD slope by genus
mean_mcwd <- apply(B_g_array[, 2, ], 1, mean)
mean_mcwd_q2pt5 <- apply(B_g_array[, 2, ], 1, quantile, .025)
mean_mcwd_q97pt5 <- apply(B_g_array[, 2, ], 1, quantile, .975)
mean_mcwd <- data.frame(genus_ID=1:nrow(B_g_array), mean_mcwd)
mean_mcwd <- cbind(mean_mcwd, q2pt5=mean_mcwd_q2pt5)
mean_mcwd <- cbind(mean_mcwd, q97pt5=mean_mcwd_q97pt5)
mean_mcwd$signif <- mean_mcwd$q97pt5 < 0
table(mean_mcwd$signif)
