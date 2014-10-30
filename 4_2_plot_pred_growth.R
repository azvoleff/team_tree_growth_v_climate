library(ggplot2)
library(runjags)
library(coda)
library(dplyr)
library(reshape2)
library(foreach)
library(doParallel)

cl <- makeCluster(4)
registerDoParallel(cl)

load("jags_fit_full_model_fixefs.RData")
load("jags_fit_full_model_ranefs_mu_B_g.RData")

load("model_data_wide.RData")
load("model_data_standardizing.RData")

n_site <- model_data$n_site
n_genus <- model_data$n_genus
n_B <- 2
n_B_g <- 5

merged <- with(model_data, data.frame(site_ID=site_ID, 
                                      genus_ID=as.integer(genus_ID),
                                      WD=WD,
                                      dbh=apply(dbh, 1, mean, na.rm=TRUE)))
d_dbh <- 5
dbhs <- seq(10, 150, d_dbh)
dbhs_std <- (dbhs - dbh_mean)/dbh_sd
dbh_class_midpoints <- diff(dbhs_std)/2 + dbhs_std[-length(dbhs_std)]
merged$dbh_class_ID <- as.numeric(cut(merged$dbh, dbhs_std, include.lowest=TRUE))

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
load("jags_fit_full_model_fixefs.RData")
# Combine the MCMC chains and convert to data.frame with coefficients in 
# columns and MCMC samples in rows
fixefs <- as.data.frame(combine.mcmc(fixefs))
n_mcmc <- nrow(fixefs)
B <- fixefs[grepl('^B\\[', names(fixefs))]
B <- B[order(names(B))]
B <- as.matrix(B)

# Random intercepts at plot, site, and period levels
load("jags_fit_full_model_ranefs.RData")
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
load("jags_fit_full_model_ranefs_B_g.RData")
B_g <- combine.mcmc(ranefs_B_g)
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
B_g <- t(matrix(t(B_g), byrow=TRUE, ncol=n_genus))
dim(B_g)
B_g <- array(B_g, dim=c(n_genus, n_B_g, n_mcmc))
dim(B_g)
stopifnot(all(B_g[1, , 1] == B_g_set1))
stopifnot(all(B_g[7, , 5] == B_g_set2))

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
B <- matrix(rep(B, each=n_genus), ncol=2)
B <- array(B, dim=c(n_genus, n_B, n_mcmc))

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
pred_growth <- function(genera_IDs=1:n_genus, dbhs=dbh_class_midpoints, mcwd=0) {
    n_genus <- length(genera_IDs)
    # Make predictions into a genus x dbh x n_mcmc array.
    growth_preds <- array(NA, dim=c(n_genus, length(dbhs), n_mcmc))
    for (n in 1:n_genus) {
        genus_ID <- genera_IDs[n]
        # Collapse betas for this genus into flat matrices (temporarily eliminating 
        # third dimension for efficient matrix multiplication)
        B_flat <- matrix(B[genus_ID, , ], nrow=dim(B)[3])
        B_g_flat <- matrix(B_g[genus_ID, , ], nrow=dim(B_g)[3], byrow=TRUE)
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
genus_mean_growth <- function(genera_IDs=1:n_genus, mcwds=0) {
    n_genus <- length(genera_IDs)
    growths <- foreach (mcwd=mcwds, .combine=rbind,
                        .export=c("pred_growth", "n_genus", "n_mcmc", "B", 
                                  "B_g", "x_all", "n_B", "n_B_g", 
                                  "dbh_class_midpoints", "dbh_sd"),
                        .packages=c("reshape2", "dplyr")) %dopar% {
        growth_preds <- pred_growth(mcwd=mcwd)
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

mcwds <- c(0, 100, 200, 500)
mcwds <- (mcwds - mcwd_mean)/mcwd_sd
growth <- genus_mean_growth(mcwds=mcwds)

genus_ID_factor_key <- read.csv("genus_ID_factor_key.csv")
growth <- merge(growth, genus_ID_factor_key, by.x="genus_ID", by.y="genus_ID_numeric")
growth$dbh <- dbh_class_midpoints[growth$dbh_class_ID]*dbh_sd + dbh_mean
growth$mcwd <- (growth$mcwd*mcwd_sd) + mcwd_mean
growth <- tbl_df(growth)
growth <- arrange(growth, genus_ID, dbh_class_ID, mcwd)

# Plot subset of growth curves
ggplot(filter(growth, genus_ID %in% c(833, 566, 784), mcwd == 0)) +
    geom_line(aes(x=dbh, y=mean, colour=genus_ID_char))+
    geom_ribbon(aes(x=dbh, ymin=q2pt5, ymax=q97pt5, fill=genus_ID_char), alpha=.25)+
    coord_cartesian(xlim=c(10, 150), ylim=c(0, 3)) +
    scale_fill_discrete("Genus") +
    scale_colour_discrete("Genus") +
    xlab("Initial diameter (cm)") +
    ylab("Growth (cm/yr)") +
    theme(legend.text=element_text(face="italic"))

# Plot ALL growth curves
ggplot(filter(growth, mcwd == 0)) +
    geom_line(aes(x=dbh, y=mean, group=genus_ID), alpha=.2)+
    #geom_ribbon(aes(x=dbh, ymin=q2pt5, ymax=q97pt5, group=genus_ID), 
    #alpha=.05)+
    coord_cartesian(xlim=c(10, 150), ylim=c(0, 3)) +
    xlab("Initial diameter (cm)") +
    ylab("Growth (cm/yr)")

# Plot specific curve versus MCWD
ggplot(filter(growth, genus_ID %in% c(784))) +
    geom_line(aes(x=dbh, y=mean, colour=factor(mcwd)))+
    geom_ribbon(aes(x=dbh, ymin=q2pt5, ymax=q97pt5, fill=factor(mcwd)), alpha=.25)+
    coord_cartesian(xlim=c(10, 150), ylim=c(0, 3)) +
    scale_fill_discrete("Genus") +
    scale_colour_discrete("Genus") +
    xlab("Initial diameter (cm)") +
    ylab("Growth (cm/yr)") +
    theme(legend.text=element_text(face="italic"))

###############################################################################
# Predict site-level mean growth based on observed frequencies of each genus at 
# each site

genus_freq <- group_by(merged, site_ID, dbh_class_ID, genus_ID) %>%
    summarize(dbh=dbh_class_midpoints[dbh_class_ID[1]],
              n=n()) %>%
    group_by(site_ID, dbh_class_ID) %>%
    mutate(freq=n/sum(n)) %>%
    group_by(site_ID, genus_ID, dbh_class_ID) %>%
    arrange(site_ID, genus_ID, dbh_class_ID)

this_site_ID <- 1
this_dbh <- dbh_class_midpoints[1]

site_means <- foreach (site_ID=c(1:n_site), .combine=rbind) %:% 
        foreach (mcwd=c(0, 100, 200, 400),
                 .packages=c("reshape2", "dplyr"),
                 .export=c("n_genus","dbh_sd"),
                 .combine=rbind) %dopar% {

    # Figure out the genera present at this site, and their frequencies
    this_genus_set <- filter(genus_freq, site_ID == this_site_ID)
    # Plot mean curves for the site
    growth_preds <- pred_growth(unique(this_genus_set$genus_ID), dbhs=this_dbh, mcwd=mcwd)
    # Calculate predictions for all genera present at this site for all MCMC 
    # samples
    this_growth_preds <- pred_growth(unique(this_genus_set$genus_ID), dbhs=dbh_class_midpoints, mcwd=mcwd)
    this_growth_preds <- this_growth_preds * dbh_sd
    # Now weight each growth prediction by the frequency of each genus within 
    # each dbh class
    weight_matrix <- matrix(0, ncol=ncol(this_growth_preds), nrow=nrow(this_growth_preds))
    this_genus_set$genus_row_num <- as.numeric(factor(this_genus_set$genus_ID)) 
    for (i in 1:nrow(this_genus_set)) {
        row_num <- this_genus_set$genus_row_num[i]
        col_num <- this_genus_set$dbh_class_ID[i]
        weight_matrix[row_num, col_num] <- this_genus_set$freq[i]
    }
    stopifnot(all(apply(weight_matrix, 2, sum) %in% c(0, 1)))

    # Repeat weight_matrix in z direction to match this_growth_preds
    weight_matrix <- array(rep(weight_matrix, dim(this_growth_preds)[3]), dim=dim(this_growth_preds))
    # Weight growth predictions by frequencies of occurrence of individuals of 
    # each genera
    this_growth_preds <- this_growth_preds*weight_matrix
    this_growth_preds <- t(apply(this_growth_preds, c(2, 3), sum))

    growth_means <- apply(this_growth_preds, 2, mean)
    growth_q2pt5 <- apply(this_growth_preds, 2, quantile, .025)
    growth_q97pt5 <- apply(this_growth_preds, 2, quantile, .975)
    data.frame(cbind(site=site_ID, mcwd=mcwd, dbh=dbh_class_midpoints, 
                     growth=growth_means, q2pt5=growth_q2pt5, 
                     q97pt5=growth_q97pt5))

}

ggplot(site_means) +
    geom_line(aes(x=dbh*dbh_sd+dbh_mean, y=growth, colour=factor(mcwd))) +
    geom_ribbon(aes(x=dbh*dbh_sd+dbh_mean, ymin=q2pt5, ymax=q97pt5, fill=factor(mcwd)), alpha=.25)+
    coord_cartesian(xlim=c(10, 150), ylim=c(0, 3)) +
    facet_wrap(~site) +
    scale_fill_discrete("MCWD") +
    scale_colour_discrete("MCWD") +
    xlab("Initial diameter (cm)") +
    ylab("Growth (cm/yr)")
