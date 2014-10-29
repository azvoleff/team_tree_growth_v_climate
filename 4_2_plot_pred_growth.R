library(ggplot2)
library(runjags)
library(coda)
library(dplyr)

load("jags_fit_full_model_fixefs.RData")
load("jags_fit_full_model_ranefs_mu_B_g.RData")

load("model_data_wide.RData")
load("model_data_standardizing.RData")


mcwd <- 0

merged <- with(model_data, data.frame(site_ID=site_ID, 
                                      genus_ID=as.integer(genus_ID),
                                      WD=WD,
                                      dbh=apply(dbh, 1, mean, na.rm=TRUE)))
d_dbh <- 5
dbhs <- seq(10, 150, d_dbh)
dbhs_std <- (dbhs - dbh_mean)/dbh_sd
dbh_class_midpoints <- diff(dbhs_std) + dbhs_std[-length(dbhs_std)]
merged$dbh_class_ID <- as.numeric(cut(merged$dbh, dbhs_std, include.lowest=TRUE))

genus_freq <- group_by(merged, site_ID, dbh_class_ID, genus_ID) %>%
    summarize(dbh=dbh_class_midpoints[dbh_class_ID[1]],
              WD=WD[1],
              mcwd=mcwd,
              n=n()) %>%
    group_by(site_ID, dbh_class_ID) %>%
    mutate(freq=n/n())

###############################################################################
# Load MCMC samples

# Fixed effects on density
load("jags_fit_full_model_fixefs.RData")
# Combine the MCMC chains and convert to data.frame with coefficients in 
# columns and MCMC samples in rows
fixefs <- as.data.frame(combine.mcmc(fixefs))
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
ranefs_B_g <- combine.mcmc(ranefs_B_g)
# Make a series of vectors to check that the array conversion is done 
# correctly.
# 
# First vector is the set of coefficients from the 1st mcmc sample for genus 1
ranefs_B_g_set1 <- ranefs_B_g[1, 1:5]
# Second vector is the set of coefficients from the 5th mcmc sample for genus 3
ranefs_B_g_set2 <- ranefs_B_g[5, 11:15]
# Now convert to a 3D array, with genera in rows, coefficients in columns, and 
# separate MCMC samples in third dimension
ranefs_B_g <- array(ranefs_B_g, dim=c(nrow(ranefs_B_g), 5, ncol(ranefs_B_g)/5))
dim(ranefs_B_g)
ranefs_B_g <- aperm(ranefs_B_g, c(3, 2, 1))
dim(ranefs_B_g)
stopifnot(all(ranefs_B_g[1, 1:5, 1] == ranefs_B_g_set1))
stopifnot(all(ranefs_B_g[3, 1:5, 5] == ranefs_B_g_set2))

###############################################################################
# Make predictions of mean genus-level growth
WD <- group_by(merged, genus_ID) %>%
    summarize(WD=WD[1])
# Construct fixed effects design matrix
x <- cbind(WD$WD, WD$WD^2)
# Construct genus-level random effects design matrix
x_g <- cbind(1, # genus-level intercept
             rep(0, nrow(WD)),
             rep(0, nrow(WD)),
             rep(dbh_class_midpoints[1], nrow(WD)),
             rep(dbh_class_midpoints[1]^2, nrow(WD)))
# Combine the two matrices:
x_all <- cbind(x, x_g)

###############################################################################
# Now make growth curves for each genus
# B <- array(B, dim=c(1, ncol(B), nrow(B)))

# dim(aperm(array(B), c(1)))
# B <- matrix(rep(t(B), nrow(WD)), ncol=ncol(B), byrow=TRUE)

# Prediction for one MCMC sample for genus 1
genus_num <- 400
mcmc_num <- 1
growths <- c()
for (n in 1:length(dbh_class_midpoints)) {
    xs <- matrix(c(x[genus_num, ], x_g[genus_num, ]), nrow=1)
    betas <- matrix(c(B[mcmc_num, ], ranefs_B_g[genus_num, , mcmc_num]))
    xs[6] <- dbh_class_midpoints[n]
    xs[7] <- dbh_class_midpoints[n]^2
    growths <- c(growths, xs %*% betas)
}
plot(dbh_class_midpoints*dbh_sd, growths*dbh_sd)
