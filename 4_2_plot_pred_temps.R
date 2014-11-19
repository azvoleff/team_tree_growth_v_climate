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

source("0_settings.R")

cl <- makeCluster(3)
registerDoParallel(cl)

model_type <- "full"
#model_type <- "testing"
precip_var <- "mcwd_run12"

temp_var <- "tmx_meanannual"
run_ID <- "vertica1.team.sdsc.edu_20141110224426_extend3" 

# temp_var <- "tmp_meanannual"
# run_ID <- "vertica1.team.sdsc.edu_20141110152032_extend1"

suffix <- paste0(model_type, '-', temp_var, '-', precip_var)

data_folder <- file.path(prefix, "TEAM", "Tree_Growth", "Data")
params_folder <- file.path(prefix, "TEAM", "Tree_Growth", "Extracted_Parameters")

n_site <- 14

###############################################################################
# Load MCMC samples

load(file.path(params_folder, paste0(suffix, "jags_fit_full_model_ranefs_B_T.RData")))

start_val <- 14000
thin_val <- 4
ranefs_B_T <- window(ranefs_B_T, start=start_val, thin=thin_val)

B_T <- combine.mcmc(ranefs_B_T)
n_mcmc <- nrow(B_T)
n_B_T <- ncol(B_T)/n_site

# Make a series of vectors to check that the array conversion is done 
# correctly.
# 
# First vector is the set of coefficients from the 1st mcmc sample for site 1
B_T_set1 <- B_T[1, c(0:(n_B_T - 1)) * n_site+ 1]
# Second vector is the set of coefficients from the 5th mcmc sample for site 7
B_T_set2 <- B_T[5, c(0:(n_B_T - 1)) * n_site + 7]
# Now convert to a 3D array, with genera in rows, coefficients in columns, and 
# separate MCMC samples in third dimension
dim(B_T)
B_T_array <- t(matrix(t(B_T), byrow=TRUE, ncol=n_site))
dim(B_T_array)
B_T_array <- array(B_T_array, dim=c(n_site, n_B_T, n_mcmc))
dim(B_T_array)
stopifnot(all(B_T_array[1, , 1] == B_T_set1))
stopifnot(all(B_T_array[7, , 5] == B_T_set2))

###############################################################################
# Now make plots of predicted temperatures at each plot for each site

site_ID_factor_key <- read.csv(file.path(data_folder,
                                         paste0("site_ID_factor_key_", suffix, 
                                                ".csv")))

# Below function predictions the temp by site for each of a set of MCMC 
# samples
pred_temp <- function(site_IDs=1:n_site, elevs=list(c(0, 100, 200, 300)), temp=0) {
    n_site <- length(site_IDs)
    # Make predictions into a n_site x n_elevs x n_mcmc array.
    temp_preds <- array(NA, dim=c(n_site, length(elevs), n_mcmc))
    for (n in 1:n_site) {
        site_ID <- site_IDs[n]
        elev <- elevs[n]
        # Collapse betas for this site into flat matrices (temporarily eliminating 
        # third dimension for efficient matrix multiplication)
        B_flat <- matrix(B_rep[site_ID, , ], nrow=dim(B_rep)[3])
        B_g_flat <- matrix(B_g_array[site_ID, , ], nrow=dim(B_g_array)[3], byrow=TRUE)
        B_all <- cbind(B_flat, B_g_flat)
        this_x_all <- matrix(rep(x_all[site_ID, ], length(elev)), ncol=n_B + n_B_g, byrow=TRUE)
        this_x_all[, 4] <- temp 
        this_x_all[, 6] <- elev 
        temp_preds_flat <- this_x_all %*% t(B_all) - rep(elev, nrow(B_all))
        temp_preds[n, , ] <- temp_preds_flat
    }
    return(temp_preds)
}

temp_means <- apply(temp_preds, c(1, 2), mean)
temp_q2pt5 <- apply(temp_preds, c(1, 2), quantile, .025)
temp_q97pt5 <- apply(temp_preds, c(1, 2), quantile, .975)

# Plot subset of temp curves
ggplot(temp, aes(elev, mean)) +
    geom_line(aes(colour=site_ID_char))+
    geom_ribbon(aes(ymin=q2pt5, ymax=q97pt5, fill=site_ID_char), alpha=.25)+
    coord_cartesian(xlim=c(10, 150), ylim=c(0, 3)) +
    scale_fill_discrete("Site") +
    scale_colour_discrete("Site") +
    xlab("Elevation (m)") +
    ylab("Temperature (degrees C)")
ggsave("temp_by_site.png", width=4, height=2, dpi=300)
