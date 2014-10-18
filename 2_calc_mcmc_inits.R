library(dplyr)
library(reshape2)
library(ggplot2)
library(lme4)

load("init_data.RData")
load("model_data_wide.RData")
load("model_data_long.RData")

#sitecode_key <- read.csv('H:/Data/TEAM/Sitecode_Key/sitecode_key.csv')
sitecode_key <- read.csv('C:/Users/azvoleff/Desktop/Sitecode_Key/sitecode_key.csv')
sitecode_key <- select(sitecode_key, sitecode, sitetype, continent)
model_data_long <- merge(model_data_long, sitecode_key, by.x="site_ID", by.y="sitecode")

img_height <- 4
img_width <- 3
img_dpi <- 300

mcwd_long <- melt(model_data$mcwd, varnames=c("tree_ID", "period_ID"), 
                  value.name="mcwd")
dbh_latent_long <- melt(init_data$dbh_latent, varnames=c("tree_ID", "period_ID"), 
                        value.name="dbh_latent_end")
dbh_latent_long <- group_by(dbh_latent_long, tree_ID) %>%
    arrange(period_ID) %>%
    mutate(dbh_latent_start=c(NA, dbh_latent_end[1:length(dbh_latent_end) - 1]))
calib_data <- merge(dbh_latent_long, mcwd_long)
calib_data <- calib_data[complete.cases(calib_data), ]

WD <- data.frame(tree_ID=seq(1, model_data$n_tree), WD=model_data$WD)
calib_data <- merge(calib_data, WD)
site_ID <- data.frame(tree_ID=seq(1, model_data$n_tree), site_ID=model_data$site_ID)
calib_data <- merge(calib_data, site_ID)
plot_ID <- data.frame(tree_ID=seq(1, model_data$n_tree), plot_ID=model_data$plot_ID)
calib_data <- merge(calib_data, plot_ID)
genus_ID <- data.frame(tree_ID=seq(1, model_data$n_tree), genus_ID=model_data$genus_ID)
calib_data <- merge(calib_data, genus_ID)

hist(calib_data$dbh_latent_start)
hist(calib_data$dbh_latent_end)
hist(calib_data$WD)
hist(calib_data$mcwd)

# obs_per_tree <- group_by(calib_data, site_ID, tree_ID) %>%
#     summarize(n=n()) %>%
#     group_by(site_ID) %>%
#     summarize(mean(n))
#
# test_model  <- lm(dbh_latent_end ~ dbh_latent_start + I(dbh_latent_start^2) +
#                    WD + I(WD^2)+ 
#                    mcwd + I(mcwd^2), data=calib_data)
# calib_model  <- lmer(dbh_latent_end ~ dbh_latent_start + I(dbh_latent_start^2) +
#                      WD + I(WD^2)+ 
#                      mcwd + I(mcwd^2) +
#                      (mcwd + I(mcwd^2) + dbh_latent_start + I(dbh_latent_start^2)|genus_ID) +
#                      (1|site_ID) + (1|plot_ID) + (1|tree_ID) + (1|period_ID), data=calib_data)
#                      control=lmerControl(optCtrl=list(maxfun=20000)))
# save(calib_model, file="calib_model.RData")
load("calib_model.RData")

get_variance <- function(model, grp, var1, var2) {
    vc <- as.data.frame(VarCorr(model))
    if (missing(var1) & missing(var2)) {
        if (grp != "Residual") {
            stop('must specify a var1 if grp is not "Residual"')
        } else {
            return(vc[vc$grp == "Residual", ]$vcov)
        }
    } else if (missing(var1) & !missing(var2)) {
        stop('must specify a var1 if var2 is specified')
    } else {
        vc <- vc[vc$grp == grp & vc$var1== var1, ]
        if (missing(var2)) {
            return(vc[is.na(vc$var2), ]$vcov)
        } else {
            return(vc$vcov[match(var2, vc$var2)])
        }
    }
}

init_data$int_ijk <- as.numeric(unlist(ranef(calib_model)$tree_ID))
init_data$int_jk <- as.numeric(unlist(ranef(calib_model)$plot_ID))
init_data$int_k <- as.numeric(unlist(ranef(calib_model)$site_ID))
init_data$int_t <- as.numeric(unlist(ranef(calib_model)$period))
init_data$sigma_int_ijk <- sqrt(get_variance(calib_model, "tree_ID", "(Intercept)"))
init_data$sigma_int_jk <- sqrt(get_variance(calib_model, "plot_ID", "(Intercept)"))
init_data$sigma_int_k <- sqrt(get_variance(calib_model, "site_ID", "(Intercept)"))
init_data$sigma_int_t <- sqrt(get_variance(calib_model, "period_ID", "(Intercept)"))
# Extract genus-level random effects
init_data$B_g_raw <- as.matrix(ranef(calib_model)$genus_ID)
# Extract variance-covariance matrix for genus-level random effects
genus_varcorr <- VarCorr(calib_model)$genus_ID
# Drop the attributes
genus_varcorr <- matrix(c(genus_varcorr), nrow=nrow(genus_varcorr))
# JAGS models the inverse variance-covariance matrix (Tau)
init_data$Tau_B_raw <- genus_varcorr^-1

save(init_data, file="init_data_with_ranefs.RData")
