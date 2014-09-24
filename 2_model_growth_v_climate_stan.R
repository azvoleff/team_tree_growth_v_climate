library(dplyr)
library(ggplot2)
library(rstan)
library(foreach)
library(doParallel)

load("growth_ctfsflagged_merged_detrended.RData")

img_height <- 4
img_width <- 3
img_dpi <- 300

growth <- tbl_df(growth)

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

n_burnin <- 200
n_chains <- 8
n_iter <- 1000

###############################################################################
#  Stan model

model_file <- "2_model_growth_v_climate_stan.stan" 

# Setup timeseries formatted dataframe (with NAs for covariates at time 1 since 
# first growth measurement can only be calculated from time 1 to time 2)
growth_wide <- data.frame(ID_tree=as.integer(factor(growth$SamplingUnitName)),
                          ID_plot=as.integer(factor(growth$plot_ID)),
                          ID_site=as.integer(factor(growth$sitecode)),
                          ID_period=as.integer(factor(growth$SamplingPeriodNumber)),
                          dbh_st=growth$diameter_start,
                          dbh_end=growth$diameter_end)
growth_ts <- arrange(growth_wide, ID_period) %>%
    group_by(ID_tree) %>%
    do(data.frame(ID_tree=c(.$ID_tree[1], .$ID_tree),
                  ID_plot=c(0, .$ID_plot),
                  ID_site=c(0, .$ID_site),
                  ID_period=c(0, .$ID_period),
                  dbh=c(.$dbh_st[1], .$dbh_end))) %>%
    arrange(ID_tree, ID_period)

# Add latent growth inits
calc_latent_growth <- function(dbh) {
    dbh_latent <- rep(NA, length(dbh))
    dbh_latent[1] <- dbh[1]
    for (i in 2:length(dbh)) {
        dbh_latent[i] <- ifelse(dbh[i] > dbh_latent[i - 1], dbh[i], dbh_latent[i - 1] + .01)
    }
    return(dbh_latent)
}
growth_ts <- arrange(growth_ts, ID_period) %>%
    group_by(ID_tree) %>%
    mutate(dbh_latent=calc_latent_growth(dbh)) %>%
    arrange(ID_tree, ID_period)

# Setup data
obs_per_tree <- group_by(growth_ts, ID_tree) %>% summarize(n=n())
n_tree <- length(unique(growth_ts$ID_tree))
n_site <- length(unique(growth_ts$ID_site))
n_plot <- length(unique(growth_ts$ID_plot))

stan_data <- list(n_dbhs=nrow(growth_ts),
                  n_tree=n_tree,
                  n_plot=n_plot,
                  n_site=n_site,
                  n_obs_per_tree=obs_per_tree$n,
                  n_growths=(sum(obs_per_tree$n) - n_tree),
                  tree_ID=as.integer(factor(growth_ts$ID_tree)),
                  plot_ID=as.integer(factor(growth_ts$ID_plot)),
                  site_ID=as.integer(factor(growth_ts$ID_site)),
                  log_dbh=log(growth_ts$dbh))

#var(log(growth_ts$dbh))
# Setup inits
stan_init <- list(list(log_dbh_latent=log(growth_ts$dbh_latent),
                       beta_0=-1.7,
                       beta_1=.008,
                       sigma_obs=.02, 
                       sigma_proc=.02, 
                       sigma_ijk=.06,
                       sigma_jk=.12,
                       sigma_k=.24,
                       b_ijk=rep(0, n_tree),
                       b_jk=rep(0, n_plot),
                       b_k=rep(0, n_site)))

## Finish adjustment for change of variables

# fit <- stan(model_file, data=stan_data, iter=n_iter, chains=n_chains,
#             init=rep(stan_init, n_chains), warmup=n_burnin)

# # Compile and run piece-by-piece
# ret <- stanc(model_file, model_name="my_model")
# ## compilation
# my_model.sm <- stan_model(stanc_ret=ret)
#
# for (s in 1:n.sims) {
#    my_model.sf <- sampling(my_model.sm, data=stan_data, iter=n_iter,
#                            chains=n_chains, init=rep(stan_init, n_chains),
#                            warmup=n_burnin)
# }

ret <- stanc(model_file, model_name="tree_growth")
tree_growth_model <- stan_model(stanc_ret=ret)

seed <- 1638
stan_fit <- stan(model_file, data=stan_data, iter=500, chains=1,
                 init=rep(stan_init, n_chains), chain_id=1)
save(stan_fit, file="stan_fit.RData")

cl <- makeCluster(n_chains)
registerDoParallel(cl)
sflist <- foreach(n=1:n_chains, .packages=c("rstan")) %dopar% {
    # Add 1 to n in order to ensure chain_id 1 is not reused
    stan(stan_fit, seed=seed, chains=1, iter=n_iter, chain_id=n+1, refresh=-1)
}
stopCluster(cl)

stan_fit_allchains <- sflist2stanfit(sflist)
save(stan_fit_allchains, file="stan_fit_allchains.RData")
