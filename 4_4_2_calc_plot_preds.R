# Forest growth model MCMC results

library(dplyr)
library(foreach)
library(reshape2)
library(stringr)
library(RPostgreSQL)

# How many mm of precip per 1 unit change? Convert precip from mm to cm in 
# order to make interpretation of results easier
mm_per_unit <- 10

prefixes <- c('D:/azvoleff/Data', # CI-TEAM
              'H:/Data', # Buffalo drive
              'O:/Data', # Blue drive
              '/localdisk/home/azvoleff/Data') # vertica1
prefix <- prefixes[match(TRUE, unlist(lapply(prefixes, function(x) file_test('-d', x))))]

base_folder <- file.path(prefix, "TEAM", "Tree_Growth")
data_folder <- file.path(base_folder, "Data")

pgsqlpwd <- as.character(read.table('~/pgsqlpwd')[[1]])
pgsqluser <- as.character(read.table('~/pgsqluser')[[1]])

###############################################################################
## Full model (interactions, uncorrelated random effects)
pg_src <- src_postgres('tree_growth', user=pgsqluser, 
                       password=pgsqlpwd)
params <- tbl(pg_src, 'interact') %>%
    # Thin chains (there are 6 chains so for 1000 samples only need every 6th 
    # sample from each chain)
    filter(iteration %% 6 == 0)

# Site-level intercepts
int_k_params <- filter(params, parameter_base == 'int_k') %>% 
    rename(site_id=row_id) %>% select(-col_id) %>% collect()

# Plot-level intercepts
int_jk_params <- filter(params, parameter_base == 'int_jk') %>% 
    rename(plot_id=row_id) %>% select(-col_id) %>% collect()

# Elevation difference between plot and CRU cell center
B_k_params <- filter(params, parameter_base == 'B_k') %>% 
    rename(site_id=row_id) %>% select(-col_id) %>% collect()

clim <- tbl(pg_src, 'climate')

load('B_g_betas_weighted_bydbhwide_byplot.RData')
B_g_betas_wide <- rename(B_g_betas, dbh_class=dbh_class_wide)

load('B_g_betas_weighted_bydbh_byplot.RData')
B_g_betas <- rbind(B_g_betas, B_g_betas_wide)
rm(B_g_betas_wide)

# Calculate dbhs to predict for from midpoints of dbh classes
dbhs <- data.frame(dbh_class=levels(B_g_betas$dbh_class))
dbhs$low <- as.numeric(gsub('^.', '', str_extract(dbhs$dbh_class, '^.[0-9]*')))
dbhs$high <- as.numeric(gsub('.$', '', str_extract(dbhs$dbh_class, '[0-9]*.$')))
dbhs$mid <- (dbhs$high - dbhs$low)/2 + dbhs$low

precips <- c(0, 7.5, 15, 22.5)
temp_diffs <- c(0:6)

preds <- foreach(this_model=c('tmn', 'tmp', 'tmx'), .combine=rbind) %do% {
    suffix <- paste0("_full-", this_model, "_meanannual-mcwd_run12")
    # Need to key plot_id to site_id
    plot_id_factor_key <- read.csv(file=file.path(data_folder,
                                                  paste0("plot_ID_factor_key", 
                                                         suffix, ".csv")))
    # Need to key plot_id to site_id
    site_id_factor_key <- read.csv(file=file.path(data_folder, 
                                                  paste0("site_ID_factor_key", 
                                                         suffix, ".csv")))
    # Need elevation data
    load(file.path(data_folder, paste0("model_data_elev_key", suffix, ".RData")))

    # Need to center the predictor variables
    load(file.path(data_folder, paste0("model_data_standardizing", suffix, 
                                       ".RData")))

    plot_id_factor_key$site_ID_char <- gsub('(VG)|([0-9]*)', '', plot_id_factor_key$plot_ID_char)

    site_key <- left_join(plot_id_factor_key, rename(elev_key, plot_ID_char=plot_ID)) %>%
        left_join(site_id_factor_key, by="site_ID_char") %>%
        rename(plot_id=plot_ID_numeric, site_id=site_ID_numeric)

    int_k_params_thismodel <- filter(int_k_params, model == this_model)
    int_jk_params_thismodel <- filter(int_jk_params, model == this_model)
    B_k_params_thismodel <- filter(B_k_params, model == this_model)

    preds_thismodel <- foreach(this_plot=c(1:82), .combine=rbind) %do% {
        print(paste(this_model, this_plot))
        this_site <- site_key$site_id[match(this_plot, site_key$plot_id)]
        this_site_char <- site_key$site_ID_char[match(this_plot, site_key$plot_id)]
        this_plot_char <- site_key$plot_ID_char[match(this_plot, site_key$plot_id)]

        this_elev_diff <- site_key$elev_diff[match(this_plot, site_key$plot_id)]

        B_g_betas_thisplot <- filter(B_g_betas, model == this_model, plot_id == this_plot)

        plot_temp_mean <- collect(filter(clim, plot_id == this_plot, variable == this_model))$mean
        # Use absolute value as in model MCWD is defined as positive
        plot_precip_mean <- abs(collect(filter(clim, plot_id == this_plot, variable == 'mcwd_run12'))$mean)

        dbhs_centered <- dbhs$mid - dbh_mean
        temps_centered <- (plot_temp_mean - temp_mean) + temp_diffs
        # Change MCWD units from cm to mm
        precips_centered <- c(plot_precip_mean, precips * mm_per_unit) - precip_mean

        # Calculate intercept
        intercept <- filter(int_k_params_thismodel, site_id == this_site, model == this_model)$value +
            filter(int_jk_params_thismodel, plot_id == this_plot, model == this_model)$value

        # Calculate elevation adjustment for this plot
        elev_adjust <- filter(B_k_params_thismodel, site_id == this_site)$value * (this_elev_diff - elev_diff_mean)

        intercept <- intercept + elev_adjust

        preds_thisplot <- foreach(this_dbh_class=unique(B_g_betas_thisplot$dbh_class), 
                               .combine=rbind) %do% {
            this_dbh <- dbhs_centered[match(this_dbh_class, dbhs$dbh_class)]

            # Convert to matrix for linear algebra
            B <- matrix(filter(B_g_betas_thisplot, dbh_class == this_dbh_class)$value, 
                        nrow=length(unique(B_g_betas_thisplot$param)), byrow=TRUE)

            X <- foreach(temp=temps_centered, .combine=rbind) %:%
                foreach(precip=precips_centered, .combine=rbind) %do% {
                    X <- matrix(1, ncol=1) # genus-level intercept
                    X <- cbind(X, precip)
                    X <- cbind(X, precip^2)
                    X <- cbind(X, temp)
                    X <- cbind(X, temp^2)
                    X <- cbind(X, this_dbh)
                    X <- cbind(X, this_dbh^2)
                    X <- cbind(X, precip*this_dbh)
                    X <- cbind(X, temp*this_dbh)
            }
            preds_thisdbh <- X %*% B 
            preds_thisdbh <- preds_thisdbh + matrix(rep(intercept, nrow(X)), 
                                                nrow=nrow(X), byrow=TRUE)
            # X[, 6] (dbh) is subtracted so that output is growth increment
            preds_thisdbh <- preds_thisdbh - matrix(rep(X[, 6], each=ncol(B)), 
                                                ncol=ncol(B), byrow=TRUE)
            medians <- apply(preds_thisdbh, 1, median)
            q2pt5 <- apply(preds_thisdbh, 1, quantile, .025)
            q97pt5 <- apply(preds_thisdbh, 1, quantile, .975)
            q5 <- apply(preds_thisdbh, 1, quantile, .05)
            q95 <- apply(preds_thisdbh, 1, quantile, .95)
            preds_thisdbh <- data.frame(model=this_model,
                                        plot_id=this_plot_char, 
                                        site_id=this_site_char,
                                        precip=(X[, 2] + precip_mean)/mm_per_unit,
                                        precip_diff=(X[, 2] + precip_mean - plot_precip_mean)/mm_per_unit,
                                        temp=X[, 4] + temp_mean,
                                        temp_diff=X[, 4] + temp_mean - plot_temp_mean,
                                        dbh_class=this_dbh_class,
                                        dbh=X[, 6] + dbh_mean,
                                        median=medians,
                                        q2pt5,
                                        q97pt5,
                                        q5,
                                        q95)

            return(preds_thisdbh)
        }
        return(preds_thisplot)
    }
    return(preds_thismodel)
}

save(preds, file='growth_predictions_byplot.RData')

con <- dbConnect(PostgreSQL(), dbname='tree_growth', user=pgsqluser, 
                 password=pgsqlpwd)
if (dbExistsTable(con, "preds_interact")) dbRemoveTable(con, "preds_interact")
dbWriteTable(con, "preds_interact", preds)
dbSendQuery(con, paste0("VACUUM ANALYZE preds_interact;"))
idx_qry <- paste0("CREATE INDEX ON preds_interact (model);",
    "CREATE INDEX ON preds_interact (plot_id);",
    "CREATE INDEX ON preds_interact (site_id);",
    "CREATE INDEX ON preds_interact (precip);",
    "CREATE INDEX ON preds_interact (temp);",
    "CREATE INDEX ON preds_interact (dbh);")
dbSendQuery(con, idx_qry)
