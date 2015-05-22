# Forest growth model MCMC results

library(dplyr)
library(foreach)
library(matrixStats) # for weightedMedian
library(doParallel)
library(reshape2)
library(RPostgreSQL)

cl <- makeCluster(20)
registerDoParallel(cl)

plot_width <- 3.5
plot_height <- 3
plot_dpi <- 300

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

# con <- dbConnect(PostgreSQL(), dbname='tree_growth', user=pgsqluser, 
# password=pgsqlpwd)
# dbSendQuery(con, paste0("DROP TABLE IF EXISTS stems"))
# dbSendQuery(con, paste0("CREATE TABLE stems (site_id integer, plot_id integer, genus_id integer)"))

###############################################################################
## Full model (interactions, uncorrelated random effects)

dbhs <- seq(10, 120, 1)
precips <- c(7.5, 15, 22.5)
temps <- c(10:32)

# this_model <- 'tmn'
# this_plot <- 1

preds <- foreach(this_model=c('tmn', 'tmp', 'tmx'), .combine=rbind) %:%
    foreach(this_plot=c(1:82), .combine=rbind,
            .packages=c('foreach', 'dplyr', 'RPostgreSQL', 'reshape2', 
                        'matrixStats'),
            .inorder=FALSE) %dopar% {
        pg_src <- src_postgres('tree_growth', user=pgsqluser, 
                               password=pgsqlpwd)
        params <- tbl(pg_src, 'interact')
        params <- filter(params, model == this_model)
        stems <- tbl(pg_src, paste0('stems_', this_model))

        # Genus-level random effects
        B_g_params <- filter(params, parameter_base == 'B_g') %>% 
            rename(genus_id=row_id, param_id=col_id)

        # Site-level intercepts
        int_k_params <- filter(params, parameter_base == 'int_k') %>% 
            rename(site_id=row_id) %>% select(-col_id) %>% collect()

        # Plot-level intercepts
        int_jk_params <- filter(params, parameter_base == 'int_jk') %>% 
            rename(plot_id=row_id) %>% select(-col_id) %>% collect()

        # Elevation difference between plot and CRU cell center
        B_k_params <- filter(params, parameter_base == 'B_k') %>% 
            rename(site_id=row_id) %>% select(-col_id) %>% collect()


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

        this_site <- site_key$site_id[match(this_plot, site_key$plot_id)]
        this_site_char <- site_key$site_ID_char[match(this_plot, site_key$plot_id)]
        this_plot_char <- site_key$plot_ID_char[match(this_plot, site_key$plot_id)]

        this_elev_diff <- site_key$elev_diff[match(this_plot, site_key$plot_id)]

        genus_weights <- filter(stems, plot_id == this_plot) %>%
            group_by(genus_id) %>%
            summarize(n=n()) %>%
            ungroup() %>%
            mutate(weight=n/sum(n)) %>%
            select(-n)

        # Calculate weighted median values of each parameter for each iteration 
        # (combining the separate chains together). Each iteration has multiple 
        # estimates of each parameter as there is one estimated parameter value 
        # for each genus.
        B_g_betas <- filter(left_join(genus_weights, B_g_params)) %>%
            mutate(param=sql("parameter_base || '_' || param_id"),
                   estimate_id=sql("chain || '_' || iteration")) %>%
            arrange(param, estimate_id, genus_id) %>%
            collect() %>%
            dcast(param + estimate_id ~ genus_id, value.var="value") %>%
            select(-estimate_id)
        B_g_betas <- data.frame(param=paste0(B_g_betas$param, '_median'),
                                value=rowWeightedMedians(as.matrix(B_g_betas[-1]), w=genus_weights$weight))

        # Convert to matrix for linear algebra
        B <- matrix(B_g_betas$value, nrow=length(unique(B_g_betas$param)), byrow=TRUE)

        # Calculate intercept
        intercept <- filter(int_k_params, site_id == this_site, model == this_model)$value +
            filter(int_jk_params, plot_id == this_plot, model == this_model)$value

        # Calculate elevation adjustment for this plot
        elev_adjust <- filter(B_k_params, site_id == this_site)$value * (this_elev_diff - elev_diff_mean)

        intercept <- intercept + elev_adjust

        dbhs_centered <- dbhs - dbh_mean
        temps_centered <- temps - temp_mean
        # Change MCWD units from cm to mm
        precips_centered <- (precips * mm_per_unit) - precip_mean

        X <- foreach(dbh=dbhs_centered, .combine=rbind) %:%
            foreach(temp=temps_centered, .combine=rbind) %:%
            foreach(precip=precips_centered, .combine=rbind) %do% {
                X <- matrix(1, ncol=1) # genus-level intercept
                X <- cbind(X, precip)
                X <- cbind(X, precip^2)
                X <- cbind(X, temp)
                X <- cbind(X, temp^2)
                X <- cbind(X, dbh)
                X <- cbind(X, dbh^2)
                X <- cbind(X, precip*dbh)
                X <- cbind(X, temp*dbh)
        }
        these_preds <- X %*% B 
        these_preds <- these_preds + matrix(rep(intercept, nrow(X)), 
                                            nrow=nrow(X), byrow=TRUE)
        # X[, 6] (dbh) is subtracted so that output is growth increment
        these_preds <- these_preds - matrix(rep(X[, 6], each=ncol(B)), 
                                            ncol=ncol(B), byrow=TRUE)
        medians <- apply(these_preds, 1, median)
        q2pt5 <- apply(these_preds, 1, quantile, .025)
        q97pt5 <- apply(these_preds, 1, quantile, .975)
        these_preds <- data.frame(model=this_model,
                                  plot_id=this_plot_char, 
                                  site_id=this_site_char,
                                  precip=(X[, 2] + precip_mean)/mm_per_unit,
                                  temp=X[, 4] + temp_mean,
                                  dbh=X[, 6] + dbh_mean,
                                  median=medians,
                                  q2pt5,
                                  q97pt5)

        head(filter(these_preds, precip == 15, temp == 20))

        return(these_preds)
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
