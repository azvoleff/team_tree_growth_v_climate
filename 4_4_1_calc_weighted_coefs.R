library(foreach)
library(dplyr)
library(reshape2)
library(matrixStats)
library(RPostgreSQL)
library(doParallel)

cl <- makeCluster(3)
registerDoParallel(cl)

pgsqlpwd <- as.character(read.table('~/pgsqlpwd')[[1]])
pgsqluser <- as.character(read.table('~/pgsqluser')[[1]])

###############################################################################
## Full model (interactions, uncorrelated random effects)
prefixes <- c('D:/azvoleff/Data', # CI-TEAM
              'H:/Data', # Buffalo drive
              'O:/Data', # Blue drive
              '/localdisk/home/azvoleff/Data') # vertica1
prefix <- prefixes[match(TRUE, unlist(lapply(prefixes, function(x) file_test('-d', x))))]
base_folder <- file.path(prefix, "TEAM", "Tree_Growth")

pg_src <- src_postgres('tree_growth', user=pgsqluser, password=pgsqlpwd)
stems <- tbl(pg_src, paste0('stems'))

load(file.path(base_folder, 'Data', 
               'model_data_standardizing_full-tmn_meanannual-mcwd_run12.RData'))

# de-standardize dbh then calculate dbh classes
stems <- mutate(stems, initial_dbh_destd=initial_dbh*dbh_sd + dbh_mean) %>%
    collect()

# Fix inaccuracies arising due to rounding error:
stems$initial_dbh_destd[abs(stems$initial_dbh_destd - 10) < .001] <- 10

stems <- mutate(stems,
                dbh_class=cut(initial_dbh_destd,
                              c(seq(10, 45, 5), seq(50, 120, 10)),
                              include.lowest=TRUE),
                dbh_class_wide=cut(initial_dbh_destd,
                                   c(10, 20, 30, 50, 80, 120),
                                   include.lowest=TRUE))

# Calculate plot-level weights, weighting within dbh classes
stems %>%
    filter(initial_dbh_destd < 120) %>%
    group_by(model, site_id, plot_id, dbh_class, genus_id) %>%
    # Count num indiv of each genus within each dbh_class within each plot  
    # within each site
    summarize(n=n()) %>%
    group_by(model, site_id, plot_id, dbh_class) %>%
    # Convert counts to plot-level weights within dbh classes
    mutate(weight=n/sum(n)) %>%
    collect() -> genus_weights_byplot_bydbh

# Calculate plot-level weights, weighting within (wide) dbh classes
stems %>%
    filter(initial_dbh_destd < 120) %>%
    group_by(model, site_id, plot_id, dbh_class_wide, genus_id) %>%
    # Count num indiv of each genus within each dbh_class within each plot  
    # within each site
    summarize(n=n()) %>%
    group_by(model, site_id, plot_id, dbh_class_wide) %>%
    # Convert counts to plot-level weights within dbh classes
    mutate(weight=n/sum(n)) %>%
    collect() -> genus_weights_byplot_bydbhwide

# Calculate site-level weights, weighting within dbh classes
stems %>%
    filter(initial_dbh_destd < 120) %>%
    group_by(model, site_id, dbh_class, genus_id) %>%
    # Count num indiv of each genus within each dbh_class within each plot  
    # within each site
    summarize(n=n()) %>%
    group_by(model, site_id, dbh_class) %>%
    # Convert counts to site-level weights within dbh classes
    mutate(weight=n/sum(n)) %>%
    collect() -> genus_weights_bysite_bydbh

# Calculate overall weights, weighting within dbh classes
group_by(genus_weights_bysite_bydbh, model, dbh_class, genus_id) %>%
    # Sum genus weights across sites
    summarise(weight=sum(weight)) %>%
    # Normalize weights overall
    mutate(weight=weight/sum(weight)) -> genus_weights_bydbh 

# Calculate overall weights, weighting all sites equally, ignoring dbh classes
stems %>%
    filter(initial_dbh_destd < 120) %>%
    group_by(model, site_id, genus_id) %>%
    # Count num indiv of each genus within each within each site
    summarize(n=n()) %>%
    group_by(model, site_id) %>%
    # Convert counts to site-level weights by genus
    mutate(weight=n/sum(n)) %>%
    collect() %>%
    # Convert site-level weights to overall weights
    group_by(model, genus_id) %>%
    # Normalize weights 
    summarise(weight=sum(weight)) %>%
    group_by(model) %>%
    # Normalize weights
    mutate(weight=weight/sum(weight)) %>%
    collect() -> genus_weights

# Check that weighting worked correctly (totals should be three as there are 
# three models:
# group_by(genus_weights, model) %>% summarise(sum(weight))
# group_by(genus_weights_bydbh, model) %>% summarise(sum(weight))
# group_by(genus_weights_bysite_bydbh, model) %>% summarise(sum(weight))

params <- tbl(src_postgres('tree_growth', user=pgsqluser, 
                           password=pgsqlpwd), 'interact')
Bs <- filter(params, parameter_base == 'B_g') %>% 
    group_by(chain) %>%
    # Thin chains (there are 6 chains so for 1000 samples only need every 
    # 6th sample from each chain)
    filter(iteration %% 6 == 0) %>%
    rename(genus_id=row_id, param_id=col_id) %>%
    collapse() %>%
    mutate(param=sql("parameter_base || '_' || param_id"),
           estimate_id=sql("chain || '_' || iteration")) %>%
    collapse() %>%
    select(model, param, estimate_id, genus_id, value) %>%
    arrange(model, param, estimate_id, genus_id) %>%
    collect() %>%
    dcast(model + param + estimate_id ~ genus_id, value.var="value") %>%
    select(-estimate_id)

###############################################################################
### Calculate coefs weighted overall
###############################################################################
B_g_betas <- foreach(this_model=c('tmn', 'tmp', 'tmx'), .combine=rbind,
                     .packages=c('dplyr', 'RPostgreSQL', 'reshape2', 
                                 'matrixStats', 'foreach'), .inorder=FALSE) %dopar% {
    these_Bs <- filter(Bs, model == this_model)

    these_weights <- filter(genus_weights, model == this_model)

    # Remember columns 1 and 2 of these_Bs are the model and parameter IDs, and 
    # columns 1 and 2 of these_weights are the model_id and genus_id
    stopifnot(names(these_Bs)[c(-1, -2)] == these_weights$genus_id)

    data.frame(model=these_Bs$model,
               param=paste0(these_Bs$param, '_median'),
               value=rowWeightedMedians(as.matrix(these_Bs[c(-1, -2)]),
                                        w=as.numeric(these_weights$weight)))
}
save(B_g_betas, file='B_g_betas_weighted.RData')

###############################################################################
### Calculate coefs weighted by dbh class (wide) and plot
###############################################################################
B_g_betas <- foreach(this_model=c('tmn', 'tmp', 'tmx'), .combine=rbind,
                     .packages=c('dplyr', 'RPostgreSQL', 'reshape2', 
                                 'matrixStats', 'foreach'), .inorder=FALSE) %dopar% {
    these_Bs <- filter(Bs, model == this_model)

    genus_weights_cast <- dcast(filter(genus_weights_byplot_bydbh, model == this_model),
                                site_id + plot_id + dbh_class_wide ~ genus_id, 
                                value.var="weight", fill=0)

    # Remember columns 1 and 2 of Bs are the model and parameter IDs, and 
    # columns 1 and 2 of genus_weights_cast are the plot_id and dbh_class_wide
    stopifnot(names(these_Bs)[c(-1, -2)] == names(genus_weights_cast[c(-1, -2, -3)]))

    these_B_g_betas <- foreach(n=1:nrow(genus_weights_cast), .combine=rbind) %do% {
        these_weights <- as.numeric(genus_weights_cast[c(-1, -2, -3)][n, ])
        data.frame(model=these_Bs$model,
                   param=paste0(these_Bs$param, '_median'),
                   site_id=genus_weights_cast$site_id[n],
                   plot_id=genus_weights_cast$plot_id[n],
                   dbh_class_wide=genus_weights_cast$dbh_class_wide[n],
                   value=rowWeightedMedians(as.matrix(these_Bs[c(-1, -2)]), 
                                            w=these_weights))
    }

}
save(B_g_betas, file='B_g_betas_weighted_bydbhwide_byplot.RData')

###############################################################################
### Calculate coefs weighted by dbh class and plot
###############################################################################
B_g_betas <- foreach(this_model=c('tmn', 'tmp', 'tmx'), .combine=rbind,
                     .packages=c('dplyr', 'RPostgreSQL', 'reshape2', 
                                 'matrixStats', 'foreach'), .inorder=FALSE) %dopar% {
    these_Bs <- filter(Bs, model == this_model)

    genus_weights_cast <- dcast(filter(genus_weights_byplot_bydbh, model == this_model),
                                site_id + plot_id + dbh_class ~ genus_id, 
                                value.var="weight", fill=0)

    # Remember columns 1 and 2 of Bs are the model and parameter IDs, and 
    # columns 1 and 2 of genus_weights_cast are the plot_id and dbh_class
    stopifnot(names(these_Bs)[c(-1, -2)] == names(genus_weights_cast[c(-1, -2, -3)]))

    these_B_g_betas <- foreach(n=1:nrow(genus_weights_cast), .combine=rbind) %do% {
        these_weights <- as.numeric(genus_weights_cast[c(-1, -2, -3)][n, ])
        data.frame(model=these_Bs$model,
                   param=paste0(these_Bs$param, '_median'),
                   site_id=genus_weights_cast$site_id[n],
                   plot_id=genus_weights_cast$plot_id[n],
                   dbh_class=genus_weights_cast$dbh_class[n],
                   value=rowWeightedMedians(as.matrix(these_Bs[c(-1, -2)]), 
                                            w=these_weights))
    }

}
save(B_g_betas, file='B_g_betas_weighted_bydbh_byplot.RData')

###############################################################################
### Calculate coefs weighted by dbh class and site
###############################################################################
B_g_betas <- foreach(this_model=c('tmn', 'tmp', 'tmx'), .combine=rbind,
                     .packages=c('dplyr', 'RPostgreSQL', 'reshape2', 
                                 'matrixStats', 'foreach'), .inorder=FALSE) %dopar% {
    these_Bs <- filter(Bs, model == this_model)

    genus_weights_cast <- dcast(filter(genus_weights_bysite_bydbh, model == this_model),
                                site_id + dbh_class ~ genus_id, value.var="weight", fill=0)

    # Remember columns 1 and 2 of Bs are the model and parameter IDs, and 
    # columns 1 and 2 of genus_weights_cast are the site_id and dbh_class
    stopifnot(names(these_Bs)[c(-1, -2)] == names(genus_weights_cast[c(-1, -2)]))

    these_B_g_betas <- foreach(n=1:nrow(genus_weights_cast), .combine=rbind) %do% {
        these_weights <- as.numeric(genus_weights_cast[c(-1, -2)][n, ])
        data.frame(model=these_Bs$model,
                   param=paste0(these_Bs$param, '_median'),
                   site_id=genus_weights_cast$site_id[n],
                   dbh_class=genus_weights_cast$dbh_class[n],
                   value=rowWeightedMedians(as.matrix(these_Bs[c(-1, -2)]), 
                                            w=these_weights))
    }

}
save(B_g_betas, file='B_g_betas_weighted_bydbh_bysite.RData')

###############################################################################
### Calculate coefs weighted by dbh class
###############################################################################
B_g_betas <- foreach(this_model=c('tmn', 'tmp', 'tmx'), .combine=rbind,
                     .packages=c('dplyr', 'RPostgreSQL', 'reshape2', 
                                 'matrixStats', 'foreach'), .inorder=FALSE) %dopar% {
    these_Bs <- filter(Bs, model == this_model)

    genus_weights_cast <- dcast(filter(genus_weights_bydbh, model == this_model),
                                dbh_class ~ genus_id, value.var="weight", fill=0)

    # Remember columns 1 and 2 of these_Bs are the model and parameter IDs, and 
    # column 1 of genus_weights_cast is the dbh_class
    stopifnot(names(these_Bs)[c(-1, -2)] == names(genus_weights_cast[-1]))

    these_B_g_Betas <- foreach(n=1:nrow(genus_weights_cast), 
                               .combine=rbind) %do% {
        these_weights <- as.numeric(genus_weights_cast[-1][n, ])
        this_dbh_class <- genus_weights_cast$dbh_class[n]
        data.frame(model=these_Bs$model,
                   param=paste0(these_Bs$param, '_median'),
                   dbh_class=this_dbh_class,
                   value=rowWeightedMedians(as.matrix(these_Bs[c(-1, -2)]), 
                                            w=these_weights))
    }
}
save(B_g_betas, file='B_g_betas_weighted_bydbh.RData')
