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

stems <- mutate(stems, dbh_class=cut(initial_dbh_destd,
                                     c(seq(10, 45, 5), seq(50, 120, 10)),
                                     include.lowest=TRUE))

# Fix error due to rounding error:

# Calculate site-level genus weights, weighting all plots equally within sites
stems %>%
    filter(initial_dbh_destd < 120) %>%
    group_by(model, site_id, dbh_class, genus_id) %>%
    # Count num indiv of each genus within each dbh_class within each plot  
    # within each site
    summarize(n=n()) %>%
    ungroup() %>%
    group_by(model, site_id, dbh_class) %>%
    # Convert counts to site-level weights within dbh classes
    mutate(weight=n/sum(n)) %>%
    collect() -> genus_weights_bysite

# Calculate overall genus weights, weighting all sites equally
genus_weights <- group_by(genus_weights_bysite, dbh_class, genus_id) %>%
    # Sum genus weights across sites
    summarise(weight=sum(weight)) %>%
    # Normalize weights overall
    mutate(weight=weight/sum(weight))

# Check that weighting worked correctly (totals should be three as there are 
# three models:
#group_by(genus_weights, site_id) %>% summarise(sum(weight))
B_g_betas <- foreach(this_model=c('tmn', 'tmp', 'tmx'), .combine=rbind,
                     .packages=c('dplyr', 'RPostgreSQL', 'reshape2', 
                                 'matrixStats'), .inorder=FALSE) %dopar% {
    params <- tbl(src_postgres('tree_growth', user=pgsqluser, 
                               password=pgsqlpwd), 'interact')
    Bs <- filter(params, parameter_base == 'B_g', model == this_model) %>% 
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

    genus_weights_cast <- dcast(genus_weights, dbh_class ~ genus_id, 
                                value.var="weight", fill=0)

    # Remember columns 1 and 2 of Bs are the model and parameter IDs, and 
    # column 1 of genus_weights_cast is the dbh_class
    stopifnot(names(Bs)[c(-1, -2)] == names(genus_weights_cast[-1]))

    these_B_g_Betas <- foreach(n=1:nrow(genus_weights_cast), 
                               .combine=rbind) %do% {
        these_weights <- as.numeric(genus_weights_cast[-1][n, ])
        this_dbh_class <- genus_weights_cast$dbh_class[n]
        data.frame(model=Bs$model,
                   param=paste0(Bs$param, '_median'),
                   dbh_class=this_dbh_class,
                   value=rowWeightedMedians(as.matrix(Bs[c(-1, -2)]), 
                                            w=these_weights))
    }
}

save(B_g_betas, file='B_g_betas_weighted_overall.RData')
