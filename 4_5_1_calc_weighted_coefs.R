library(foreach)
library(dplyr)
library(matrixStats)
library(RPostgreSQL)
library(doParallel)

cl <- makeCluster(3)
registerDoParallel(cl)

pgsqlpwd <- as.character(read.table('~/pgsqlpwd')[[1]])
pgsqluser <- as.character(read.table('~/pgsqluser')[[1]])

###############################################################################
## Full model (interactions, uncorrelated random effects)

pg_src <- src_postgres('tree_growth', user=pgsqluser, password=pgsqlpwd)

stems <- tbl(pg_src, paste0('stems'))

# Calculate site-level genus weights, weighting all plots equally within sites
stems %>%
    group_by(model, site_id, plot_id, genus_id) %>%
    # Count num indiv of each genus within plot
    summarize(n=n()) %>%
    ungroup() %>%
    group_by(model, site_id, plot_id) %>%
    # Convert counts to plot-level weights
    mutate(weight=n/sum(n)) %>%
    collect() %>%
    group_by(model, site_id, genus_id) %>%
    # Convert plot-level weights to site-level weights for each genus
    summarise(weight=sum(weight)) %>%
    group_by(model, site_id) %>%
    # Normalize weights by site
    mutate(weight=weight/sum(weight)) -> genus_weights_bysite

# Calculate overall genus weights, weighting all sites equally
genus_weights <- group_by(genus_weights_bysite, genus_id) %>%
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
    # Remember columns 1 and 2 are the model and parameter IDs
    stopifnot(names(Bs)[c(-1, -2)] == genus_weights$genus_id)
    data.frame(model=Bs$model,
               param=paste0(Bs$param, '_median'),
               value=rowWeightedMedians(as.matrix(Bs[c(-1, -2)]), 
                                        w=genus_weights$weight))
}

save(B_g_betas, 'B_g_betas_weighted_overall.RData')
