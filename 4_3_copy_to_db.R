library(dplyr)
library(foreach)
library(RPostgreSQL)
library(doParallel)

cl <- makeCluster(2)
registerDoParallel(cl)

pgsqlpwd <- as.character(read.table('~/pgsqlpwd')[[1]])
pgsqluser <- as.character(read.table('~/pgsqluser')[[1]])

prefixes <- c('D:/azvoleff/Data', # CI-TEAM
              'H:/Data', # Buffalo drive
              'O:/Data', # Blue drive
              '/localdisk/ci_share/azvoleff/Data', # vertica1
              '/localdisk/home/azvoleff/Data') # vertica1
prefix <- prefixes[match(TRUE, unlist(lapply(prefixes, function(x) file_test('-d', x))))]

base_folder <- file.path(prefix, "TEAM", "Tree_Growth")

# load(file.path(base_folder, "Extracted_Parameters", "parameter_estimates_simple.RData"))
# write.csv(params, file.path(base_folder, "Extracted_Parameters", "parameter_estimates_simple.csv"), row.names=FALSE)
#
# load(file.path(base_folder, "Extracted_Parameters", "parameter_estimates_interact.RData"))
# write.csv(params, file.path(base_folder, "Extracted_Parameters", "parameter_estimates_interact.csv"), row.names=FALSE)
#
# load(file.path(base_folder, "Extracted_Parameters", "parameter_estimates_correlated.RData"))
# write.csv(params, file.path(base_folder, "Extracted_Parameters", "parameter_estimates_correlated.csv"), row.names=FALSE)

###############################################################################
### Add results data
###############################################################################
foreach (model_type=c('simple', 'interact', 'correlated'),
         .packages=c('RPostgreSQL')) %dopar% {
    con <- dbConnect(PostgreSQL(), dbname='tree_growth', user=pgsqluser, 
                     password=pgsqlpwd)

    message(paste("Copying", model_type))
    dbSendQuery(con, paste0("DROP TABLE IF EXISTS ", model_type))
    create_qry <- paste0("CREATE TABLE ", model_type, "(
        iteration integer,
        chain integer,
        parameter text,
        value double precision,
        parameter_base text,
        row_ID integer,
        col_ID integer,
        model text
    )")
    dbSendQuery(con, create_qry)

    csv_file <- file.path(base_folder, "Extracted_Parameters", 
                          paste0("parameter_estimates_", model_type, ".csv"))

    dbSendQuery(con, paste0("COPY ", model_type, " FROM '", csv_file,
                            "' WITH (FORMAT CSV, HEADER);"))
    dbSendQuery(con, paste0("VACUUM ANALYZE ", model_type, ";"))

    idx_qry <- paste0("CREATE INDEX ON ", model_type, " (parameter);",
        "CREATE INDEX ON ", model_type, " (parameter_base);",
        "CREATE INDEX ON ", model_type, " (row_id);",
        "CREATE INDEX ON ", model_type, " (col_id);",
        "CREATE INDEX ON ", model_type, " (model);")
    dbSendQuery(con, idx_qry)
}

###############################################################################
### Add stems data
###############################################################################
prefixes <- c('D:/azvoleff/Data', # CI-TEAM
              'H:/Data', # Buffalo drive
              'O:/Data', # Blue drive
              '/localdisk/home/azvoleff/Data') # vertica1
prefix <- prefixes[match(TRUE, unlist(lapply(prefixes, function(x) file_test('-d', x))))]
base_folder <- file.path(prefix, "TEAM", "Tree_Growth")
stems <- foreach (temp_var=c('tmn', 'tmp', 'tmx'), .combine=rbind) %dopar% {
    load(file.path(base_folder, 'Data', paste0("model_data_wide_full-", 
                                               temp_var, 
                                               "_meanannual-mcwd_run12.RData")))

    # Use linear indexing to calculate index of first observation for each stem
    first_obs_row <- match(c(1:nrow(model_data$dbh)), model_data$obs_indices[, 1])
    # First climate obs is associated with second dbh observation
    second_obs_row <- match(c(1:nrow(model_data$dbh)), model_data$obs_indices[, 1]) + 1
    obs_indices_linear <- (model_data$obs_indices[, 2] - 1) * nrow(model_data$dbh) + model_data$obs_indices[, 1]
    initial_dbh <- model_data$dbh[obs_indices_linear[first_obs_row]]
    initial_precip <- model_data$precip[obs_indices_linear[second_obs_row]]
    initial_temp <- model_data$temp[obs_indices_linear[second_obs_row]]

    data.frame(model=temp_var,
               site_id=model_data$site_ID, 
               plot_id=model_data$plot_ID, 
               genus_id=model_data$genus_ID,
               initial_dbh=initial_dbh,
               initial_precip=initial_precip,
               initial_temp=initial_temp)
}

con <- dbConnect(PostgreSQL(), dbname='tree_growth', user=pgsqluser, password=pgsqlpwd)
if (dbExistsTable(con, 'stems')) dbRemoveTable(con, 'stems')
dbWriteTable(con, 'stems', stems)
dbSendQuery(con, paste0("VACUUM ANALYZE stems;"))
idx_qry <- paste0("CREATE INDEX ON stems (model);",
    "CREATE INDEX ON stems (site_id);",
    "CREATE INDEX ON stems (plot_id);",
    "CREATE INDEX ON stems (genus_id);")
dbSendQuery(con, idx_qry)
