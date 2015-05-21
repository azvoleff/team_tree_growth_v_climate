library(foreach)
library(RPostgreSQL)

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

con <- dbConnect(PostgreSQL(), dbname='tree_growth', user=pgsqluser, 
                 password=pgsqlpwd)

make_table <- function(con, model_type) {
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
}

create_indices <- function(con, model_type) {
    idx_qry <- paste0("CREATE INDEX ON ", model_type, " (parameter);",
        "CREATE INDEX ON ", model_type, " (parameter_base);",
        "CREATE INDEX ON ", model_type, " (row_id);",
        "CREATE INDEX ON ", model_type, " (col_id);",
        "CREATE INDEX ON ", model_type, " (model);")
    dbSendQuery(con, idx_qry)
}

copy_csv <- function(con, csv_file, model_type) {
    dbSendQuery(con, paste0("COPY ", model_type, " FROM '", csv_file,
                            "' WITH (FORMAT CSV, HEADER);"))
    dbSendQuery(con, paste0("VACUUM ANALYZE ", model_type, ";"))
}

copy_stems <- function(con, temp_var, stems) {
    stem_table <- paste('stems', temp_var, sep='_')
    if (dbExistsTable(con, stem_table)) dbRemoveTable(con, stem_table)
    dbWriteTable(con, stem_table, stems)
    dbSendQuery(con, paste0("VACUUM ANALYZE ", stem_table, ";"))
    idx_qry <- paste0("CREATE INDEX ON ", stem_table, " (site_id);",
        "CREATE INDEX ON ", stem_table, " (plot_id);",
        "CREATE INDEX ON ", stem_table, " (genus_id);")
    dbSendQuery(con, idx_qry)
}

foreach (model_type=c('simple', 'interact', 'correlated')) %do% {
    message(paste("Copying", model_type))
    make_table(con, model_type)
    csv_file <- file.path(base_folder, "Extracted_Parameters", 
                          paste0("parameter_estimates_", model_type, ".csv"))
    copy_csv(con, csv_file, model_type)
    create_indices(con, model_type)
}

foreach (temp_var=c('tmn', 'tmp', 'tmx')) %do% {
    load(file.path(data_folder, paste0("model_data_wide_full-", this_model, 
                                       "_meanannual-mcwd_run12.RData")))
    stems <- data.frame(site_id=model_data$site_ID, 
                        plot_id=model_data$plot_ID, 
                        genus_id=model_data$genus_ID)
    copy_stems(con, temp_var, stems)
}
