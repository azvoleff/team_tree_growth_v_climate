library(foreach)
library(RPostgreSQL)

pgsqlpwd <- as.character(read.table('~/pgsqlpwd')[[1]])

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

con <- dbConnect(PostgreSQL(), dbname='tree_growth', user='cistaff', 
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

foreach (model_type=c('simple', 'interact', 'correlated')) %do% {
    message(paste("Copying", model_type))
    make_table(con, model_type)
    csv_file <- file.path(base_folder, "Extracted_Parameters", 
                          paste0("parameter_estimates_", model_type, ".csv"))
    copy_csv(con, csv_file, model_type)
    create_indices(con, model_type)
}

# Below is dplyr command
library(dplyr)
pgsrc <- src_postgres('tree_growth', user='azvoleff', password='DEBUGent')
simple_data <- tbl(pgsrc, 'interact')
