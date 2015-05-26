# Forest growth model MCMC results
library(RColorBrewer)
library(ggmcmc)
library(gridExtra)
library(foreach)
library(dplyr)

prefixes <- c('D:/azvoleff/Data', # CI-TEAM
              'H:/Data', # Buffalo drive
              'O:/Data', # Blue drive
              '/localdisk/ci_share/azvoleff/Data', # vertica1
              '/localdisk/home/azvoleff/Data') # vertica1
prefix <- prefixes[match(TRUE, unlist(lapply(prefixes, function(x) file_test('-d', x))))]

base_folder <- file.path(prefix, "TEAM", "Tree_Growth")

plot_width <- 3.5
plot_height <- 3
plot_dpi <- 300

# How many mm of precip per 1 unit change? Convert precip from mm to cm in 
# order to make interpretation of results easier
mm_per_unit <- 10

# What thinning was used? This does not change the thinning - it only is used 
# to properly align the plots.
nThin <- 100

# What burn-in was used? This does not change the thinning - it only is used to 
# properly align the plots.
nBurnin  <- 100000

caterpillar <- function(mods, labels=NULL) {
    cis <- group_by(mods, model) %>%
        rename(Parameter=param) %>%
        do(ci(.)) %>%
        rename(param=Parameter)
    # cis$model <- ordered(cis$model, levels=rev(unique(cis$model)))
    # # Flip parameters so first parameters show up on top
    # cis$param <- ordered(cis$param, levels=rev(unique(cis$param)))
    p <- ggplot(cis, aes(x=reorder(param, rev(1:length(param))), 
                         y=median, colour=model, fill=model)) +
        theme_bw(base_size=8) +
        geom_point(position=position_dodge(width=.4), size=1.25) +
        geom_linerange(aes(ymin=Low, ymax=High), size=.75, position=position_dodge(width=.4)) +
        geom_linerange(aes(ymin=low, ymax=high), size=.25, position=position_dodge(width=.4)) +
        xlab('Parameter') +
        ylab('Effect (cm)') +
        scale_colour_manual(values=c('#882255', '#117733', '#88CCEE'))
    if (!is.null(labels)) {
        p <- p + scale_x_discrete(labels=labels) +
            theme(axis.text.x=element_text(vjust=0))
    }
    p <- p + coord_flip()
    p
}

load('B_g_betas_weighted_overall.RData')

B_g_labels <- c('mu_B_g[1]'='int.',
                'mu_B_g[2]'='P',
                'mu_B_g[3]'=expression(P^2),
                'mu_B_g[4]'='T',
                'mu_B_g[5]'=expression(T^2),
                'mu_B_g[6]'='D',
                'mu_B_g[7]'=expression(D^2),
                'mu_B_g[8]'=expression(D%*%P),
                'mu_B_g[9]'=expression(D%*%T))

# Change MCWD units from mm to cm
B_g_betas$value[B_g_betas$param == "B_g_2_mean"] <- B_g_betas$value[B_g_betas$param == "B_g_2_mean"] * mm_per_unit
B_g_betas$value[B_g_betas$param == "B_g_3_mean"] <- B_g_betas$value[B_g_betas$param == "B_g_3_mean"] * mm_per_unit^2
B_g_betas$value[B_g_betas$param == "B_g_8_mean"] <- B_g_betas$value[B_g_betas$param == "B_g_8_mean"] * mm_per_unit

p <- caterpillar(filter(B_g_betas, param %in% paste0('B_g_', c(2:9), '_median')),
                 labels=c('B_g_2_median'='P',
                          'B_g_3_median'=expression(P^2),
                          'B_g_4_median'='T',
                          'B_g_5_median'=expression(T^2),
                          'B_g_6_median'='D',
                          'B_g_7_median'=expression(D^2),
                          'B_g_8_median'=expression(D%*%P),
                          'B_g_9_median'=expression(D%*%T)))
ggsave('caterpillar_climate_interact_weighted.png', p, width=plot_width*1.5, 
       height=plot_height, dpi=plot_dpi)

p <- caterpillar(filter(B_g_betas, param %in% paste0('B_g_', c(2:5, 8:9), '_median')),
                 labels=c('B_g_2_median'='P',
                          'B_g_3_median'=expression(P^2),
                          'B_g_4_median'='T',
                          'B_g_5_median'=expression(T^2),
                          'B_g_8_median'=expression(D%*%P),
                          'B_g_9_median'=expression(D%*%T)))
ggsave('caterpillar_climate_interact_weighted_noD.png', p, width=plot_width*1.5, 
       height=plot_height, dpi=plot_dpi)

# Function to make a growth prediction as a function of a particular climate 
# variable
make_preds <- function(B, temp, precip, dbhs) {
    X <- matrix(rep(1, length(dbhs)), ncol=1)
    X <- cbind(X, precip)
    X <- cbind(X, precip^2)
    X <- cbind(X, temp)
    X <- cbind(X, temp^2)
    X <- cbind(X, dbhs)
    X <- cbind(X, dbhs^2)
    X <- cbind(X, precip*dbhs)
    X <- cbind(X, temp*dbhs)
    # dbhs is subtracted so that output is growth increment
    these_preds <- X %*% B - rep(dbhs, ncol(B))
    medians <- apply(these_preds, 1, median)
    q2pt5 <- apply(these_preds, 1, quantile, .025)
    q97pt5 <- apply(these_preds, 1, quantile, .975)
    data.frame(median=medians, q2pt5, q97pt5)
}

dbhs <- seq(10, 120, 1)
preds <- foreach(this_model=unique(B_g_betas$model), .combine=rbind) %do% {
    # Need to center the predictor variables
    load(file.path(base_folder, 'Data',
                   paste0("model_data_standardizing_full-", this_model, 
                          "_meanannual-mcwd_run12.RData")))

    # Model effects of temp variation:
    B <- filter(B_g_betas, model == this_model) %>%
        ungroup() %>%
        arrange(param) %>%
        select(param, value)
    # Convert to matrix for linear algebra
    B <- matrix(B$value, nrow=length(unique(B$param)), byrow=TRUE)

    dbhs_centered <- dbhs - dbh_mean

    # Remember temp is standardized, so round it then center it just so the 
    # units look pretty
    temps <- round(temp_mean) - temp_mean
    temps <- c(temps - 2, temps, temps + 2)
    temp_preds <- foreach(temp=temps, .combine=rbind) %do% {
        temp_preds <- make_preds(B, temp, (150 - precip_mean)/mm_per_unit, dbhs_centered)
        temp_preds <- cbind(Panel="Temperature", clim=temp + temp_mean, 
                            dbh=dbhs,
                            temp_preds)
        return(temp_preds)
    }

    # Center, and set units
    precips <- (c(75, 150, 150 + 75) - precip_mean)/mm_per_unit
    precip_preds <- foreach(precip=precips, .combine=rbind) %do% {
        precip_preds <- make_preds(B, 0, precip, dbhs_centered)
        precip_preds <- cbind(Panel="MCWD",
                              clim=precip + precip_mean/mm_per_unit,
                              dbh=dbhs,
                              precip_preds)
        return(precip_preds)
    }

    preds <- rbind(temp_preds, precip_preds)
    preds <- cbind(model=this_model, preds)

    return(preds)
}







dbhs <- seq(10, 120, 1)
precips <- c(7.5, 15, 22.5)
temps <- c(10:32)

# this_model <- 'tmn'
# this_plot <- 1
preds <- foreach(this_model=c('tmn', 'tmp', 'tmx'), .combine=rbind,
                 .packages=c('foreach', 'dplyr', 'RPostgreSQL', 'reshape2', 
                             'matrixStats'), .inorder=FALSE) %dopar% {
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


















preds$clim <- factor(sprintf('%.02f', preds$clim))
preds$model <- factor(preds$model, levels=c('tmn', 'tmp', 'tmx'), 
                      labels=c('Min. Temp. Model',
                               'Mean Temp. Model',
                               'Max. Temp. Model'))
ann_text <- group_by(preds, Panel, model) %>%
    summarise(lab=paste0('Lower=', unique(clim)[1], '\nUpper=', unique(clim)[2]))
ann_text$dbh <- 90
ann_text$mean <- 5

# ggplot(preds, aes(x=dbh, y=mean)) +
#     geom_line(aes(colour=clim)) +
#     geom_ribbon(aes(ymin=q2pt5, ymax=q97pt5, fill=clim), alpha=.2) +
#     facet_grid(Panel~model) +
#     ylab('Increment (cm)') +
#     geom_text(data=ann_text, aes(label=lab))
#
# head(preds)

ps <- foreach(this_panel=unique(preds$Panel), .combine=c) %:%
    foreach(this_model=unique(preds$model)) %do% {
    p <- ggplot(filter(preds, this_model == model, this_panel == Panel), aes(x=dbh, y=median)) +
        geom_line(aes(colour=clim)) +
        geom_ribbon(aes(ymin=q2pt5, ymax=q97pt5, fill=clim), alpha=.2) +
        facet_wrap(~model) +
        xlab('Initial size (cm)') +
        ylab('Growth increment (cm)') +
        coord_cartesian(ylim=c(0, .7)) +
        theme(legend.position=c(.20, .2)) +
        guides(fill=guide_legend(this_panel),
               colour=guide_legend(this_panel))
    return(p)
}

p <- arrangeGrob(ps[[1]], ps[[2]], ps[[3]],
                 ps[[4]], ps[[5]], ps[[6]], ncol=3, nrow=2)
ggsave(paste0('predicted_growth_increments.png'), width=plot_width*2.5,
       height=plot_height*2.5, dpi=plot_dpi, plot=p)
ggsave(paste0('predicted_growth_increments.svg'), width=plot_width*2.5,
       height=plot_height*2.5, dpi=plot_dpi, plot=p)

###############################################################################
## Full model (correlated random effects)

# load(file.path(base_folder, "Extracted_Parameters", "parameter_estimates_correlated.RData"))
