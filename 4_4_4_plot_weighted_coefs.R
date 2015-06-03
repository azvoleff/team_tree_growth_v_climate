# Forest growth model MCMC results
library(RColorBrewer)
library(ggplot2)
library(ggmcmc)
library(gridExtra)
library(foreach)
library(dplyr)
library(reshape2)

prefixes <- c('D:/azvoleff/Data', # CI-TEAM
              'H:/Data', # Buffalo drive
              'O:/Data', # Blue drive
              '/localdisk/ci_share/azvoleff/Data', # vertica1
              '/localdisk/home/azvoleff/Data') # vertica1
prefix <- prefixes[match(TRUE, unlist(lapply(prefixes, function(x) file_test('-d', x))))]

base_folder <- file.path(prefix, "TEAM", "Tree_Growth")

plot_width <- 3.5
plot_height <- 3
plot_dpi <- 600

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
        geom_point(aes(shape=model), position=position_dodge(width=.4), size=1.25) +
        geom_linerange(aes(ymin=Low, ymax=High), size=.5, position=position_dodge(width=.4)) +
        geom_linerange(aes(ymin=low, ymax=high), size=.25, position=position_dodge(width=.4)) +
        xlab('Parameter') +
        ylab('Effect (cm)') +
        scale_shape_manual("Model", values=c(23, 4, 17)) + 
        scale_colour_brewer("Model", palette = "Dark2") + 
        scale_fill_brewer("Model", palette = "Dark2")
    if (!is.null(labels)) {
        p <- p + scale_x_discrete(labels=labels) +
            theme(axis.text.x=element_text(vjust=0))
    }
    p <- p + coord_flip()
    p
}

load('B_g_betas_weighted.RData')

B_g_betas$model <- factor(B_g_betas$model, levels=c('tmn', 'tmp', 'tmx'), 
                          labels=c('Min. Temp.',
                                   'Mean Temp.',
                                   'Max. Temp.'))

# Change MCWD units from mm to cm
B_g_betas$value[B_g_betas$param == "B_g_2_median"] <- B_g_betas$value[B_g_betas$param == "B_g_2_median"] * mm_per_unit
B_g_betas$value[B_g_betas$param == "B_g_3_median"] <- B_g_betas$value[B_g_betas$param == "B_g_3_median"] * mm_per_unit^2
B_g_betas$value[B_g_betas$param == "B_g_8_median"] <- B_g_betas$value[B_g_betas$param == "B_g_8_median"] * mm_per_unit

p <- caterpillar(filter(B_g_betas, param %in% paste0('B_g_', c(1:9), '_median')),
                 labels=c('B_g_1_median'='Intercept',
                          'B_g_2_median'='P',
                          'B_g_3_median'=expression(P^2),
                          'B_g_4_median'='T',
                          'B_g_5_median'=expression(T^2),
                          'B_g_6_median'='D',
                          'B_g_7_median'=expression(D^2),
                          'B_g_8_median'=expression(D%*%P),
                          'B_g_9_median'=expression(D%*%T)))
ggsave('caterpillar_climate_interact_weighted.png', p, width=plot_width*1.5, 
       height=plot_height, dpi=plot_dpi)
ggsave('caterpillar_climate_interact_weighted.pdf', p, width=plot_width*1.5, 
       height=plot_height, dpi=plot_dpi)

p <- caterpillar(filter(B_g_betas, param %in% paste0('B_g_', c(1:9), '_median')),
                 labels=c('B_g_1_median'='Intercept',
                          'B_g_2_median'='P',
                          'B_g_3_median'=expression(P^2),
                          'B_g_4_median'='T',
                          'B_g_5_median'=expression(T^2),
                          'B_g_6_median'='D',
                          'B_g_7_median'=expression(D^2),
                          'B_g_8_median'=expression(D%*%P),
                          'B_g_9_median'=expression(D%*%T)))
ggsave('caterpillar_climate_interact_weighted.png', p, width=plot_width*1.5, 
       height=plot_height, dpi=plot_dpi)
ggsave('caterpillar_climate_interact_weighted.pdf', p, width=plot_width*1.5, 
       height=plot_height, dpi=plot_dpi)


p <- caterpillar(filter(B_g_betas, param %in% paste0('B_g_', c(1:5, 7:9), '_median')),
                 labels=c('B_g_1_median'='Intercept',
                          'B_g_2_median'='P',
                          'B_g_3_median'=expression(P^2),
                          'B_g_4_median'='T',
                          'B_g_5_median'=expression(T^2),
                          'B_g_7_median'=expression(D^2),
                          'B_g_8_median'=expression(D%*%P),
                          'B_g_9_median'=expression(D%*%T)))
ggsave('caterpillar_climate_interact_weighted_noD.png', p, width=plot_width*1.5, 
       height=plot_height, dpi=plot_dpi)
ggsave('caterpillar_climate_interact_weighted_noD.pdf', p, width=plot_width*1.5, 
       height=plot_height, dpi=plot_dpi)
ggsave('caterpillar_climate_interact_weighted_noD.svg', p, width=plot_width*1.5, 
       height=plot_height, dpi=plot_dpi)

p <- caterpillar(filter(B_g_betas, param %in% paste0('B_g_6_median')),
                 labels=c('B_g_6_median'='D'))
ggsave('caterpillar_climate_interact_weighted_onlyD.png', p, width=plot_width*1.5, 
       height=plot_height/2, dpi=plot_dpi)
ggsave('caterpillar_climate_interact_weighted_onlyD.pdf', p, width=plot_width*1.5, 
       height=plot_height/2, dpi=plot_dpi)
ggsave('caterpillar_climate_interact_weighted_onlyD.svg', p, width=plot_width*1.5, 
       height=plot_height/2, dpi=plot_dpi)



# Plot standardized coefficients
load(file.path(base_folder, 'Data',
               paste0("model_data_standardizing_full-tmn_meanannual-mcwd_run12.RData")))
B_g_betas_std <- filter(B_g_betas, param %in% paste0('B_g_', c(1:9), '_median'))
B_g_betas_std$value[B_g_betas_std$param == "B_g_2_median"] <- B_g_betas_std$value[B_g_betas_std$param == "B_g_2_median"] / mm_per_unit
B_g_betas_std$value[B_g_betas_std$param == "B_g_3_median"] <- B_g_betas_std$value[B_g_betas_std$param == "B_g_3_median"] / mm_per_unit^2
B_g_betas_std$value[B_g_betas_std$param == "B_g_8_median"] <- B_g_betas_std$value[B_g_betas_std$param == "B_g_8_median"] / mm_per_unit
B_g_betas_std <- filter(B_g_betas, param %in% paste0('B_g_', c(1:9), '_median'))
B_g_betas_std$value[B_g_betas_std$param == 'B_g_1_median'] <- B_g_betas_std$value[B_g_betas_std$param == 'B_g_1_median'] / dbh_sd
B_g_betas_std$value[B_g_betas_std$param == 'B_g_2_median'] <- B_g_betas_std$value[B_g_betas_std$param == 'B_g_2_median'] / (dbh_sd/precip_sd)
B_g_betas_std$value[B_g_betas_std$param == 'B_g_3_median'] <- B_g_betas_std$value[B_g_betas_std$param == 'B_g_3_median'] / (dbh_sd/precip_sd^2)
B_g_betas_std$value[B_g_betas_std$param == 'B_g_4_median'] <- B_g_betas_std$value[B_g_betas_std$param == 'B_g_4_median'] / (dbh_sd/temp_sd)
B_g_betas_std$value[B_g_betas_std$param == 'B_g_5_median'] <- B_g_betas_std$value[B_g_betas_std$param == 'B_g_5_median'] / (dbh_sd/temp_sd^2)
B_g_betas_std$value[B_g_betas_std$param == 'B_g_6_median'] <- B_g_betas_std$value[B_g_betas_std$param == 'B_g_6_median'] / 1
B_g_betas_std$value[B_g_betas_std$param == 'B_g_7_median'] <- B_g_betas_std$value[B_g_betas_std$param == 'B_g_7_median'] / (1/dbh_sd)
B_g_betas_std$value[B_g_betas_std$param == 'B_g_8_median'] <- B_g_betas_std$value[B_g_betas_std$param == 'B_g_8_median'] / (1/precip_sd)
B_g_betas_std$value[B_g_betas_std$param == 'B_g_9_median'] <- B_g_betas_std$value[B_g_betas_std$param == 'B_g_9_median'] / (1/temp_sd)

p <- caterpillar(filter(B_g_betas_std, param %in% paste0('B_g_', c(1:9), '_median')),
                 labels=c('B_g_1_median'='Intercept',
                          'B_g_2_median'='P',
                          'B_g_3_median'=expression(P^2),
                          'B_g_4_median'='T',
                          'B_g_5_median'=expression(T^2),
                          'B_g_6_median'='D',
                          'B_g_7_median'=expression(D^2),
                          'B_g_8_median'=expression(D%*%P),
                          'B_g_9_median'=expression(D%*%T)))
ggsave('caterpillar_climate_interact_weighted_std.png', p, 
       width=plot_width*1.5, height=plot_height, dpi=plot_dpi)
ggsave('caterpillar_climate_interact_weighted_std.pdf', p, 
       width=plot_width*1.5, height=plot_height, dpi=plot_dpi)

p <- caterpillar(filter(B_g_betas_std, param %in% paste0('B_g_', c(1:5, 7:9), '_median')),
                 labels=c('B_g_1_median'='Intercept',
                          'B_g_2_median'='P',
                          'B_g_3_median'=expression(P^2),
                          'B_g_4_median'='T',
                          'B_g_5_median'=expression(T^2),
                          'B_g_7_median'=expression(D^2),
                          'B_g_8_median'=expression(D%*%P),
                          'B_g_9_median'=expression(D%*%T)))
ggsave('caterpillar_climate_interact_weighted_std_noD.png', p, 
       width=plot_width*1.5, height=plot_height, dpi=plot_dpi)
ggsave('caterpillar_climate_interact_weighted_std_noD.pdf', p, 
       width=plot_width*1.5, height=plot_height, dpi=plot_dpi)

# Save table of parameter estimates for use in Word
cis <- group_by(B_g_betas, model) %>%
    rename(Parameter=param) %>%
    do(ci(.)) %>%
    rename(param=Parameter)

# Make a string that I can use in the table in Word
cis$estimate <- paste0(signif(cis$median, 4), ' (',
                       signif(cis$Low, 4), ', ',
                       signif(cis$High, 4), ')')

coef_table <- dcast(cis, param ~ model, value.var='estimate')
write.csv(coef_table, file='parameter_estimates.csv', row.names=FALSE)

###############################################################################
# Make combined plot of coefs for publication
labels <- c('B_g_1_median'='Intercept',
            'B_g_2_median'='P',
            'B_g_3_median'=expression(P^2),
            'B_g_4_median'='T',
            'B_g_5_median'=expression(T^2),
            'B_g_6_median'='D',
            'B_g_7_median'=expression(D^2),
            'B_g_8_median'=expression(D%*%P),
            'B_g_9_median'=expression(D%*%T))
mods <- filter(B_g_betas, param %in% paste0('B_g_', c(1:9), '_median'))
cis <- group_by(mods, model) %>%
    rename(Parameter=param) %>%
    do(ci(.)) %>%
    rename(param=Parameter)
cis$which_facet <- 'main'
cis$which_facet[cis$param == 'B_g_1_median'] <- 'intercept'
cis$which_facet[cis$param == 'B_g_6_median'] <- 'D'
cis$which_facet <- ordered(cis$which_facet, levels=c('main', 'intercept', 'D'))
p <- ggplot(cis, aes(x=param,
                     y=median, colour=model, fill=model)) +
    theme_bw(base_size=8) +
    geom_point(aes(shape=model), position=position_dodge(width=.4), size=1.25) +
    geom_linerange(aes(ymin=Low, ymax=High), size=.5, position=position_dodge(width=.4)) +
    geom_linerange(aes(ymin=low, ymax=high), size=.25, position=position_dodge(width=.4)) +
    xlab('Parameter') +
    ylab('Effect (cm)') +
    scale_shape_manual("Model", values=c(23, 4, 17)) + 
    scale_colour_brewer("Model", palette = "Dark2") + 
    scale_fill_brewer("Model", palette = "Dark2") +
    scale_x_discrete(labels=labels) +
    theme(axis.text.x=element_text(vjust=0),
          strip.background=element_blank(),
          strip.text=element_blank(),
          panel.grid.minor.x=element_blank(),
          panel.grid.major.x=element_blank()) +
    facet_wrap(~which_facet, scales='free', drop=TRUE)
ggsave('caterpillar_climate_interact_weighted_forpub.png', p, width=plot_width*2, 
       height=plot_height, dpi=plot_dpi)
ggsave('caterpillar_climate_interact_weighted_forpub.pdf', p, width=plot_width*2, 
       height=plot_height, dpi=plot_dpi)
ggsave('caterpillar_climate_interact_weighted_forpub.svg', p, width=plot_width*2, 
       height=plot_height, dpi=plot_dpi)

