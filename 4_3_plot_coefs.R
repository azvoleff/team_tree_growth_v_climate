# Forest growth model MCMC results

prefixes <- c('D:/azvoleff/Data', # CI-TEAM
              'H:/Data', # Buffalo drive
              'O:/Data', # Blue drive
              '/localdisk/home/azvoleff/Data') # vertica1
prefix <- prefixes[match(TRUE, unlist(lapply(prefixes, function(x) file_test('-d', x))))]

data_folder <- file.path(prefix, "TEAM", "Tree_Growth", "Data")

library(ggmcmc)
library(foreach)
library(ggthemes)

plot_width <- 3.5
plot_height <- 3
plot_dpi <- 300

multimodel_caterpillar <- function(mods, labels=NULL) {
    cis <- group_by(mods, Model) %>%
        do(calc_cis(.))
    p <- ggplot(cis, aes(x=median, y=Parameter, colour=Model)) +
        theme_bw(base_size=8) +
        geom_point(position=position_dodge(width=.4), size=1.25) +
        geom_linerange(aes(ymin=Low, ymax=High), size=.75, position=position_dodge(width=.4)) +
        geom_linerange(aes(ymin=low, ymax=high), size=.25, position=position_dodge(width=.4)) +
        ylab('Effect (cm)')

    if (!is.null(labels)) {
        p <- p + scale_x_discrete(labels=labels) +
            theme(axis.text.x=element_text(vjust=0))
    }
    p
}

# Function to calculated weighted coefficients from an array of MCMC results.  
# Weights should be a 2 column data.frame with IDs in the first column, and 
# weights in the second. D should be a ggs object with parameter IDs added 
# using the jags_param_ids function
weight_coef <- function(d, w) {
    d <- left_join(d, w)
    d <- group_by(d, Model, Chain, Iteration, param_ID) %>%
        summarise(Parameter=paste(Parameter_Base[1], param_ID[1], 'mean', sep='_'),
                  value=mean(weight * value))
    return(d)
}

###############################################################################
## Simple model (no random effects)
##

simple_labels_B <- data.frame(Parameter=paste0('B[', 1:10, ']'),
    Label=c("int",
            "P" ,
            "P^2",
            "T",
            "T^2",
            "D",
            "D^2"))

load("smp_mods.RData")

p <- multimodel_caterpillar(filter(simple_B, Parameter %in% c('P', 'P^2', 'T', 'T^2')),
                            labels=c('P'='P',
                                     'P^2'=expression(P^2), 
                                     'T'='T',
                                     'T^2'=expression(T^2)))
ggsave('simple_caterpillar_climate.png', p, width=plot_width, 
       dpi=plot_dpi)

###############################################################################
## Full model (interactions, uncorrelated random effects)

interact_labels <- data.frame(Label=c("int",
            "P",
            "P^2",
            "T",
            "T^2",
            "D",
            "D^2",
            "D*P",
            "D*T"))
interact_labels$Parameter <- paste0('mu_B_g[', 1:nrow(interact_labels), ']')

load("int_mods.RData")

# Calculate weights for each genus ID (doesn't matter which temp_var is used 
# since the genus IDs and frequencies are the same across all simulations).
load(file.path(data_folder, paste0("model_data_wide_full-tmn_meanannual-mcwd_run12.RData")))
merged <- tbl_df(data.frame(site_ID=model_data$site_ID, genus_ID=as.integer(model_data$genus_ID)))
genus_weights <- group_by(merged, genus_ID) %>%
    summarize(n=n()) %>%
    ungroup() %>%
    mutate(weight=n/sum(n)) %>%
    select(-n) %>%
    arrange(desc(weight))

g_betas <- weight_coef(filter(int_mods, Parameter_Base == 'B_g'), genus_weights)

p <- multimodel_caterpillar(clim_betas,
                            labels=c('mu_B_g[2]'='P',
                                     'mu_B_g[3]'=expression(P^2), 
                                     'mu_B_g[4]'='T',
                                     'mu_B_g[5]'=expression(T^2)))
ggsave('interact_caterpillar_climate.png', p, width=plot_width, height=plot_height, 
       dpi=plot_dpi)

# mu_B_g_vals <- ggs(as.mcmc.list(jags_fit, "mu_B_g"), par_labels=interact_labels_mu_B_g)
# ggs_caterpillar(mu_B_g_vals)
# ggs_caterpillar(mu_B_g_vals) + theme_tufte()
# ggs_caterpillar(mu_B_g_vals) + theme_solarized_2()

###############################################################################
## Full model (correlated random effects)

correlated_labels_mu_B_g <- data.frame(Parameter=paste0('mu_B_g[', 1:8, ']'),
    Label=c("int",
            "P",
            "P^2",
            "T",
            "T^2",
            "D",
            "D^2",
            "E"))

load("cor_mods.RData")
