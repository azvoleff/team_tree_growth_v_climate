# Forest growth model MCMC results

prefixes <- c('D:/azvoleff/Data', # CI-TEAM
              'H:/Data', # Buffalo drive
              'O:/Data', # Blue drive
              '/localdisk/home/azvoleff/Data') # vertica1
prefix <- prefixes[match(TRUE, unlist(lapply(prefixes, function(x) file_test('-d', x))))]

library(mcmcplots)
library(runjags)
library(ggmcmc)
library(foreach)
library(ggthemes)
library(coda)

plot_width <- 3.5
plot_height <- 3
plot_dpi <- 300

start_val <- 30000
thin_val <- 50

plot_mcmc <- function(fit, params, title) {
    mcmc_results <- as.mcmc.list(fit, c(params))
    mcmc_results <- window(mcmc_results, start=start_val, thin=thin_val)
    mcmcplot(mcmc_results, title=paste0(title))
    #caterplot(mcmc_results)
}

multimodel_caterpillar <- function(mods, labels=NULL) {
    cis <- foreach (mod=mods, name=names(mods), .combine=rbind) %do% {
        ret <- ci(mod)
        ret$Model <- name
        ret
    }
    p <- ggplot(cis, aes(x=Parameter, y=median, colour=Model)) +
        theme_bw(base_size=8) +
        geom_point(position=position_dodge(width=.4), size=1.25) +
        geom_linerange(aes(ymin=Low, ymax=High), size=.75, position=position_dodge(width=.4)) +
        geom_linerange(aes(ymin=low, ymax=high), size=.25, position=position_dodge(width=.4)) +
        ylab('Value (standardized)')
    if (!is.null(labels)) {
        p <- p + scale_x_discrete(labels=labels) +
            theme(axis.text.x=element_text(vjust=0))
    }
    p
}

###############################################################################
## Simple model (no random effects)
##

plot_simple_model <- function(title) {
    plot_mcmc(jags_fit, c("B", "sigma_", "int_"), title)
}

simple_labels_<- data.frame(Parameter=paste0('B[', 1:10, ']'),
    Label=c("int",
            "W",
            "W^2",
            "P",
            "P^2",
            "T",
            "T^2",
            "D",
            "D^2",
            "E"))

### Minimum temperature
load(file.path(prefix, "TEAM", "Tree_Growth", "MCMC_Chains", "simple", 
               "simple_tmn.RData"))
plot_simple_model('simple min temp')
simple_tmn_B_vals <- ggs(as.mcmc.list(jags_fit, "B"), par_labels=simple_labels_B)

### Mean temperature
load(file.path(prefix, "TEAM", "Tree_Growth", "MCMC_Chains", "simple", 
               "simple_tmp.RData"))
plot_simple_model('simple mean temp')
simple_tmp_B_vals <- ggs(as.mcmc.list(jags_fit, "B"), par_labels=simple_labels_B)

### Maximum temperature
load(file.path(prefix, "TEAM", "Tree_Growth", "MCMC_Chains", "simple", 
               "simple_tmx.RData"))
plot_simple_model('simple max temp')
simple_tmx_B_vals <- ggs(as.mcmc.list(jags_fit, "B"), par_labels=simple_labels_B)

simple_climate <- list("Min temp"=filter(simple_tmn_B_vals, Parameter %in% c('P', 'P^2', 'T', 'T^2')),
                       "Mean temp"=filter(simple_tmp_B_vals, Parameter %in% c('P', 'P^2', 'T', 'T^2')),
                       "Max temp"=filter(simple_tmx_B_vals, Parameter %in% c('P', 'P^2', 'T', 'T^2')))
p <- multimodel_caterpillar(simple_climate,
                            labels=c('P'='P',
                                     'P^2'=expression(P^2), 
                                     'T'='T',
                                     'T^2'=expression(T^2)))
ggsave('simple_caterpillar_climate.png', p, width=plot_width, height=plot_height, 
       dpi=plot_dpi)

###############################################################################
## Model with only sites with large elevation gradients (UDZ, BIF, VB only)

# plot_small_model <- function(title) {
#     plot_mcmc(jags_fit, c("B[1]", "B[2]"), title)
#     plot_mcmc(jags_fit, "mu_B_g", title)
#     plot_mcmc(jags_fit, "sigma_", title)
# }
#
# ### Mean temperature
# load(file.path(prefix, "TEAM", "Tree_Growth", "MCMC_Chains", "highelev", 
#                "highelev_tmp.RData"))
# plot_small_model("small mean temp")
#
# ### Minimum temperature
# load(file.path(prefix, "TEAM", "Tree_Growth", "MCMC_Chains", "highelev", 
#                "highelev_tmn.RData"))
# plot_small_model("small min temp")
#
# ### Maximum temperature
# load(file.path(prefix, "TEAM", "Tree_Growth", "MCMC_Chains", "highelev", 
#                "highelev_tmx.RData"))
# plot_small_model("small max temp")

###############################################################################
## Full model (interactions, uncorrelated random effects)

plot_full_interact_model <- function(title) {
    plot_mcmc(jags_fit, c("B[1]", "B[2]", "mu_B_g", "sigma_", "int_"), title)
}

interact_labels_mu_B_g <- data.frame(Parameter=paste0('mu_B_g[', 1:12, ']'),
    Label=c("int",
            "P",
            "P^2",
            "T",
            "T^2",
            "D",
            "D^2",
            "W*P",
            "D*P",
            "W*T",
            "D*T",
            "E"))

### Min temperature
load(file.path(prefix, "TEAM", "Tree_Growth", "MCMC_Chains", "interact", 
               "interact_tmn.RData"))
plot_full_interact_model("interact min temp")
interact_tmn_B_g_vals <- ggs(as.mcmc.list(jags_fit, "mu_B_g"), par_labels=interact_labels_mu_B_g)

### Mean temperature
load(file.path(prefix, "TEAM", "Tree_Growth", "MCMC_Chains", "interact", 
               "interact_tmp.RData"))
plot_full_interact_model("interact mean temp")
interact_tmp_B_g_vals <- ggs(as.mcmc.list(jags_fit, "mu_B_g"), par_labels=interact_labels_mu_B_g)

### Maximum temperature
load(file.path(prefix, "TEAM", "Tree_Growth", "MCMC_Chains", "interact", 
               "interact_tmx.RData"))
plot_full_interact_model("interact max temp")
interact_tmx_B_g_vals <- ggs(as.mcmc.list(jags_fit, "mu_B_g"), par_labels=interact_labels_mu_B_g)

interact_climate <- list("Min temp"=filter(interact_tmn_B_g_vals, Parameter %in% c('P', 'P^2', 'T', 'T^2')),
                       "Mean temp"=filter(interact_tmp_B_g_vals, Parameter %in% c('P', 'P^2', 'T', 'T^2')),
                       "Max temp"=filter(interact_tmx_B_g_vals, Parameter %in% c('P', 'P^2', 'T', 'T^2')))
p <- multimodel_caterpillar(interact_climate,
                            labels=c('P'='P',
                                     'P^2'=expression(P^2), 
                                     'T'='T',
                                     'T^2'=expression(T^2)))
ggsave('interact_caterpillar_climate.png', p, width=plot_width, height=plot_height, 
       dpi=plot_dpi)

# mu_B_g_vals <- ggs(as.mcmc.list(jags_fit, "mu_B_g"), par_labels=interact_labels_mu_B_g)
# ggs_caterpillar(mu_B_g_vals)
# ggs_caterpillar(mu_B_g_vals) + theme_tufte()
# ggs_caterpillar(mu_B_g_vals) + theme_solarized_2()

###############################################################################
## Full model (correlated random effects)
start_val <- 20000
thin_val <- 10
plot_full_correlated_model <- function(title) {
    plot_mcmc(jags_fit, c("B[1]", "B[2]", "mu_B_g", "sigma_", "int_"), title)
}

correlated_labels_mu_B_g <- data.frame(Parameter=paste0('mu_B_g[', 1:8, ']'),
    Label=c("int",
            "P",
            "P^2",
            "T",
            "T^2",
            "D",
            "D^2",
            "E"))

### Min temperature
load(file.path(prefix, "TEAM", "Tree_Growth", "MCMC_Chains", "correlated", 
               "correlated_tmn.RData"))
plot_full_correlated_model("correlated min temp")
correlated_tmn_B_g_vals <- ggs(as.mcmc.list(jags_fit, "mu_B_g"), par_labels=correlated_labels_mu_B_g)

### Mean temperature
load(file.path(prefix, "TEAM", "Tree_Growth", "MCMC_Chains", "correlated", 
               "correlated_tmp.RData"))
plot_full_correlated_model("correlated mean temp")
correlated_tmp_B_g_vals <- ggs(as.mcmc.list(jags_fit, "mu_B_g"), par_labels=correlated_labels_mu_B_g)

### Maximum temperature
load(file.path(prefix, "TEAM", "Tree_Growth", "MCMC_Chains", "correlated", 
               "correlated_tmx.RData"))
plot_full_correlated_model("correlated max temp")
correlated_tmx_B_g_vals <- ggs(as.mcmc.list(jags_fit, "mu_B_g"), par_labels=correlated_labels_mu_B_g)

correlated_climate <- list("Min temp"=filter(correlated_tmn_B_g_vals, Parameter %in% c('P', 'P^2', 'T', 'T^2')),
                           "Mean temp"=filter(correlated_tmp_B_g_vals, Parameter %in% c('P', 'P^2', 'T', 'T^2')),
                           "Max temp"=filter(correlated_tmx_B_g_vals, Parameter %in% c('P', 'P^2', 'T', 'T^2')))
p <- multimodel_caterpillar(correlated_climate,
                            labels=c('P'='P',
                                     'P^2'=expression(P^2), 
                                     'T'='T',
                                     'T^2'=expression(T^2)))
ggsave('correlated_caterpillar_climate.png', p, width=plot_width, height=plot_height, 
       dpi=plot_dpi)

