
# Plot Rhats for each model
foreach(this_model=unique(params$model)) %do% {
    pars <- filter(params, model == this_model)
    attributes(pars)$nChains <- max(params$chain)
    attributes(pars)$nIterations <- max(params$iteration)
    attributes(pars)$nThin <- nThin
    attributes(pars)$nBurnin <- nBurnin

    ggs_Rhat(pars, 'B_k')
    ggsave(paste0('Rhat_interact_B_k_', this_model, '.png'), width=plot_width*2, 
           height=plot_height*2, dpi=plot_dpi)

    ggs_traceplot(pars, 'mu_B_g')
    ggsave(paste0('traceplot_interact_mu_B_g_', this_model, '.png'), width=plot_width*2, 
           height=plot_height*4, dpi=plot_dpi)
    ggs_Rhat(pars, 'mu_B_g')
    ggsave(paste0('Rhat_interact_mu_B_g_', this_model, '.png'), width=plot_width*2, 
           height=plot_height*2, dpi=plot_dpi)

    ggs_traceplot(pars, 'sigma_')
    ggsave(paste0('traceplot_interact_sigma_', this_model, '.png'), width=plot_width*2, 
           height=plot_height*6, dpi=plot_dpi)
    ggs_Rhat(pars, 'sigma_')
    ggsave(paste0('Rhat_interact_sigma_', this_model, '.png'), width=plot_width*2, 
           height=plot_height*2, dpi=plot_dpi)

    # ggs_traceplot(pars, 'int_')
    # ggsave(paste0('traceplot_interact_int_', this_model, '.png'), width=plot_width*2, 
    #        height=plot_height*3, dpi=plot_dpi)
    ggs_Rhat(pars, 'int_')
    ggsave(paste0('Rhat_interact_int_', this_model, '.png'), width=plot_width*2, 
           height=plot_height*3, dpi=plot_dpi)
}

