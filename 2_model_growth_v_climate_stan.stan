data {
    int<lower=0> n_obs;
    int<lower=0> n_tree;
    int<lower=0> n_plot;
    int<lower=0> n_site;
    int<lower=0, upper=13> n_per_tree[n_tree];
    int<lower=0> n_growths;
    //int<lower=0> n_period;
    int<lower=0, upper=n_tree> tree_ID[n_obs];
    int<lower=0, upper=n_plot> plot_ID[n_obs];
    int<lower=0, upper=n_site> site_ID[n_obs];
    //int<lower=1, upper=n_period> period_ID;
    vector[n_obs] log_dbh;
}

parameters {
    real beta_0;
    real beta_1;
    real<lower=0> sigma;
    real<lower=0> sigma_ijk;
    real<lower=0> sigma_jk;
    real<lower=0> sigma_k;
    real log_dbh_latent[n_obs];
    real<lower=0> w;
    vector[n_tree] b_ijk;
    vector[n_plot] b_jk;
    vector[n_site] b_k;
    //vector[n_period] b_t;
}

transformed parameters {
    vector[n_growths] log_dbh_latent_st;
    vector[n_growths] growth;
    vector[n_growths] growth_hat;
    real<lower=0> diff; // constrain growth to be positive
    /* real a7; */
    /* real a8; */
    /* a7 <- 10 * n_obs; // Clark (2007), pg. 1947 */
    /* a8 <- .001 * (a7 - 1); // Clark (2007), pg. 1947 */
    // to use obs_index and growths_index below, need to set them up to be 
    // local variables, by enclosing w/in brackets
    {
        int obs_index;
        // Below tracks position in the growths vector (shorter than the 
        // observations vector since growth can't be calculated for the first 
        // time period).
        int growths_index;
        obs_index <- 1;
        growths_index <- 1;
        for (i in 1:n_tree) {
            obs_index <- obs_index + 1;
            for (j in 2:n_per_tree[i]) {
                diff <- exp(log_dbh_latent[obs_index]) - exp(log_dbh_latent[obs_index - 1]);
                growth[growths_index] <- log_diff_exp(log_dbh_latent[obs_index], log_dbh_latent[obs_index - 1]);
                // save log_dbh_latent_st since we need a vector of starting 
                // dbh values for use in the likelihood that is of the same 
                // length as growths
                log_dbh_latent_st[growths_index] <- log_dbh_latent[obs_index - 1];
                obs_index <- obs_index + 1;
                growths_index <- growths_index + 1;
            }
        }
    }

    // Also need to account for varying period length with a predictor, and 
    // include random effect for period
    
    for (i in 1:n_growths)
        growth_hat[i] <- beta_0 + beta_1 * exp(log_dbh_latent_st[i]) + b_ijk[tree_ID[i]] + b_jk[plot_ID[i]] + b_k[site_ID[i]];

    /* print("growth=", growth[1], */
    /*       ", growth_hat=", growth_hat[1], */
    /*       ", dbh_latent_st=", log_dbh_latent[2 - 1], */
    /*       ", dbh_latent_end=", log_dbh_latent[2], */
    /*       ", log_dbh_st=", log_dbh[2 - 1], */
    /*       ", log_dbh_end=", log_dbh[2]) */
}

model {
    log_dbh ~ normal(log_dbh_latent, w);
    //w ~ inv_gamma(a7, a8);
    w ~ inv_gamma(.001, .001);
    // TODO: Fix below increment_log_prob statement
    growth ~ normal(growth_hat, sigma);
    // Account for log transform of latent growth variables:
    //increment_log_prob(-log(growth));
    //beta_0 ~ normal(0, 100);
    //beta_1 ~ normal(0, 100);
    b_ijk ~ normal(0, sigma_ijk);
    b_jk ~ normal(0, sigma_jk);
    b_k ~ normal(0, sigma_k);
    sigma ~ inv_gamma(.001, .001);
    sigma_ijk ~ inv_gamma(.001, .001);
    sigma_jk ~ inv_gamma(.001, .001);
    sigma_k ~ inv_gamma(.001, .001);
}
