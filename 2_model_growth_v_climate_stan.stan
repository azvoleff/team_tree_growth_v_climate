/*
    TODO:

    1) account for varying period length with a predictor
    2) include random effect for period
    3) include random effect for genus
*/

data {
    int<lower=0> n_dbhs;
    int<lower=0> n_tree;
    int<lower=0> n_plot;
    int<lower=0> n_site;
    int<lower=0, upper=13> n_obs_per_tree[n_tree];
    int<lower=0> n_growths;
    //int<lower=0> n_period;
    //int<lower=0> n_genus;
    int<lower=0, upper=n_tree> tree_ID[n_dbhs];
    int<lower=0, upper=n_plot> plot_ID[n_dbhs];
    int<lower=0, upper=n_site> site_ID[n_dbhs];
    //int<lower=1, upper=n_period> period_ID;
    //int<lower=1, upper=n_genus> genus_ID;
    vector[n_dbhs] log_dbh;
}

parameters {
    real<lower=-3, upper=3> beta_0;
    real<lower=-3, upper=3> beta_1;
    real<lower=0> sigma_obs;
    real<lower=0> sigma_proc;
    real<lower=0> sigma_ijk;
    real<lower=0> sigma_jk;
    real<lower=0> sigma_k;
    real log_dbh_latent[n_dbhs];
    vector[n_tree] b_ijk;
    vector[n_plot] b_jk;
    vector[n_site] b_k;
    //vector[n_period] b_t;
    //vector[n_genus] b_g;
}

transformed parameters {
    vector[n_growths] log_dbh_latent_st;
    vector[n_growths] growth;
    vector[n_growths] log_growth_hat;
    // to use obs_index and growths_index below, need to set them up to be 
    // local variables, by enclosing w/in brackets
    {
        // obs_index tracks the current observation
        int obs_index;
        obs_index <- 1;
        // growths_index tracks position in the growths vector (shorter than 
        // the observations vector since growth can't be calculated for the 
        // first time period).
        int growths_index;
        growths_index <- 1;
        for (tree_num in 1:n_tree) {
            obs_index <- obs_index + 1;
            for (obs_num in 2:n_obs_per_tree[tree_num]) {
                growth[growths_index] <- exp(log_dbh_latent[obs_index]) - exp(log_dbh_latent[obs_index - 1]);
                // save log_dbh_latent_st since we need a vector of starting 
                // dbh values for use in the likelihood that is of the same 
                // length as growths
                log_dbh_latent_st[growths_index] <- log_dbh_latent[obs_index - 1];
                obs_index <- obs_index + 1;
                growths_index <- growths_index + 1;
            }
        }
    }
    
    for (i in 1:n_growths)
        log_growth_hat[i] <- beta_0 + beta_1 * exp(log_dbh_latent_st[i]) + b_ijk[tree_ID[i]] + b_jk[plot_ID[i]] + b_k[site_ID[i]];

    /* print("growth=", growth[1], */
    /*       ", log_growth_hat=", log_growth_hat[1], */
    /*       ", dbh_latent_st=", log_dbh_latent[2 - 1], */
    /*       ", dbh_latent_end=", log_dbh_latent[2], */
    /*       ", log_dbh_st=", log_dbh[2 - 1], */
    /*       ", log_dbh_end=", log_dbh[2]) */
}

model {
    log_dbh ~ normal(log_dbh_latent, sigma_obs);
    sigma_obs ~ cauchy(0, 1);
    log(growth) ~ normal(log_growth_hat, sigma_proc);
    sigma_proc ~ cauchy(0, 1);
    // Jacobian adjustment to account for log transform of latent growth 
    // variables (log absolute determinant of transform):
    for (i in 1:n_growths)
        increment_log_prob(-log(fabs(growth[i])));
    b_ijk ~ normal(0, sigma_ijk);
    sigma_ijk ~ cauchy(0, 1);
    b_jk ~ normal(0, sigma_jk);
    sigma_jk ~ cauchy(0, 1);
    b_k ~ normal(0, sigma_k);
    sigma_k ~ cauchy(0, 1);
}
