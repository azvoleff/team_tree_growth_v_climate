/*
    TODO:

    1) account for varying period length with a predictor
    2) include random effect for period
    3) include random effect for genus
*/

data {
    int<lower=0> n_tree;
    int<lower=0> n_plot;
    int<lower=0> n_site;
    int<lower=0> max_obs_per_tree;
    int<lower=0> first_obs_period[n_tree];
    int<lower=0> last_obs_period[n_tree];
    //int<lower=0> n_period;
    //int<lower=0> n_genus;
    int<lower=0, upper=n_plot> plot_ID[n_tree];
    int<lower=0, upper=n_site> site_ID[n_tree];
    //int<lower=1, upper=n_period> period_ID;
    //int<lower=1, upper=n_genus> genus_ID;
    matrix[n_tree, max_obs_per_tree] log_dbh;
    vector[n_tree] WD;
    matrix[n_tree, max_obs_per_tree] spi;
}

parameters {
    matrix<lower=2, upper=5.6>[n_tree, max_obs_per_tree] log_dbh_latent;
    real<lower=-3, upper=3> inter;
    real<lower=-3, upper=3> slp_dbh;
    real<lower=-3, upper=3> slp_WD;
    real<lower=-3, upper=3> slp_WD_sq;
    real<lower=-3, upper=3> slp_spi;
    real<lower=-3, upper=3> slp_spi_sq;
    real<lower=0, upper=100> sigma_obs;
    real<lower=0, upper=100> sigma_proc;
    real<lower=0, upper=100> sigma_ijk;
    real<lower=0, upper=100> sigma_jk;
    real<lower=0, upper=100> sigma_k;
    vector[n_tree] b_ijk_std; // standardized ranef for tree
    vector[n_plot] b_jk_std; // standardized ranef for plot
    vector[n_site] b_k_std; // standardized ranef for site
    //vector[n_period] b_t;
    //vector[n_genus] b_g;
}

transformed parameters {
    vector[n_tree] b_ijk;
    vector[n_plot] b_jk;
    vector[n_site] b_k;

    // Matt trick- see http://bit.ly/1qz4NC6
    b_ijk <- sigma_ijk * b_ijk_std; // b_ijk ~ normal(0, sigma_ijk)
    b_jk <- sigma_jk * b_jk_std; // b_jk ~ normal(0, sigma_jk)
    b_k <- sigma_k * b_k_std; // b_k ~ normal(0, sigma_k)
    
    /* print("dbh_latent_st=", log_dbh_latent[2 - 1], */
    /*       ", dbh_latent_end=", log_dbh_latent[2], */
    /*       ", log_dbh_st=", log_dbh[2 - 1], */
    /*       ", log_dbh_end=", log_dbh[2]) */
}

model {
    matrix[n_tree, max_obs_per_tree] log_growth_hat;
    matrix<lower=0>[n_tree, max_obs_per_tree] growth;
    //growth <- rep_matrix(0, n_tree, max_obs_per_tree);
    //log_growth_hat <- rep_matrix(0, n_tree, max_obs_per_tree);
    for (tree_num in 1:n_tree) {
        // Make sure first latent dbh is also modeled
        log_dbh[tree_num, first_obs_period[tree_num]] ~ normal(log_dbh_latent[tree_num, first_obs_period[tree_num]], sigma_obs);
        for (obs_num in (first_obs_period[tree_num] + 1):last_obs_period[tree_num]) {
            // Distribution of dbh, with  mean latent dbh
            log_dbh[tree_num, obs_num] ~ normal(log_dbh_latent[tree_num, obs_num], sigma_obs);

            // Growth model
            growth[tree_num, obs_num] <- exp(log_dbh_latent[tree_num, obs_num]) -
                exp(log_dbh_latent[tree_num, obs_num - 1]);
            log_growth_hat[tree_num, obs_num] <- inter +
                slp_dbh * exp(log_dbh_latent[tree_num, obs_num - 1]) +
                slp_WD * WD[tree_num] +
                slp_WD_sq * pow(WD[tree_num], 2) +
                slp_spi * spi[tree_num, obs_num] +
                slp_spi_sq * pow(spi[tree_num, obs_num], 2) +
                b_ijk[tree_num] + b_jk[plot_ID[tree_num]] + b_k[site_ID[tree_num]];
            log(growth[tree_num, obs_num]) ~ normal(log_growth_hat[tree_num, obs_num], sigma_proc);

            // Jacobian adjustment to account for log transform of latent growth 
            // variables (log absolute determinant of transform):
            increment_log_prob(-log(fabs(growth[tree_num, obs_num])));
        }
    }

    //sigma_obs ~ cauchy(0, 1);
    //sigma_proc ~ cauchy(0, 1);

    // Random effects
    b_ijk_std ~ normal(0, 1); // Matt trick
    //sigma_ijk ~ cauchy(0, 1);
    b_jk_std ~ normal(0, 1); // Matt trick
    //sigma_jk ~ cauchy(0, 1);
    b_k_std ~ normal(0, 1); // Matt trick
    //sigma_k ~ cauchy(0, 1);
}
